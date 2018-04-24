/************************************************************************
*                                                                       *
* TAPSim - an atom probe data simulation program                        *
*                                                                       *
* Copyright (C) 2011 Christian Oberdorfer                               *
*                                                                       *
* This program is free software: you can redistribute it and/or modify  *
* it under the terms of the GNU General Public License as published by  *
* the Free Software Foundation, either version 3 of the License, or any *
* any later version.                                                    *
*                                                                       *
* This program is distributed in the hope that it will be useful,       *
* but WITHOUT ANY WARRANTY; without even the implied warranty of        *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
* GNU General Public License for more details.                          *
*                                                                       *
* You should have received a copy of the GNU General Public License     *
* along with this program.  If not, see 'http://www.gnu.org/licenses'   *
*                                                                       *
************************************************************************/ 

#include "grid_3d.h"

#include <limits>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include <algorithm>
#include <deque>
#include <stack>

#include <ctime>

#include <unistd.h>

#include "../vector/vector.h"

#include "debug.h"

#include "utils.h"

//#define DEBUG_MODE
//#undef DEBUG_MODE

#ifdef DEBUG_MODE
#include "info.h"
#endif

// ***** ===== GRID_3D::TABLE::NODEDATA ===== *****

Grid_3d::Node::Node()
	:
	  _id(),
	  _number(),
	  _charge(0.0f),
	  _properties(0),
	  _numNeighbours(0),
	  _neighbours(0),
	  _couplingValues(0)
{
	for (int index = 0; index < 2; index++)
		_phi[index] = 0.0f;
}

Grid_3d::Node::~Node()
{
	if (_neighbours != 0) delete[] _neighbours;
	if (_couplingValues != 0) delete[] _couplingValues;
}

void Grid_3d::Node::setBoundary(const bool value)
{
	if (value)
		_properties |= BOUNDARY;
	else
		_properties &= ~BOUNDARY;
}

void Grid_3d::Node::setNeumannBoundary(const bool value)
{
	if (value)
		_properties |= NEUMANN;
	else
		_properties &= ~NEUMANN;
}

void Grid_3d::Node::setDirichletBoundary(const bool value)
{
	if (value)
		_properties |= DIRICHLET;
	else
		_properties &= ~DIRICHLET;
}

void Grid_3d::Node::setNumNeighbours(const int value)
{
	if (value > 0)
	{
		if (value != _numNeighbours)
		{
			_numNeighbours = value;

			if (_neighbours != 0) delete[] _neighbours;
			_neighbours = new int[_numNeighbours];

			if (_couplingValues != 0) delete[] _couplingValues;
			_couplingValues = new float[_numNeighbours];
		}

		for (int index = 0; index < _numNeighbours; index++)
		{
			_neighbours[index] = -1;
			_couplingValues[index] = 1.0f / _numNeighbours;
		}
	}
	else
	{
		_numNeighbours = 0;

		if (_neighbours != 0)
		{
			delete[] _neighbours;
			_neighbours = 0;
		}

		if (_couplingValues != 0)
		{
			delete[] _couplingValues;
			_couplingValues = 0;
		}
	}
}

// ***** ===== GRID_3D::TABLE ===== *****

const std::string Grid_3d::Table::_binaryVersion("1.0");

Grid_3d::Table::Table(const unsigned int threadNum)
	: _threadNum(threadNum),
	  _threadedRelax_ids(),
	  _threadedRelax_params(this),
	  _threadedLocalRelax_ids(),
	  _threadedLocalRelax_params(this),
	  _phiSwitch(0),
	  _phiLocalSwitch(0)
{
	if (_threadNum == 0) _threadNum = sysconf(_SC_NPROCESSORS_ONLN);

	// ***

	_threadedRelax_ids = new pthread_t[_threadNum];

	for (unsigned int i = 0; i < _threadNum; i++)
		pthread_create(&_threadedRelax_ids[i],0,cwrapper_threadedRelax,&_threadedRelax_params);

	// ***

	_threadedLocalRelax_ids = new pthread_t[_threadNum];

	for (unsigned int i = 0; i < _threadNum; i++)
		pthread_create(&_threadedLocalRelax_ids[i],0,cwrapper_threadedLocalRelax,&_threadedLocalRelax_params);
}

Grid_3d::Table::~Table()
{
	_threadedRelax_params.lock();
	_threadedRelax_params.keepRunning = false;

	pthread_cond_broadcast(&_threadedRelax_params.runCondition);

	_threadedRelax_params.unlock();

	for (unsigned int i = 0; i < _threadNum; i++)
		pthread_join(_threadedRelax_ids[i],0);

	delete[] _threadedRelax_ids;

	// ***

	_threadedLocalRelax_params.lock();
	_threadedLocalRelax_params.keepRunning = false;

	pthread_cond_broadcast(&_threadedLocalRelax_params.runCondition);

	_threadedLocalRelax_params.unlock();

	for (unsigned int i = 0; i < _threadNum; i++)
		pthread_join(_threadedLocalRelax_ids[i],0);

	delete[] _threadedLocalRelax_ids;
}

void Grid_3d::Table::setThreadNum(const unsigned int value)
{
	// *** stop running threads
	
	_threadedRelax_params.lock();
	_threadedRelax_params.keepRunning = false;

	pthread_cond_broadcast(&_threadedRelax_params.runCondition);

	_threadedRelax_params.unlock();

	// ***
	
	for (unsigned int i = 0; i < _threadNum; i++)
		pthread_join(_threadedRelax_ids[i],0);

	delete[] _threadedRelax_ids;

	// ***

	_threadedLocalRelax_params.lock();
	_threadedLocalRelax_params.keepRunning = false;

	pthread_cond_broadcast(&_threadedLocalRelax_params.runCondition);

	_threadedLocalRelax_params.unlock();

	for (unsigned int i = 0; i < _threadNum; i++)
		pthread_join(_threadedLocalRelax_ids[i],0);

	delete[] _threadedLocalRelax_ids;
	
	// *** restart threads

	_threadNum = (value == 0 ? sysconf(_SC_NPROCESSORS_ONLN) : value); 

	// ***
	
	_threadedRelax_ids = new pthread_t[_threadNum];
	_threadedRelax_params.keepRunning = true;

	for (unsigned int i = 0; i < _threadNum; i++)
		pthread_create(&_threadedRelax_ids[i],0,cwrapper_threadedRelax,&_threadedRelax_params);

	// ***
	
	_threadedLocalRelax_ids = new pthread_t[_threadNum];
	_threadedLocalRelax_params.keepRunning = true;

	for (unsigned int i = 0; i < _threadNum; i++)
		pthread_create(&_threadedLocalRelax_ids[i],0,cwrapper_threadedLocalRelax,&_threadedLocalRelax_params);
}


void Grid_3d::Table::operator<<(std::ifstream& stream)
{
	std::string version;
	{
		char c;
		stream.read(reinterpret_cast<char*>(&c),sizeof(char));


		while (c != '\0')
		{
			version.push_back(c);
			stream.read(reinterpret_cast<char*>(&c),sizeof(char));
		}
	}

	if (version != _binaryVersion) throw std::runtime_error("Grid_3d::Table<<(): wrong version!");

	// ***

	float infinityFactor;
	stream.read(reinterpret_cast<char*>(&infinityFactor),sizeof(float));
	if (infinityFactor != _infinityFactor) throw std::runtime_error("Grid_3d::Table::operator<<()");

	stream.read(reinterpret_cast<char*>(&_relaxationFactor),sizeof(double));
	stream.read(reinterpret_cast<char*>(&_phiSwitch),sizeof(int));

	// ***

	{
		int nodeSize;
		stream.read(reinterpret_cast<char*>(&nodeSize),sizeof(int));
		_nodes.allocate(nodeSize);
	}

	for (unsigned int i = 0; i < _nodes.size(); i++)
	{
		Configuration::NodeId id;
		stream.read(reinterpret_cast<char*>(&id),sizeof(Configuration::NodeId));
		_nodes[i].setId(id);

		Configuration::NodeNumber number;
		stream.read(reinterpret_cast<char*>(&number),sizeof(Configuration::NodeNumber));
		_nodes[i].setNumber(number);

		float phiOne;
		stream.read(reinterpret_cast<char*>(&phiOne),sizeof(float));
		_nodes[i].setPhi(0,phiOne);

		float phiTwo;
		stream.read(reinterpret_cast<char*>(&phiTwo),sizeof(float));
		_nodes[i].setPhi(1,phiTwo);

		float charge;
		stream.read(reinterpret_cast<char*>(&charge),sizeof(float));
		_nodes[i].setCharge(charge);

		bool boundary;
		stream.read(reinterpret_cast<char*>(&boundary),sizeof(bool));
		_nodes[i].setBoundary(boundary);

		bool dirichletBoundary;
		stream.read(reinterpret_cast<char*>(&dirichletBoundary),sizeof(bool));
		_nodes[i].setDirichletBoundary(dirichletBoundary);

		bool neumannBoundary;
		stream.read(reinterpret_cast<char*>(&neumannBoundary),sizeof(bool));
		_nodes[i].setNeumannBoundary(neumannBoundary);

		int neighSize;
		stream.read(reinterpret_cast<char*>(&neighSize),sizeof(int));
		_nodes[i].setNumNeighbours(neighSize);

		for (int j = 0; j < _nodes[i].numNeighbours(); j++)
		{
			int neighbour;
			stream.read(reinterpret_cast<char*>(&neighbour),sizeof(int));
			_nodes[i].setNeighbour(j,neighbour);

			float coupling;
			stream.read(reinterpret_cast<char*>(&coupling),sizeof(float));
			_nodes[i].setCoupling(j,coupling);
		}
	}
}

void Grid_3d::Table::operator>>(std::ofstream& stream) const
{
	const char* version = _binaryVersion.c_str();
	stream.write(version,_binaryVersion.size()+1);

	// ***

	const float infinityFactor = _infinityFactor;
	stream.write(reinterpret_cast<const char*>(&infinityFactor),sizeof(float));

	stream.write(reinterpret_cast<const char*>(&_relaxationFactor),sizeof(double));
	stream.write(reinterpret_cast<const char*>(&_phiSwitch),sizeof(int));

	// ***

	const unsigned int nodeSize = _nodes.size();
	stream.write(reinterpret_cast<const char*>(&nodeSize),sizeof(unsigned int));

	for (unsigned int i = 0; i < _nodes.size(); i++)
	{
		const Configuration::NodeId id = _nodes[i].id();
		stream.write(reinterpret_cast<const char*>(&id),sizeof(Configuration::NodeId));

		const Configuration::NodeNumber number = _nodes[i].number();
		stream.write(reinterpret_cast<const char*>(&number),sizeof(Configuration::NodeNumber));

		const float phiOne = _nodes[i].phi(0);
		stream.write(reinterpret_cast<const char*>(&phiOne),sizeof(float));

		const float phiTwo = _nodes[i].phi(1);
		stream.write(reinterpret_cast<const char*>(&phiTwo),sizeof(float));

		const float charge = _nodes[i].charge();
		stream.write(reinterpret_cast<const char*>(&charge),sizeof(float));

		const bool boundary = _nodes[i].boundary();
		stream.write(reinterpret_cast<const char*>(&boundary),sizeof(bool));

		const bool dirichletBoundary = _nodes[i].dirichletBoundary();
		stream.write(reinterpret_cast<const char*>(&dirichletBoundary),sizeof(bool));

		const bool neumannBoundary = _nodes[i].neumannBoundary();
		stream.write(reinterpret_cast<const char*>(&neumannBoundary),sizeof(bool));

		const int numNeighbours = _nodes[i].numNeighbours();
		stream.write(reinterpret_cast<const char*>(&numNeighbours),sizeof(int));

		for (int j = 0; j < _nodes[i].numNeighbours(); j++)
		{
			const int neighbour = _nodes[i].neighbour(j);
			stream.write(reinterpret_cast<const char*>(&neighbour),sizeof(int));

			const float coupling = _nodes[i].coupling(j);
			stream.write(reinterpret_cast<const char*>(&coupling),sizeof(float));
		}
	}
}

void Grid_3d::Table::sync(const Geometry_3d::Table& geomTable, const Configuration::Table& configTable)
{
	for (unsigned int index = 0; index < _nodes.size(); index++)
	{
		//const float chargeValue = computeCharge(index,geomTable,configTable);
		//_nodes[index].setCharge(chargeValue);
	
		const bool boundary = geomTable.node(index).boundary;
		_nodes[index].setBoundary(boundary);
	
		const bool dirichlet = configTable[_nodes[index].id()].dirichletBoundary();
		_nodes[index].setDirichletBoundary(dirichlet);
	
		const bool neumann = configTable[_nodes[index].id()].neumannBoundary();
		_nodes[index].setNeumannBoundary(neumann);
	
		// ***
	
		{
			std::set<int> delaunayConnections;
			geomTable.adjacentNodes(index,&delaunayConnections);
	
			_nodes[index].setNumNeighbours(delaunayConnections.size());
	
			int j = 0;
			std::set<int>::const_iterator k = delaunayConnections.begin();
	
			while (delaunayConnections.end() != k)
			{
				_nodes[index].setNeighbour(j,*k);
	
				const float couplingValue = computeCoupling(index,*k,geomTable,configTable);
				_nodes[index].setCoupling(j,couplingValue);
	
				j++, k++;
			}
		}

		// ***

		const float chargeValue = computeCharge(index,geomTable,configTable);
		_nodes[index].setCharge(chargeValue);
	}
}

void Grid_3d::Table::resync(const int nodeIndex, const Geometry_3d::Table& geomTable, const Configuration::Table& configTable)
{
	std::stack<int> indices;

	indices.push(nodeIndex);
	for (int i = 0; i < _nodes[nodeIndex].numNeighbours(); i++)
		indices.push(_nodes[nodeIndex].neighbour(i));

	// ***

	while(!indices.empty());
	{
		//const float chargeValue = computeCharge(indices.top(),geomTable,configTable);
		//_nodes[indices.top()].setCharge(chargeValue);
	
		const bool boundary = geomTable.node(indices.top()).boundary;
		_nodes[indices.top()].setBoundary(boundary);
	
		const bool dirichlet = configTable[_nodes[indices.top()].id()].dirichletBoundary();
		_nodes[indices.top()].setDirichletBoundary(dirichlet);
	
		const bool neumann = configTable[_nodes[indices.top()].id()].neumannBoundary();
		_nodes[indices.top()].setNeumannBoundary(neumann);
	
		// ***
	
		{
			std::set<int> delaunayConnections;
			geomTable.adjacentNodes(indices.top(),&delaunayConnections);
	
			_nodes[indices.top()].setNumNeighbours(delaunayConnections.size());
	
			int j = 0;
			std::set<int>::const_iterator k = delaunayConnections.begin();
	
			while (delaunayConnections.end() != k)
			{
				_nodes[indices.top()].setNeighbour(j,*k);
	
				const float couplingValue = computeCoupling(indices.top(),*k,geomTable,configTable);
				_nodes[indices.top()].setCoupling(j,couplingValue);
	
				j++, k++;
			}
		}

		// ***

		const float chargeValue = computeCharge(indices.top(),geomTable,configTable);
		_nodes[indices.top()].setCharge(chargeValue);

		indices.pop();
	}
}

void Grid_3d::Table::fastSync(const Geometry_3d::Table& geomTable, const Configuration::Table& configTable)
{
	for (unsigned int nodeIndex = 0; nodeIndex < _nodes.size(); nodeIndex++)
	{
		//const float chargeValue = computeCharge(nodeIndex,geomTable,configTable);
		//_nodes[nodeIndex].setCharge(chargeValue);
	
		const bool boundary = geomTable.node(nodeIndex).boundary;
		_nodes[nodeIndex].setBoundary(boundary);
	
		const bool dirichlet = configTable[_nodes[nodeIndex].id()].dirichletBoundary();
		_nodes[nodeIndex].setDirichletBoundary(dirichlet);
	
		const bool neumann = configTable[_nodes[nodeIndex].id()].neumannBoundary();
		_nodes[nodeIndex].setNeumannBoundary(neumann);

		std::set<int> delaunayConnections;
		geomTable.adjacentNodes(nodeIndex,&delaunayConnections);

		_nodes[nodeIndex].setNumNeighbours(delaunayConnections.size());

		int j = 0;
		std::set<int>::const_iterator k = delaunayConnections.begin();

		while (delaunayConnections.end() != k)
		{
			_nodes[nodeIndex].setNeighbour(j,*k);
			_nodes[nodeIndex].setCoupling(j,std::numeric_limits<float>::quiet_NaN());

			j++, k++;
		}

		const float chargeValue = computeCharge(nodeIndex,geomTable,configTable);
		_nodes[nodeIndex].setCharge(chargeValue);
	}

	// *** initialize coupling factors

	for (int nodeIndex = 0; nodeIndex < (int) _nodes.size(); nodeIndex++)
	{
		for (int i = 0; i < _nodes[nodeIndex].numNeighbours(); i++)
		{
			if (_nodes[nodeIndex].coupling(i) != _nodes[nodeIndex].coupling(i))
			{
				const int neighIndex = _nodes[nodeIndex].neighbour(i);

				int j(0);
				while (_nodes[neighIndex].neighbour(j) != nodeIndex)
					j++;

				if (_nodes[neighIndex].coupling(j) != _nodes[neighIndex].coupling(j))
				{
					const float couplingValue = computeCoupling(nodeIndex,neighIndex,geomTable,configTable);

					_nodes[nodeIndex].setCoupling(i,couplingValue);
					_nodes[neighIndex].setCoupling(j,couplingValue);
				}
				else
					_nodes[nodeIndex].setCoupling(i,_nodes[neighIndex].coupling(j));
			}
		}
	}
}

void Grid_3d::Table::fastResync(const int nodeIndex, const Geometry_3d::Table& geomTable, const Configuration::Table& configTable)
{
	//const float chargeValue = computeCharge(nodeIndex,geomTable,configTable);
	//_nodes[nodeIndex].setCharge(chargeValue);

	const bool boundary = geomTable.node(nodeIndex).boundary;
	_nodes[nodeIndex].setBoundary(boundary);

	const bool dirichlet = configTable[_nodes[nodeIndex].id()].dirichletBoundary();
	_nodes[nodeIndex].setDirichletBoundary(dirichlet);

	const bool neumann = configTable[_nodes[nodeIndex].id()].neumannBoundary();
	_nodes[nodeIndex].setNeumannBoundary(neumann);

	// ***

	for (int i = 0; i < _nodes[nodeIndex].numNeighbours(); i++)
	{
		const int neighIndex = _nodes[nodeIndex].neighbour(i);

		const float couplingValue = computeCoupling(nodeIndex,neighIndex,geomTable,configTable);

		int j(0);
		while (_nodes[neighIndex].neighbour(j) != nodeIndex)
			j++;

		_nodes[nodeIndex].setCoupling(i,couplingValue);
		_nodes[neighIndex].setCoupling(j,couplingValue);
	}

	// ***

        const float chargeValue = computeCharge(nodeIndex,geomTable,configTable);
	_nodes[nodeIndex].setCharge(chargeValue);
}

double Grid_3d::Table::relax(const unsigned int cycleSize)
{
	doRelax(cycleSize);
	
	// *** compute deviation

	double deviation(0.0);

	for (unsigned int i = 0; i < _nodes.size(); i++)
	{
		double localDeviation = _nodes[i].phi(0);
		localDeviation -= _nodes[i].phi(1);
		localDeviation = std::fabs(localDeviation);

		if (localDeviation > deviation) deviation = localDeviation;
	}

	return deviation;
}

signed long Grid_3d::Table::relax(const double threshold, const unsigned int cycleSize, const unsigned int queueSize)
{
	double deviation;
	signed long iteration(0);

	std::deque<double> queue;

	do
	{
		deviation = relax(cycleSize);

		// ***

		queue.push_back(deviation);
		while (queue.size() > queueSize)
			queue.pop_front();

		if (queueSize > 0 && queue.size() == queueSize)
		{
			const double xMean = (queueSize - 1) / 2.0f;
			
			double yMean(0.0f);
			for (std::deque<double>::const_iterator i = queue.begin(); i != queue.end(); i++)
				yMean += *i;
			yMean /= queueSize;

			double nominator(0.0f);

			int x;
			std::deque<double>::const_iterator j;
			for (x = 0, j = queue.begin(); j != queue.end(); x++, j++)
				nominator += (x - xMean) * (*j - yMean);

			double denominator = queueSize;
			denominator *= queueSize;
			denominator -= 1;
			denominator *= queueSize;
			denominator *= cycleSize;
			denominator /= 12.0f;

			const double slope = nominator / denominator;

			if (slope == 0.0f)  throw RelaxationFault("Grid_3d::Table::relax()",iteration,deviation,slope);
		}

		// ***

		iteration += cycleSize;

		#ifdef DEBUG_MODE
		info::begin() << "Grid_3d::Table::adaptiveRelax(): ";
		info::out() << threshold / deviation * 100.0f << "% (deviation = " << deviation << ", ";
		info::out() << " iteration #" << iteration << ")" << std::string(10,' ') << '\r';
		info::out().flush();
		#endif

	}
	while (deviation > threshold);

	return iteration;
}

double Grid_3d::Table::localRelax(const int index, const unsigned int order, const unsigned int cycleSize)
{
	std::set<int> nodes;
	find_localNodes(index,order,&nodes);

	{
		std::set<int>::iterator i = nodes.begin();
		while (nodes.end() != i)
		{
			if (_nodes[*i].dirichletBoundary())
				nodes.erase(i++);
			else
				i++;
		}
	}

	doLocalRelax(cycleSize,nodes);

	// *** compute deviation

	double deviation(0.0);

	for (std::set<int>::const_iterator i = nodes.begin(); i != nodes.end(); i++)
	{
		double localDeviation = _nodes[*i].phi(0);
		localDeviation -= _nodes[*i].phi(1);
		localDeviation = std::fabs(localDeviation);

		if (localDeviation > deviation) deviation = localDeviation;
	}

	return deviation;
}

signed long Grid_3d::Table::localRelax(const int index, const unsigned int order, const double threshold, const unsigned int cycleSize, const unsigned int queueSize)
{
	std::set<int> nodes;
	find_localNodes(index,order,&nodes);

	{
		std::set<int>::iterator i = nodes.begin();
		while (nodes.end() != i)
		{
			if (_nodes[*i].dirichletBoundary())
				nodes.erase(i++);
			else
				i++;
		}
	}

	// ***

	double deviation;
	signed long iteration(0);

	std::deque<double> queue;

	do
	{
		doLocalRelax(cycleSize,nodes);

		// *** compute deviation
	
		deviation = 0.0f;
		for (std::set<int>::const_iterator i = nodes.begin(); i != nodes.end(); i++)
		{
			double localDeviation = _nodes[*i].phi(0);
			localDeviation -= _nodes[*i].phi(1);
			localDeviation = std::fabs(localDeviation);
	
			if (localDeviation > deviation) deviation = localDeviation;
		}

		// ***

		queue.push_back(deviation);
		while (queue.size() > queueSize)
			queue.pop_front();

		if (queueSize > 0 && queue.size() == queueSize)
		{
			const double xMean = (queueSize - 1) / 2.0f;
			
			double yMean(0.0f);
			for (std::deque<double>::const_iterator i = queue.begin(); i != queue.end(); i++)
				yMean += *i;
			yMean /= queueSize;

			double nominator(0.0f);

			int x;
			std::deque<double>::const_iterator j;
			for (x = 0, j = queue.begin(); j != queue.end(); x++, j++)
				nominator += (x - xMean) * (*j - yMean);

			double denominator = queueSize;
			denominator *= queueSize;
			denominator -= 1;
			denominator *= queueSize;
			denominator *= cycleSize;
			denominator /= 12.0f;

			const double slope = nominator / denominator;

			if (slope == 0.0f)  throw RelaxationFault("Grid_3d::Table::relax()",iteration,deviation,slope);
		}

		// ***

		iteration += cycleSize;

		#ifdef DEBUG_MODE
		info::begin() << "Grid_3d::Table::adaptiveRelax(): ";
		info::out() << threshold / deviation * 100.0f << "% (deviation = " << deviation << ", ";
		info::out() << " iteration #" << iteration << ")" << std::string(10,' ') << '\r';
		info::out().flush();
		#endif

	}
	while (deviation > threshold);

	return iteration;
}

void Grid_3d::Table::reset(const Configuration::Table& configTable)
{
	for (unsigned int index = 0; index < _nodes.size(); index++)
		reset(index,configTable);
}

void Grid_3d::Table::reset(const int index, const Configuration::Table& configTable)
{
	const float phi = configTable[_nodes[index].id()].phi();

	for (int i = 0; i < 2; i++)
		_nodes[index].setPhi(i,phi);
}

void Grid_3d::Table::randomReset(const Configuration::Table& configTable)
{
	const std::list<Configuration::NodeId> idList = configTable.ids();

	float minPhi = configTable[idList.front()].phi();
	float maxPhi = configTable[idList.front()].phi();

	for (std::list<Configuration::NodeId>::const_iterator i = idList.begin(); i != idList.end(); i++)
	{
		if (configTable[*i].phi() < minPhi)
			minPhi = configTable[*i].phi();
		else
		{
			if (configTable[*i].phi() > maxPhi)
				maxPhi = configTable[*i].phi();
		}
	}

	// ***

	for (unsigned int index = 0; index < _nodes.size(); index++)
	{
		float phi;

		if (_nodes[index].dirichletBoundary())
			phi = configTable[_nodes[index].id()].phi();
		else
		{
			phi = static_cast<float>(std::rand());
			phi /= RAND_MAX;
			phi *= maxPhi - minPhi;
			phi += minPhi;
		}

		for (int i = 0; i < 2; i++)
			_nodes[index].setPhi(i,phi);
	}
}

void Grid_3d::Table::randomReset(const int index, const Configuration::Table& configTable)
{
	const std::list<Configuration::NodeId> idList = configTable.ids();

	float minPhi = configTable[idList.front()].phi();
	float maxPhi = configTable[idList.front()].phi();

	for (std::list<Configuration::NodeId>::const_iterator i = idList.begin(); i != idList.end(); i++)
	{
		if (configTable[*i].phi() < minPhi)
			minPhi = configTable[*i].phi();
		else
		{
			if (configTable[*i].phi() > maxPhi)
				maxPhi = configTable[*i].phi();
		}
	}

	// ***

	float phi;

	if (_nodes[index].dirichletBoundary())
		phi = configTable[_nodes[index].id()].phi();
	else
	{
		phi = static_cast<float>(std::rand());
		phi /= RAND_MAX;
		phi *= maxPhi - minPhi;
		phi += minPhi;
	}

	for (int i = 0; i < 2; i++)
		_nodes[index].setPhi(i,phi);
}

float Grid_3d::Table::flux(const int index, const Geometry_3d::Table& geomTable) const
{
	float value(0.0f);

	for (int i = 0; i < _nodes[index].numNeighbours(); i++)
	{
		const int neighIndex = _nodes[index].neighbour(i);

		// ***

		float localFlux = potential(neighIndex);
		localFlux -= potential(index);
		localFlux /= -1.0f * geomTable.distance(index,neighIndex);
		localFlux *= geomTable.voronoiArea(index,neighIndex);

		value += localFlux;
	}

	if (value != value) throw std::runtime_error("Grid_3d::Table::flux()");

	return value;
}

MathVector3d<float> Grid_3d::Table::field_o1(const int index, const Geometry_3d::Table& geomTable) const
{
	if (index < 0) throw std::runtime_error("Grid_3d::Table::field_o1()");

	if (_nodes[index].boundary()) return MathVector3d<float>(0.0f);

	std::vector<int> nodeIndizes;
	nodeIndizes.reserve(_nodes[index].numNeighbours());

	for (int i = 0; i < _nodes[index].numNeighbours(); i++)
		if (_nodes[index].coupling(i) != 0.0f) nodeIndizes.push_back( _nodes[index].neighbour(i));
		

	std::vector<Geometry_3d::Point> coords;
	for (std::vector<int>::const_iterator i = nodeIndizes.begin(); i != nodeIndizes.end(); i++)
		coords.push_back(geomTable.nodeCoords(*i) - geomTable.nodeCoords(index));

	// ***

	float a[3*3];

	a[0*3+0] = 0.0;
	for (unsigned int k = 0; k < nodeIndizes.size(); k++)
		a[0*3+0] += coords[k].x() * coords[k].x();

	a[0*3+1] = 0.0;
	for (unsigned int k = 0; k < nodeIndizes.size(); k++)
		a[0*3+1] += coords[k].x() * coords[k].y();

	a[0*3+2] = 0.0;
	for (unsigned int k = 0; k < nodeIndizes.size(); k++)
		a[0*3+2] += coords[k].x() * coords[k].z();

	a[1*3+1] = 0.0;
	for (unsigned int k = 0; k < nodeIndizes.size(); k++)
		a[1*3+1] += coords[k].y() * coords[k].y();

	a[1*3+2] = 0.0;
	for (unsigned int k = 0; k < nodeIndizes.size(); k++)
		a[1*3+2] += coords[k].y() * coords[k].z();

	a[2*3+2] = 0.0;
	for (unsigned int k = 0; k < nodeIndizes.size(); k++)
		a[2*3+2] += coords[k].z() * coords[k].z();

	// ***

	a[1*3+0] = a[0*3+1];
	a[2*3+0] = a[0*3+2];
	a[2*3+1] = a[1*3+2];

	// ***

	float b[3];

	for (int i = 0; i < 3; i++)
	{
		b[i] = 0.0;
		for (unsigned int j = 0; j < nodeIndizes.size(); j++)
			b[i] += (potential(nodeIndizes[j]) - potential(index)) * coords[j][i];
	}

	// ***

	float d;
	int perIndex[3];
		
	if (!lu_decmp<float,3>(a,perIndex,&d)) throw std::runtime_error("Grid_3d::Table::field_o1()");
	lu_solve<float,3>(a,perIndex,b);

	// ***

	MathVector3d<float> value;

	for (int i = 0; i < 3; i++)
	{
		if (b[i] != b[i]) throw std::runtime_error("Grid_3d::Table::field_o1()");

		value[i] = -1.0f * b[i];
	}

	return value;
}

MathVector3d<float> Grid_3d::Table::field_o2(const int index, const Geometry_3d::Table& geomTable, float* charge) const
{
	if (index < 0) throw std::runtime_error("Grid_3d::Table::field_o2()");

	if (_nodes[index].boundary())
	{
		if (charge != 0) *charge = 0.0f;
		return MathVector3d<float>(0.0f);
	}

	std::vector<int> nodeIndizes;
	nodeIndizes.reserve(_nodes[index].numNeighbours());

	for (int i = 0; i < _nodes[index].numNeighbours(); i++)
		nodeIndizes.push_back( _nodes[index].neighbour(i));

	std::vector<Geometry_3d::Point> coords;
	coords.reserve(nodeIndizes.size());

	std::vector<float> potentials;
	potentials.reserve(nodeIndizes.size());
	
	for (std::vector<int>::const_iterator i = nodeIndizes.begin(); i != nodeIndizes.end(); i++)
	{
		coords.push_back(geomTable.nodeCoords(*i) - geomTable.nodeCoords(index));
		potentials.push_back(potential(*i) - potential(index));
	}
	
	// ***

	double a[9*9];
	
	// => neglecting symmetry issues of the solution matrix !!!

	for (int i = 0; i < 9; i++)
	{
		for (int j = 0; j < 9; j++)
			a[i*9+j] = 0.0;
	}
	
	for (std::vector<Geometry_3d::Point>::const_iterator n = coords.begin(); n != coords.end(); n++)
	{
		a[0*9+0] += n->x() * n->x();
		a[0*9+1] += n->x() * n->y();
		a[0*9+2] += n->x() * n->z();
		a[0*9+3] += n->x() * n->x() * n->x() * 0.5;
		a[0*9+4] += n->x() * n->x() * n->y(); 
		a[0*9+5] += n->x() * n->x() * n->z();
		a[0*9+6] += n->x() * n->y() * n->y() * 0.5;
		a[0*9+7] += n->x() * n->y() * n->z();
		a[0*9+8] += n->x() * n->z() * n->z() * 0.5;

		//a[1*9+0] += n->y() * n->x();
		a[1*9+1] += n->y() * n->y();
		a[1*9+2] += n->y() * n->z();
		a[1*9+3] += n->y() * n->x() * n->x() * 0.5;
		a[1*9+4] += n->y() * n->x() * n->y(); 
		a[1*9+5] += n->y() * n->x() * n->z();
		a[1*9+6] += n->y() * n->y() * n->y() * 0.5;
		a[1*9+7] += n->y() * n->y() * n->z();
		a[1*9+8] += n->y() * n->z() * n->z() * 0.5;

		//a[2*9+0] += n->z() * n->x();
		//a[2*9+1] += n->z() * n->y();
		a[2*9+2] += n->z() * n->z();
		a[2*9+3] += n->z() * n->x() * n->x() * 0.5;
		a[2*9+4] += n->z() * n->x() * n->y(); 
		a[2*9+5] += n->z() * n->x() * n->z();
		a[2*9+6] += n->z() * n->y() * n->y() * 0.5;
		a[2*9+7] += n->z() * n->y() * n->z();
		a[2*9+8] += n->z() * n->z() * n->z() * 0.5;

		//a[3*9+0] += n->x() * n->x() * 0.5 * n->x();
		//a[3*9+1] += n->x() * n->x() * 0.5 * n->y();
		//a[3*9+2] += n->x() * n->x() * 0.5 * n->z();
		a[3*9+3] += n->x() * n->x() * 0.5 * n->x() * n->x() * 0.5;
		a[3*9+4] += n->x() * n->x() * 0.5 * n->x() * n->y();
		a[3*9+5] += n->x() * n->x() * 0.5 * n->x() * n->z();
		a[3*9+6] += n->x() * n->x() * 0.5 * n->y() * n->y() * 0.5;
		a[3*9+7] += n->x() * n->x() * 0.5 * n->y() * n->z();
		a[3*9+8] += n->x() * n->x() * 0.5 * n->z() * n->z() * 0.5;

		//a[4*9+0] += n->x() * n->y() * n->x();
		//a[4*9+1] += n->x() * n->y() * n->y();
		//a[4*9+2] += n->x() * n->y() * n->z();
		//a[4*9+3] += n->x() * n->y() * n->x() * n->x() * 0.5;
		a[4*9+4] += n->x() * n->y() * n->x() * n->y();
		a[4*9+5] += n->x() * n->y() * n->x() * n->z();
		a[4*9+6] += n->x() * n->y() * n->y() * n->y() * 0.5;
		a[4*9+7] += n->x() * n->y() * n->y() * n->z();
		a[4*9+8] += n->x() * n->y() * n->z() * n->z() * 0.5;

		//a[5*9+0] += n->x() * n->z() * n->x();
		//a[5*9+1] += n->x() * n->z() * n->y();
		//a[5*9+2] += n->x() * n->z() * n->z();
		//a[5*9+3] += n->x() * n->z() * n->x() * n->x() * 0.5;
		//a[5*9+4] += n->x() * n->z() * n->x() * n->y();
		a[5*9+5] += n->x() * n->z() * n->x() * n->z();
		a[5*9+6] += n->x() * n->z() * n->y() * n->y() * 0.5;
		a[5*9+7] += n->x() * n->z() * n->y() * n->z();
		a[5*9+8] += n->x() * n->z() * n->z() * n->z() * 0.5;

		//a[6*9+0] += n->y() * n->y() * 0.5 * n->x();
		//a[6*9+1] += n->y() * n->y() * 0.5 * n->y();
		//a[6*9+2] += n->y() * n->y() * 0.5 * n->z();
		//a[6*9+3] += n->y() * n->y() * 0.5 * n->x() * n->x() * 0.5;
		//a[6*9+4] += n->y() * n->y() * 0.5 * n->x() * n->y();
		//a[6*9+5] += n->y() * n->y() * 0.5 * n->x() * n->z();
		a[6*9+6] += n->y() * n->y() * 0.5 * n->y() * n->y() * 0.5;
		a[6*9+7] += n->y() * n->y() * 0.5 * n->y() * n->z();
		a[6*9+8] += n->y() * n->y() * 0.5 * n->z() * n->z() * 0.5;

		//a[7*9+0] += n->y() * n->z() * n->x();
		//a[7*9+1] += n->y() * n->z() * n->y();
		//a[7*9+2] += n->y() * n->z() * n->z();
		//a[7*9+3] += n->y() * n->z() * n->x() * n->x() * 0.5;
		//a[7*9+4] += n->y() * n->z() * n->x() * n->y();
		//a[7*9+5] += n->y() * n->z() * n->x() * n->z();
		//a[7*9+6] += n->y() * n->z() * n->y() * n->y() * 0.5;
		a[7*9+7] += n->y() * n->z() * n->y() * n->z();
		a[7*9+8] += n->y() * n->z() * n->z() * n->z() * 0.5;

		//a[8*9+0] += n->z() * n->z() * 0.5 * n->x();
		//a[8*9+1] += n->z() * n->z() * 0.5 * n->y();
		//a[8*9+2] += n->z() * n->z() * 0.5 * n->z();
		//a[8*9+3] += n->z() * n->z() * 0.5 * n->x() * n->x() * 0.5;
		//a[8*9+4] += n->z() * n->z() * 0.5 * n->x() * n->y();
		//a[8*9+5] += n->z() * n->z() * 0.5 * n->x() * n->z();
		//a[8*9+6] += n->z() * n->z() * 0.5 * n->y() * n->y() * 0.5;
		//a[8*9+7] += n->z() * n->z() * 0.5 * n->y() * n->z();
		a[8*9+8] += n->z() * n->z() * 0.5 * n->z() * n->z() * 0.5;
	}

	// set additional values due to symmetry of the matrix

	for (int i = 0; i < 9; i++)
	{
		for (int j = 0; j < i; j++)
			a[i*9+j] = a[j*9+i];
	}
	
	// ***
	
	double b[9];
	
	for (int i = 0; i < 9; i++)
		b[i] = 0.0;
	
	for (unsigned int n  = 0; n < coords.size(); n++)
	{
		double tmp1[3];
		for (int i = 0; i < 3; i++)
		{
			tmp1[i] = potentials[n];
			tmp1[i] *= coords[n][i];
			
			b[i] += tmp1[i];
		}

		int rowIndex(3);
		for (int i = 0; i < 3; i++)
		{
			for (int j = i; j < 3; j++)
			{
				float tmp2 = tmp1[i];
				tmp2 *= coords[n][j];
				if (i == j) tmp2 *= 0.5;
				
				b[rowIndex++] += tmp2;
			}
		}
	}

	// ***

	double d;
	int perIndex[9];

	if (!lu_decmp<double,9>(a,perIndex,&d)) throw std::runtime_error("Grid_3d::Table::field_o2()"); 
	lu_solve<double,9>(a,perIndex,b);

	// ***

	MathVector3d<float> value;

	for (int i = 0; i < 3; i++)
	{
		if (b[i] != b[i]) throw std::runtime_error("Grid_3d::Table::field_o2()");

		value[i] = -1.0f * b[i];
	}

	if (charge != 0)
	{
		*charge = b[3];
		*charge += b[6];
		*charge += b[8];
		*charge *= -1.0f;
		*charge *= _epsilon0;
	}

	#ifdef DEMO_MESH_MODE
	value.y() = 0.0f;
	#endif

	return value;
}

MathVector3d<float> Grid_3d::Table::force(const int index, const Geometry_3d::Table& geomTable, const Configuration::Table& configTable) const
{
	if (index < 0) throw std::runtime_error("Grid_3d::Table::force()");

	const float iPotential = potential(index);
	const float iEpsilon = 1.0f / configTable[_nodes[index].id()].epsilon();
	const MathVector3d<float> iCoords = geomTable.nodeCoords(index);
	
	MathVector3d<float> value(0.0f);
	for (int j = 0; j < _nodes[index].numNeighbours(); j++)
	{
		const int neighIndex = _nodes[index].neighbour(j);

		float magnitude = potential(neighIndex);
		magnitude -= iPotential;
		magnitude *= magnitude;
		magnitude *= _epsilon0;
		magnitude *= _nodes[index].coupling(j);

		const float jEpsilon = 1.0f / configTable[_nodes[neighIndex].id()].epsilon();

		magnitude *= jEpsilon - iEpsilon;

		MathVector3d<float> localValue = geomTable.nodeCoords(neighIndex);
		localValue -= iCoords;
		localValue /= (localValue*localValue);
		localValue *= magnitude;

		value += localValue;
	}

	for (int i = 0; i < 3; i++)
		if (value[i] != value[i]) throw std::runtime_error("Grid_3d::Table::force()");

	return value;
}

void Grid_3d::Table::doRelax(const unsigned int n)
{
	_threadedRelax_params.lock();
	_threadedRelax_params.phiSwitch = _phiSwitch;

	for (unsigned int i = 0; i < n; i++)
	{
		_threadedRelax_params.cycleNum = 1;
		
		_threadedRelax_params.index = 0;
		_threadedRelax_params.delta = _nodes.size() / _threadNum;
	
		pthread_cond_broadcast(&_threadedRelax_params.runCondition);
		
		while (_threadedRelax_params.cycleNum > 0 || _threadedRelax_params.workCnt > 0)
			pthread_cond_wait(&_threadedRelax_params.workCondition,&_threadedRelax_params.mutex);
	}

	_phiSwitch = _threadedRelax_params.phiSwitch;
	_threadedRelax_params.unlock();
}

void Grid_3d::Table::doLocalRelax(const unsigned int n, const std::set<int>& nodes)
{
	_threadedLocalRelax_params.lock();
	_threadedLocalRelax_params.phiSwitch = _phiLocalSwitch;

	_threadedLocalRelax_params.localNodes = &nodes;

	for (unsigned int i = 0; i < n; i++)
	{
		_threadedLocalRelax_params.cycleNum = 1;
	
		_threadedLocalRelax_params.index = nodes.begin();
		_threadedLocalRelax_params.delta = nodes.size() / _threadNum;
	
		pthread_cond_broadcast(&_threadedLocalRelax_params.runCondition);
		
		while (_threadedLocalRelax_params.cycleNum > 0 || _threadedLocalRelax_params.workCnt > 0)
			pthread_cond_wait(&_threadedLocalRelax_params.workCondition,&_threadedLocalRelax_params.mutex);
	}

	_phiLocalSwitch = _threadedLocalRelax_params.phiSwitch;
	_threadedLocalRelax_params.unlock();
}

void Grid_3d::Table::find_localNodes(const int primaryNode, const int limit, std::set<int>* secondaryNodes) const
{
	if (primaryNode < 0) throw std::runtime_error("Grid_3d::Table::find_localNodes()");

	// ***

	std::list<int> untestedNodes;
	untestedNodes.push_back(primaryNode);

	for (int i = 0; i < limit; i++)
	{
		std::list<int> nextCandidates;

		while (!untestedNodes.empty())
		{
			for (int j = 0; j < _nodes[untestedNodes.front()].numNeighbours(); j++)
			{
				const int node = _nodes[untestedNodes.front()].neighbour(j);

				if (secondaryNodes->find(node) == secondaryNodes->end())
				{
					secondaryNodes->insert(node);
					nextCandidates.push_back(node);
				}
			}

			untestedNodes.pop_front();
		}

		untestedNodes.splice(untestedNodes.end(),nextCandidates);
	}
}

double Grid_3d::Table::computePotential(const int index) const
{
	double localValue(0.0);
	double weightSum(0.0);

	const Node& localNode = _nodes[index];

	if (localNode.dirichletBoundary())
		return localNode.phi(_phiSwitch);
	else if (localNode.neumannBoundary())
	{
		for (int i = 0; i < localNode.numNeighbours(); i++)
		{
			const Node& neighNode = _nodes[localNode.neighbour(i)];

			// ***

			const double coupling = localNode.coupling(i);

			weightSum += coupling;

			double adjacentValue = neighNode.phi(_phiSwitch);
			adjacentValue *=  coupling;

			localValue += adjacentValue;
		}
	}
	else
	{
		for (int i = 0; i < localNode.numNeighbours(); i++)
		{
			const Node& neighNode = _nodes[localNode.neighbour(i)];

			if (neighNode.neumannBoundary()) continue;

			// ***

			const double coupling = localNode.coupling(i);

			weightSum += coupling;

			double adjacentValue = neighNode.phi(_phiSwitch);
			adjacentValue *=  coupling;

			localValue += adjacentValue;
		}
	}

	localValue -= localNode.charge();
	localValue /= weightSum;

	if (localValue != localValue) throw std::runtime_error("Grid_3d::Table::computePotential()");

	return localValue;
}

double Grid_3d::Table::computeCharge(const int index, const Geometry_3d::Table& geomTable, const Configuration::Table& configTable) const
{
	float meanDistance(0.0f);
	for (int i = 0; i < _nodes[index].numNeighbours(); i++)
		meanDistance += geomTable.distance(_nodes[index].neighbour(i),index);
	meanDistance /= _nodes[index].numNeighbours();
	meanDistance /= 2.0f;

	// force finite volume at the boundary
	double	value = geomTable.voronoiVolume(index,meanDistance);

	value *= configTable[_nodes[index].id()].chargeDensity();
	value /= _epsilon0;

	if (value != value) throw std::runtime_error("Grid_3d::Table::computeCharge()");

	return value;
}

double Grid_3d::Table::computeCoupling(const int n1, const int n2, const Geometry_3d::Table& geomTable, const Configuration::Table& configTable) const
{
	if (!_nodes[n1].neumannBoundary() && _nodes[n2].neumannBoundary()) return 0.0;

	double value(2.0);
	value *= configTable[_nodes[n1].id()].epsilon() * configTable[_nodes[n2].id()].epsilon();
	value /= configTable[_nodes[n1].id()].epsilon() + configTable[_nodes[n2].id()].epsilon();

	// force finite facet at the boundary
	const float distance = 0.5 * geomTable.distance(n1,n2);
	value *= geomTable.voronoiArea(n1,n2,distance);

	value /= geomTable.distance(n1,n2);

	if (value != value) throw std::runtime_error("Grid_3d::Table::computeCoupling()");

	return value;
}

// ***** ===== GRID_3D::TABLE::THREADED_RELAX_PARAMS ===== *****

Grid_3d::Table::ThreadedRelax_Params::ThreadedRelax_Params(Grid_3d::Table* ptr)
	: obj(ptr),
	  keepRunning(true),
	  cycleNum(0),
	  index(),
	  delta(),
	  workCnt(0)
{
	pthread_mutex_init(&mutex,0);
	pthread_cond_init(&runCondition,0);
	pthread_cond_init(&workCondition,0);
}

Grid_3d::Table::ThreadedRelax_Params::~ThreadedRelax_Params()
{
	pthread_mutex_destroy(&mutex);
	pthread_cond_destroy(&runCondition);
	pthread_cond_destroy(&workCondition);
}

// ***** ===== GRID_3D::TABLE::THREADED_LOCAL_RELAX_PARAMS ===== *****

Grid_3d::Table::ThreadedLocalRelax_Params::ThreadedLocalRelax_Params(Grid_3d::Table* ptr)
	: obj(ptr),
	  keepRunning(true),
	  cycleNum(0),
	  localNodes(0),
	  index(),
	  delta(),
	  workCnt(0)
{
	pthread_mutex_init(&mutex,0);
	pthread_cond_init(&runCondition,0);
	pthread_cond_init(&workCondition,0);
}

Grid_3d::Table::ThreadedLocalRelax_Params::~ThreadedLocalRelax_Params()
{
	pthread_mutex_destroy(&mutex);
	pthread_cond_destroy(&runCondition);
	pthread_cond_destroy(&workCondition);
}

// ***** ===== GRID_3D::CWRAPPER-FUNCTIONS ===== *****

void* Grid_3d::cwrapper_threadedRelax(void* params)
{
	Grid_3d::Table::ThreadedRelax_Params* data = static_cast<Grid_3d::Table::ThreadedRelax_Params*>(params);

	while (true)
	{
		data->lock();

		while (data->cycleNum == 0)
		{
			if (!data->keepRunning)
			{
				data->unlock();
				return 0;
			}

			pthread_cond_wait(&data->runCondition,&data->mutex);
		}

		data->workCnt += 1;

		int min = data->index;
		data->index += data->delta;

		int max = std::min(data->index,data->obj->numNodes());

		int localSwitch = data->phiSwitch;
		localSwitch += 1;
		localSwitch %= 2;

		if (data->index >= data->obj->numNodes())
		{
			data->index = 0;

			data->phiSwitch += 1;
			data->phiSwitch %= 2;

			data->cycleNum -= 1;
		}

		data->unlock();

		for (int i = min; i < max; i++)
		{
			const double phi = data->obj->computePotential(i);
			data->obj->_nodes[i].setPhi(localSwitch,phi);
		}

		data->lock();

		data->workCnt -= 1;

		if (data->cycleNum == 0 && data->workCnt == 0) pthread_cond_signal(&data->workCondition);

		data->unlock();
	}
}

void* Grid_3d::cwrapper_threadedLocalRelax(void* params)
{
	Grid_3d::Table::ThreadedLocalRelax_Params* data = static_cast<Grid_3d::Table::ThreadedLocalRelax_Params*>(params);

	while (true)
	{
		data->lock();

		while (data->cycleNum == 0)
		{
			if (!data->keepRunning)
			{
				data->unlock();
				return 0;
			}

			pthread_cond_wait(&data->runCondition,&data->mutex);
		}

		data->workCnt += 1;

		std::set<int>::const_iterator min = data->index;

		int loopCnt(data->delta);
		while (data->localNodes->end() != data->index && loopCnt-- > 0)
			data->index++;

		std::set<int>::const_iterator max = data->index;

		int localSwitch = data->phiSwitch;
		localSwitch += 1;
		localSwitch %= 2;

		if (data->localNodes->end() == data->index)
		{
			data->index = data->localNodes->begin();

			data->phiSwitch += 1;
			data->phiSwitch %= 2;

			data->cycleNum -= 1;
		}

		data->unlock();

		for (std::set<int>::const_iterator i = min; i != max; i++)
		{
			const double phi = data->obj->computePotential(*i);
			data->obj->_nodes[*i].setPhi(localSwitch,phi);
		}

		data->lock();

		data->workCnt -= 1;

		if (data->cycleNum == 0 && data->workCnt == 0) pthread_cond_signal(&data->workCondition);

		data->unlock();
	}
}

// ***** ===== GRID_3D::POTENTIAL() ===== *****

float Grid_3d::potential(const Geometry_3d::Table& geomTable, const Grid_3d::Table& gridTable, const Geometry_3d::Point& pos)
{
	const int guessIndex = geomTable.findTetrahedron(pos);
	return potential(geomTable,gridTable,pos,guessIndex);
}

float Grid_3d::potential(const Geometry_3d::Table& geomTable, const Grid_3d::Table& gridTable, const Geometry_3d::Point& pos, const int guessIndex)
{
	const int tetIndex = geomTable.findTetrahedron(pos,guessIndex);
	if (tetIndex < 0) return 0.0f;

	// ***

	Quadruple<float> weights;
	geomTable.tetWeights(tetIndex,pos,&weights);

	float value(0.0f);
	for (int index = 0; index < 4; index++)
		value += weights[index] * gridTable.potential(geomTable.tetrahedron(tetIndex).vertices[index]);

	if (value != value) throw std::runtime_error("Grid_3d::potential()");

	return value;
}

// ***** ===== GRID_3D::FIELD() ===== *****

MathVector3d<float> Grid_3d::field(const Geometry_3d::Table& geomTable, const Grid_3d::Table& gridTable, const Geometry_3d::Point& pos)
{
	const int guessIndex = geomTable.findTetrahedron(pos);
	return field(geomTable,gridTable,pos,guessIndex);
}

MathVector3d<float> Grid_3d::field(const Geometry_3d::Table& geomTable, const Grid_3d::Table& gridTable, const Geometry_3d::Point& pos, const int guessIndex)
{
	const int tetIndex = geomTable.findTetrahedron(pos,guessIndex);
	if (tetIndex < 0) return MathVector3d<float>(0.0f);

	// ***

	Quadruple<float> weights;
	geomTable.tetWeights(tetIndex,pos,&weights);

	MathVector3d<float> value(0.0f);
	for (int index = 0; index < 4; index++)
		value += weights[index] * gridTable.field_o1(geomTable.tetrahedron(tetIndex).vertices[index],geomTable);

	for (int index = 0; index < 3; index++)
		if (value[index] != value[index]) throw std::runtime_error("Grid_3d::field()");

	return value;
}

// ***** ===== GRID_3D::FASTFIELD ===== *****

Grid_3d::FastField::FastField(const Table* gridTable, const Geometry_3d::Table* geomTable, const unsigned int maxSize)
	: _gridTable(gridTable),
	  _geomTable(geomTable),
	  _maxSize(maxSize),
	  _buffer()
{}

void Grid_3d::FastField::setGrid(const Table* obj)
{
	if (obj == 0) throw std::runtime_error("Grid_3d::FastField::setGrid()");

	 _gridTable = obj;
	_buffer.clear();
}

void Grid_3d::FastField::setGeometry(const Geometry_3d::Table* obj)
{
	if (obj == 0) throw std::runtime_error("Grid_3d::FastField::setGeometry()");

	_geomTable = obj;
	_buffer.clear();
}

MathVector3d<float> Grid_3d::FastField::compute(const Geometry_3d::Point& pos)
{
	const int guessIndex = _geomTable->findTetrahedron(pos);
	return compute(pos,guessIndex);
}

MathVector3d<float> Grid_3d::FastField::compute(const Geometry_3d::Point& pos, const int guessIndex)
{
	const int tetIndex = _geomTable->findTetrahedron(pos,guessIndex);
	if (tetIndex < 0) return MathVector3d<float>(0.0f);

	// ***

	Quadruple<field_t> nodeFields;

	for (int index = 0; index < 4; index++)
	{
		const int nodeIndex = _geomTable->tetrahedron(tetIndex).vertices[index];

		// ***

		std::map<int,field_t>::const_iterator j = _buffer.find(nodeIndex);

		if (_buffer.end() == j)
		{
			nodeFields[index] = _gridTable->field_o1(nodeIndex,*_geomTable);

			if (_maxSize > 0 && _buffer.size() == _maxSize)
			{
				while (_buffer.size() > _maxSize * 0.75f)
				{
					unsigned int offset = std::rand();
					offset %= _buffer.size() - 1;

					std::map<int,field_t>::iterator k = _buffer.begin();
					advance(k,offset);

					_buffer.erase(k);
				}
			}

			_buffer[nodeIndex] = nodeFields[index];
		}
		else
			nodeFields[index] = j->second;
	}

	// ***

	Quadruple<float> weights;
	_geomTable->tetWeights(tetIndex,pos,&weights);

	MathVector3d<float> value(0.0f);
	for (int index = 0; index < 4; index++)
		value += weights[index] * nodeFields[index];

	for (int index = 0; index < 3; index++)
		if (value[index] != value[index]) throw std::runtime_error("Grid_3d::FastField::compute()");

	#ifdef DEMO_MESH_MODE
		value.y() = 0.0f;
	#endif
		
	return value;
}

#undef DEBUG_MODE
