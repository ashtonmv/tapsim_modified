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

#include "logSurface.h"

#include <iostream>
#include <fstream>
#include <cstdio>
#include <ctime>

#include <pthread.h>
#include <unistd.h>

#include "file_util.h"

#include "trajectory_3d.h"
#include "logTrajectory.h"

const std::string File_Io::LogSurface::_fileVersion("1.0");

File_Io::LogSurface::LogSurface(const char* name, const System_3d* system, const int ioMode, const unsigned int delay, const unsigned int bufferSize, const bool neighbours)
	: _filename(),
	  _system(0),
	  _keepRunning(false),
	  _queue()
{
	pthread_mutex_init(&_queueMutex,0);
	pthread_cond_init(&_queueCondition,0);

	if (name != 0) init(name,system,ioMode,delay,bufferSize,neighbours);
}

File_Io::LogSurface::~LogSurface()
{
	suspend();

	pthread_mutex_destroy(&_queueMutex);
	pthread_cond_destroy(&_queueCondition);

	// ***

	while (!_queue.empty())
	{
		delete _queue.front();
		_queue.pop();
	}
}

void File_Io::LogSurface::init(const char* name, const System_3d* system, const int ioMode, const unsigned int delay, const unsigned int bufferSize, const bool neighbours)
{
	if (_keepRunning) suspend();

	// ***

	_filename = std::string(name);

	if (system == 0)
	{
		_filename.clear();
		throw std::runtime_error("File_Io::LogSurface::init()");
	}
	else
		_system = system;

	_ioMode = ioMode;

	_delay = delay;
	_bufferSize = bufferSize;

	_withNeighbours = neighbours;

	// ***
	
	enable();
}

void File_Io::LogSurface::setDelay(const int delay)
{
	_delay = delay;

	pthread_mutex_lock(&_queueMutex);
	pthread_cond_signal(&_queueCondition);
	pthread_mutex_unlock(&_queueMutex);
}

void File_Io::LogSurface::enable()
{
	if (_keepRunning) return;

	_keepRunning = true;

	pthread_create(&_threadId,0,thread_writeSurface,static_cast<void*>(this));
}

void File_Io::LogSurface::suspend()
{
	if (!_keepRunning) return;

	_keepRunning = false;

	pthread_mutex_lock(&_queueMutex);
	pthread_cond_signal(&_queueCondition);
	pthread_mutex_unlock(&_queueMutex);

	pthread_join(_threadId,0);
}

void File_Io::LogSurface::push(const int id, const Surface_3d::Table& surfaceTable)
{
	if (!_keepRunning)
	{
		if (_filename.empty())
			throw std::runtime_error("File_Io::LogSurface::push(): Bad initialization!");
		else
			if (_bufferSize == 0) throw std::runtime_error("File_Io::LogSurface::push(): Cannot process data!");
	}

	// *** prepare data for caching

	Dataset* obj = new Dataset;

	prepare(id,surfaceTable,*_system,obj,_withNeighbours);

	// *** store data

	pthread_mutex_lock(&_queueMutex);

	if (_bufferSize > 0 && _queue.size() == _bufferSize)
	{
		pthread_mutex_unlock(&_queueMutex);
		throw std::overflow_error("File_Io::LogSurface::push(): Buffer overflow!");
	}

	_queue.push(obj);

	if (_delay == 0) pthread_cond_signal(&_queueCondition);

	pthread_mutex_unlock(&_queueMutex);
}

void File_Io::LogSurface::prepare(const int id, const Surface_3d::Table& surfaceTable, const System_3d& system, Dataset* obj, const bool withNeighbours)
{
	obj->id = id;
	
	for (Surface_3d::Nodeset::const_iterator i = surfaceTable.nodes().begin(); i != surfaceTable.nodes().end(); i++)
	{
		CellData item;
		item.index = i->index();
		item.type = system.gridTable.id(i->index()).toValue();
		item.number = system.gridTable.number(i->index()).toValue();
		item.probability = i->probability();

		obj->cells[item.index] = item;

		// ***
		
		// 1) create entry in obj->nodes (also needed if there wouldn´t be any data)
		// 2) additionally set node data

		obj->nodes[item.index].normal = surfaceTable.normal(i->index(),system);

		// ***

		for (int j = 0; j < system.gridTable.node(item.index).numNeighbours(); j++)
		{
			const int neighIndex = system.gridTable.node(item.index).neighbour(j);

			// ***

			if (withNeighbours) obj->cells[item.index].neighbours.push_back(neighIndex);

			obj->nodes[neighIndex]; // 3) create entry (possibly empty: no normal data)
		}
	}

	for (std::map<int,NodeData>::iterator i = obj->nodes.begin(); i != obj->nodes.end(); i++)
	{
		i->second.phi = system.gridTable.potential(i->first);
		i->second.field = system.gridTable.field_o1(i->first,system.geomTable);
	}	
}

void File_Io::LogSurface::writer(const char* filename, const int ioMode, const Dataset& obj, const bool withNeighbours)
{
	std::ofstream stream(filename,std::ofstream::out|std::ofstream::trunc);
	if (!stream.good()) throw std::runtime_error("File_Io::LogSurface::writer(): Error writing data!");

	// *** header

	stream << "TAPSIM SURFACE DATA" << '\n';
	stream << "VERSION " << _fileVersion << '\n';
	
	if (ioMode == BINARY)
		stream << "BINARY" << '\n';
	else if (ioMode == ASCII)
		stream << "ASCII" <<  '\n';
	else
		throw std::runtime_error("File_Io::LogSurface::writer(): Invalid mode!");
	
	stream << '\n';

	// *** body

	if (ioMode == BINARY)
	{
		stream << "CELLS " << obj.cells.size() << '\n';
		for (std::map<int,CellData>::const_iterator i = obj.cells.begin(); i != obj.cells.end(); i++)
		{
			File_Utils::binWrite<int>(stream,i->second.index);
			File_Utils::binWrite<short>(stream,i->second.type);
			File_Utils::binWrite<unsigned int>(stream,i->second.number);
			File_Utils::binWrite<float>(stream,i->second.probability);

			if (withNeighbours)
			{
				File_Utils::binWrite<unsigned int>(stream,i->second.neighbours.size());
				for (std::list<int>::const_iterator j = i->second.neighbours.begin(); j != i->second.neighbours.end(); j++)
					File_Utils::binWrite<int>(stream,*j);
			}
			else
				File_Utils::binWrite<unsigned int>(stream,0);
		}

		stream << '\n';

		stream << "NODES " << obj.nodes.size() << '\n';
		for (std::map<int,NodeData>::const_iterator i = obj.nodes.begin(); i != obj.nodes.end(); i++)
		{
			File_Utils::binWrite<int>(stream,i->first);
			File_Utils::binWrite<float>(stream,i->second.phi);
			File_Utils::binWrite(stream,i->second.field);
			File_Utils::binWrite(stream,i->second.normal);
		}
	}
	else
	{
		const char spacer = '\t';

		stream.setf(std::ios::scientific);
		
		stream << "CELLS " << obj.cells.size() << '\n';
		for (std::map<int,CellData>::const_iterator i = obj.cells.begin(); i != obj.cells.end(); i++)
		{
			stream << i->second.index << spacer;
			stream << i->second.type << spacer;
			stream << i->second.number << spacer;
			stream << i->second.probability;

			if (withNeighbours)
			{
				stream << spacer << i->second.neighbours.size();
				for (std::list<int>::const_iterator j = i->second.neighbours.begin(); j != i->second.neighbours.end(); j++)
					stream << spacer << *j;
			}
			
			stream << '\n';
		}

		stream << '\n';

		stream << "NODES " << obj.nodes.size() << '\n';
		for (std::map<int,NodeData>::const_iterator i = obj.nodes.begin(); i != obj.nodes.end(); i++)
		{
			stream << i->first << spacer;
			stream << i->second.phi << spacer;
			stream << i->second.field.x() << spacer;
			stream << i->second.field.y() << spacer;
			stream << i->second.field.z() << spacer;
			stream << i->second.normal.x() << spacer;
			stream << i->second.normal.y() << spacer;
			stream << i->second.normal.z() << '\n';
		}
	}

	// ***

	stream.close();
}

void* File_Io::LogSurface::thread_writeSurface(void* params)
{
	LogSurface* root = static_cast<LogSurface*>(params);

	pthread_mutex_lock(&root->_queueMutex);

	while (root->_keepRunning)
	{
		if (root->_queue.empty())
		{

			struct timespec timeout;

			if (root->_delay > 0)
			{
				timeout.tv_sec = std::time(0) + root->_delay;
				timeout.tv_nsec = 0;

				pthread_cond_timedwait(&root->_queueCondition,&root->_queueMutex, &timeout);
			}
			else
				pthread_cond_wait(&root->_queueCondition,&root->_queueMutex);
		}

		while (!root->_queue.empty())
		{
			const Dataset* obj = root->_queue.front();
			root->_queue.pop();

			pthread_mutex_unlock(&root->_queueMutex);

			// ***

			char fileSuffix[25];
			if (std::sprintf(fileSuffix,"%08u",obj->id) < 0)
				throw std::runtime_error("File_Io::LogSurface::thread_writeResults(): sprintf()");

			std::string filename(root->_filename);
			filename += '.';
			filename += fileSuffix;
			
			try
			{
				writer(filename.c_str(),root->_ioMode,*obj,root->_withNeighbours);
			}
			catch   (std::runtime_error&)
			{
				throw std::runtime_error("File_Io::LogSurface::thread_writeResults(): Error writing data!");
			}

			// ***

			delete obj;
			
			pthread_mutex_lock(&root->_queueMutex);
		}
	}

	pthread_mutex_unlock(&root->_queueMutex);

	return 0;
}

// ***

namespace
{
	struct ComputeTrajectory_Data
	{
		unsigned int index;
		unsigned int number;
		Trajectory_3d trajectory;
	};

	struct ComputeTrajectory_Params
	{
		bool terminate;

		float delta;
		float accuracy;

		pthread_cond_t queue_empty;
		pthread_cond_t queue_refilled;

		pthread_mutex_t queue_mutex;
		std::queue<ComputeTrajectory_Data> queue;

		pthread_mutex_t statistics_mutex;
		std::map<int,int> statistics;

		File_Io::LogTrajectory logger;
	};
	
	void* computeTrajectory(void* data)
	{
		ComputeTrajectory_Params* params = static_cast<ComputeTrajectory_Params*>(data);

		pthread_mutex_lock(&params->queue_mutex);

		while (!params->terminate)
		{
			while (!params->queue.empty())
			{
				ComputeTrajectory_Data obj = params->queue.front();
				params->queue.pop();

				pthread_mutex_unlock(&params->queue_mutex);

				obj.trajectory.integrate(params->delta,params->accuracy);

				pthread_mutex_lock(&params->statistics_mutex);
				params->statistics[obj.trajectory.status()]++;
				pthread_mutex_unlock(&params->statistics_mutex);

				try
				{
					params->logger.push(obj.index,obj.number,1.0f,obj.trajectory);
				}
				catch (...)
				{
					std::cout << "push()" << std::endl;
				}

				pthread_mutex_lock(&params->queue_mutex);
			}

			pthread_cond_signal(&params->queue_empty);
			if (!params->terminate) pthread_cond_wait(&params->queue_refilled,&params->queue_mutex);
		}

		pthread_mutex_unlock(&params->queue_mutex);

		return 0;
	}
}

void File_Io::write_surfaceNodes(const char* filename, const Surface_3d::Table& surfaceTable, const System_3d& system, int ioMode)
{
	std::cout << "*** Writing surface nodes to file:" << std::endl;

	std::cout << "\t-> filename: " << filename;

	std::ofstream stream(filename,std::ofstream::out|std::ofstream::trunc);
	if (stream.fail()) throw std::runtime_error("File_Io::write_surfaceNodes(): Error opening file!");

	// *** header information

	if (ioMode == BINARY)
	{
		stream << "# node [1] (int)\n";
		stream << "# number [1] (unsigned int)\n";
		stream << "# id [1] (short)\n";
		stream << "# point [m] (float[3])\n";
		stream << "# field [V/m] (float[3])\n";
		stream << "# normal [1] (float[3])\n";
		stream << "# probability [1] (float)\n";

		stream << "BINARY " << surfaceTable.nodes().size() << '\n';

		std::cout << " (binary mode)" << std::endl;
	}
	else
	{
		stream.setf(std::ios::scientific);

		stream << "# node [1]\n";
		stream << "# number [1]\n";
		stream << "# id [1]\n";
		stream << "# point [m]\n";
		stream << "# field [V/m]\n";
		stream << "# normal [1]\n";
		stream << "# probability [1]\n";

		stream << "ASCII " << surfaceTable.nodes().size() << '\n';

		std::cout << " (ascii mode)" << std::endl;
	}

	// *** data content

	unsigned int progressCnt(0);

	if (ioMode == BINARY)
	{
		for (Surface_3d::Nodeset::const_iterator i = surfaceTable.nodes().begin(); i != surfaceTable.nodes().end(); i++)
		{
			if (++progressCnt % (1+surfaceTable.nodes().size()/100) == 0)
			{
				std::cout << "\t-> writing data: " << progressCnt << " / " << surfaceTable.nodes().size();
				std::cout << std::string(10,' ') << "\r";
				std::cout.flush();
			}

			const int index = i->index();
			stream.write(reinterpret_cast<const char*>(&index),sizeof(int));

			const unsigned int number = system.gridTable.node(i->index()).number().toValue();
			stream.write(reinterpret_cast<const char*>(&number),sizeof(unsigned int));

			const short id = system.gridTable.node(i->index()).id().toValue();
			stream.write(reinterpret_cast<const char*>(&id),sizeof(short));

			const Geometry_3d::Point& point = system.geomTable.nodeCoords(i->index());
			stream.write(reinterpret_cast<const char*>(&point),sizeof(Geometry_3d::Point));

			const MathVector3d<float> field = system.gridTable.field_o1(i->index(),system.geomTable);
			stream.write(reinterpret_cast<const char*>(&field),sizeof(MathVector3d<float>));

			const MathVector3d<float> normal = surfaceTable.normal(i->index(),system);
			stream.write(reinterpret_cast<const char*>(&normal),sizeof(MathVector3d<float>));
			
			const float& probability = i->probability();
			stream.write(reinterpret_cast<const char*>(&probability),sizeof(float));

		}
	}
	else
	{
		const char spacer = '\t';

		for (Surface_3d::Nodeset::const_iterator i = surfaceTable.nodes().begin(); i != surfaceTable.nodes().end(); i++)
		{
			if (++progressCnt % (1+surfaceTable.nodes().size()/100) == 0)
			{
				std::cout << "\t-> writing data: " << progressCnt << " / " << surfaceTable.nodes().size();
				std::cout << std::string(10,' ') << "\r";
				std::cout.flush();
			}

			stream << i->index() << spacer;

			stream << system.gridTable.node(i->index()).number().toValue() << spacer;

			stream << system.gridTable.node(i->index()).id().toValue() << spacer;

			const Geometry_3d::Point& point = system.geomTable.nodeCoords(i->index());
			for (int j = 0; j < 3; j++)
				stream << point[j] << spacer;

			const MathVector3d<float> field = system.gridTable.field_o1(i->index(),system.geomTable);
			for (int j = 0; j < 3; j++)
				stream << field[j] << spacer;

			const MathVector3d<float> normal = surfaceTable.normal(i->index(),system);
			for (int j = 0; j < 3; j++)
				stream << normal[j] << spacer;

			stream << i->probability() << '\n';
		}
	}

	// ***

	stream.close();

	std::cout << "\t-> writing data: finished" << std::string(25,' ') << std::endl;
}

void File_Io::write_surfaceCells(const char* filename, const Surface_3d::Table& surfaceTable, const System_3d& system, const int ioMode, const bool withNeighbours)
{
	std::cout << "*** Writing surface cells to file:" << std::endl;

	// ***
	
	std::cout << "\t-> filename: \"" << filename << "\" ";

	if (ioMode == ASCII)
		std::cout << "(ascii mode)";
	else if (ioMode == BINARY)
		std::cout << "(binary mode)";
	else
		throw std::runtime_error("File_Io::write_surfaceCells(): Invalid mode!");

	std::cout << std::endl;
	
	// ***
	
	std::cout << "\t-> preparing data: ";

	LogSurface::Dataset obj;
	LogSurface::prepare(0,surfaceTable,system,&obj,withNeighbours);

	std::cout << obj.cells.size() << " cells, " << obj.nodes.size() << " nodes" << std::endl;

	// ***
	
	std::cout << "\t-> writing data: ";

	LogSurface::writer(filename,ioMode,obj,withNeighbours);

	std::cout << "finished" << std::endl;
}

void File_Io::write_surfaceTrajectories(const char* filename, const Surface_3d::Table& surfaceTable, const System_3d& system, const float accuracy, int ioMode)
{
	std::cout << "*** Writing surface trajectories to file:" << std::endl;

	std::cout << "\t-> filename: \"" << filename << "\" ";

	if (ioMode == BINARY)
		std::cout << "(binary mode)" << std::endl;
	else
		std::cout << "(ascii mode)" << std::endl;

	// *** initialize threading

	const int threadNum = sysconf(_SC_NPROCESSORS_ONLN);
	std::cout << "\t-> number of computing threads: " << threadNum << std::endl;

	std::cout << "\t-> accuracy: " << accuracy << std::endl;

	ComputeTrajectory_Params params;

	{
		params.terminate = false;

		params.delta = 1e-14f;
		params.accuracy = accuracy;

		pthread_cond_init(&params.queue_empty,0);
		pthread_cond_init(&params.queue_refilled,0);

		pthread_mutex_init(&params.queue_mutex,0);

		pthread_mutex_init(&params.statistics_mutex,0);

		params.logger.init(filename,ioMode,60);
	}

	pthread_t* threadIds = new pthread_t[threadNum];

	for (int i = 0; i < threadNum; i++)
		pthread_create(&threadIds[i],0,computeTrajectory,reinterpret_cast<void*>(&params));

	// *** compute trajectories, refilling the working queue continuously

	const unsigned int queueSize(256);

	int progressCnt(0);
	Surface_3d::Nodeset::const_iterator current = surfaceTable.nodes().begin();

	pthread_mutex_lock(&params.queue_mutex);

	while (surfaceTable.nodes().end() != current)
	{
		while (surfaceTable.nodes().end() != current && params.queue.size() < queueSize)
		{
			if (progressCnt % (1+surfaceTable.nodes().size()/100) == 0)
			{
				std::cout << "\t-> processing data: " << progressCnt << " / " << surfaceTable.nodes().size();
				std::cout << std::string(20,' ') << '\r';
				std::cout.flush();
			}
			
			const Trajectory_3d::Vector3d velocity = Trajectory_3d::Vector3d(0.0f);
	
			const float charge = 3.0f * Trajectory_3d::eCharge;
			const float mass = 63.4f * Trajectory_3d::amu2kg;
	
			ComputeTrajectory_Data obj;
			obj.trajectory.setGeometry(&system.geomTable);
			obj.trajectory.setGrid(&system.gridTable);
			obj.trajectory.setIntegratorType(Trajectory_3d::O5_INTEGRATOR);
			obj.trajectory.setStepperType(Trajectory_3d::ERROR_RESTRICTED);
			obj.trajectory.setTimeStepLimit(1e-16);

			obj.index = ++progressCnt;
			obj.number = system.gridTable.node(current->index()).number().toValue();
			obj.trajectory.init(system.geomTable.nodeCoords(current->index()),velocity,charge,mass);

			params.queue.push(obj);
			
			current++;
		}

		pthread_cond_broadcast(&params.queue_refilled);
		pthread_cond_wait(&params.queue_empty,&params.queue_mutex);
	}
	
	std::cout << "\t-> processing data: " << surfaceTable.nodes().size() << " / " << surfaceTable.nodes().size();
	std::cout << std::string(20,' ') << '\r';
	std::cout.flush();
	
	// ***
	
	params.terminate = true;
	pthread_cond_broadcast(&params.queue_refilled);

	pthread_mutex_unlock(&params.queue_mutex);

	// *** clean up

	for (int i = 0; i < threadNum; i++)
		pthread_join(threadIds[i],0);
	
	delete[] threadIds;

	pthread_cond_destroy(&params.queue_empty);
	pthread_cond_destroy(&params.queue_refilled);
	pthread_mutex_destroy(&params.queue_mutex);
	pthread_mutex_destroy(&params.statistics_mutex);

	std::cout << "\t-> processing data: finished" << std::string(25,' ') << std::endl;

	// ***

	std::cout << "\t-> status statistics: " << std::endl;
	for (std::map<int,int>::const_iterator i = params.statistics.begin(); i != params.statistics.end(); i++)
	{
		std::cout << "\t\t=> trajectories with status code #(" << i->first << "): " << i->second;
		std::cout << " (\"" << Trajectory_3d::status_str(i->first) << "\")" << std::endl;
	}
}
