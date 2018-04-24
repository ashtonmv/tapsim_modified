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

#include "logGrid.h"

#include <fstream>
#include <stdexcept>
#include <sstream>
#include <ctime>
#include <cstring>

#include <iostream>

#include "file_util.h"

const std::string File_Io::LogGrid::_fileVersion("1.0");

File_Io::LogGrid::LogGrid(const char* name, const Geometry_3d::Table* geomTable, const Configuration::Table* configTable, int ioMode, const unsigned int delay, const unsigned int bufferSize)
	: _filename(),
	  _keepRunning(false),
	  _queue()
{
	pthread_mutex_init(&_queueMutex,0);
	pthread_cond_init(&_queueCondition,0);

	if (name != 0) init(name,geomTable,configTable,ioMode,delay,bufferSize);
}

File_Io::LogGrid::~LogGrid()
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

void File_Io::LogGrid::init(const char* name, const Geometry_3d::Table* geomTable, const Configuration::Table* configTable, const int ioMode, const unsigned int delay, const unsigned int bufferSize)
{
	if (_keepRunning) suspend();

	_filename = std::string(name);

	if (geomTable == 0 || configTable == 0)
	{
		_filename.clear();
		throw std::runtime_error("File_Io::LogGrid::init()");
	}
	else
	{
		_geomTable = geomTable;
		_configTable = configTable;
	}

	_ioMode = ioMode;

	_delay = delay;
	_bufferSize = bufferSize;

	enable();
}

void File_Io::LogGrid::setDelay(const int delay)
{
	_delay = delay;

	pthread_mutex_lock(&_queueMutex);
	pthread_cond_signal(&_queueCondition);
	pthread_mutex_unlock(&_queueMutex);
}

void File_Io::LogGrid::enable()
{
	if (_keepRunning) return;

	_keepRunning = true;

	pthread_create(&_threadId,0,thread_writeGrid,static_cast<void*>(this));
}

void File_Io::LogGrid::suspend()
{
	if (!_keepRunning) return;

	_keepRunning = false;

	pthread_mutex_lock(&_queueMutex);
	pthread_cond_signal(&_queueCondition);
	pthread_mutex_unlock(&_queueMutex);

	pthread_join(_threadId,0);
}

void File_Io::LogGrid::push(const int index, const Grid_3d::Table& gridTable)
{
	if (!_keepRunning)
	{
		if (_filename.empty())
			throw std::runtime_error("File_Io::LogGrid::push(): Bad initialization!");
		else
			if (_bufferSize == 0) throw std::runtime_error("File_Io::LogGrid::push(): Cannot process data!");
	}

	Dataset* obj = new Dataset;

	obj->index = index;
	for (int i = 0; i < gridTable.numNodes(); i++)
	{
		const Grid_3d::Node& node = gridTable.node(i);

		NodeData item;

		item.id = node.id().toValue();
		item.number = node.number().toValue();

		item.phi = gridTable.potential(i);
		item.field = gridTable.field_o1(i,*_geomTable);
		item.force = gridTable.force(i,*_geomTable,*_configTable);
		
		obj->data.push_back(item);
	}

	pthread_mutex_lock(&_queueMutex);

	if (_bufferSize > 0 && _queue.size() == _bufferSize)
	{
		pthread_mutex_unlock(&_queueMutex);
		throw std::overflow_error("File_Io::LogGrid::push(): Buffer overflow!");
	}

	_queue.push(obj);

	if (_delay == 0) pthread_cond_signal(&_queueCondition);

	pthread_mutex_unlock(&_queueMutex);
}

void* File_Io::LogGrid::thread_writeGrid(void* params)
{
	LogGrid* root = static_cast<LogGrid*>(params);

	std::ofstream stream;
	stream.setf(std::ios::scientific);

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
			if (std::sprintf(fileSuffix,"%08u",obj->index) < 0)
				throw std::runtime_error("File_Io::thread_writeGrid(): sprintf()");

			std::string filename(root->_filename);
			filename += '.';
			filename += fileSuffix;

			std::ofstream file(filename.c_str(),std::ofstream::out);
			if (!file.good()) throw std::runtime_error("File_Io::thread_writeGrid()");

			
			if (root->_ioMode == File_Io::ASCII)
			{
				file.setf(std::ios_base::scientific|std::ios_base::showpoint);
				
				const char spacer = '\t';
				const char endline = '\n';

				file << "TAPSIM GRID DATA" << endline;
				file << "VERSION " << _fileVersion << endline;
				file << "ASCII " << obj->data.size() << endline;
				for (std::vector<NodeData>::const_iterator i = obj->data.begin(); i != obj->data.end(); i++)
				{
					file << i->id;
					file << spacer << i->number;
					file << spacer << i->phi;
					for (int j = 0; j < 3; j++)
						file << spacer << i->field[j];
					for (int j = 0; j< 3; j++)
						file << spacer << i->force[j];
					file << endline;
				}
			}
			else
			{
				file << "TAPSIM GRID DATA" << '\n';
				file << "VERSION " << _fileVersion << '\n';
				file << "BINARY " << obj->data.size() << '\n';
				for (std::vector<NodeData>::const_iterator i = obj->data.begin(); i != obj->data.end(); i++)
				{
					File_Utils::binWrite<short>(file,i->id);
					File_Utils::binWrite<unsigned int>(file,i->number);
					File_Utils::binWrite<float>(file,i->phi);
					for (int j = 0; j < 3; j++)
						File_Utils::binWrite<float>(file,i->field[j]);
					for (int j = 0; j < 3; j++)
						File_Utils::binWrite<float>(file,i->force[j]);
				}
			}

			// ***

			delete obj;

			pthread_mutex_lock(&root->_queueMutex);
		}
	}

	pthread_mutex_unlock(&root->_queueMutex);

	if (stream.is_open()) stream.close();

	return 0;
}

void File_Io::writeGrid(const char* filename, const Grid_3d::Table& gridTable, const int ioMode)
{
	if (ioMode == ASCII)
	{
		std::ofstream file(filename,std::ofstream::out|std::ofstream::trunc);
		if (!file.good()) throw std::runtime_error("File_Io::writeGrid()");
	
		file.setf(std::ios_base::scientific|std::ios_base::showpoint);
		file.precision(5);
	
		// ***

		const char spacer = '\t';
		const char endline = '\n';

		std::multimap<int,int> connections;
		std::multimap<int,int>::iterator hint = connections.end();

		file << "ASCII" << endline;
		file << endline;
		
		// ***

		file << "NODE_DATA " << gridTable.numNodes() << endline;
		for (int i = 0; i < gridTable.numNodes(); i++)
		{
			const Grid_3d::Node& node = gridTable.node(i);
			
			file << node.id().toValue() << spacer;
			file << node.number().toValue() << spacer;

			file << gridTable.potential(i) << spacer;

			file << node.charge() << spacer;

			file << (node.dirichletBoundary() ? 1 : 0) << spacer;

			file << (node.neumannBoundary() ? 1 : 0) << endline;
			
			// ***
			
			for (int j = 0; j < node.numNeighbours(); j++)
			{
				if (i < node.neighbour(j))
				{
					const std::pair<int,int> obj(i,node.neighbour(j));
					hint = connections.insert(hint,obj);
				}
			}
		}
		file << endline;

		// ***

		file << "COUPLING_DATA " << connections.size() << endline;
		for (std::multimap<int,int>::const_iterator i = connections.begin(); i != connections.end(); i++)
		{
			int index(0);
			while (gridTable.node(i->first).neighbour(index) != i->second)
				index++;
			
			file << i->first << spacer;
			file << i->second << spacer;
			file << gridTable.node(i->first).coupling(index) << endline;
		}
		
		// ***

		file.close();
	}
	else if (ioMode == BINARY)
	{
		std::ofstream file(filename,std::ofstream::out|std::ofstream::binary|std::ofstream::trunc);
		if (!file.good()) throw std::runtime_error("File_Io::writeGrid()");

		const char endline = '\n';

		std::multimap<int,int> connections;
		std::multimap<int,int>::iterator hint = connections.end();

		file << "BINARY" << endline;
		
		// ***

		file << "NODE_DATA " << gridTable.numNodes() << endline;
		for (int i = 0; i < gridTable.numNodes(); i++)
		{
			const Grid_3d::Node& node = gridTable.node(i);

			File_Utils::binWrite<short>(file,node.id().toValue());
			File_Utils::binWrite<unsigned int>(file,node.number().toValue());

			File_Utils::binWrite<float>(file,gridTable.potential(i));

			File_Utils::binWrite<float>(file,node.charge());
			
			const char dirichletBoundary = (node.dirichletBoundary() ? 1 : 0);
			File_Utils::binWrite<char>(file,dirichletBoundary);
			
			char neumannBoundary = (node.neumannBoundary() ? 1 : 0);
			File_Utils::binWrite<char>(file,neumannBoundary);
			
			// ***
			
			for (int j = 0; j < node.numNeighbours(); j++)
			{
				if (i < node.neighbour(j))
				{
					const std::pair<int,int> obj(i,node.neighbour(j));
					hint = connections.insert(hint,obj);
				}
			}
		}
		file << endline;

		// ***

		file << "COUPLING_DATA " << connections.size() << endline;
		for (std::multimap<int,int>::const_iterator i = connections.begin(); i != connections.end(); i++)
		{
			int index(0);
			while (gridTable.node(i->first).neighbour(index) != i->second)
				index++;

			File_Utils::binWrite<int>(file,i->first);
			File_Utils::binWrite<int>(file,i->second);
			File_Utils::binWrite<float>(file,gridTable.node(i->first).coupling(index));
		}
		file << endline;
		
		// ***

		file.close();
	}
	else
		throw std::runtime_error("File_Io::writeGrid()");
}
