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

#include "logTrajectory.h"

#include <fstream>
#include <stdexcept>
#include <ctime>

#include "file_util.h"

const std::string File_Io::LogTrajectory::_fileVersion("1.0");

File_Io::LogTrajectory::LogTrajectory(const char* name, const int ioMode, const unsigned int delay, const unsigned int chunkSize, const unsigned int bufferSize)
	: _filename(),
	  _keepRunning(false),
	  _queue()
{
	pthread_mutex_init(&_queueMutex,0);
	pthread_cond_init(&_queueCondition,0);

	if (name != 0) init(name,ioMode,delay,chunkSize,bufferSize);
}

File_Io::LogTrajectory::~LogTrajectory()
{
	suspend();

	pthread_mutex_destroy(&_queueMutex);
	pthread_cond_destroy(&_queueCondition);

	// ***

	while(!_queue.empty())
	{
		delete _queue.front();
		_queue.pop();
	}
}

void File_Io::LogTrajectory::init(const char* name, const int ioMode, const unsigned int delay, const unsigned int chunkSize, const unsigned int bufferSize)
{
	if (_keepRunning) suspend();

	_filename = std::string(name);
	_ioMode = ioMode;

	_delay = delay;
	_chunkSize = chunkSize;
	_bufferSize = bufferSize;

	enable();
}

void File_Io::LogTrajectory::setDelay(const int delay)
{
	_delay = delay;

	pthread_mutex_lock(&_queueMutex);
	pthread_cond_signal(&_queueCondition);
	pthread_mutex_unlock(&_queueMutex);
}

void File_Io::LogTrajectory::enable()
{
	if (_keepRunning) return;

	_keepRunning = true;

	pthread_create(&_threadId,0,thread_writeTrajectory,static_cast<void*>(this));
}

void File_Io::LogTrajectory::suspend()
{
	if (!_keepRunning) return;

	_keepRunning = false;

	pthread_mutex_lock(&_queueMutex);
	pthread_cond_signal(&_queueCondition);
	pthread_mutex_unlock(&_queueMutex);

	pthread_join(_threadId,0);
}

void File_Io::LogTrajectory::push(const unsigned int index, const unsigned int number, const float scale, const Trajectory_3d& data)
{
	if (!_keepRunning)
	{
		if (_filename.empty() == 0)
			throw std::runtime_error("File_Io::LogTrajectory::push(): Bad initialization!");
		else
			if (_bufferSize == 0) throw std::runtime_error("File_Io::LogTrajectory::push(): Cannot process data!");
	}

	pthread_mutex_lock(&_queueMutex);

	if (_bufferSize > 0 && _queue.size() == _bufferSize)
	{
		pthread_mutex_unlock(&_queueMutex);
		throw std::overflow_error("File_Io::LogTrajectory::push(): Buffer overflow!");
	}

	Dataset* item = new Dataset;
	{
		item->eventIndex = index;
		item->nodeNumber = number;
		item->timeScale = scale;
		item->trajectory = data;
	}

	_queue.push(item);

	if (_delay == 0) pthread_cond_signal(&_queueCondition);
	
	pthread_mutex_unlock(&_queueMutex);
}

void File_Io::LogTrajectory::writeHeader(std::ofstream& stream, const int ioMode)
{
	stream << "TAPSIM TRAJECTORY DATA" << '\n';
	stream << "VERSION " << _fileVersion << '\n';

	if (ioMode == BINARY)
	{
		stream << "# time [s] (float)\n";
		stream << "# positionX [m] (float)\n";
		stream << "# positionY [m] (float)\n";
		stream << "# positionZ [m] (float)\n";
		stream << "# velocityX [m/s] (float)\n";
		stream << "# velocityY [m/s] (float)\n";
		stream << "# velocityZ [m/s] (float)\n";
		stream << "# tetrahedron [1] (int)\n";
	}
	else if (ioMode == ASCII)
	{
		stream << "# time [s]\n";
		stream << "# positionX [m]\n";
		stream << "# positionY [m]\n";
		stream << "# positionZ [m]\n";
		stream << "# velocityX [m/s]\n";
		stream << "# velocityY [m/s]\n";
		stream << "# velocityZ [m/s]\n";
		stream << "# tetrahedron [1]\n";
	}
	else
		throw std::runtime_error("File_Io::LogTrajectory::header(): Invalid output mode!");

	stream << '\n';
}

void* File_Io::LogTrajectory::thread_writeTrajectory(void* params)
{
	LogTrajectory* root = static_cast<LogTrajectory*>(params);

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
		
		if (!stream.good()) throw std::runtime_error("File_Io::LogTrajectory::thread_writeTrajectory(): Error writing data!");

		while (!root->_queue.empty())
		{
			const Dataset* obj = root->_queue.front();
			root->_queue.pop();

			pthread_mutex_unlock(&root->_queueMutex);

			// ***
			
			if (stream.is_open())
			{
				if (root->_chunkSize > 0 && obj->eventIndex % root->_chunkSize == 0)
				{
					stream.close();

					char fileSuffix[25];
					if (std::sprintf(fileSuffix,"%08u",obj->eventIndex) < 0)
						throw std::runtime_error("File_Io::LogTrajectory::thread_writeTrajectory(): sprintf()");

					std::string filename(root->_filename);
					filename += '.';
					filename += fileSuffix;
	
					stream.open(filename.c_str(),std::ofstream::out|std::ofstream::trunc);
	
					writeHeader(stream,root->_ioMode);
				}
			}
			else
			{
				if (root->_chunkSize == 0)
					stream.open(root->_filename.c_str(),std::ofstream::out|std::ofstream::trunc);
				else
				{
					char fileSuffix[25];
					if (std::sprintf(fileSuffix,"%08u",obj->eventIndex) < 0)
						throw std::runtime_error("File_Io::LogTrajectory::thread_writeTrajectory(): sprintf()");

					std::string filename(root->_filename);
					filename += '.';
					filename += fileSuffix;
					
					stream.open(filename.c_str(),std::ofstream::out|std::ofstream::trunc);
				}
				
				writeHeader(stream,root->_ioMode);
			}
			
			// ***

			stream << "# index = " << obj->eventIndex << '\n';
			stream << "# unique number = " << obj->nodeNumber << '\n';
			stream << "# time scale = " << obj->timeScale << '\n';
			stream << "# integrator status = " << obj->trajectory.status() << " => " << Trajectory_3d::status_str(obj->trajectory.status()) << '\n';
			stream << "# charge [C] / mass [kg] = " << obj->trajectory.charge() << " / " << obj->trajectory.mass() << '\n';

			if (root->_ioMode == BINARY)
			{
				stream << "BINARY " << obj->trajectory.data().size() << '\n';
				
				std::vector<Trajectory_3d::phaseVector>::const_iterator i = obj->trajectory.data().begin();
				while (obj->trajectory.data().end() != i)
				{
					File_Utils::binWrite<float>(stream,i->time()*obj->timeScale);
					File_Utils::binWrite(stream,i->position());
					File_Utils::binWrite(stream,i->velocity());
					File_Utils::binWrite<int>(stream,i->tetIndex());

					i++;
				}

				stream << '\n';
			}
			else
			{
				stream << "ASCII " << obj->trajectory.data().size() << '\n';
				
				std::vector<Trajectory_3d::phaseVector>::const_iterator i = obj->trajectory.data().begin();
				while (obj->trajectory.data().end() != i)
				{
					stream << i->time() * obj->timeScale << '\t';
					for (int j = 0; j < 3; j++)
						stream << i->position(j) << '\t';
					for (int j = 0; j < 3; j++)
						stream << i->velocity(j) << '\t';
					stream << i->tetIndex() << '\n';
					
					i++;
				}
			} 
		
			stream << "# final error estimate (position) = (";
			stream << obj->trajectory.error_estimate().position().x() << "/";
			stream << obj->trajectory.error_estimate().position().y() << "/";
			stream << obj->trajectory.error_estimate().position().z() << ")";
			stream << '\n';
		
			stream << "# final error estimate (velocity) = (";
			stream << obj->trajectory.error_estimate().velocity().x() << "/";
			stream << obj->trajectory.error_estimate().velocity().y() << "/";
			stream << obj->trajectory.error_estimate().velocity().z() << ")";
			stream << '\n';
		
			stream << '\n';

			// ***

			delete obj;

			pthread_mutex_lock(&root->_queueMutex);
		}
	}

	pthread_mutex_unlock(&root->_queueMutex);

	if (stream.is_open()) stream.close();

	return 0;
}
