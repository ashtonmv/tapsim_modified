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

#include "logResults.h"

#include <fstream>
#include <stdexcept>
#include <sstream>
#include <ctime>
#include <cstring>

#include <iostream>

#include "file_util.h"

const std::string File_Io::LogResults::_fileVersion("1.0");

File_Io::LogResults::LogResults(const char* name,const int ioMode, const unsigned int delay, const unsigned int chunkSize, const unsigned int bufferSize)
	: _filename(),
	  _keepRunning(false),
	  _queue()
{
	pthread_mutex_init(&_queueMutex,0);
	pthread_cond_init(&_queueCondition,0);

	if (name != 0) init(name,ioMode,delay,chunkSize,bufferSize);
}

File_Io::LogResults::~LogResults()
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

void File_Io::LogResults::init(const char* name, const int ioMode, const unsigned int delay, const unsigned int chunkSize, const unsigned int bufferSize)
{
	if (_keepRunning) suspend();

	_filename = std::string(name);
	_ioMode = ioMode;

	_delay = delay;
	_chunkSize = chunkSize;
	_bufferSize = bufferSize;

	enable();
}

void File_Io::LogResults::setDelay(const int delay)
{
	_delay = delay;

	pthread_mutex_lock(&_queueMutex);
	pthread_cond_signal(&_queueCondition);
	pthread_mutex_unlock(&_queueMutex);
}

void File_Io::LogResults::enable()
{
	if (_keepRunning) return;

	_keepRunning = true;

	pthread_create(&_threadId,0,thread_writeResults,static_cast<void*>(this));
}

void File_Io::LogResults::suspend()
{
	if (!_keepRunning) return;

	_keepRunning = false;

	pthread_mutex_lock(&_queueMutex);
	pthread_cond_signal(&_queueCondition);
	pthread_mutex_unlock(&_queueMutex);

	pthread_join(_threadId,0);
}

void File_Io::LogResults::push(const Dataset& data)
{
	if (!_keepRunning)
	{
		if (_filename.empty())
			throw std::runtime_error("File_Io::LogResults::push(): Bad initialization!");
		else
			if (_bufferSize == 0) throw std::runtime_error("File_Io::LogResults::push(): Cannot process data!");
	}

	Dataset* obj = new Dataset(data);

	pthread_mutex_lock(&_queueMutex);

	if (_bufferSize > 0 && _queue.size() == _bufferSize)
	{
		pthread_mutex_unlock(&_queueMutex);
		throw std::overflow_error("File_Io::LogResults::push(): Buffer overflow!");
	}

	_queue.push(obj);

	if (_delay == 0) pthread_cond_signal(&_queueCondition);

	pthread_mutex_unlock(&_queueMutex);
}

void File_Io::LogResults::writeHeader(std::ofstream& stream,const int ioMode, const std::string& header)
{
	stream << "TAPSIM RESULTS DATA" << '\n';
	stream << "VERSION " << _fileVersion << '\n';

	if (!header.empty()) stream << header;
	
	if (ioMode == BINARY)
	{
		stream << "# index [1] (int)\n";

		stream << "# id [1] (int)\n";
		stream << "# number [1] (unsigned int)\n";
		
		stream << "# voltage [V] (float)\n";
		
		stream << "# startX [m] (float)\n";
		stream << "# startY [m] (float)\n";
		stream << "# startZ [m] (float)\n";
		stream << "# stopX [m] (float)\n";
		stream << "# stopY [m] (float)\n";
		stream << "# stopZ [m] (float)\n";

		stream << "# tof [s] (float)\n";

		stream << "# probability [1] (float)\n";
		
		stream << "# potentialBefore [1] (float)\n";
		stream << "# fieldBeforeX [V/m] (float)\n";
		stream << "# fieldBeforeY [V/m] (float)\n";
		stream << "# fieldBeforeZ [V/m] (float)\n";

		stream << "# potentialAfter [1] (float)\n";
		stream << "# fieldAfterX [V/m] (float)\n";
		stream << "# fieldAfterY [V/m] (float)\n";
		stream << "# fieldAfterZ [V/m] (float)\n";

		stream << "# normalX [m] (float)\n";
		stream << "# normalY [m] (float)\n";
		stream << "# normalZ [m] (float)\n";

		stream << "# apexX [m] (float)\n";
		stream << "# apexY [m] (float)\n";
		stream << "# apexZ [m] (float)\n";
		
		stream << "BINARY\n";
	}
	else if (ioMode == ASCII)
	{
		stream << "# index [1]\n";

		stream << "# id [1]\n";
		stream << "# number [1]\n";

		stream << "# voltage [V]\n";

		stream << "# startX [m]\n";
		stream << "# startY [m]\n";
		stream << "# startZ [m]\n";
		stream << "# stopX [m]\n";
		stream << "# stopY [m]\n";
		stream << "# stopZ [m]\n";

		stream << "# tof [s]\n";

		stream << "# probability [1]\n";

		stream << "# potentialBefore [V]\n";
		stream << "# fieldBeforeX [V/m]\n";
		stream << "# fieldBeforeY [V/m]\n";
		stream << "# fieldBeforeZ [V/m]\n";

		stream << "# potentialAfter [V]\n";
		stream << "# fieldAfterX [V/m]\n";
		stream << "# fieldAfterY [V/m]\n";
		stream << "# fieldAfterZ [V/m]\n";

		stream << "# normalX [m]\n";
		stream << "# normalY [m]\n";
		stream << "# normalZ [m]\n";

		stream << "# apexX [m]\n";
		stream << "# apexY [m]\n";
		stream << "# apexZ [m]\n";

		stream << "ASCII\n";
	}
	else
		throw std::runtime_error("File_Io::LogResults::header(): Invalid output mode!");
}


void* File_Io::LogResults::thread_writeResults(void* params)
{
	LogResults* root = static_cast<LogResults*>(params);

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

			if (stream.is_open())
			{
				if (root->_chunkSize > 0 && obj->eventIndex % root->_chunkSize == 0)
				{
					stream.close();

					char fileSuffix[25];
					if (std::sprintf(fileSuffix,"%08u",obj->eventIndex) < 0)
						throw std::runtime_error("File_Io::LogResults::thread_writeResults");

					std::string filename(root->_filename);
					filename += '.';
					filename += fileSuffix;
					
					stream.open(filename.c_str(),std::ofstream::out|std::ofstream::trunc);
					
					writeHeader(stream,root->_ioMode,std::string());
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
						throw std::runtime_error("File_Io::LogResults::thread_writeResults(): sprintf()");

					std::string filename(root->_filename);
					filename += '.';
					filename += fileSuffix;
					
					stream.open(filename.c_str(),std::ofstream::out|std::ofstream::trunc);
				}
	
				writeHeader(stream,root->_ioMode,root->_header);
			}

			// ***

			if (!stream.good()) throw std::runtime_error("File_Io::LogResults::thread_writeResults(): Error writing data!");

			if (root->_ioMode == BINARY)
			{
				File_Utils::binWrite<unsigned int>(stream,obj->eventIndex);

				File_Utils::binWrite<short>(stream,obj->nodeId.toValue());
				File_Utils::binWrite<unsigned int>(stream,obj->nodeNumber.toValue());

				File_Utils::binWrite<float>(stream,obj->voltage);
				
				File_Utils::binWrite(stream,obj->trajectory.data().front().position());
				File_Utils::binWrite(stream,obj->trajectory.data().back().position());
				
				File_Utils::binWrite<float>(stream,obj->trajectory.data().back().time()*obj->timeScale);
				
				File_Utils::binWrite<float>(stream,obj->probability);

				File_Utils::binWrite<float>(stream,obj->potentialBefore);
				File_Utils::binWrite(stream,obj->fieldBefore);
				
				File_Utils::binWrite<float>(stream,obj->potentialAfter);
				File_Utils::binWrite(stream,obj->fieldAfter);
				
				File_Utils::binWrite(stream,obj->normal);
				File_Utils::binWrite(stream,obj->apex);
			}
			else
			{
				const char spacer = '\t';
				const char endline = '\n';

				stream << obj->eventIndex << spacer;

				stream << obj->nodeId.toValue() << spacer;
				stream << obj->nodeNumber.toValue() << spacer;

				stream << obj->voltage << spacer;

				for (short i = 0; i < 3; i++)
					stream << obj->trajectory.data().front().position(i) << spacer;
				for (short i = 0; i < 3; i++)
					stream << obj->trajectory.data().back().position(i) << spacer;

				stream << obj->trajectory.data().back().time() * obj->timeScale << spacer;

				stream << obj->probability << spacer;

				stream << obj->potentialBefore << spacer;
				for (short i = 0; i < 3; i++)
					stream << obj->fieldBefore[i] << spacer;
				
				stream << obj->potentialAfter << spacer;
				for (short i = 0; i < 3; i++)
					stream << obj->fieldAfter[i] << spacer;

				for (short i = 0; i < 3; i++)
					stream << obj->normal[i] << spacer;

				for (short i = 0; i < 2; i++)
					stream << obj->apex[i] << spacer;

				stream << obj->apex[2] << endline;
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
