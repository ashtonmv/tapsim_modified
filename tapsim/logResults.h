#ifndef TAPSIM_LOG_RESULTS_H
#define TAPSIM_LOG_RESULTS_H

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

#include <queue>
#include <string>
#include <fstream>

#include <pthread.h>

#include "file_io.h"
#include "geometry_3d.h"
#include "configuration.h"
#include "trajectory_3d.h"

namespace File_Io
{
	class LogResults
	{
		public:
			struct Dataset
			{
				float simTime;
				
				unsigned int eventIndex;
	
				int nodeIndex;
				Configuration::NodeId nodeId;
				Configuration::NodeNumber nodeNumber;

				float probability;
				
				float potentialBefore;
				MathVector3d<float> fieldBefore;
				
				float potentialAfter;
				MathVector3d<float> fieldAfter;
				
				float voltage;
				
				Geometry_3d::Point apex;
				Geometry_3d::Point normal;

				float timeScale;
				Trajectory_3d trajectory;
			};
	
			LogResults(const char* name =0, const int =ASCII, const unsigned int =0, const unsigned int =0, const unsigned int =0);
			~LogResults();

			void init(const char* name, const int =ASCII, const unsigned int =0, const unsigned int =0, const unsigned int =0);

			const char* filename() const { return _filename.c_str(); }

			void setHeader(const char* value) { _header = value; }
			const std::string& header() const { return _header; }

			void setDelay(const int);
			unsigned int delay() const { return _delay; }

			void setChunk(const unsigned int value) { _chunkSize = value; }
			unsigned int chunk() const { return _chunkSize; }

			void setBuffer(const int value) { _bufferSize = value; }
			unsigned int buffer() const { return _bufferSize; }

			void enable();
			void suspend();

			bool isEnabled() const { return _keepRunning; }
			
			void push(const Dataset& data);
	
		private:
			LogResults(const LogResults&);
			LogResults& operator=(const LogResults&);

			static void writeHeader(std::ofstream&, const int, const std::string&);
			static void* thread_writeResults(void*);

			static const std::string _fileVersion;

			std::string _filename;

			std::string _header;

			int _ioMode;
			unsigned int _delay;
			unsigned int _chunkSize;
			unsigned int _bufferSize;
	
			bool _keepRunning;
	
			pthread_t _threadId;
			pthread_mutex_t _queueMutex;
			pthread_cond_t _queueCondition;

			std::queue<Dataset*> _queue;
	};
	
}

#endif
