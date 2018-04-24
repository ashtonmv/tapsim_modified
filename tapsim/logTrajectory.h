#ifndef TAPSIM_LOG_TRAJECTORY_H
#define TAPSIM_LOG_TRAJECTORY_H

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

#include <string>
#include <queue>

#include <pthread.h>

#include "file_io.h"
#include "trajectory_3d.h"

namespace File_Io
{
	class LogTrajectory
	{
		public:
			
			LogTrajectory(const char* =0, const int =ASCII, const unsigned int =0, const unsigned int =0, const unsigned int =0);
			~LogTrajectory();

			void init(const char*, const int =ASCII, const unsigned int =0, const unsigned int =0, const unsigned int =0);

			const char* filename() const { return _filename.c_str(); }
			
			void setDelay(const int);
			unsigned int delay() const { return _delay; }
			
			void setBuffer(const int value) { _bufferSize = value; }
			unsigned int buffer() const { return _bufferSize; }

			void setChunk(const unsigned int value) { _chunkSize = value; }
			unsigned int chunk() const { return _chunkSize; }

			void enable();
			void suspend();

			bool isEnabled() const { return _keepRunning; }
			
			void push(const unsigned int, const unsigned int, const float, const Trajectory_3d&);
			

		private:
			struct Dataset
			{
				unsigned int eventIndex;
				unsigned int nodeNumber;
				
				float timeScale;
				Trajectory_3d trajectory;
			};

			LogTrajectory(const LogTrajectory&);
			LogTrajectory& operator=(const LogTrajectory&);

			static void writeHeader(std::ofstream&, const int);
			static void* thread_writeTrajectory(void*);

			static const std::string _fileVersion;

			std::string _filename;

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
