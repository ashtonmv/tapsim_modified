#ifndef TAPSIM_LOG_GRID_H
#define TAPSIM_LOG_GRID_H

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

#include "file_io.h"
#include "grid_3d.h"
#include "configuration.h"

#include <pthread.h>

#include <string>
#include <vector>
#include <queue>

namespace File_Io
{
	class LogGrid
	{
		public:
			LogGrid(const char* name =0, const Geometry_3d::Table* =0, const Configuration::Table* =0, int =ASCII, const unsigned int =0, const unsigned int =0);
			~LogGrid();

			void init(const char* name, const Geometry_3d::Table*, const Configuration::Table*, const int =ASCII, const unsigned int =0, const unsigned int =0);

			const char* filename() const { return _filename.c_str(); }

			void setDelay(const int);
			unsigned int delay() const { return _delay; }

			void setBuffer(const int value) { _bufferSize = value; }
			unsigned int buffer() const { return _bufferSize; }

			void enable();
			void suspend();

			bool isEnabled() const { return _keepRunning; }

			void push(const int, const Grid_3d::Table&);

		private:
			struct NodeData
			{
				short id;
				unsigned int number;

				float phi;
				MathVector3d<float> field;
				MathVector3d<float> force;
			};

			struct Dataset
			{
				int index;
				std::vector<NodeData> data;
			};
			
			// ***

			LogGrid(const LogGrid&);
			LogGrid& operator=(const LogGrid&);

			static void* thread_writeGrid(void*);

			static const std::string _fileVersion;

			std::string _filename;

			const Geometry_3d::Table* _geomTable;
			const Configuration::Table* _configTable;

			int _ioMode;
			unsigned int _delay;
			unsigned int _bufferSize;

			bool _keepRunning;

			pthread_t _threadId;
			pthread_mutex_t _queueMutex;
			pthread_cond_t _queueCondition;

			std::queue<Dataset*> _queue;
	};

	void writeGrid(const char*, const Grid_3d::Table&, const int =ASCII);
}

#endif
