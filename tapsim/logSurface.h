#ifndef TAPSIM_LOG_SURFACE_H
#define TAPSIM_LOG_SURFACE_H

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
#include "system_3d.h"
#include "surface_3d.h"

#include <pthread.h>

#include <map>
#include <queue>
#include <string>

namespace File_Io
{
	class LogSurface
	{
		friend void write_surfaceCells(const char*, const Surface_3d::Table&, const System_3d&, const int, const bool);
		
		public:
			LogSurface(const char* name =0, const System_3d* =0, const int =ASCII, const unsigned int =0, const unsigned int =0, const bool =false);
			~LogSurface();

			void init(const char* name, const System_3d* =0, const int =ASCII, const unsigned int =0, const unsigned int =0, const bool =false);

			const char* filename() const { return _filename.c_str(); }

			void setDelay(const int);
			unsigned int delay() const { return _delay; }

			void setBuffer(const int value) { _bufferSize = value; }
			unsigned int buffer() const { return _bufferSize; }

			void setWithNeighbours(const bool value) { _withNeighbours = value; }
			bool withNeighbours() const { return _withNeighbours; }

			void enable();
			void suspend();

			bool isEnabled() const { return _keepRunning; }

			void push(const int, const Surface_3d::Table&);

		private:
			struct NodeData
			{
				float phi;
				MathVector3d<float> field;
				MathVector3d<float> normal;
			};

			struct CellData
			{
				int index;
				
				short type;
				unsigned int number;
				
				float probability;

				std::list<int> neighbours;
			};

			struct Dataset
			{
				int id;
				std::map<int,CellData> cells;
				std::map<int,NodeData> nodes;
				
			};

			static void prepare(const int, const Surface_3d::Table&, const System_3d&, Dataset*, const bool);
			static void writer(const char*, const int, const Dataset&, const bool);
			
			// ***

			LogSurface(const LogSurface&);
			LogSurface& operator=(const LogSurface&);

			static void* thread_writeSurface(void*);

			static const std::string _fileVersion;

			std::string _filename;

			const System_3d* _system;
			
			int _ioMode;
			unsigned int _delay;
			unsigned int _bufferSize;

			bool _withNeighbours;

			bool _keepRunning;

			pthread_t _threadId;
			pthread_mutex_t _queueMutex;
			pthread_cond_t _queueCondition;

			std::queue<Dataset*> _queue;
	};
	
	void write_surfaceNodes(const char*, const Surface_3d::Table&, const System_3d&, const int =ASCII);
	void write_surfaceCells(const char*, const Surface_3d::Table&, const System_3d&, const int =ASCII, const bool =false);

	void write_surfaceTrajectories(const char*, const Surface_3d::Table&, const System_3d&, const float =1e-6f, int =ASCII);
}

#endif
