#ifndef TAPSIM_INFO_H
#define TAPSIM_INFO_H

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

#include <ostream>

struct cmdLineOptions;

class cellData;
class cellList;
class cellConfiguration;

class gridSystem;

namespace info
{
	std::ostream& stream();
	void setStream(std::ostream* stream);

	void open(const char* tag, const char* attributes =0);
	void close(const char* tag);

	std::ostream& begin();
	std::ostream& out();

	void log_commandLine(const cmdLineOptions&);

	//void log_cellData(const cellData& cell);
	//void log_cellList(const cellList& list);
	//void log_cellConfiguration(const cellConfiguration& config);

	//void log_gridSystem(const gridSystem& grid);
}

#endif
