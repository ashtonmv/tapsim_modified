#ifndef TAPSIM_FILE_IO_H
#define TAPSIM_FILE_IO_H

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

#include <map>
#include <list>
#include <string>

#include "system_3d.h"

namespace File_Io
{
	enum { ASCII, BINARY };
	
	enum { NODES =0x01, TETRAHEDRA =0x02, VORONOI_CELLS =0x04, VORONOI_FACETS =0x08, VORONOI_VERTICES =0x10 };

	typedef std::pair<std::string,std::string> KeyValue;

	void writeInitialization(const char*, const std::list< File_Io::KeyValue >&, const char* =" = ");
	void readInitialization(const char*, std::map<std::string,std::string>*, const char* =" = ");

	void writeConfig(const char*, const Configuration::Table&, const int =ASCII);
	void readConfig(const char*, Configuration::Table*); 

	void writeGeometry(const char*, const Geometry_3d::Table&, const unsigned char =0xff, const int =ASCII);
	
	void writeGeomGrid(const char*, const Geometry_3d::Table&, const Grid_3d::Table&, const int =ASCII, const int withNumbers =1, const int withPotentials =1);
	void readGeomGrid(const char*, Geometry_3d::Table*, Grid_3d::Table*, int* withNumbers =0, int* withPotentials =0);

	void writeSystem(const char* filename, const System_3d& system, const unsigned int index);
	void readSystem(const char*, System_3d*, unsigned int* =0);
}

#endif
