#ifndef TAPSIM_SYSTEM_3D_H
#define TAPSIM_SYSTEM_3D_H

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

#include "configuration.h"
#include "geometry_3d.h"
#include "grid_3d.h"

struct System_3d
{
	System_3d() 
		: configTable(),
		  geomTable(),
		  gridTable()
	{}
	
	// ***
	
	Configuration::Table configTable;
	Geometry_3d::Table geomTable;
	Grid_3d::Table gridTable;
};

void preset_sphericalPotential(System_3d*, const Geometry_3d::Point&, const float, const float);

#endif