#ifndef TAPSIM_DEBUG_H
#define TAPSIM_DEBUG_H

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

#include "system_3d.h"
#include "surface_3d.h"
#include "process.h"

namespace Debug
{
	void printNode(std::ostream*, const int, const System_3d&);
	void printTetrahedron(std::ostream*, const int, const System_3d&);

	void printVoronoiFace(std::ostream*, const Geometry_3d::VoronoiFace&);
	void printVoronoiCell(std::ostream*, const Geometry_3d::VoronoiCell&);
	
	// ***
	
	void printEvapOptions(std::ostream*, const Process::EvaporationOptions&);
	
	// ***

	void writeData_PARAVIEW(const char*, const System_3d&);
	void writeData_PARAVIEW(const char*, const System_3d&, const Geometry_3d::VoronoiFace&);
	void writeData_PARAVIEW(const char*, const System_3d&, const Geometry_3d::VoronoiCell&);

	void writeData_PARAVIEW_potentials(const char*, const System_3d&);
	void writeData_PARAVIEW_voronoiCells(const char*, const System_3d&);
	//void writeData_PARAVIEW_constrainedVoronoiCells(const char*, const System_3d&);
	//void writeData_PARAVIEW_delaunayTetrahedra(const char*, const System_3d&);

	// ***

	void writeData_PARAVIEW_surfaceNodes(const char*, const Surface_3d::Table&, const System_3d&);
	void writeData_PARAVIEW_surfaceCells(const char*, const Surface_3d::Table&, const System_3d&);

	void writeData_PARAVIEW_surfaceFieldlines(const char*, const Surface_3d::Table&, const System_3d&, const float, const float =1.0f);
	void writeData_PARAVIEW_surfaceTrajectories(const char*, const Surface_3d::Table&, const System_3d&, const float, const float =1e-3f, const float =1.0f);

	// ***

	void writeData(const char*, const System_3d&, const char spacer =',', const char endline ='\n');

	void writeData_surfaceTrajectories(const char*, const Surface_3d::Table&, const System_3d&, const float =1e-6f, const float =1.0f);

	// ***

	void testCoupling(std::ostream&, const Grid_3d::Table&);
	void testNeighbourDistance(const char*, const Grid_3d::Table&);

	// ***
	
	void surfaceSize_statistics(const System_3d&);
}

#endif
