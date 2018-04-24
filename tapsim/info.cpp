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

#include "info.h"

#include <iostream>

#include <vector.h>

static std::ostream* logStream = &std::cout;

static int indentLevel = 0;

namespace
{
	void indent()
	{
		for (int i = 0; i < indentLevel; i++)
			info::out() << '\t';
	}
}

std::ostream& info::stream()
{
	return *logStream;
}

void info::setStream(std::ostream* stream)
{
	logStream = stream;
}

void info::open(const char* tag, const char* attributes)
{
	indent();
	out() << "<" << tag;
	if (attributes != 0) out() << " " << attributes;
	out() << ">" << std::endl;

	indentLevel++;
}

void info::close(const char* tag)
{
	indentLevel--;

	indent();
	out() << "</" << tag << ">";
	out() << std::endl;
}

std::ostream& info::begin()
{
	indent();
	return *logStream;
}

std::ostream& info::out()
{
	return *logStream;
}

// void info::log_commandLine(const cmdLineOptions& options)
// {
// 	open("log_entry","type=\"*** command line options ***\"");
// 
// 	begin() << "threadEnable = " << options.threadEnable << std::endl;
// 	begin() << "threadNum = " << options.threadNum << std::endl;
// 	begin() << "logDelay = " << options.logDelay << "s" << std::endl;
// 	begin() << "logInterval = " << options.logInterval << std::endl;
// 
// 	begin() << "logFilename = ";
// 	if (options.logFilename.empty())
// 		out() << "*** standard out ***";
// 	else
// 		out() << options.logFilename;
// 	out() << std::endl;
// 
// 	begin() << "relaxationThreshold = " << options.relaxationThreshold << std::endl;
// 	begin() << "SOR-Factor = " << options.sor_factor << std::endl;
// 	begin() << "trajectoryErrorThreshold = " << options.trajectoryErrorThreshold << std::endl;
// 	begin() << "evaporation-limit = " << options.evaporationLimit << std::endl;
// 	begin() << "layer-shrink-limit = " << options.shrinkLimit << std::endl;
// 
// 	begin() << "probabbility estimation = ";
// 	switch (options.probabilityMethod)
// 	{
// 		case cmdLineOptions::LINEAR:
// 			out() << "LINEAR";
// 			break;
// 		case cmdLineOptions::BOLTZMANN:
// 			out() << "BOLTZMANN";
// 			break;
// 		default:
// 			out() << "UNKNOWN / ERROR";
// 	}
// 	out() << std::endl;
// 
// 	begin() << "monteCarlo = " << options.monteCarlo << std::endl;
// 	begin() << "fixed initial position = " << options.fixedInitialPosition << std::endl;
// 	begin() << "initial velocity = " << options.initialVelocity << std::endl;
// 	begin() << "temperature = " << options.temperature << std::endl;
// 
// 	close("log_entry");
// }

// void info::log_cellData(const cellData& cell)
// {
// 	open("log_entry","type=\"*** cellData ***\"");
// 
// 	begin() << "id = " << cell.id().id() << std::endl;
// 	begin() << "weight = " << cell.weight() << std::endl;
// 	begin() << "name = \"" << cell.name() << "\"" << std::endl;
// 	begin() << "charge = " << cell.charge() << " C" << std::endl;
// 	begin() << "epsilon = " << cell.epsilon() << std::endl;
// 
// 	begin() << "properties = ";
// 	switch (cell.properties())
// 	{
// 		case gridSystem::potData::DEFAULT:
// 			out() << "DEFAULT";
// 			break;
// 		case gridSystem::potData::DIRICHLET:
// 			out() << "DIRICHLET";
// 			break;
// 		case gridSystem::potData::NEUMANN:
// 			out() << "NEUMANN";
// 			break;
// 	}
// 	out() << std::endl;
// 
// 	begin() << "phi = " << cell.phi() << " V" << std::endl;
// 	begin() << "removable = " << cell.removable() << std::endl;
// 	begin() << "mass = " << cell.mass() << " amu" << std::endl;
// 	begin() << "charge state = " << cell.evapCharge() << std::endl;
// 	begin() << "evaporation field strength = " << cell.evapField() << " V/m" << std::endl;
// 	begin() << "activation energy = " << cell.evapEnergy() << " eV" << std::endl;
// 
// 	close("log_entry");
// }

// void info::log_cellList(const cellList& list)
// {
// 	open("log_entry","type=\"*** cellList ***\"");
// 
// 	const std::list<cellId> ids = list.ids();
// 	
// 	for (std::list<cellId>::const_iterator i = ids.begin(); i != ids.end(); i++)
// 		log_cellData(list[*i]);
// 
// 	close("log_entry");
// }

// void info::log_cellConfiguration(const cellConfiguration& config)
// {
// 	open("log_entry","type=\"*** cellConfiguration ***\"");
// 
// 	begin() << "size (cell units) = ";
// 	out() << config.size(Vector3d<unsigned long>::X) << " / ";
// 	out() << config.size(Vector3d<unsigned long>::Y) << " / ";
// 	out() << config.size(Vector3d<unsigned long>::Z);
// 	out() << std::endl;
// 
// 	begin() << "extents (real units) = ";
// 	out() << config.extents(Vector3d<unsigned long>::X) << " / ";
// 	out() << config.extents(Vector3d<unsigned long>::Y) << " / ";
// 	out() << config.extents(Vector3d<unsigned long>::Z) << " m³";
// 	out() << std::endl;
// 
// 	begin() << "temperature = " << config.temperature() << " K" << std::endl;
// 
// 	begin() << "center-of-projection (cell units) = ";
// 	out() << config.centerOfProjection(Vector3d<unsigned long>::X) << " / ";
// 	out() << config.centerOfProjection(Vector3d<unsigned long>::Y) << " / ";
// 	out() << config.centerOfProjection(Vector3d<unsigned long>::Z);
// 	out() << std::endl;
// 
// 	begin() << "flight distance (cell units) = " << config.flightDistance() << std::endl;
// 
// 	begin() << "single cell size (computing-grid) = " << config.cellSize() << std::endl;
// 	begin() << "single cell extents (real units) = " << config.cellExtents() << " m" << std::endl;
// 
// 	log_cellList(config.description());
// 
// 		open("statistics");
// 	
// 		const std::list<cellId> ids = config.description().ids();
// 	
// 		for (std::list<cellId>::const_iterator i = ids.begin(); i != ids.end(); i++)
// 		{
// 			begin() << "Number of cells with id \"" << i->id() << "\"";
// 			out() << " / name \"" << config.description(*i).name() << "\": ";
// 			out() << config.countIds(*i);
// 			out() << std::endl;
// 		}
// 	
// 		begin() << "Number of cells with invalid id: " << config.countIds(cellId::invalid()) << std::endl;
// 	
// 		{
// 			unsigned long totalCells = config.size(Vector3d<unsigned long>::X);
// 			totalCells *= config.size(Vector3d<unsigned long>::Y);
// 			totalCells *= config.size(Vector3d<unsigned long>::Z);
// 
// 			begin() << "Total number of cell sites: " << totalCells << std::endl;
// 		}
// 	
// 		{
// 			unsigned long cellCnt(0);
// 			unsigned long removableCnt(0);
// 		
// 			Vector3d<unsigned long> cellPos;
// 			for (cellPos.x() = 0; cellPos.x() < config.size(Vector3d<unsigned long>::X); cellPos.x()++)
// 			{
// 				for (cellPos.y() = 0; cellPos.y() < config.size(Vector3d<unsigned long>::Y); cellPos.y()++)
// 				{
// 					for (cellPos.z() = 0; cellPos.z() < config.size(Vector3d<unsigned long>::Z); cellPos.z()++)
// 					{
// 						if (config.id(cellPos) == cellId::invalid()) continue;
// 	
// 						if (config.description(cellPos).removable() == true) removableCnt++;
// 	
// 						cellCnt++;
// 					}
// 				}
// 			}
// 	
// 			begin() << "Occupied cells counted: " << cellCnt << std::endl;
// 			begin() << "Removable occupied cells counted: " << removableCnt << std::endl;
// 		}
// 	
// 		close("statistics");
// 
// 	close("log_entry");
// }

// void info::log_gridSystem(const gridSystem& grid)
// {
// 	open("log_entry","type=\"*** gridSystem ***\"");
// 
// 	begin() << "SOR-Factor = " << grid.relaxationFactor() << std::endl;
// 	begin() << "Deviation threshold for relaxation (relative value) =  " << grid.deviationLimit() << std::endl;
// 	begin() << "Iteration limit = " << grid.iterationLimit() << std::endl;
// 
// 	begin() << "Threads enabled = ";
// 	if (grid.threadEnabled())
// 		out() << "true (" << grid.threadNumber() << ")";
// 	else
// 		out() << "false";
// 	out() << std::endl;
// 
// 	begin() << "Slice axis = " << grid.sliceAxis() << std::endl;
// 	begin() << "Slice size =  " << grid.sliceSize() << std::endl;
// 
// 	begin() << "Knot size = ";
// 	out() << grid.knotSize(Vector3d<unsigned long>::X) << " / ";
// 	out() << grid.knotSize(Vector3d<unsigned long>::Y) << " / ";
// 	out() << grid.knotSize(Vector3d<unsigned long>::Z);
// 	out() << std::endl;
// 
// 	begin() << "Vol-cell size = ";
// 	out() << grid.volSize(Vector3d<unsigned long>::X) << " / ";
// 	out() << grid.volSize(Vector3d<unsigned long>::Y) << " / ";
// 	out() << grid.volSize(Vector3d<unsigned long>::Z);
// 	out() << std::endl;
// 
// 	begin() << "Mesh extents (real coordinates) = ";
// 	out() << grid.meshExtents(Vector3d<double>::X) << " / ";
// 	out() << grid.meshExtents(Vector3d<double>::Y) << " / ";
// 	out() << grid.meshExtents(Vector3d<double>::Z) << " m³";
// 	out() << std::endl;
// 
// 	close("log_entry");
// }
