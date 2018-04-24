#ifndef TAPSIM_COMMANDS_H
#define TAPSIM_COMMANDS_H

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

#include <list>
#include <string>

#include "process.h"
#include "file_io.h"

struct MakeIniCommands
{
	MakeIniCommands();

	static void setDefaults(MakeIniCommands*);

	std::string iniFilename;
};

struct RelaxationCommands
{
	RelaxationCommands();
	
	static void setDefaults(RelaxationCommands*);
	
	std::string iConfigFile;
	std::string oConfigFile;
	int oConfigMode;

	std::string iNodeFile;
	std::string oNodeFile;
	int oNodeMode;

	std::string iDumpFile;
	std::string oDumpFile;

	float threshold;
	unsigned int cycleSize;
	unsigned int queueSize;

	int resetPotentials;
};

struct EvaporationCommands
{
	EvaporationCommands();

	static void setDefaults(EvaporationCommands*, const std::string&);

	std::string iDumpFile;
	
	std::string iConfigFile;
	std::string iNodeFile;

	// ***
	
	int initRelax;
	
	float initThreshold;
	unsigned int initCycleSize;
	unsigned int initQueueSize;
	
	int resetPotentials;

	// ***

	Process::EvaporationOptions evapParams;
};

struct ResumptionCommands
{
	ResumptionCommands();

	static void setDefaults(ResumptionCommands*, const std::string&);

	std::string iDumpFile;

	// ***

	Process::EvaporationOptions evapParams;
};


struct CommandLineOptions
{
	enum { MAKE_INI, RELAXATION, EVAPORATION, RESUMPTION };

	CommandLineOptions();

	std::string outputHeader;
	
	int threadNum;

	int mode;

	std::string iniFilename;

	// *** special 'make-ini' mode parameters
	MakeIniCommands makeIniParams;

	// *** special 'relaxation' mode parameters
	RelaxationCommands relaxationParams;
	
	// *** special 'evaporation' mode parameters
	EvaporationCommands evaporationParams;
	
	// *** special 'resumption' mode parameters
	ResumptionCommands resumptionParams;
};

void parseCommandLine(int, char**, CommandLineOptions*);

void makeIniList(std::list<File_Io::KeyValue>*);

#endif
