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

#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <list>
#include <stdexcept>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include "version.h"

#include "info.h"

#include "commands.h"

#include "surface_3d.h"

#include "fieldline_3d.h"
#include "trajectory_3d.h"

#include "system_3d.h"

#include "file_io.h"
#include "logResults.h"
#include "logTrajectory.h"
#include "logSurface.h"

#include "debug.h"

const char* tapSimVersion = TAPSIM_VERSION;
const char* tapSimRevision = TAPSIM_REVISION;

const char* tapSim_distVersion = TAPSIM_DIST_VERSION;
const char* tapSim_distDate = TAPSIM_DIST_DATE;
const char* tapSim_distTime = TAPSIM_DIST_TIME;

const char* uniqueDate = __DATE__;
const char* uniqueTime = __TIME__;

const char* defaultIniFilename = "tapsim.ini";

void doMakeIni(const CommandLineOptions& options)
{
	info::open("make-ini-mode");

	std::string filename;
	if (options.makeIniParams.iniFilename.empty())
		filename = defaultIniFilename;
	else
		filename = options.makeIniParams.iniFilename;

	info::begin() << "Writing internal default data into initialization file! Filenname is \"";
	info::out() << filename << "\"" << std::endl;

	std::list<File_Io::KeyValue> iniParams;
	makeIniList(&iniParams);

	File_Io::writeInitialization(filename.c_str(),iniParams);

	info::close("make-ini-mode");
}

void doRelaxation(const CommandLineOptions& options)
{
	const RelaxationCommands& localOptions = options.relaxationParams;

	// ***

	info::open("relaxation-mode");
	
	System_3d system;
	system.gridTable.setThreadNum(options.threadNum);
	
	unsigned int dumpIndex(0);

	if (!localOptions.iDumpFile.empty())
	{
		info::begin() << "Reading dump from file \"" << localOptions.iDumpFile << "\"" << std::endl;
		File_Io::readSystem(localOptions.iDumpFile.c_str(),&system,&dumpIndex);
	}
	
	if (!localOptions.iConfigFile.empty() && !localOptions.iNodeFile.empty())
		Process::initialization(localOptions.iConfigFile.c_str(),localOptions.iNodeFile.c_str(),&system);

	if (localOptions.resetPotentials)
	{
		system.gridTable.randomReset(system.configTable);
		info::begin() << "Potential values are reset." << std::endl;
	}

	Process::relaxation(localOptions.threshold,localOptions.cycleSize,localOptions.queueSize,&system);

	if (!localOptions.oConfigFile.empty())
	{
		info::begin() << "Saving configuration in file \"" << localOptions.oConfigFile << "\" ";

		if (localOptions.oConfigMode == File_Io::ASCII)
			info::out() << "(ascii mode)" << std::endl;
		else
			info::out() << "(binary mode)" << std::endl;

		File_Io::writeConfig(localOptions.oConfigFile.c_str(),system.configTable,localOptions.oConfigMode);
	}
	
	if (!localOptions.oNodeFile.empty())
	{
		info::begin() << "Saving relaxed nodes in file \"" << localOptions.oNodeFile << "\" ";

		if (localOptions.oNodeMode == File_Io::ASCII)
			info::out() << "(ascii mode)" << std::endl;
		else
			info::out() << "(binary mode)" << std::endl;

		File_Io::writeGeomGrid(localOptions.oNodeFile.c_str(),system.geomTable,system.gridTable,localOptions.oNodeMode);
	}

	if (!localOptions.oDumpFile.empty())
	{
		info::begin() << "Saving dump in file \"" << localOptions.oDumpFile << "\"" << std::endl;
		File_Io::writeSystem(localOptions.oDumpFile.c_str(),system,dumpIndex);
	}

	info::close("relaxation-mode");
}

void doEvaporation(CommandLineOptions& options)
{
	const EvaporationCommands& localOptions = options.evaporationParams;

	// ***

	info::open("evaporation-mode");
	
	System_3d system;
	system.gridTable.setThreadNum(options.threadNum);

	if (!localOptions.iDumpFile.empty() && !localOptions.iConfigFile.empty())
	{
		info::begin() << "Reading dump from file \"" << localOptions.iDumpFile << "\"" << std::endl;
		File_Io::readSystem(localOptions.iDumpFile.c_str(),&system,0);
		
		info::begin() << "Reading configuration from file \"" << localOptions.iConfigFile << "\"" << std::endl;
		File_Io::readConfig(localOptions.iConfigFile.c_str(),&system.configTable);

		time_t startTime = std::time(0);

		system.gridTable.fastSync(system.geomTable,system.configTable);

		info::begin() << "Syncronizing data structures ";
		info::out() << "(elapsed time = " << std::difftime(std::time(0),startTime) / 60.0f;
		info::out() << " min)" << std::endl;
	}
	else if (!localOptions.iConfigFile.empty() && !localOptions.iNodeFile.empty())
		Process::initialization(localOptions.iConfigFile.c_str(),localOptions.iNodeFile.c_str(),&system);
	else
		throw std::runtime_error("doEvaporation()");

	if (localOptions.resetPotentials)
	{
		system.gridTable.randomReset(system.configTable);
		info::begin() << "Potential values are reset." << std::endl;
	}

	if (localOptions.initRelax) Process::relaxation(localOptions.initThreshold,localOptions.initCycleSize,localOptions.initQueueSize,&system);

	Process::evaporation(localOptions.evapParams,options.outputHeader,&system);

	info::close("evaporation-mode");
}

void doResumption(CommandLineOptions& options)
{
	const ResumptionCommands& localOptions = options.resumptionParams;

	// ***

	info::open("resumption-mode");

	System_3d system;
	system.gridTable.setThreadNum(options.threadNum);

	if (!localOptions.iDumpFile.empty())
	{
		unsigned int dumpIndex;
		
		info::begin() << "Reading dump from file \"" << localOptions.iDumpFile << "\"" << std::endl;
		File_Io::readSystem(options.resumptionParams.iDumpFile.c_str(),&system,&dumpIndex);

		// ***

		info::begin() << "Initial event index is set to ";

		if (std::numeric_limits<unsigned int>::max() == options.resumptionParams.evapParams.initEventCnt)
		{
			options.resumptionParams.evapParams.initEventCnt = dumpIndex;
			info::out() << "\"" << dumpIndex << "\" (read from dump)";
		}
		else
			info::out() << "\"" << options.resumptionParams.evapParams.initEventCnt << "\"";

		info::out() << std::endl;
	}

	Process::evaporation(localOptions.evapParams,options.outputHeader,&system);
	
	info::close("resumption-mode");
}

int main(const int argc, char** argv)
{
	std::srand(99); // initialize random number generator

	// ***

	std::cout << "TAPSim " << tapSimVersion << std::endl;
	std::cout << std::endl;
	std::cout << "Copyright (C) 2012 Christian Oberdorfer" << std::endl;
	std::cout << std::endl;

	// ***

	CommandLineOptions cmdOptions;
	
	try
	{
		parseCommandLine(argc,argv,&cmdOptions);
	}
	catch (int& status)
	{
		return status;
	}
	catch (std::runtime_error& cmdError)
	{
		throw cmdError;
	}
	catch (...)
	{
		std::cout << "An error occurred while parsing the command line!" << std::endl;
		std::cout << std::endl;
		return -1;
	}

	// ***

	info::begin() << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
	info::open("tapsim-log-output");

	info::open("command-line-options");

	info::begin() << "number of threads: " << cmdOptions.threadNum << std::endl;
	info::begin() << "operation-mode: " << cmdOptions.mode << std::endl;
	
	if (cmdOptions.iniFilename.empty())
		info::begin() << "initialization file: none (internal defaults)" << std::endl;
	else if (cmdOptions.iniFilename == defaultIniFilename)
		info::begin() << "initialization file: '" << defaultIniFilename << "' (external defaults)" << std::endl;
	else
		info::begin() << "initialization file: '" << cmdOptions.iniFilename << "'" << std::endl;

	switch (cmdOptions.mode)
	{
		case CommandLineOptions::RELAXATION:
		{
			info::open("relaxation-mode-options");

			info::begin() << "input configuration: \"" << cmdOptions.relaxationParams.iConfigFile << "\"" << std::endl;
			info::begin() << "output configuration: \"" << cmdOptions.relaxationParams.oConfigFile << "\"" << std::endl;
			info::begin() << "output configuration mode: \"" << cmdOptions.relaxationParams.oConfigMode << "\"" << std::endl;

			info::begin() << "input nodes: \"" << cmdOptions.relaxationParams.iNodeFile << "\"" << std::endl;
			info::begin() << "output nodes: \"" << cmdOptions.relaxationParams.oNodeFile << "\"" << std::endl;
			info::begin() << "output nodes mode: \"" << cmdOptions.relaxationParams.oNodeMode << "\"" << std::endl;

			info::begin() << "input raw file: \"" << cmdOptions.relaxationParams.iDumpFile << "\"" << std::endl;
			info::begin() << "output raw file: \"" << cmdOptions.relaxationParams.oDumpFile << "\"" << std::endl;
		
			info::begin() << "relaxation threshold: " << cmdOptions.relaxationParams.threshold << std::endl;
			info::begin() << "iteration cycle-size: " << cmdOptions.relaxationParams.cycleSize << std::endl;
			info::begin() << "deviation queue-size: " << cmdOptions.relaxationParams.queueSize << std::endl;
			
			info::begin() << "reset-potentials: " << cmdOptions.relaxationParams.resetPotentials << std::endl;

			info::close("relaxation-mode-options");

			break;
		}
		case CommandLineOptions::EVAPORATION:
		{
			info::open("evaporation-mode-options");

			info::begin() << "input dump file: \"" << cmdOptions.evaporationParams.iDumpFile << "\"" << std::endl;

			info::begin() << "input configuration: \"" << cmdOptions.evaporationParams.iConfigFile << "\"" << std::endl;
			info::begin() << "input nodes: \"" << cmdOptions.evaporationParams.iNodeFile << "\"" << std::endl;

			info::begin() << "do initial relaxation: " << cmdOptions.evaporationParams.initRelax << std::endl;
			info::begin() << "initial relaxation threshold: " << cmdOptions.evaporationParams.initThreshold << std::endl;
			info::begin() << "initial realxation iteration cycle-size: " << cmdOptions.evaporationParams.initCycleSize << std::endl;
			info::begin() << "initial relaxation deviation queue-size: " << cmdOptions.evaporationParams.initQueueSize << std::endl;

			info::begin() << "reset-potentials: " << cmdOptions.evaporationParams.resetPotentials << std::endl;

			info::open("evap-options");
			Debug::printEvapOptions(&info::out(),cmdOptions.evaporationParams.evapParams);
			info::close("evap-options");

			info::close("evaporation-mode-options");

			break;
		}
		case CommandLineOptions::RESUMPTION:
		{
			info::open("resumption-mode-options");

			info::begin() << "input dump file: \"" << cmdOptions.resumptionParams.iDumpFile << "\"" << std::endl;

			info::open("evap-options");
			Debug::printEvapOptions(&info::out(),cmdOptions.resumptionParams.evapParams);
			info::close("evap-options");

			info::close("resumption-mode-options");

			break;
		}
	}
	
	info::close("command-line-options");

	// ***

	#ifdef DEMO_MESH_MODE
	info::begin() << "***** !!!!! >>>>> DEMO MESH MODE <<<<< !!!!! *****" << std::endl;
	#endif

	// ***
	
	switch (cmdOptions.mode)
	{
		case CommandLineOptions::MAKE_INI:
			doMakeIni(cmdOptions);
			break;
		case CommandLineOptions::RELAXATION:
			doRelaxation(cmdOptions);
			break;
		case CommandLineOptions::EVAPORATION:
			doEvaporation(cmdOptions);
			break;
		case CommandLineOptions::RESUMPTION:
			doResumption(cmdOptions);
		default:
			break;
	}
	
	info::close("tapsim");

	return 0;
}
