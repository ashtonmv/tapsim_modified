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

#include "commands.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <limits>

#include <cstdio>
#include <cstring>
#include <cstdlib>

#include <getopt.h>

#include "file_io.h"
#include "surface_3d.h"
#include "trajectory_3d.h"

extern const char* tapSimVersion;
extern const char* tapSimRevision;

extern const char* tapSim_distVersion;
extern const char* tapSim_distDate;
extern const char* tapSim_distTime;

extern const char* uniqueDate;
extern const char* uniqueTime;

extern const char* defaultIniFilename;

namespace
{
	void parseMakeIniCommands(int argc, char** argv, MakeIniCommands* options)
	{
		/* ------------------------------------------------------------------
		*
		* Calling convention:
		*
		*	./tapsim [...] make-ini
		*	./tapsim [...] make-ini filename
		*
		* -------------------------------------------------------------------
		*/

		struct option localOptions[] =
		{
			// invalid global options in this context:
			{ "no-file-output", no_argument, 0, -1 },
			{ "write-ascii", no_argument, 0, -1 },
			{ "write-binary", no_argument,0, -1 },
			{ "threads", required_argument, 0, -1 },
			{ "ini-file", required_argument, 0, -1 },

			{ "help", no_argument, 0, 1 },

			{ 0, 0, 0, 0 }
		};

		MakeIniCommands::setDefaults(options);

		// ***

		optind = 1; // reset option index
		opterr = 1; // enable output of error messages by getopt()

		while (true)
		{
			int key;

			do
				key = getopt_long(argc,argv,"",localOptions,0);
			while (key == 0 || key == '?');

			if (key == -1) break;

			switch (key)
			{
				case 1: // --help
					std::cout << "'Make-ini'-mode syntax:" << std::endl;
					std::cout << std::endl;
					std::cout << "./tapsim make-ini" << std::endl;
					std::cout << "./tapsim make-ini filename" << std::endl;
					std::cout << std::endl;
					std::cout << "List of additional parameters in 'make-ini'-mode:" << std::endl;
					std::cout << std::endl;
					std::cout << "\t<< none >>" << std::endl;
					std::cout << std::endl;
					std::cout << "=> if no filename is supplied an initialization file with a default is created" << std::endl;
					std::cout << std::endl;
					std::cout << "The default initialization file  will be used for operation,  if it exists. In" << std::endl;
					std::cout << "the case  it doesn't exist and use of  a different file is not forced from the" << std::endl;
					std::cout << "command line, a fall-back on internal defaults occurs." << std::endl;
					std::cout << std::endl;
					throw int(0);
				default:
					std::string errorString("Unknown or invalid argument in command line: \"");
					errorString += optarg;
					errorString += "\"!";

					throw std::runtime_error(errorString);
			}
		}

		// *** process additional arguments
		
		if (argc - 1 - optind >= 0 && std::strcmp(argv[optind],"make-ini") != 0)
			throw std::runtime_error("parseMakeIniCommands()");

		if (argc - 1 - optind == 0)
			return;
		else if (argc - 1 - optind == 1)
		{
			options->iniFilename = argv[optind+1];
			return;
		}
		else
			throw std::runtime_error("Missing arguments for processing in 'make-ini'-mode!");
	}
	
	void parseRelaxationCommands(int argc, char** argv, RelaxationCommands* options)
	{
		/* ---------------------------------------------------------------------------
		* 
		* Calling convention:
		* 
		*	1) ./tapsim [...] relaxation iconfig inode
		*	2) ./tapsim [...] relaxation iconfig inode onode
		*	3) ./tapsim [...] relaxation --read-dump=? oconfig onode
		* 
		*       => optionally the "--write-dump=?" paramater may be added in each case
		* 
		* ----------------------------------------------------------------------------
		*/
		
		struct option localOptions[] =
		{
			// global options: recognized but skipped
			{ "no-file-output", no_argument, 0, 0 },
			{ "write-ascii", no_argument, 0, 0 },
			{ "write-binary", no_argument,0, 0 },
			{ "threads", required_argument, 0, 0 },

			// invalid global options in this context:
			{ "ini-file", required_argument, 0, -1 },

			// local options:
			{ "init-potentials", no_argument, &options->resetPotentials, 1},
			{ "reset-potentials", no_argument, &options->resetPotentials, 1 },

			{ "threshold", required_argument, 0, 1 },
			{ "cycle-size", required_argument, 0, 2 },
			{ "queue-size", required_argument, 0, 3 },
			{ "read-dump", required_argument, 0, 4 },
			{ "write-dump", required_argument, 0, 5 },

			{ "help", no_argument, 0, 6 },

			{ 0, 0, 0, 0 }
		};

		RelaxationCommands::setDefaults(options);
		
		// ***

		optind = 1; // reset option index 
		opterr = 1; // enable output of error messages by getopt()

		while (true)
		{	
			int key;
		
			do
				key = getopt_long(argc,argv,"",localOptions,0);
			while (key == 0 || key == '?');
		
			if (key == -1) break;

			switch (key)
			{
				case 1: // --threshold
					if (optarg == 0 || std::sscanf(optarg,"%f",&options->threshold) != 1)
						throw std::runtime_error("Error parsing 'threshold' argument!");
					break;
				case 2: // --cycle-size
					if (optarg == 0 || std::sscanf(optarg,"%u",&options->cycleSize) != 1)
						throw std::runtime_error("Error parsing 'cyclce-size' argument!");
					break;
				case 3: // --queue-size
					if (optarg == 0 || std::sscanf(optarg,"%u",&options->queueSize) != 1)
						throw std::runtime_error("Error parsing 'queue-size' argument!");
					break;
				case 4: // --read-dump
					if (optarg == 0) throw std::runtime_error("Error parsing 'read-dump' argument!");
					options->iDumpFile = optarg;
					break;
				case 5: // --write-dump
					if (optarg == 0) throw std::runtime_error("Error parsing 'write-dump' argument!");
					options->oDumpFile = optarg;
					break;
				case 6: // --help
					std::cout << "Relaxation mode syntax:" << std::endl;
					std::cout << std::endl;
					std::cout << "./tapsim relaxation --read-dump=? oconfig onode" << std::endl;
					std::cout << "./tapsim relaxation iconfig inode [onode]" << std::endl;
					std::cout << std::endl;
					std::cout << "List of additional parameters in relaxation mode:" << std::endl;
					std::cout << "\t--no-file-output, --write-binary, --write-ascii" << std::endl;
					std::cout << "\t--threads=?" << std::endl;
					std::cout << "\t--init-potentials, --reset-potentials" << std::endl;
					std::cout << "\t--threshold=?," << std::endl;
					std::cout << "\t--cycle-size=?," << std::endl;
					std::cout << "\t--queue-size=?," << std::endl;
					std::cout << "\t--write-dump=?," << std::endl;
					std::cout << std::endl;
					throw int(0);
				default:
					std::string errorString("Unknown or invalid argument in command line: \"");
					errorString += optarg;
					errorString += "\"!";
					
					throw std::runtime_error(errorString);
			}
		}

		// *** process additional arguments

		if (std::strcmp(argv[optind],"relaxation") != 0) throw std::runtime_error("parseRelaxationCommands()");

		switch (argc - 1 - optind)
		{
			case 2:
				if (options->iDumpFile.empty())
				{
					// *** ./tapsim [...] relaxation iconfig inode
					options->iConfigFile = argv[optind+1];
					options->iNodeFile = argv[optind+2];
					options->oNodeFile = argv[optind+2];
				}
				else
				{
					// *** ./tapsim [...] relaxation --read-dump [--write-dump] oconfig onode
					options->oConfigFile = argv[optind+1];
					options->oNodeFile = argv[optind+2];
				}
				break;
			case 3:
				if (options->iDumpFile.empty())
				{
					// *** ./tapsim [...] relaxation iconfig inode onode
					options->iConfigFile = argv[optind+1];
					options->iNodeFile = argv[optind+2];
					options->oNodeFile = argv[optind+3];
				}
				break;
			default:
				throw std::runtime_error("Missing arguments for processing in 'relaxation'-mode!");
		}
	}

	void parseEvaporationCommands(int argc, char** argv, EvaporationCommands* options, const std::string& iniFile)
	{
		/* -----------------------------------------------------------------------------
		* 
		* Calling convention:
		* 
		*	./tapsim [...] evaporation iconfig inode [results]
		*	./tapsim [...] evaporation --input-dump=? iconfig [results]
		* 
		* ------------------------------------------------------------------------------
		*/

		struct option localOptions[] =
		{
			// global options: recognized but skipped
			{ "no-file-output", no_argument, 0, 0 },
			{ "write-ascii", no_argument, 0, 0 },
			{ "write-binary", no_argument,0, 0 },
			{ "threads", required_argument, 0, 0 },
			{ "ini-file", required_argument, 0, 0 },

			// local options:
			{ "skip-relaxation", no_argument, &options->initRelax, 0 },
			{ "init-potentials", no_argument, &options->resetPotentials, 1},
			{ "reset-potentials", no_argument, &options->resetPotentials, 1 },

			{ "threshold", required_argument, 0, 1 },
			{ "cycle-size", required_argument, 0, 2 },
			{ "queue-size", required_argument, 0, 3 },

			{ "no-trajectories", no_argument, 0, 4 },
			{ "write-trajectories", required_argument, 0, 4 },

			{ "no-grids", no_argument, 0, 5 },
			{ "write-grids", required_argument, 0, 5 },

			{ "no-surfaces", no_argument, 0, 6 },
			{ "write-surfaces", required_argument, 0 , 6 },

			{ "no-chunks", no_argument, 0, 7 },
			{ "chunk-size", required_argument, 0, 7 },

			{ "snapshot-interval", required_argument, 0, 8 },

			{ "input-dump", required_argument, 0, 9 },

			{ "no-dump", no_argument, 0, 10 },
			{ "write-dump", required_argument, 0, 10 },
			{ "dump-interval", required_argument, 0, 11 },
			
			{ "vacuum-identifier", required_argument, 0, 12 },
			
			{ "trajectory-error-threshold", required_argument, 0, 13 },

			{ "shrink-limit", required_argument, 0, 14 },
			{ "event-limit", required_argument, 0, 15 },

			{ "voltage-queue-limit", required_argument, 0, 16},
			
			{ "evap-threshold", required_argument, 0, 17 },
			{ "evap-shell-size", required_argument, 0, 18 },
			{ "evap-cycle-size", required_argument, 0, 19 },
			{ "evap-queue-size", required_argument, 0, 20 },
			{ "evap-global-cycles", required_argument, 0, 21 },
			
			{ "refresh-interval", required_argument, 0, 22 },
			{ "refresh-threshold", required_argument, 0, 23 },
			{ "refresh-cycle-size", required_argument, 0, 24 },
			{ "refresh-queue-size", required_argument, 0, 25 },

			{ "help", no_argument, 0, 26 },

			{ 0, 0, 0, 0 }
		};

		EvaporationCommands::setDefaults(options,iniFile);
		
		// ***

		optind = 1; // reset option index 
		opterr = 1; // enable output of error messages by getopt()

		while (true)
		{	
			int key;
		
			do
				key = getopt_long(argc,argv,"",localOptions,0);
			while (key == 0 || key == '?');
		
			if (key == -1) break;

			switch (key)
			{
				case 1: // --threshold
					if (optarg == 0 || std::sscanf(optarg,"%f",&options->initThreshold) != 1)
						throw std::runtime_error("Error parsing 'threshold' argument!");
					break;
				case 2: // --cycle-size
					if (optarg == 0 || std::sscanf(optarg,"%u",&options->initCycleSize) != 1)
						throw std::runtime_error("Error parsing 'cycle-size' argument!");
					break;
				case 3: // --queue-size
					if (optarg == 0 || std::sscanf(optarg,"%u",&options->initQueueSize) != 1)
						throw std::runtime_error("Error parsing 'queue-size' argument!");
					break;
				case 4: // --no-trajectories, --write-trajectories
					if (optarg == 0)
						options->evapParams.trajectoryFile.clear();
					else
						options->evapParams.trajectoryFile = optarg;
					break;
				case 5: // --no-grids, --write-grids
					if (optarg == 0)
						options->evapParams.gridFile.clear();
					else
						options->evapParams.gridFile = optarg;
					break;
				case 6: // --no-surfaces, --write-surfaces
					if (optarg == 0)
						options->evapParams.surfaceFile.clear();
					else
						options->evapParams.surfaceFile = optarg;
					break;
				case 7: // --no-chunks, --chunk-size
					if (optarg == 0)
					{
						options->evapParams.resultsChunkSize = 0;
						options->evapParams.trajectoryChunkSize = 0;
					}
					else
					{
						unsigned int value;
						if (std::sscanf(optarg,"%u",&value) != 1) 
							throw std::runtime_error("Error parsing 'chunk-size' argument!");
						
						options->evapParams.resultsChunkSize = value;
						options->evapParams.trajectoryChunkSize = value;
					}
					break;
				case 8: // --snapshot-interval
					{
						unsigned int value;
						if (optarg == 0 || std::sscanf(optarg,"%u",&value) != 1)
							throw std::runtime_error("Error parsing 'snapshot-interval' argument!");
						
						options->evapParams.gridInterval = value;
						options->evapParams.surfaceFile = value;
					}
					break;
				case 9: // --input--dump
					if (optarg == 0) throw std::runtime_error("Error parsing 'input-dump' argument!");
					options->iDumpFile = optarg;
					break;
				case 10: // --no-dump, --write-dump
					if (optarg == 0)
						options->evapParams.dumpFile.clear();
					else
						options->evapParams.dumpFile = optarg;
					break;
				case 11: // --dump-interval
					if (optarg == 0 || std::sscanf(optarg,"%u",&options->evapParams.dumpInterval) != 1)
						throw std::runtime_error("Error parsing 'dump-interval' argument!");
					break;
				case 12: // --vacuum-identifier
					if (optarg == 0) throw std::runtime_error("Error parsing 'vacuum-identifier' argument!");
					options->evapParams.vacuumName = optarg;
					break;
				case 13: // --trajectory-error-threshold
					if (optarg == 0 || std::sscanf(optarg,"%f",&options->evapParams.trajectoryErrorThreshold) != 1)
						throw std::runtime_error("Error parsing 'trajectory-error-threshold' argument!");
					break;
				case 14: // --shrink-limit
					if (optarg == 0 || std::sscanf(optarg,"%f",&options->evapParams.shrinkLimit) != 1)
						throw std::runtime_error("Error parsing 'shrink-limit' argument!");
					break;
				case 15: // --event-limit
					if (optarg == 0 || std::sscanf(optarg,"%u",&options->evapParams.eventCntLimit) != 1)
						throw std::runtime_error("Error parsing 'event-limit' argument!");
					break;
				case 16: // --voltage-queue-size
					if (optarg == 0 || std::sscanf(optarg,"%u",&options->evapParams.voltageQueueSize) != 1)
						throw std::runtime_error("Error parsing 'voltage-queue-limit' argument!");
					break;
				case 17: // --evap-threshold
					if (optarg == 0 || std::sscanf(optarg,"%f",&options->evapParams.relaxThreshold) != 1)
						throw std::runtime_error("Error parsing 'evap-threshold' argument!");
					break;
				case 18: // --evap-shell-size
					if (optarg == 0 || std::sscanf(optarg,"%u",&options->evapParams.relaxShellSize) != 1)
						throw std::runtime_error("Error parsing 'evap-shell-size' argument!");
					break;
				case 19: // --evap-cycle-size
					if (optarg == 0 || std::sscanf(optarg,"%u",&options->evapParams.relaxCycleSize) != 1)
						throw std::runtime_error("Error parsing 'evap-cycle-size' argument!");
					break;
				case 20: // --evap-queue-size
					if (optarg == 0 || std::sscanf(optarg,"%u",&options->evapParams.relaxQueueSize) != 1)
						throw std::runtime_error("Error parsing 'evap-queue-size' argument!");
				case 21: // --evap-global-cycles
					if (optarg == 0 || std::sscanf(optarg,"%u",&options->evapParams.relaxGlobalCycles) != 1)
						throw std::runtime_error("Error parsing 'evap-global-cycles' argument!");
					break;
				case 22: // --refresh-interval
					if (optarg == 0 || std::sscanf(optarg,"%u",&options->evapParams.refreshInterval) != 1)
						throw std::runtime_error("Error parsing 'refresh-interval' argument!");
					break;
				case 23: // --refresh-threshold
					if (optarg == 0 || std::sscanf(optarg,"%f",&options->evapParams.refreshThreshold) != 1)
						throw std::runtime_error("Error parsing 'refresh-threshold' argument!");
					break;
				case 24: // --refresh-cycle-size
					if (optarg == 0 || std::sscanf(optarg,"%u",&options->evapParams.refreshCycleSize) != 1)
						throw std::runtime_error("Error parsing 'refresh-cycle-size' argument!");
					break;
				case 25: // --refresh-queue-size
					if (optarg == 0 || std::sscanf(optarg,"%u",&options->evapParams.refreshQueueSize) != 1)
						throw std::runtime_error("Error parsing 'refresh-queue-size' argument!");
					break;
				case 26: // --help
					std::cout << "Evaporation mode syntax:" << std::endl;
					std::cout << std::endl;
					std::cout << "./tapsim evaporation --input-dump=? iconfig [results]" << std::endl;
					std::cout << "./tapsim evaporation iconfig inode [results]" << std::endl;
					std::cout << std::endl;
					std::cout << "List of additional parameters in evaporation mode:" << std::endl;
					std::cout << "\t--no-file-output, --write-binary, --write-ascii" << std::endl;
					std::cout << "\t--threads=?" << std::endl;
					std::cout << "\t--ini-file=?" << std::endl;
					std::cout << "\t--skip-relaxation" << std::endl;
					std::cout << "\t--init-potentials, --reset-potentials" << std::endl;
					std::cout << "\t--threshold=?," << std::endl;
					std::cout << "\t--cycle-size=?," << std::endl;
					std::cout << "\t--queue-size=?," << std::endl;
					std::cout << "\t--no-trajectories, --write-trajectories=?," << std::endl;
					std::cout << "\t--no-grids, --write-grids=?," << std::endl;
					std::cout << "\t--no-surfaces, --write-surfaces=?" << std::endl;
					std::cout << "\t--no-chunks, --chunk-size=?" << std::endl;
					std::cout << "\t--snapshot-interval=?" << std::endl;
					std::cout << "\t--no-dump, --write-dump=?" << std::endl;
					std::cout << "\t--dump-interval=?" << std::endl;
					std::cout << "\t--vacuum-identifier=?" << std::endl;
					std::cout << "\t--trajectory-error-threshold=?" << std::endl;
					std::cout << "\t--shrink-limit=?" << std::endl;
					std::cout << "\t--event-limit=?" << std::endl;
					std::cout << "\t--voltage-queue-size=?" << std::endl;
					std::cout << "\t--evap-threshold=?" << std::endl;
					std::cout << "\t--evap-cycle-size=?" << std::endl;
					std::cout << "\t--evap-shell-size=?" << std::endl;
					std::cout << "\t--evap-queue-size=?" << std::endl;
					std::cout << "\t--evap-global-cycles=?" << std::endl;
					std::cout << "\t--refresh-interval=?" << std::endl;
					std::cout << "\t--refresh-threshold=?" << std::endl;
					std::cout << "\t--refresh-cycle-size=?" << std::endl;
					std::cout << "\t--refresh-queue-size=?" << std::endl;
					std::cout << std::endl;
					throw int(0);
				default:
					std::string errorString("Unknown or invalid argument in command line: \"");
					errorString += optarg;
					errorString += "\"!";
					
					throw std::runtime_error(errorString);
			}
		}

		// *** process additional arguments
		
		if (std::strcmp(argv[optind],"evaporation") == 0)
		{
			switch (argc - 1 - optind)
			{
				case 1: // **** ./tapsim [...] evaporation --input-dump=? iconfig
					options->iConfigFile = argv[optind+1];
					break;
				case 2:
					if (options->iDumpFile.empty())
					{
						// **** ./tapsim [...] evaporation iconfig inode
						options->iConfigFile = argv[optind+1];
						options->iNodeFile = argv[optind+2];
					}
					else
					{
						// **** ./tapsim [...] evaporation --input-dump=? iconfig results
						options->iConfigFile = argv[optind+1];
						options->evapParams.resultsFile = argv[optind+2];
					}
					break;
				case 3:
					// **** ./tapsim [...] evaporation iconfig inode results
					options->iConfigFile = argv[optind+1];
					options->iNodeFile = argv[optind+2];
					options->evapParams.resultsFile = argv[optind+3];
					break;
				default:
					throw std::runtime_error("Missing arguments for processing in 'evaporation'-mode!");
			}
		}
		else
			throw std::runtime_error("parseEvaporationCommands()");
	}

	void parseResumptionCommands(int argc, char** argv, ResumptionCommands* options, const std::string& iniFile)
	{
		/* -----------------------------------------------------------------------------
		* 
		* Calling convention:
		* 
		*	./tapsim [...] resumption idump [--last-event=?] [results]
		* 
		* ------------------------------------------------------------------------------
		*/

		struct option localOptions[] =
		{
			// global options: recognized but skipped
			{ "no-file-output", no_argument, 0, 0 },
			{ "write-ascii", no_argument, 0, 0 },
			{ "write-binary", no_argument,0, 0 },
			{ "threads", required_argument, 0, 0 },
			{ "ini-file", required_argument, 0, 0 },

			// local options:
			{ "no-trajectories", no_argument, 0, 1 },
			{ "write-trajectories", required_argument, 0, 1 },

			{ "no-grids", no_argument, 0, 2 },
			{ "write-grids", required_argument, 0, 2 },

			{ "no-surfaces", no_argument, 0, 3 },
			{ "write-surfaces", required_argument, 0 , 3 },

			{ "no-chunks", no_argument, 0, 4 },
			{ "chunk-size", required_argument, 0, 4 },

			{ "snapshot-interval", required_argument, 0, 5 },

			{ "no-dump", no_argument, 0, 6 },
			{ "write-dump", required_argument, 0, 6 },

			{ "dump-interval", required_argument, 0, 7 },

			{ "vacuum-identifier", required_argument, 0, 8 },

			{ "trajectory-error-threshold", required_argument, 0, 9 },

			{ "shrink-limit", required_argument, 0, 10 },
			{ "last-event", required_argument, 0, 11 },
			{ "event-limit", required_argument, 0, 12 },
			
			{ "voltage-queue-limit", required_argument, 0, 13},
			
			{ "evap-threshold", required_argument, 0, 14 },
			{ "evap-shell-size", required_argument, 0, 15 },
			{ "evap-cycle-size", required_argument, 0, 16 },
			{ "evap-queue-size", required_argument, 0, 17 },
			{ "evap-global-cycles", required_argument, 0, 18 },

			{ "refresh-interval", required_argument, 0, 19 },
			{ "refresh-threshold", required_argument, 0, 20 },
			{ "refresh-cycle-size", required_argument, 0, 21 },
			{ "refresh-queue-size", required_argument, 0, 22 },

			{ "help", no_argument, 0, 23 },

			{ 0, 0, 0, 0 }
		};

		// ***
		
		ResumptionCommands::setDefaults(options,iniFile);

		// ***

		optind = 1; // reset option index
		opterr = 1; // enable output of error messages by getopt()

		while (true)
		{
			int key;

			do
				key = getopt_long(argc,argv,"",localOptions,0);
			while (key == 0 || key == '?');

			if (key == -1) break;

			switch (key)
			{
				case 1: // --no-trajectories, --write-trajectories
					if (optarg == 0)
						options->evapParams.trajectoryFile.clear();
					else
						options->evapParams.trajectoryFile = optarg;
					break;
				case 2: // -no-grids, --write-grids
					if (optarg == 0)
						options->evapParams.gridFile.clear();
					else
						options->evapParams.gridFile = optarg;
					break;
				case 3: // -no-surfaces, --write-surfaces
					if (optarg == 0)
						options->evapParams.surfaceFile.clear();
					else
						options->evapParams.surfaceFile = optarg;
					break;
				case 4: // --no-chunks, --chunk-size
					if (optarg == 0)
					{
						options->evapParams.resultsChunkSize = 0;
						options->evapParams.trajectoryChunkSize = 0;
					}
					else
					{
						unsigned int value;
						if (std::sscanf(optarg,"%u",&value) != 1) 
							throw std::runtime_error("Error parsing 'chunk-size' argument!");
						
						options->evapParams.resultsChunkSize = value;
						options->evapParams.trajectoryChunkSize = value;
					}
					break;
				case 5: // --snapshot-interval
					{
						unsigned int value;
						if (optarg == 0 || std::sscanf(optarg,"%u",&value) != 1)
							throw std::runtime_error("Error parsing 'snapshot-interval' argument!");
						
						options->evapParams.gridInterval = value;
						options->evapParams.surfaceFile = value;
					}
					break;
				case 6: // --no-dump, --write-dump
					if (optarg == 0)
						options->evapParams.dumpFile.clear();
					else
						options->evapParams.dumpFile = optarg;
					break;
				case 7: // --dump-interval
					if (optarg == 0 || std::sscanf(optarg,"%u",&options->evapParams.dumpInterval) != 1)
						throw std::runtime_error("Error parsing 'dump-interval' argument!");
					break;
				case 8: // --vacuum-identifier
					if (optarg == 0) throw std::runtime_error("Error parsing 'vacuum-identifier' argument!");
					options->evapParams.vacuumName = optarg;
					break;
				case 9: // --trajectory-error-threshold
					if (optarg == 0 || std::sscanf(optarg,"%f",&options->evapParams.trajectoryErrorThreshold) != 1)
						throw std::runtime_error("Error parsing 'trajectory-error-threshold' argument!");
					break;
				case 10: // --shrink-limit
					if (optarg == 0 || std::sscanf(optarg,"%f",&options->evapParams.shrinkLimit) != 1)
						throw std::runtime_error("Error parsing 'shrink-limit' argument!");
					break;
				case 11: // --last-event
					if (optarg == 0 || std::sscanf(optarg,"%u",&options->evapParams.initEventCnt) != 1)
						throw std::runtime_error("Error parsing 'last-event' argument!");
					break;
				case 12:// --event-limit
					if (optarg == 0 || std::sscanf(optarg,"%u",&options->evapParams.eventCntLimit) != 1)
						throw std::runtime_error("Error parsing 'event-limit' argument!");
					break;
				case 13: // --voltage-queue-size
					if (optarg == 0 || std::sscanf(optarg,"%u",&options->evapParams.voltageQueueSize) != 1)
						throw std::runtime_error("Error parsing 'voltage-queue-limit' argument!");
					break;
				case 14: // --evap-threshold
					if (optarg == 0 || std::sscanf(optarg,"%f",&options->evapParams.relaxThreshold) != 1)
						throw std::runtime_error("Error parsing 'evap-threshold' argument!");
					break;
				case 15: // --evap-shell-size
					if (optarg == 0 || std::sscanf(optarg,"%u",&options->evapParams.relaxShellSize) != 1)
						throw std::runtime_error("Error parsing 'evap-shell-size' argument!");
					break;
				case 16: // --evap-cycle-size
					if (optarg == 0 || std::sscanf(optarg,"%u",&options->evapParams.relaxCycleSize) != 1)
						throw std::runtime_error("Error parsing 'evap-cycle-size' argument!");
					break;
				case 17: // --evap-queue-size
					if (optarg == 0 || std::sscanf(optarg,"%u",&options->evapParams.relaxQueueSize) != 1)
						throw std::runtime_error("Error parsing 'evap-queue-size' argument!");
				case 18: // --evap-global-cycles
					if (optarg == 0 || std::sscanf(optarg,"%u",&options->evapParams.relaxGlobalCycles) != 1)
						throw std::runtime_error("Error parsing 'evap-global-cycles' argument!");
					break;
				case 19: // --refresh-interval
					if (optarg == 0 || std::sscanf(optarg,"%u",&options->evapParams.refreshInterval) != 1)
						throw std::runtime_error("Error parsing 'refresh-interval' argument!");
					break;
				case 20: // --refresh-threshold
					if (optarg == 0 || std::sscanf(optarg,"%f",&options->evapParams.refreshThreshold) != 1)
						throw std::runtime_error("Error parsing 'refresh-threshold' argument!");
					break;
				case 21: // --refresh-cycle-size
					if (optarg == 0 || std::sscanf(optarg,"%u",&options->evapParams.refreshCycleSize) != 1)
						throw std::runtime_error("Error parsing 'refresh-cycle-size' argument!");
					break;
				case 22: // --refresh-queue-size
					if (optarg == 0 || std::sscanf(optarg,"%u",&options->evapParams.refreshQueueSize) != 1)
						throw std::runtime_error("Error parsing 'refresh-queue-size' argument!");
					break;
				case 23: // --help
					std::cout << "Resumption mode syntax:" << std::endl;
					std::cout << std::endl;
					std::cout << "./tapsim resumption idump [--last-event=?] [results]" << std::endl;
					std::cout << std::endl;
					std::cout << "List of additional parameters in resumption mode:" << std::endl;
					std::cout << "\t--no-file-output, --write-binary, --write-ascii" << std::endl;
					std::cout << "\t--threads=?" << std::endl;
					std::cout << "\t--ini-file=?" << std::endl;
					std::cout << "\t--no-trajectories, --write-trajectories=?," << std::endl;
					std::cout << "\t--no-grids, --write-grids=?," << std::endl;
					std::cout << "\t--no-surfaces, --write-surfaces=?" << std::endl;
					std::cout << "\t--no-chunks, --chunk-size=?" << std::endl;
					std::cout << "\t--snapshot-interval=?" << std::endl;
					std::cout << "\t--no-dump, --write-dump=?" << std::endl;
					std::cout << "\t--dump-interval=?" << std::endl;
					std::cout << "\t--vacuum-identifier=?" << std::endl;
					std::cout << "\t--trajectory-error-threshold=?" << std::endl;
					std::cout << "\t--shrink-limit=?" << std::endl;
					std::cout << "\t--event-limit=?" << std::endl;
					std::cout << "\t--voltage-queue-size=?" << std::endl;
					std::cout << "\t--evap-threshold=?" << std::endl;
					std::cout << "\t--evap-cycle-size=?" << std::endl;
					std::cout << "\t--evap-shell-size=?" << std::endl;
					std::cout << "\t--evap-queue-size=?" << std::endl;
					std::cout << "\t--evap-global-cycles=?" << std::endl;
					std::cout << "\t--refresh-interval=?" << std::endl;
					std::cout << "\t--refresh-threshold=?" << std::endl;
					std::cout << "\t--refresh-cycle-size=?" << std::endl;
					std::cout << "\t--refresh-queue-size=?" << std::endl;
					std::cout << std::endl;
					throw int(0);
				default:
					std::string errorString("Unknown or invalid argument in command line: \"");
					errorString += optarg;
					errorString += "\"!";

					throw std::runtime_error(errorString);
			}
		}

		// *** process additional arguments

		if (std::strcmp(argv[optind],"resumption") == 0)
		{
			switch (argc - 1 - optind)
			{
				case 1:
					// **** ./tapsim [...] resumption idump [--last-event=?]
					options->iDumpFile = argv[optind+1];
					break;
				case 2:
					// **** ./tapsim [...] resumption idump [--last-event=?] results
					options->iDumpFile = argv[optind+1];
					options->evapParams.resultsFile = argv[optind+2];
					break;
				default:
					throw std::runtime_error("Missing arguments for processing in 'resumption'-mode!");
			}
		}
		else
			throw std::runtime_error("parseResumptionCommands()");
	}
}

MakeIniCommands::MakeIniCommands()
	: iniFilename()
{}

void MakeIniCommands::setDefaults(MakeIniCommands* obj)
{
	obj->iniFilename = defaultIniFilename;
}


// ***

RelaxationCommands::RelaxationCommands()
	: iConfigFile(),
	  oConfigFile(),
	  oConfigMode(0),
	  iNodeFile(),
	  oNodeFile(),
	  oNodeMode(0),
	  iDumpFile(),
	  oDumpFile(),
	  threshold(0.0f),
	  cycleSize(0),
	  queueSize(0),
	  resetPotentials(0)  
{}

void RelaxationCommands::setDefaults(RelaxationCommands* obj)
{
	obj->iConfigFile = std::string();

	obj->oConfigFile = std::string();
	obj->oConfigMode = File_Io::ASCII;
	
	obj->iNodeFile = std::string();

	obj->oNodeFile = std::string();
	obj->oNodeFile = File_Io::BINARY;
	

	obj->iDumpFile = std::string();
	obj->oDumpFile = std::string();

	obj->threshold = 0.0f;
	obj->cycleSize = 100;
	obj->queueSize = 25;

	obj->resetPotentials = 0;
}

// ***

EvaporationCommands::EvaporationCommands()
	: iDumpFile(),
	  iConfigFile(),
	  iNodeFile(),
	  evapParams()
{}

void EvaporationCommands::setDefaults(EvaporationCommands* obj, const std::string& iniFile)
{
	obj->iDumpFile = std::string();
	
	obj->iConfigFile = std::string();
	obj->iNodeFile = std::string();

	// ***

	obj->initRelax = 1;

	//obj->initThreshold = 0.0f;
	//obj->initCycleSize = 100;
	//obj->initQueueSize = 25;

	obj->resetPotentials = 0;
	
	// ***

	Process::EvaporationOptions::setDefaults(iniFile.c_str(),&obj->evapParams);
	
	// *** set the refresh relaxation parameters as the default for the initial relaxation
	
	obj->initThreshold = obj->evapParams.refreshThreshold;
	obj->initCycleSize = obj->evapParams.refreshCycleSize;
	obj->initQueueSize = obj->evapParams.refreshQueueSize;
}

// ***

ResumptionCommands::ResumptionCommands()
	: iDumpFile(),
	  evapParams()
{}

void ResumptionCommands::setDefaults(ResumptionCommands* obj, const std::string& iniFile)
{
	obj->iDumpFile = std::string();

	// ***

	Process::EvaporationOptions::setDefaults(iniFile.c_str(),&obj->evapParams);
	
	obj->evapParams.initEventCnt = std::numeric_limits<unsigned int>::max();
}

// ***

CommandLineOptions::CommandLineOptions()
	: outputHeader(),
	  threadNum(0),
	  mode(0),
	  iniFilename(),
	  relaxationParams(),
	  evaporationParams()
{}

// ***

void parseCommandLine(int argc, char** argv, CommandLineOptions* options)
{
	if (argc == 1)
	{
		std::cout << "For details on the program type './tapsim --about'." << std::endl;
		std::cout << std::endl;
		std::cout << "No arguments! Please try './tapsim --help'!" << std::endl;
		std::cout << std::endl;
		
		throw int(0);
	}
	
	// ***

	int helpFlag(0);
	int aboutFlag(0);
	int versionFlag(0);

	int enableFileOutput(1);
	int asciiFileOutputMode(0);
	int binaryFileOutputMode(0);
	
	struct option globalOptions[] =
	{
		{ "help", no_argument, &helpFlag, 1 },
		{ "about", no_argument, &aboutFlag, 1 },
		{ "version", no_argument, &versionFlag, 1 },

		{ "no-file-output", no_argument, &enableFileOutput, 0 },
		{ "write-ascii", no_argument, &asciiFileOutputMode, 1 },
		{ "write-binary", no_argument, &binaryFileOutputMode, 1 },

		{ "threads", required_argument, 0, 1 },
		{ "ini-file", required_argument, 0, 2 },

		
		{ 0, 0, 0, 0 }
	};

	opterr = 0; // disable printing of error messages
	
	int key(0);
	int index(0);

	while (key >= 0)
	{	
		do
			key = getopt_long(argc,argv,"",globalOptions,&index);
		while (key == 0 || key == '?');
		
		if (key == -1) break;

		switch (key)
		{
			case 1: // --threadNum
				if (std::sscanf(optarg,"%d",&options->threadNum) != 1)
					throw std::runtime_error("Error parsing 'threadNum' argument!");
				break;
			case 2: // --ini-file
				if (optarg == 0) throw std::runtime_error("Error parsing 'ini-file' argument!");
				options->iniFilename = optarg;
				break;
			default:
				std::string errorString("Unknown or invalid argument in command line: \"");
				errorString += optarg;
				errorString += "\"!";
				
				throw std::runtime_error(errorString);
		}
	}

	// *** set default ini file if it exists
	
	if (options->iniFilename.empty())
	{
		std::ifstream file(defaultIniFilename,std::ifstream::in);
		if (file.good()) options->iniFilename = defaultIniFilename;
	}
	
	// ***

	options->outputHeader.append("# Generated by TAPSim ");
	options->outputHeader.append(tapSimVersion);
	options->outputHeader.append(" (");
	options->outputHeader.append(tapSimRevision);
	options->outputHeader.append(")\n");
	options->outputHeader.append("#\n");

	options->outputHeader.append("# Compilation time: ");
	options->outputHeader.append(uniqueDate);
	options->outputHeader.append(", ");
	options->outputHeader.append(uniqueTime);
	options->outputHeader.append("\n");
	
	options->outputHeader.append("# Initialization data: ");

	if (options->iniFilename.empty())
		options->outputHeader.append("internal defaults");
	else
	{
		options->outputHeader.append("'");
		options->outputHeader.append(options->iniFilename);
		options->outputHeader.append("'");
	}
	options->outputHeader.append("\n");

	options->outputHeader.append("#\n");
	options->outputHeader.append("# Invoked command line parameters:\n#\n");;

	for (int i = 1; i < argc; i++)
	{
		options->outputHeader.append("#\t");
		options->outputHeader.append(argv[i]);
		options->outputHeader.append("\n");
	}

	options->outputHeader.append("#\n");

	// ***

	if (versionFlag)
	{
		std::cout << "Distribution: \"" << tapSim_distVersion << "\" (" << tapSim_distDate << ", " << tapSim_distTime << ")" << std::endl;
		std::cout << std::endl;
		std::cout << "Compilation time-stamp: " << uniqueDate << ", " << uniqueTime << std::endl;
		std::cout << std::endl;
		throw int(0);
	}
	else if (aboutFlag)
	{
		std::cout << "TAPSim is an atom probe (3DAP) data simulation program written by Christian" << std::endl;
		std::cout << "Oberdorfer.  Its source code is licensed under the terms of the GNU General" << std::endl;
		std::cout << "Public License Version 3 or later ( GNU GPL v3+ )  as published by the Free" << std::endl;
		std::cout << "Software Foundation. This  program is distributed  in the hope that it will" << std::endl;
		std::cout << "be useful,  but WITHOUT ANY WARRANTY;  without even the implied warranty of" << std::endl;
		std::cout << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GPL for  more" << std::endl;
		std::cout << "details." << std::endl;
		
		std::cout << "You should have receeived a copy of the GNU GPL along with this program. If" << std::endl;
		std::cout << "not, look at 'http://www.gnu.org/licenses/gpl.txt'." << std::endl;
		std::cout << std::endl;
		std::cout << "More information about the program and  its application is available in the" << std::endl;
		std::cout << "internet, see 'http://www.uni-muenster.de/physik.mp/schmitz/tapsim'."   << std::endl;
		std::cout << std::endl;

		std::cout << "In the case of problems or bugs, please contact:" << std::endl;
		std::cout << std::endl;
		std::cout << "Dipl. Phys. Christian Oberdorfer" << std::endl;
		std::cout << std::endl;
		std::cout << "E-Mail: oberdorc@uni-muenster.de" << std::endl;
		std::cout << std::endl;
		std::cout << "Affiliation:" << std::endl;
		std::cout << "Westfälische Wilhelms-Universität Münster" << std::endl;
		std::cout << "Institut für Materialphysik" << std::endl;
		std::cout << "Wilhelm-Klemm-Str. 10" << std::endl;
		std::cout << "48149 Münster, Germany" << std::endl;
		std::cout << std::endl;
		throw int(0);
	}
	else if (!helpFlag)
	{
		std::cout << "This is free software, and you are welcome to redistribute" << std::endl;
		std::cout << "it under certain conditions." << std::endl;
		std::cout << std::endl;
		std::cout << "For further details on the program type './tapsim --about'." << std::endl;
		std::cout << std::endl;
	}
	
	// ***

	if (optind < argc)
	{
		if (std::strcmp(argv[optind],"make-ini") == 0)
		{
			options->mode = CommandLineOptions::MAKE_INI;
			parseMakeIniCommands(argc,argv,&options->makeIniParams);
		}
		else if (std::strcmp(argv[optind],"relaxation") == 0)
		{
			options->mode = CommandLineOptions::RELAXATION;
			parseRelaxationCommands(argc,argv,&options->relaxationParams);
		}
		else if (std::strcmp(argv[optind],"evaporation") == 0)
		{
			options->mode = CommandLineOptions::EVAPORATION;
			parseEvaporationCommands(argc,argv,&options->evaporationParams,options->iniFilename);

			// *** set global options which are not parsed locally
			// << empty >>
		}
		else if (std::strcmp(argv[optind],"resumption") == 0)
		{
			options->mode = CommandLineOptions::RESUMPTION;
			parseResumptionCommands(argc,argv,&options->resumptionParams,options->iniFilename);

			// *** set global options which are not parsed locally
			// << empty >>
		}
		else
			throw std::runtime_error("parseCommandLine(): Invalid operation mode!");
	}
	else if (helpFlag)
	{
		std::cout << "Call: ./tapsim mode [file(s)] [parameter(s)]" << std::endl;
		std::cout << std::endl;
		std::cout << "Valid modes: " << std::endl;
		std::cout << "\t1) make-ini" << std::endl;
		std::cout << "\t2) relaxation" << std::endl;
		std::cout << "\t3) evaporation" << std::endl;
		std::cout << "\t4) resumption" << std::endl;
		std::cout << std::endl;
		std::cout << "Neccessary files depend on chosen mode!" << std::endl;
		std::cout << "Possible parameters depend on mode!" << std::endl;
		std::cout << std::endl;
		std::cout << "Try './tapsim mode --help' for detailed information!" << std::endl;
		std::cout << std::endl;
		throw int(0);
	}
	else
		throw std::runtime_error("parseCommandLine()");

	// *** override ini default settings with options from 
	// *** the command line - global parameters only

	if (!enableFileOutput)
	{
		options->relaxationParams.oConfigFile.clear();
		options->relaxationParams.oNodeFile.clear();
		options->relaxationParams.oDumpFile.clear();

		options->evaporationParams.evapParams.resultsFile.clear();
		options->evaporationParams.evapParams.trajectoryFile.clear();
		options->evaporationParams.evapParams.gridFile.clear();
		options->evaporationParams.evapParams.surfaceFile.clear();
		options->evaporationParams.evapParams.geometryFile.clear();
		options->evaporationParams.evapParams.dumpFile.clear();

		options->resumptionParams.evapParams.resultsFile.clear();
		options->resumptionParams.evapParams.trajectoryFile.clear();
		options->resumptionParams.evapParams.gridFile.clear();
		options->resumptionParams.evapParams.surfaceFile.clear();
		options->resumptionParams.evapParams.geometryFile.clear();
		options->resumptionParams.evapParams.dumpFile.clear();
	}
	
	if (asciiFileOutputMode)
	{
		options->relaxationParams.oConfigMode = File_Io::ASCII;
		options->relaxationParams.oNodeMode = File_Io::ASCII;
		
		options->evaporationParams.evapParams.resultsMode = File_Io::ASCII;
		options->evaporationParams.evapParams.trajectoryMode = File_Io::ASCII;
		options->evaporationParams.evapParams.surfaceMode = File_Io::ASCII;
		options->evaporationParams.evapParams.gridMode = File_Io::ASCII;
		options->evaporationParams.evapParams.geometryMode = File_Io::ASCII;

		options->resumptionParams.evapParams.resultsMode = File_Io::ASCII;
		options->resumptionParams.evapParams.trajectoryMode = File_Io::ASCII;
		options->resumptionParams.evapParams.surfaceMode = File_Io::ASCII;
		options->resumptionParams.evapParams.gridMode = File_Io::ASCII;
		options->resumptionParams.evapParams.geometryMode = File_Io::ASCII;
	}
	
	if (binaryFileOutputMode)
	{
		options->relaxationParams.oConfigMode = File_Io::BINARY;
		options->relaxationParams.oNodeMode = File_Io::BINARY;
		
		options->evaporationParams.evapParams.resultsMode = File_Io::BINARY;
		options->evaporationParams.evapParams.trajectoryMode = File_Io::BINARY;
		options->evaporationParams.evapParams.surfaceMode = File_Io::BINARY;
		options->evaporationParams.evapParams.gridMode = File_Io::BINARY;
		options->evaporationParams.evapParams.geometryMode = File_Io::BINARY;

		options->resumptionParams.evapParams.resultsMode = File_Io::BINARY;
		options->resumptionParams.evapParams.trajectoryMode = File_Io::BINARY;
		options->resumptionParams.evapParams.surfaceMode = File_Io::BINARY;
		options->resumptionParams.evapParams.gridMode = File_Io::BINARY;
		options->resumptionParams.evapParams.geometryMode = File_Io::BINARY;
	}
}

void makeIniList(std::list<File_Io::KeyValue>* obj)
{
	if (obj == 0) throw std::runtime_error("makeIniList()");

	obj->clear();

	obj->push_back(File_Io::KeyValue("RESULTS_FILENAME","results_data"));
	obj->push_back(File_Io::KeyValue("RESULTS_BINARY_OUTPUT","0"));
	obj->push_back(File_Io::KeyValue("RESULTS_CHUNK_SIZE","10000"));
	obj->push_back(File_Io::KeyValue());

	obj->push_back(File_Io::KeyValue("TRAJECTORY_FILENAME","trajectory_data"));
	obj->push_back(File_Io::KeyValue("TRAJECTORY_BINARY_OUTPUT","1"));
	obj->push_back(File_Io::KeyValue("TRAJECTORY_CHUNK_SIZE","10000"));
	obj->push_back(File_Io::KeyValue());

	obj->push_back(File_Io::KeyValue("GRID_FILENAME","grid_data"));
	obj->push_back(File_Io::KeyValue("GRID_BINARY_OUTPUT","1"));
	{
		std::ostringstream gridInterval;
		gridInterval << std::numeric_limits<int>::max();
		obj->push_back(File_Io::KeyValue("GRID_INTERVAL",gridInterval.str()));
	}
	obj->push_back(File_Io::KeyValue());

	obj->push_back(File_Io::KeyValue("SURFACE_FILENAME","surface_data"));
	obj->push_back(File_Io::KeyValue("SURFACE_BINARY_OUTPUT","1"));
	{
		std::ostringstream surfaceInterval;
		surfaceInterval << std::numeric_limits<int>::max();
		obj->push_back(File_Io::KeyValue("SURFACE_INTERVAL",surfaceInterval.str()));
	}
	obj->push_back(File_Io::KeyValue());

	obj->push_back(File_Io::KeyValue("GEOMETRY_FILENAME","geometry.dat"));
	obj->push_back(File_Io::KeyValue("GEOMETRY_BINARY_OUTPUT","1"));
	obj->push_back(File_Io::KeyValue("GEOMETRY_CONTENTS","0xff"));
	obj->push_back(File_Io::KeyValue());

	obj->push_back(File_Io::KeyValue("DUMP_FILENAME","dump"));
	obj->push_back(File_Io::KeyValue("DUMP_INTERVAL","10000"));
	obj->push_back(File_Io::KeyValue());

	obj->push_back(File_Io::KeyValue("OUTPUT_DELAY_TIME","60"));
	obj->push_back(File_Io::KeyValue());

	obj->push_back(File_Io::KeyValue("METHOD_EVAPORATION_PROBABILITY","0"));
	obj->push_back(File_Io::KeyValue("METHOD_EVAPORATION_CHOICE","0"));
	obj->push_back(File_Io::KeyValue());

	obj->push_back(File_Io::KeyValue("VACUUM_CELL_IDENTIFIER","Vacuum"));
	obj->push_back(File_Io::KeyValue());

	obj->push_back(File_Io::KeyValue("TRAJECTORY_INTEGRATOR_TYPE","1"));
	obj->push_back(File_Io::KeyValue("TRAJECTORY_STEPPER_TYPE","1"));
	obj->push_back(File_Io::KeyValue("TRAJECTORY_INITIAL_TIME_STEP","1e-14"));
	obj->push_back(File_Io::KeyValue("TRAJECTORY_MINIMUM_TIME_STEP","1e-16"));
	obj->push_back(File_Io::KeyValue("TRAJECTORY_ERROR_THRESHOLD","1e-6"));
	obj->push_back(File_Io::KeyValue("TRAJECTORY_NON_RANDOM_START_POSITION","1"));
	obj->push_back(File_Io::KeyValue("TRAJECTORY_NO_INITIAL_VELOCITY","1"));
	obj->push_back(File_Io::KeyValue());

	obj->push_back(File_Io::KeyValue("VOLTAGE_QUEUE_SIZE","250"));
	obj->push_back(File_Io::KeyValue());

	obj->push_back(File_Io::KeyValue("LOCAL_RELAXATION_ERROR_THRESHOLD","1e-2"));
	obj->push_back(File_Io::KeyValue("LOCAL_RELAXATION_SHELL_SIZE","3"));
	obj->push_back(File_Io::KeyValue("LOCAL_RELAXATION_CYCLE_SIZE","5"));
	obj->push_back(File_Io::KeyValue("LOCAL_RELAXATION_QUEUE_SIZE","5"));
	obj->push_back(File_Io::KeyValue("GLOBAL_RELAXATION_STEPS","1"));
	obj->push_back(File_Io::KeyValue());

	obj->push_back(File_Io::KeyValue("REFRESH_RELAXATION_INTERVAL","0"));
	obj->push_back(File_Io::KeyValue("REFRESH_RELAXATION_ERROR_THRESHOLD","0.0"));
	obj->push_back(File_Io::KeyValue("REFRESH_RELAXATION_CYCLE_SIZE","100"));
	obj->push_back(File_Io::KeyValue("REFRESH_RELAXATION_QUEUE_SIZE","25"));
	obj->push_back(File_Io::KeyValue());
}