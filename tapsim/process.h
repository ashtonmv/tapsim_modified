#ifndef TAPSIM_PROCESS_H
#define TAPSIM_PROCESS_H

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

#include <string>

#include "system_3d.h"

namespace Process
{
	void initialization(const char*, const char*, System_3d*);

	// ***

	void relaxation(const float, const unsigned int, const unsigned int, System_3d*);

	// ***

	struct EvaporationOptions
	{
		EvaporationOptions();

		static void setDefaults(const char* filename, Process::EvaporationOptions* obj);

		std::string resultsFile;
		int resultsMode;
		unsigned int resultsChunkSize;

		std::string trajectoryFile;
		int trajectoryMode;
		unsigned int trajectoryChunkSize;

		std::string gridFile;
		int gridMode;
		unsigned int gridInterval;

		std::string surfaceFile;
		int surfaceMode;
		unsigned int surfaceInterval;

		// ***

		std::string geometryFile;
		unsigned char geometryParams;
		int geometryMode;

		std::string dumpFile;
		unsigned int dumpInterval;

		// ***

		unsigned int delayTime;

		// ***

		int probMode;
		int evapMode;

		std::string vacuumName;

		// ***

		int trajectoryIntegrator;
		int trajectoryStepper;

		float trajectoryTimeStep;
		float trajectoryTimeStepLimit;
		float trajectoryErrorThreshold;

		int fixedInitialPosition;
		int noInitialVelocity;

		// ***

		float initShrinkage;
		float shrinkLimit;

		unsigned int initEventCnt;
		unsigned int eventCntLimit;

		// ***

		unsigned int voltageQueueSize;

		// ***

		float relaxThreshold;

		unsigned int relaxShellSize;
		unsigned int relaxCycleSize;
		unsigned int relaxQueueSize;

		unsigned int relaxGlobalCycles;

		// ***

		int refreshInterval;

		float refreshThreshold;
		unsigned int refreshCycleSize;
		unsigned int refreshQueueSize;
	};

	void evaporation(const EvaporationOptions&, const std::string&, System_3d*);
}

#endif
