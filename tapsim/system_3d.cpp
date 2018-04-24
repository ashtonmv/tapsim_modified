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

#include "system_3d.h"

void preset_sphericalPotential(System_3d* system, const Geometry_3d::Point& center, const float min, const float max)
{
	float phiSlope;
	{
		const float minDistance = (system->geomTable.min() - center).length();
		const float maxDistance = (system->geomTable.min() - center).length();

		phiSlope = max;
		phiSlope -= min;

		if (minDistance > maxDistance)
			phiSlope /= minDistance;
		else
			phiSlope /= maxDistance;
	}

	// ***

	for (int index = 0; index < system->gridTable.numNodes(); index++)
	{
		float phi;

		if (system->gridTable.node(index).dirichletBoundary())
			phi = system->configTable[system->gridTable.node(index).id()].phi();
		else
		{
			phi = max;
			phi -= phiSlope * (system->geomTable.nodeCoords(index) - center).length();
		}

		for (int i = 0; i < 2; i++)
			system->gridTable.node(index).setPhi(i,phi);
	}
}
