/************************************************************************
*                                                                       *
* MeshGen - a mesh generator utility for TapSim                         *
*                                                                       *
* Copyright (C) 2012 Christian Oberdorfer                               *
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
#include <vector>
#include <list>
#include <map>
#include <set>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits>

#include <cmath>
#include <cstring>

#include <getopt.h>

#include "../tetgen/dist_tetgen.h"

enum { VACUUM =0, TOP =1, BOTTOM =2, SIDES =3 };

const double pi = 3.141592653589793238462643383279502884197169399375105820974944;

std::string trim(const std::string& source)
{
	const char* whiteValues = " \t\n";

	const size_t start = source.find_first_not_of(whiteValues,0);
	if (start == std::string::npos) return std::string();

	const size_t end = source.find_last_not_of(whiteValues,std::string::npos);
	if (end == std::string::npos) return std::string();

	return source.substr(start,end-start+1);
}

void readInitialization(const char* filename, std::map<std::string,std::string>* values, const char* delimiter)
{
	if (values == 0) throw std::runtime_error("readInitialization()");

	values->clear();

	// ***

	std::ifstream file(filename,std::ifstream::in);
	if (!file.good()) throw std::runtime_error("readInitialization(): cannot open file!");


	while (file.good())
	{
		std::string buffer;

		while (file.good() && (buffer.empty() || buffer.at(0) == '#'))
			std::getline(file,buffer);

		if (!file.good()) break;

		std::pair<std::string,std::string> dataEntry;

		const size_t pos = buffer.find(delimiter);
		if (pos == std::string::npos) throw std::runtime_error("readInitialization()");

		dataEntry.first = trim(buffer.substr(0,pos));
		dataEntry.second = trim(buffer.substr(pos+std::strlen(delimiter)));

		const std::pair<std::map<std::string,std::string>::iterator,bool> status = values->insert(dataEntry);
		if (!status.second) throw std::runtime_error("readInitialization()");
	}

	file.close();
}

void writeInitialization(const char* filename, const char* delimiter )
{
	std::ofstream stream(filename,std::ofstream::out);
	if (!stream.good()) throw std::runtime_error("writeInitialization()");

	// ***
	
	const char endline = '\n';
	
	stream << "LEVEL_LIMIT" << delimiter << 0 << endline;
	stream << endline;
	
	// ***
	
	stream << "INPUT_XY_EXPANSION" << delimiter << 1.5 << endline;
	stream << "INPUT_Z_EXPANSION" << delimiter << 1.1 << endline;
	
	stream << "INPUT_PLC_CONTOUR_RADIUS" << delimiter << 100 << endline;
	stream << "INPUT_PLC_CONTOUR_PHI" << delimiter << 16 << endline;
	stream << "INPUT_PLC_CONTOUR_HEIGHT" << delimiter << 100 << endline;

	stream << "INPUT_PLC_RECOVERY_RADIUS" << delimiter << 50 << endline;
	stream << "INPUT_PLC_RECOVERY_PHI" << delimiter << 32 << endline;
	stream << "INPUT_PLC_RECOVERY_HEIGHT" << delimiter << 1e-8 << endline;

	stream << "INPUT_CYLINDER_XY_EXPANSION" << delimiter << 7.5 << endline;
	stream << "INPUT_CYLINDER_Z_EXPANSION" << delimiter << 7.5 << endline;
	stream << "INPUT_CYLINDER_RESOLUTION" << delimiter << 32 << endline;

	stream << endline;

	// ***

	stream << "FIRST_LEVEL_TETGEN_PARAMS" << delimiter << "zpqQa" << endline;
	stream << "FIRST_LEVEL_VOLUME_CONSTRAINTS" << delimiter << 8e-27/6.0 << endline;
	stream << endline;

	// ***

	stream << "SECOND_LEVEL_TETGEN_PARAMS" << delimiter << "zpqQ" << endline;
	stream << "SECOND_LEVEL_MAXIMUM_CONSTRAINTS" << delimiter << std::numeric_limits<REAL>::max() << endline;
	stream << endline;
	
	// ***

	stream << "THIRD_LEVEL_CYLINDER_RADIUS" << delimiter << 4e-6 << endline;
	stream << "THIRD_LEVEL_CYLINDER_HEIGHT" << delimiter << 1e-6 << endline;
	stream << "THIRD_LEVEL_INNER_CYLINDER_RESOLUTION" << delimiter << 64 << endline;
	stream << "THIRD_LEVEL_OUTWARD_CYLINDER_RESOLUTION" << delimiter << 64 << endline;
	stream << "THIRD_LEVEL_MAXIMUM_CONSTRAINTS" << delimiter << std::numeric_limits<REAL>::max() << endline;
	stream << endline;

	// ***

	stream << "FOURTH_LEVEL_CYLINDER_RADIUS" << delimiter << 4e-1 << endline;
	stream << "FOURTH_LEVEL_CYLINDER_HEIGHT" << delimiter << 1e-1 << endline;
	stream << "FOURTH_LEVEL_INNER_CYLINDER_RESOLUTION" << delimiter << 32 << endline;
	stream << "FOURTH_LEVEL_OUTWARD_CYLINDER_RESOLUTION" << delimiter << 512 << endline;
	stream << "FOURTH_LEVEL_MAXIMUM_CONSTRAINTS" << delimiter << (1e-6/6.0) << endline;
	stream << endline;

	// ***

	stream << "MINIMUM_VOLUME" << delimiter << (1.25e-26/6.0) << endline;
	stream << "VOLUME_DECAY_EXPONENT" << delimiter << 3.0 << endline;
	stream << endline;
	
	stream << "DISTANCE_MEASURE_RADIUS_RESOLUTION" << delimiter << 1000 << endline;
	stream << "DISTANCE_MEASURE_THETA_RESOLUTION" << delimiter << 32 << endline;
	stream << "DISTANCE_MEASURE_PHI_RESOLUTION" << delimiter << 64 << endline;

	stream << endline;

	// ***

	stream << "MAXIMUM_REFINE_CYCLES" << delimiter << 15 << endline;
	stream << "REFINE_THRESHOLD" << delimiter << 0.0 << endline;
	stream << endline;

	// ***

	stream << "MAXIMUM_EDGE_FILTER_CYCLES" << delimiter << 15 << endline;
	stream << "EDGE_FILTER_LENGTH" << delimiter << 5e-10 << endline;
	stream << endline;

	// ***

	stream << "VACUUM_MARKER" << delimiter << VACUUM << endline;
	stream << "TOP_MARKER" << delimiter << TOP << endline;
	stream << "BOTTOM_MARKER" << delimiter << BOTTOM << endline;
	stream << "SIDES_MARKER" << delimiter << SIDES << endline;
	stream << endline;

	// ***

	stream << "ASCII_OUTPUT_MODE" << delimiter << 0 << endline;
	stream << endline;

	// ***

	stream.close();
}

void readSample(const char* filename, TetGen::tetgenio* obj)
{
	std::ifstream input(filename,std::ifstream::in);
	if (!input.good()) throw std::runtime_error("readSample()");

	std::string firstLine;
	std::getline(input,firstLine);
	
	bool asciiMode;
	if (firstLine.substr(0,5) == "ASCII")
	{
		asciiMode = true;
		firstLine = firstLine.substr(5,std::string::npos);
	}
	else if (firstLine.substr(0,6) == "BINARY")
	{
		asciiMode = false;
		firstLine = firstLine.substr(6,std::string::npos);
	}
	else
		throw std::runtime_error("readSample()");
	
	int size;

	const bool withIds(true);
	
	bool withNumbers;
	bool withPotentials;
	
	{
		std::istringstream iParams(firstLine);
		
		iParams >> size;
		iParams >> std::noboolalpha >> withNumbers;
		iParams >> std::noboolalpha >> withPotentials;
		
		if (iParams.fail()) throw std::runtime_error("readSample()");
	}
	
	// *** initialize tetgen object

	obj->initialize();
	
	obj->firstnumber = 0;
	obj->numberofpoints = size;
	obj->pointlist = new REAL[obj->numberofpoints*3];

	if (withIds)
		obj->pointmarkerlist = new int[obj->numberofpoints];

	if (withNumbers)
	{
		obj->numberofpointattributes = 1;
		obj->pointattributelist = new REAL[obj->numberofpoints];
	}

	if (asciiMode)
	{
		for (int i =0 ; i < obj->numberofpoints; i++)
		{
			const char spacer = '\t';

			std::string buffer;
			std::getline(input,buffer);

			std::list<std::string> values;
			std::istringstream line(buffer);

			do
			{
				std::string tmp;
				std::getline(line,tmp,spacer);
				values.push_back(tmp);
			}
			while (line.good());

			// ***

			std::list<std::string>::const_iterator k = values.begin();

			for (int j = 0; j < 3; j++)
			{
				if (values.end() == k) throw std::runtime_error("readSample()");

				float tmpCoordinate;
				std::istringstream tmpStream(*k++);
				tmpStream >> tmpCoordinate;
				if (tmpStream.fail()) throw std::runtime_error("readSample()");

				obj->pointlist[i*3+j] = static_cast<REAL>(tmpCoordinate);
			}

			if (withIds)
			{
				if (values.end() == k) throw std::runtime_error("readSample()");

				short tmpId;
				std::istringstream tmpStream(*k++);
				tmpStream >> tmpId;
				if (tmpStream.fail()) throw std::runtime_error("readSample()");

				obj->pointmarkerlist[i] = tmpId;
			}

			if (withNumbers)
			{
				if (values.end() == k) throw std::runtime_error("readSample()");
	
				unsigned int tmpNumber;
				std::istringstream tmpStream(*k++);
				tmpStream >> tmpNumber;
				if (tmpStream.fail()) throw std::runtime_error("readSample()");

				obj->pointattributelist[i] = static_cast<REAL>(tmpNumber);
			}

			if (withPotentials)
			{
				if (values.end() == k) throw std::runtime_error("readSample()");
	
				float dummyValue;
				std::istringstream tmpStream(*k++);
				tmpStream >> dummyValue;
				if (tmpStream.fail()) throw std::runtime_error("readSample()");
			}
			
			if (values.end() != k) throw std::runtime_error("readSample()"); // check wether all data is read
		}
	}
	else
	{
		for (int i = 0; i < obj->numberofpoints; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				float tmpCoordinate;
				input.read(reinterpret_cast<char*>(&tmpCoordinate),sizeof(float));
				obj->pointlist[i*3+j] = static_cast<REAL>(tmpCoordinate);
			}

			if (withIds)
			{
				short tmpId;
				input.read(reinterpret_cast<char*>(&tmpId),sizeof(short));
				obj->pointmarkerlist[i] = tmpId;
			}

			if (withNumbers)
			{
				unsigned int tmpNumber;
				input.read(reinterpret_cast<char*>(&tmpNumber),sizeof(unsigned int));
				obj->pointattributelist[i] = static_cast<REAL>(tmpNumber);
			}
			
			if (withPotentials)
			{
				float dummyValue;
				input.read(reinterpret_cast<char*>(&dummyValue),sizeof(float));
			}
		}
	}
	
	input.close();
}

void writeTapSim_Nodefile(const char* filename, const TetGen::tetgenio& obj, const int asciiMode =1)
{
	const char endline = '\n';

	std::ofstream output(filename,std::ifstream::out);
	if (!output.good()) throw std::runtime_error("writeTapSim_Nodefile()");

	if (asciiMode)
		output << "ASCII ";
	else
		output << "BINARY ";

	output << obj.numberofpoints << " ";
	
	if (obj.pointattributelist != 0)
		// write numbers but no potentials
		output << "1 0" << endline;
	else
		// neither write numbers nor potentials
		output << "0 0" << endline;

	if (obj.numberofpointattributes > 1) throw std::runtime_error("writeTapSim_Nodefile(): too many point attributes!");

	if (asciiMode == 1)
	{
		const char spacer = '\t';

		output.setf(std::ios_base::showbase);
		output.setf(std::ios_base::scientific);
		output.precision(10);

		for (int i = 0; i < obj.numberofpoints; i++)
		{
			for (int j = 0; j < 3; j++)
				output << obj.pointlist[i*3+j] << spacer;

			if (obj.pointmarkerlist != 0)
				output << obj.pointmarkerlist[i] << spacer;
			else
				output << short(-1) << spacer;

			if (obj.pointattributelist != 0)
				output << static_cast<unsigned int>(obj.pointattributelist[i]) << endline;
			else
				output << static_cast<unsigned int>(0) << endline;
		}
	}
	else if (asciiMode == 0)
	{
		for (int i = 0; i < obj.numberofpoints; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				const float tmpCoordinate = static_cast<float>(obj.pointlist[i*3+j]);
				output.write(reinterpret_cast<const char*>(&tmpCoordinate),sizeof(float));
			}

			if (obj.pointmarkerlist != 0)
			{
				const short tmpId(obj.pointmarkerlist[i]);
				output.write(reinterpret_cast<const char*>(&tmpId),sizeof(short));
			}
			else
			{
				//const short tmpId(-1); // converts to tapsim invalid id value
				//output.write(reinterpret_cast<const char*>(&tmpId),sizeof(short));
			}

			if (obj.pointattributelist != 0)
			{
				const unsigned int tmpNumber = static_cast<unsigned int>(obj.pointattributelist[i]);
				output.write(reinterpret_cast<const char*>(&tmpNumber),sizeof(unsigned int));
			}
			else
			{
				//const unsigned int tmpNumber(0); // converts to tapsim invalid number value
				//output.write(reinterpret_cast<const char*>(&tmpNumber),sizeof(unsigned int));
			}
		}
	}
	else
		throw std::runtime_error("writeTapSim_Nodefile()");

	output.close();
}

void writeCSV(const char* filename, const TetGen::tetgenio& obj)
{
	const char spacer = ',';
	const char endline = '\n';

	std::ofstream output(filename,std::ifstream::out);
	if (!output.good()) throw std::runtime_error("writeTapSim_Nodefile()");

	if (obj.numberofpointattributes > 1) throw std::runtime_error("writeTapSim_Nodefile(): too many point attributes!");

	output.setf(std::ios_base::showbase);
	output.setf(std::ios_base::scientific);
	output.precision(10);

	for (int i = 0; i < obj.numberofpoints; i++)
	{
		for (int j = 0; j < 3; j++)
			output << obj.pointlist[i*3+j] << spacer;

		if (obj.pointmarkerlist != 0)
			output << obj.pointmarkerlist[i] << spacer;
		else
			output << short(-1) << spacer;

		if (obj.pointattributelist != 0)
			output << static_cast<unsigned int>(obj.pointattributelist[i]) << endline;
		else
			output << static_cast<unsigned int>(0) << endline;
	}

	output.close();
}

struct ConfigTemplate_Params
{
	int vacuumId;
	int topId;
	int bottomId;
	int sidesId;
	
	std::set<int> userIds;
};

void writeConfig(const char* filename, const ConfigTemplate_Params& obj)
{
	const float defaultTemperature = 0.0f;
	const float defaultSamplePotential = 1000.0f;

	std::ofstream stream(filename,std::ofstream::out);
	if (stream.fail()) throw std::runtime_error("writeConfig()");

	const char* spacer = " = ";
	const char endline = '\n';

	stream.setf(std::ios_base::showbase);
	stream.setf(std::ios_base::scientific);
	stream.precision(5);

	stream << "ASCII" << endline;
	stream << endline;

	stream << "TEMPERATURE" << spacer << defaultTemperature << endline;
	stream << endline;

	{
		stream << "ID" << spacer << obj.vacuumId << endline;
		stream << "NAME" << spacer << "Vacuum" << endline;
		stream << "CHARGE_DENSITY" << spacer << 0.0f << endline;
		stream << "DIELECTRICITY" << spacer << 1.0f << endline;
		stream << "REMOVABLE" << spacer << 0 << endline;
		stream << "NEUMANN_BOUNDARY" << spacer << 0 << endline;
		stream << "DIRICHLET_BOUNDARY" <<  spacer << 0 << endline;
		stream << "POTENTIAL" << spacer << 0.0f << endline;
		stream << "MASS" << spacer << 0.0f << endline;
		stream << "EVAPORATION_CHARGE_STATE" << spacer << 0 << endline;
		stream << "EVAPORATION_FIELD_STRENGTH" << spacer << 0.0f << endline;
		stream << "EVAPORATION_ACTIVATION_ENERGY" << spacer << 0.0f << endline;
		stream << endline;
	}

	{
		stream << "ID" << spacer << obj.topId << endline;
		stream << "NAME" << spacer << "Top" << endline;
		stream << "CHARGE_DENSITY" << spacer << 0.0f << endline;
		stream << "DIELECTRICITY" << spacer << 1.0f << endline;
		stream << "REMOVABLE" << spacer << 0 << endline;
		stream << "NEUMANN_BOUNDARY" << spacer << 0 << endline;
		stream << "DIRICHLET_BOUNDARY" <<  spacer << 1 << endline;
		stream << "POTENTIAL" << spacer << 0.0f << endline;
		stream << "MASS" << spacer << 0.0f << endline;
		stream << "EVAPORATION_CHARGE_STATE" << spacer << 0 << endline;
		stream << "EVAPORATION_FIELD_STRENGTH" << spacer << 0.0f << endline;
		stream << "EVAPORATION_ACTIVATION_ENERGY" << spacer << 0.0f << endline;
		stream << endline;
	}

	{
		stream << "ID" << spacer << obj.bottomId << endline;
		stream << "NAME" << spacer << "Bottom" << endline;
		stream << "CHARGE_DENSITY" << spacer << 0.0f << endline;
		stream << "DIELECTRICITY" << spacer << 1.0f << endline;
		stream << "REMOVABLE" << spacer << 0 << endline;
		stream << "NEUMANN_BOUNDARY" << spacer << 1 << endline;
		stream << "DIRICHLET_BOUNDARY" <<  spacer << 0 << endline;
		stream << "POTENTIAL" << spacer << defaultSamplePotential << endline;
		stream << "MASS" << spacer << 0.0f << endline;
		stream << "EVAPORATION_CHARGE_STATE" << spacer << 0 << endline;
		stream << "EVAPORATION_FIELD_STRENGTH" << spacer << 0.0f << endline;
		stream << "EVAPORATION_ACTIVATION_ENERGY" << spacer << 0.0f << endline;
		stream << endline;
	}

	{
		stream << "ID" << spacer << obj.sidesId << endline;
		stream << "NAME" << spacer << "Sides" << endline;
		stream << "CHARGE_DENSITY" << spacer << 0.0f << endline;
		stream << "DIELECTRICITY" << spacer << 1.0f << endline;
		stream << "REMOVABLE" << spacer << 0 << endline;
		stream << "NEUMANN_BOUNDARY" << spacer << 0 << endline;
		stream << "DIRICHLET_BOUNDARY" <<  spacer << 1 << endline;
		stream << "POTENTIAL" << spacer << 0.0f << endline;
		stream << "MASS" << spacer << 0.0f << endline;
		stream << "EVAPORATION_CHARGE_STATE" << spacer << 0 << endline;
		stream << "EVAPORATION_FIELD_STRENGTH" << spacer << 0.0f << endline;
		stream << "EVAPORATION_ACTIVATION_ENERGY" << spacer << 0.0f << endline;
		stream << endline;
	}

	for (std::set<int>::const_iterator i = obj.userIds.begin(); i != obj.userIds.end(); i++)
	{
		stream << "ID" << spacer << *i << endline;
		stream << "NAME" << spacer << "***_SET_NAME_HERE_***" << endline;
		stream << "CHARGE_DENSITY" << spacer << 0.0f << endline;
		stream << "DIELECTRICITY" << spacer << 1.0f << endline;
		stream << "REMOVABLE" << spacer << 1 << endline;
		stream << "NEUMANN_BOUNDARY" << spacer << 0 << endline;
		stream << "DIRICHLET_BOUNDARY" <<  spacer << 1 << endline;
		stream << "POTENTIAL" << spacer << defaultSamplePotential << endline;
		stream << "MASS" << spacer << "***_SET_MASS_HERE_***" << endline;
		stream << "EVAPORATION_CHARGE_STATE" << spacer << 1 << endline;
		stream << "EVAPORATION_FIELD_STRENGTH" << spacer << "***_SET_EVAPORATION_FIELD_STRENGTH_HERE_***" << endline;
		stream << "EVAPORATION_ACTIVATION_ENERGY" << spacer << 1.0f << endline;
		stream << endline;
	}

	stream.close();
}

int testLocation(const TetGen::tetgenio& obj, const REAL* point, const int in =-1, const int out =1, const int boundary =0, const REAL* min =0, const REAL* max =0)
{
	if (min != 0 && max != 0)
	{
		for (int j = 0; j < 3; j++)
			if (point[j] < min[j] || point[j] > max[j]) return out;
	}

	bool isBoundary(false);
	for (int i = 0; i < obj.numberoftrifaces; i++)
	{
		REAL* aVertex = &obj.pointlist[obj.trifacelist[i*3]*3];
		REAL* bVertex = &obj.pointlist[obj.trifacelist[i*3+1]*3];
		REAL* cVertex = &obj.pointlist[obj.trifacelist[i*3+2]*3];
		
		const REAL orient = TetGen::orient3d(aVertex,bVertex,cVertex,const_cast<REAL*>(point));
		
		if (orient > 0.0f) 
			return out;
		else
			if (orient == 0.0) isBoundary = true;
	}

	if (isBoundary) 
		return boundary;
	else
		return in;
}

bool testLocationWithin(const TetGen::tetgenio& obj, const REAL* point, const REAL* min =0, const REAL* max =0)
{
	// *** location on the boundary counts as within
	// *** => result of orient3d(): '>' outside, '<=' inside

	const int result = testLocation(obj,point,-1,1,0,min,max);
	
	if (result == 1)
		return false;
	else
		return true;
}

void markPoints(const TetGen::tetgenio& obj, TetGen::tetgenio* pointSource, const int* marker, int* cnts =0)
{
	// *** Sets the pointmarker values of the intput 'pointSource' to 'insideMarker' for points within the
	// *** convex hull of 'obj'and to 'outsideMarker' for points without or on the convex hull of 'obj'.

	// marker[0] => insideMarker, marker[1] => outsideMarker, marker[2] => boundaryMarker
	// cnts[0] => insideCnts, cnts[1] => outsideCnts, cnts[2] => boundaryCnts

	if (marker == 0) throw std::runtime_error("markPoints(): invalid pointer to marker");

	REAL min[3], max[3];

	if (obj.numberoftrifaces == 0 || obj.trifacelist == 0)
		throw std::runtime_error("markPoints(): no trifaces!");

	for (int i = 0; i < 3; i++)
	{
		min[i] = obj.pointlist[obj.trifacelist[0]*3+i];
		max[i] = obj.pointlist[obj.trifacelist[0]*3+i];
	}

	for (int i = 0; i < obj.numberoftrifaces*3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (obj.pointlist[obj.trifacelist[i]*3+j] < min[j])
				min[j] =obj.pointlist[obj.trifacelist[i]*3+j];
			else
				if (obj.pointlist[obj.trifacelist[i]*3+j] > max[j]) max[j] = obj.pointlist[obj.trifacelist[i]*3+j];
		}
	}

	// ***

	if (pointSource->pointmarkerlist == 0)
		pointSource->pointmarkerlist = new int[pointSource->numberofpoints];

	int withinCnt(0);
	int boundaryCnt(0);
	for (int i = 0; i < pointSource->numberofpoints; i++)
	{
		bool isWithin(true);
		bool isBoundary(false);

		for (int j = 0; j < 3; j++)
		{
			if (pointSource->pointlist[i*3+j] < min[j] || pointSource->pointlist[i*3+j] > max[j])
			{
				isWithin = false;
				break;
			}
		}

		if (isWithin)
		{
			for (int j = 0; j < obj.numberoftrifaces; j++)
			{
				REAL* aVertex = &obj.pointlist[obj.trifacelist[j*3]*3];
				REAL* bVertex = &obj.pointlist[obj.trifacelist[j*3+1]*3];
				REAL* cVertex = &obj.pointlist[obj.trifacelist[j*3+2]*3];

				const REAL orient = TetGen::orient3d(aVertex,bVertex,cVertex,&pointSource->pointlist[i*3]);
				
				if (orient > 0.0)
				{
					isWithin = false;
					break;
				}
				else 
					if (orient == 0.0) isBoundary = true;
			}
		}

		if (isWithin)
		{
			if (isBoundary)
			{
				boundaryCnt++;
				pointSource->pointmarkerlist[i] = marker[2];
			}
			else
			{
				withinCnt++;
				pointSource->pointmarkerlist[i] = marker[0];
			}
		}
		else
			pointSource->pointmarkerlist[i] = marker[1];
	}

	if (cnts != 0)
	{
		cnts[0] = withinCnt;
		cnts[1] = pointSource->numberofpoints - withinCnt - boundaryCnt;
		cnts[2] = boundaryCnt;
	}
}

int markPointsWithin(const TetGen::tetgenio& obj, TetGen::tetgenio* pointSource, const int insideMarker, const int outsideMarker)
{
	// *** location on the boundary counts as within
	// *** => result of orient3d(): '>' outside, '<=' inside
	
	int marker[3] = { insideMarker, outsideMarker, insideMarker };
	int cnts[3];
	
	markPoints(obj,pointSource,marker,cnts);
	
	return cnts[0]+cnts[2];
}

int removePoints(const TetGen::tetgenio& in, TetGen::tetgenio* unequal, TetGen::tetgenio* equal, const int marker)
{
	if (in.pointmarkerlist == 0)
		throw std::runtime_error("removePoints(): no pointmarkerlist!");
	
	int equalCnt(0);
	for (int i = 0; i < in.numberofpoints; i++)
		if (in.pointmarkerlist[i] == marker) equalCnt++;

	if (unequal == 0 && equal == 0) return equalCnt;

	// ***

	int unequalNumber = in.numberofpoints - equalCnt;

	REAL* unequalPointlist = new REAL[unequalNumber*3];

	int* unequalPointmarkerlist = new int[unequalNumber];

	REAL* unequalAttributes;
	if (in.pointattributelist != 0 && in.numberofpointattributes > 0)
		unequalAttributes = new REAL[unequalNumber*in.numberofpointattributes];
	else
		unequalAttributes = 0;

	// ***

	int equalNumber = equalCnt;
	
	REAL* equalPointlist = new REAL[equalNumber*3];

	int* equalPointmarkerlist = new int[equalNumber];

	REAL* equalAttributes;
	if (in.pointattributelist != 0 && in.numberofpointattributes > 0)
		equalAttributes = new REAL[equalNumber*in.numberofpointattributes];
	else
		equalAttributes = 0;

	// ***

	int unequalIndex(0);
	int equalIndex(0);

	for (int i = 0; i < in.numberofpoints; i++)
	{
		if (in.pointmarkerlist[i] != marker)
		{
			for (int j = 0; j < 3; j++)
				unequalPointlist[unequalIndex*3+j] = in.pointlist[i*3+j];

			unequalPointmarkerlist[unequalIndex] = in.pointmarkerlist[i];

			if (in.pointattributelist != 0 && in.numberofpointattributes > 0)
			{
				for (int j = 0; j < in.numberofpointattributes; j++)
					unequalAttributes[unequalIndex*in.numberofpointattributes+j] = in.pointattributelist[i*in.numberofpointattributes+j];
			}

			unequalIndex++;
		}
		else
		{
			for (int j = 0; j < 3; j++)
				equalPointlist[equalIndex*3+j] = in.pointlist[i*3+j];

			equalPointmarkerlist[equalIndex] = in.pointmarkerlist[i];

			if (in.pointattributelist != 0 && in.numberofpointattributes > 0)
			{
				for (int j = 0; j < in.numberofpointattributes; j++)
					equalAttributes[equalIndex*in.numberofpointattributes+j] = in.pointattributelist[i*in.numberofpointattributes+j];
			}

			equalIndex++;
		}
	}

	// ***

	if (unequal != 0)
	{
		unequal->initialize();
		unequal->numberofpoints = unequalNumber;
		unequal->pointlist = unequalPointlist;
		unequal->pointmarkerlist = unequalPointmarkerlist;
		unequal->numberofpointattributes = in.numberofpointattributes;
		unequal->pointattributelist = unequalAttributes;
	}
	else
	{
		delete[] unequalPointlist;
		delete[] unequalAttributes;
		delete[] unequalPointmarkerlist;
	}
	
	if (equal != 0)
	{
		equal->initialize();
		equal->numberofpoints = equalNumber;
		equal->pointlist = equalPointlist;
		equal->pointmarkerlist = equalPointmarkerlist;
		equal->numberofpointattributes = in.numberofpointattributes;
		equal->pointattributelist = equalAttributes;
	}
	else
	{
		delete[] equalPointlist;
		delete[] equalAttributes;
		delete[] equalPointmarkerlist;
	}

	return equalCnt;
}

void addPoints(TetGen::tetgenio* destination, const TetGen::tetgenio& source, const int dummyMarker =0, const REAL dummyAttribute =0.0)
{
	if (destination == 0)
		throw std::runtime_error("addPoints(): no output object!");

	// ***

	int numberofpoints = destination->numberofpoints + source.numberofpoints;
	REAL* pointlist = new REAL[numberofpoints*3];

	int* pointmarkerlist;
	if (destination->pointmarkerlist != 0 || source.pointmarkerlist != 0)
		pointmarkerlist = new int[numberofpoints];
	else
		pointmarkerlist = 0;

	int numberofpointattributes;
	if (destination->numberofpointattributes > 0 || source.numberofpointattributes > 0)
		numberofpointattributes = std::max(destination->numberofpointattributes,source.numberofpointattributes);
	else
		numberofpointattributes = 0;

	REAL* pointattributelist;
	if (destination->pointattributelist != 0 || source.pointattributelist != 0)
		pointattributelist = new REAL[numberofpoints*numberofpointattributes];
	else
		pointattributelist = 0;

	// ***
	
	int index(0);
	
	for (int i = 0; i < destination->numberofpoints; i++)
	{
		for (int j = 0; j < 3; j++)
			pointlist[index*3+j] = destination->pointlist[i*3+j];

		if (pointmarkerlist != 0)
		{
			if (destination->pointmarkerlist != 0)
				pointmarkerlist[index] = destination->pointmarkerlist[i];
			else
				pointmarkerlist[index] = dummyMarker;
		}

		if (pointattributelist != 0)
		{
			if (destination->pointattributelist != 0)
			{
				for (int j = 0; j < destination->numberofpointattributes; j++)
					pointattributelist[index*numberofpointattributes+j] = destination->pointattributelist[i*destination->numberofpointattributes+j];

				for (int j = destination->numberofpointattributes; j < numberofpointattributes; j++)
					pointattributelist[index*numberofpointattributes+j] = dummyAttribute;
			}
			else
			{
				for (int j = 0; j < numberofpointattributes; j++)
					pointattributelist[index*numberofpointattributes+j] = dummyAttribute;
			}
		}

		index++;
	}

	for (int i = 0; i < source.numberofpoints; i++)
	{
		for (int j = 0; j < 3; j++)
			pointlist[index*3+j] = source.pointlist[i*3+j];

		if (pointmarkerlist != 0)
		{
			if (source.pointmarkerlist != 0)
				pointmarkerlist[index] = source.pointmarkerlist[i];
			else
				pointmarkerlist[index] = dummyMarker;
		}

		if (pointattributelist != 0)
		{
			if (destination->pointattributelist != 0)
			{
				for (int j = 0; j < source.numberofpointattributes; j++)
					pointattributelist[index*numberofpointattributes+j] = source.pointattributelist[i*source.numberofpointattributes+j];

				for (int j = source.numberofpointattributes; j < numberofpointattributes; j++)
					pointattributelist[index*numberofpointattributes+j] = dummyAttribute;
			}
			else
			{
				for (int j = 0; j < numberofpointattributes; j++)
					pointattributelist[index*numberofpointattributes+j] = dummyAttribute;
			}
		}

		index++;
	}

	// ***
	
	destination->deinitialize();
	destination->initialize();

	destination->numberofpoints = numberofpoints;
	destination->pointlist = pointlist;
	destination->pointmarkerlist = pointmarkerlist;
	destination->numberofpointattributes = numberofpointattributes;
	destination->pointattributelist = pointattributelist;
}

void setPointMarker(TetGen::tetgenio* obj, const int marker)
{
	if (obj == 0) throw std::runtime_error("setPointMarker()");
	
	if (obj->pointmarkerlist == 0)
	{
		if (obj->numberofpoints > 0)
			obj->pointmarkerlist = new int[obj->numberofpoints];
		else
			return;
	}

	for (int i = 0; i < obj->numberofpoints; i++)
		obj->pointmarkerlist[i] = marker;
}

class EdgeList
{
	public:
		EdgeList()
			: _size(0),
			  _edges(0)
		{}
		
		EdgeList(const TetGen::tetgenio& obj)
			: _size(0),
			  _edges(0)
		{
			init(obj);
		}
		
		~EdgeList()
		{
			if (_edges != 0) delete[] _edges;
		}

		void init(const TetGen::tetgenio& obj)
		{
			if (_size != 0) _size = 0;

			if (_edges != 0) 
			{
				delete[] _edges;
				_edges = 0;
			}
			
			// ***

			std::map<int,IntSet> edgeTree;

			for (int i = 0; i < obj.numberoftetrahedra; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					for (int k = j+1; k < 4; k++)
					{
						int n1 = obj.tetrahedronlist[i*4+j];
						int n2 = obj.tetrahedronlist[i*4+k];

						if (n2 > n1) std::swap<int>(n1,n2);

						edgeTree.operator[](n1).insert(n2);
					}
				}
			}

			// ***

			for (std::map<int,IntSet>::const_iterator i = edgeTree.begin(); i != edgeTree.end(); i++)
				_size += i->second.size();

			_edges = new int[_size*2];

			int index(0);
			for (std::map<int,IntSet>::const_iterator i = edgeTree.begin(); i != edgeTree.end(); i++)
			{
				for (IntSet::const_iterator j = i->second.begin(); j != i->second.end(); j++)
				{
					_edges[index*2] = i->first;
					_edges[index*2+1] = *j;
					
					index++;
				}
			}
		}

		// ***

		const int* edges() const { return _edges; }

		int e1(const int index) const { return _edges[index*2]; }
		int e2(const int index) const { return _edges[index*2+1]; }

		int size() const { return _size; }

	private:
		typedef std::set<int> IntSet;

		int _size;
		int* _edges;
};

int filterEdges(TetGen::tetgenio* obj, REAL edgeConstraint, TetGen::tetgenio* removedPoints)
{
	if (obj->tetrahedronlist == 0 || obj->numberoftetrahedra == 0)
		throw std::runtime_error("filterEdges()");
	
	edgeConstraint = std::pow(edgeConstraint,2.0);
	
	// ***
	
	EdgeList edgelist(*obj);
	
	// process edges - first run

	int* markerlist = new int[obj->numberofpoints];

	for (int j = 0; j < obj->numberofpoints; j++)
		markerlist[j] = 0;
	
	int* edgemarkerlist = new int[edgelist.size()];
	
	for (int j = 0; j < edgelist.size(); j++)
		edgemarkerlist[j] = 0;
	
	for (int j = 0; j < edgelist.size(); j++)
	{
		const int e1 = edgelist.e1(j);
		const int e2 = edgelist.e2(j);
		
		REAL distance(0.0);
		for (int k = 0; k < 3; k++)
			distance += std::pow(obj->pointlist[e1*3+k]-obj->pointlist[e2*3+k],2.0);
		
		if (edgeConstraint > distance)
		{
			edgemarkerlist[j] = 1;
			
			markerlist[e1]++;
			markerlist[e2]++;
		}
	}

	// process edges - second run

	for (int j = 0; j < obj->numberofpoints; j++)
		obj->pointmarkerlist[j] = 0;

	int removeCnt(0);
	for (int j = 0; j < edgelist.size(); j++)
	{
		if (edgemarkerlist[j] == 0) continue;
		
		const int e1 = edgelist.e1(j);
		const int e2 = edgelist.e2(j);
		
		if (markerlist[e1] > markerlist[e2])
			obj->pointmarkerlist[e1] = 1;
		else
			obj->pointmarkerlist[e2] = 1;
		
		removeCnt++;
	}

	if (removedPoints != 0)
		removePoints(*obj,obj,removedPoints,1);
	else
		removePoints(*obj,obj,0,1);

	// clean-up

	delete[] markerlist;
	delete[] edgemarkerlist;

	return removeCnt;
}

class plc_from_triface_hull
{
	public:
		plc_from_triface_hull()
			: _outwardRadius(),
			  _outwardHeight(),
			  _outwardPhiSegments(32),
			  _outwardTopMarker(-1),
			  _outwardSidesMarker(-1),
			  _outwardBottomMarker(-1),
			  _rContourResolution(100),
			  _phiContourResolution(16),
			  _zContourResolution(100),
			  _inner_rResolution(50),
			  _inner_phiResolution(32),
			  _inner_zResolution(),
			  _innerMarker(-1),
			  _point_marker(-1)
		{}

		void defineOutwardCylinder(const REAL radius, const REAL height)
		{
			_outwardRadius = radius;
			_outwardHeight = height;
		}

		void setOutwardCylinderResolution(const unsigned int phi)
		{
			_outwardPhiSegments = phi;
		}

		void setMarker_at_top(const int marker)
		{
			_outwardTopMarker = marker;
		}

		void setMarker_at_bottom(const int marker)
		{
			_outwardBottomMarker = marker;
		}

		void setMarker_at_sides(const int marker)
		{
			_outwardSidesMarker = marker;
		}

		void defineInnerContour(const unsigned int radius, const unsigned phi, const unsigned int height)
		{
			_rContourResolution = radius;
			_phiContourResolution = phi;
			_zContourResolution = height;
		}

		void defineInnerMesh(const unsigned int radius, const unsigned int phi, const REAL z)
		{
			_inner_rResolution = radius;
			_inner_phiResolution = phi;
			_inner_zResolution = z;
		}

		void setInnerMarker(const int marker)
		{
			_innerMarker = marker;
		}

		void setPointMarker(const int marker)
		{
			_point_marker = marker;
		}
		
		void generateUnbounded(TetGen::tetgenio* obj, const TetGen::tetgenio& input) const
		{
			const unsigned int zSize = trifacePLC(obj,input);

			const bool withFacetMarker = _innerMarker >= 0 && _outwardBottomMarker >= 0 && _outwardTopMarker >= 0 && _outwardSidesMarker >= 0;
			
			// *** construct facets
			
			obj->numberoffacets = 1;
			obj->numberoffacets += 2*(zSize-1)*_inner_phiResolution;
			obj->numberoffacets += _inner_phiResolution;

			obj->facetlist = new TetGen::tetgenio::facet[obj->numberoffacets];
			obj->numberofholes = 0;
			obj->holelist = 0;

			if (withFacetMarker) obj->facetmarkerlist = new int[obj->numberoffacets];

			TetGen::tetgenio::facet* fPtr;
			TetGen::tetgenio::polygon* pPtr;

			// bottom facet

			{
				fPtr = &obj->facetlist[0];
				fPtr->numberofpolygons = 1;
				fPtr->polygonlist = new TetGen::tetgenio::polygon[fPtr->numberofpolygons];
				fPtr->numberofholes = 0;
				fPtr->holelist = 0;

				pPtr = &fPtr->polygonlist[0];
				pPtr->numberofvertices = _inner_phiResolution;
				pPtr->vertexlist = new int[pPtr->numberofvertices];

				for (int i = 0; i < pPtr->numberofvertices; i++)
					pPtr->vertexlist[i] = i;

				if (withFacetMarker) obj->facetmarkerlist[0] = _outwardBottomMarker;
			}

			// side facets

			for (unsigned int i = 1; i < zSize; i++)
			{
				for (unsigned int j = 0; j < _inner_phiResolution; j++)
				{
					fPtr = &obj->facetlist[1+2*(_inner_phiResolution*(i-1)+j)];
					fPtr->numberofpolygons = 1;
					fPtr->polygonlist = new TetGen::tetgenio::polygon[fPtr->numberofpolygons];
					fPtr->numberofholes = 0;
					fPtr->holelist = 0;

					pPtr = &fPtr->polygonlist[0];
					pPtr->numberofvertices = 3;
					pPtr->vertexlist = new int[pPtr->numberofvertices];

					pPtr->vertexlist[0] = (i-1)*_inner_phiResolution+j;
					pPtr->vertexlist[1] = (i-1)*_inner_phiResolution+((j+1)%_inner_phiResolution);
					pPtr->vertexlist[2] = i*_inner_phiResolution+j;

					if (withFacetMarker) obj->facetmarkerlist[1+2*(_inner_phiResolution*(i-1)+j)] = 1+2*(_inner_phiResolution*(i-1)+j);

					// ***

					fPtr = &obj->facetlist[1+2*(_inner_phiResolution*(i-1)+j)+1];
					fPtr->numberofpolygons = 1;
					fPtr->polygonlist = new TetGen::tetgenio::polygon[fPtr->numberofpolygons];
					fPtr->numberofholes = 0;
					fPtr->holelist = 0;

					pPtr = &fPtr->polygonlist[0];
					pPtr->numberofvertices = 3;
					pPtr->vertexlist = new int[pPtr->numberofvertices];

					pPtr->vertexlist[0] = (i-1)*_inner_phiResolution+((j+1)%_inner_phiResolution);
					pPtr->vertexlist[1] = i*_inner_phiResolution+((j+1)%_inner_phiResolution);
					pPtr->vertexlist[2] = i*_inner_phiResolution+j;

					if (withFacetMarker) obj->facetmarkerlist[1+2*(_inner_phiResolution*(i-1)+j)+1] = _innerMarker;
				}
			}

			// top facets

			for (unsigned int i = 0; i < _inner_phiResolution; i++)
			{
				fPtr = &obj->facetlist[1+2*(zSize-1)*_inner_phiResolution+i];
				fPtr->numberofpolygons = 1;
				fPtr->polygonlist = new TetGen::tetgenio::polygon[fPtr->numberofpolygons];
				fPtr->numberofholes = 0;
				fPtr->holelist = 0;

				pPtr = &fPtr->polygonlist[0];
				pPtr->numberofvertices = 3;
				pPtr->vertexlist = new int[pPtr->numberofvertices];

				pPtr->vertexlist[0] = (zSize-1)*_inner_phiResolution+i;
				pPtr->vertexlist[1] = (zSize-1)*_inner_phiResolution+((i+1)%_inner_phiResolution);
				pPtr->vertexlist[2] = zSize*_inner_phiResolution;

				if (withFacetMarker) obj->facetmarkerlist[1+2*(zSize-1)*_inner_phiResolution+i] = _innerMarker;
			}
		}
		
		void generateBounded(TetGen::tetgenio* obj, const TetGen::tetgenio& input) const
		{
			TetGen::tetgenio tmpObj;
			generateUnbounded(&tmpObj,input);
			
			// *** check geometry

			for (int i = 0; i < tmpObj.numberofpoints; i++)
			{
				REAL tmpRadius(0.0);
				for (int j = 0; j < 2; j++)
					tmpRadius += std::pow(tmpObj.pointlist[i*3+j],2.0);

				tmpRadius = std::sqrt(tmpRadius);

				if (tmpRadius > _outwardRadius || tmpObj.pointlist[i*3+2] > _outwardHeight)
					throw std::runtime_error("plc_from_triface_hull::generateBounded(): invalid geometry!");
			}

			// ***
			
			const int outwardPoints = 2 * _outwardPhiSegments;

			obj->initialize();
			obj->firstnumber = 0;

			obj->numberofpoints = tmpObj.numberofpoints;
			obj->numberofpoints += outwardPoints;
			obj->pointlist = new REAL[obj->numberofpoints*3];

			if (_point_marker >= 0) obj->pointmarkerlist = new int[obj->numberofpoints];

			// *** copy existing data

			if (tmpObj.pointlist == 0)
				throw std::runtime_error("plc_from_triface_hull::generateBounded()");
			else
			{
				for (int i = 0; i < tmpObj.numberofpoints*3; i++)
					obj->pointlist[i] = tmpObj.pointlist[i];
			}

			if (_point_marker >= 0)
			{
				if (tmpObj.pointmarkerlist != 0)
				{
					for (int i = 0; i < tmpObj.numberofpoints; i++)
						obj->pointmarkerlist[i] = tmpObj.pointmarkerlist[i];
				}
				else
				{
					for (int i = 0; i< tmpObj.numberofpoints; i++)
						obj->pointmarkerlist[i] = 0;
				}
			}

			// *** add outward cylinder points

			{
				const REAL deltaPhi = REAL(2.0) * pi / _outwardPhiSegments;

				for (unsigned int i = 0; i < _outwardPhiSegments; i++)
				{
					// bottom points
					obj->pointlist[(tmpObj.numberofpoints+i)*3] = _outwardRadius * std::cos(deltaPhi*i);
					obj->pointlist[(tmpObj.numberofpoints+i)*3+1] = _outwardRadius * std::sin(deltaPhi*i);
					obj->pointlist[(tmpObj.numberofpoints+i)*3+2] = 0.0f;

					// top points
					obj->pointlist[(tmpObj.numberofpoints+_outwardPhiSegments+i)*3] = _outwardRadius * std::cos(deltaPhi*i);
					obj->pointlist[(tmpObj.numberofpoints+_outwardPhiSegments+i)*3+1] = _outwardRadius * std::sin(deltaPhi*i);
					obj->pointlist[(tmpObj.numberofpoints+_outwardPhiSegments+i)*3+2] = _outwardHeight;
				}

				if (_point_marker >= 0)
				{
					for (int i = tmpObj.numberofpoints; i < obj->numberofpoints; i++)
						obj->pointmarkerlist[i] = _point_marker;
				}
			}

			//*** construct facets

			const bool withFacetMarker = _innerMarker >= 0 && _outwardBottomMarker >= 0 && _outwardTopMarker >= 0 && _outwardSidesMarker >= 0;

			const int outwardFaces = 1 + _outwardPhiSegments;

			obj->numberoffacets = tmpObj.numberoffacets;
			obj->numberoffacets += outwardFaces;
			obj->facetlist = new TetGen::tetgenio::facet[obj->numberoffacets];

			if (withFacetMarker) obj->facetmarkerlist = new int[obj->numberoffacets];

			// *** copy existing facet data

			if (tmpObj.facetlist == 0)
				throw std::runtime_error("plc_from_triface_hull::generateBounded()");
			else
			{
				TetGen::tetgenio::facet* src_fPtr;
				TetGen::tetgenio::polygon* src_pPtr;

				TetGen::tetgenio::facet* dest_fPtr;
				TetGen::tetgenio::polygon* dest_pPtr;
				
				for (int i = 0; i < tmpObj.numberoffacets; i++)
				{
					if (i == 0) continue; // skip bottom facet data

					dest_fPtr = &obj->facetlist[i];
					src_fPtr = &tmpObj.facetlist[i];

					dest_fPtr->numberofpolygons = src_fPtr->numberofpolygons;

					if (dest_fPtr->numberofpolygons == 0)
						dest_fPtr->polygonlist = 0;
					else
					{
						dest_fPtr->polygonlist = new TetGen::tetgenio::polygon[dest_fPtr->numberofpolygons];

						for (int j = 0; j < dest_fPtr->numberofpolygons; j++)
						{
							dest_pPtr = &dest_fPtr->polygonlist[j];
							src_pPtr = &src_fPtr->polygonlist[j];

							dest_pPtr->numberofvertices = src_pPtr->numberofvertices;

							if (dest_pPtr->numberofvertices == 0)
								dest_pPtr->vertexlist = 0;
							else
							{
								dest_pPtr->vertexlist = new int[dest_pPtr->numberofvertices];

								for (int k = 0; k < dest_pPtr->numberofvertices; k++)
									dest_pPtr->vertexlist[k] = src_pPtr->vertexlist[k];
							}
						}
					}

					dest_fPtr->numberofholes = src_fPtr->numberofholes;

					if (dest_fPtr->numberofholes == 0)
						dest_fPtr->holelist = 0;
					else
					{
						dest_fPtr->holelist = new REAL[dest_fPtr->numberofholes*3];

						for (int j = 0; j < dest_fPtr->numberofholes*3; j++)
							dest_fPtr->holelist[j] = src_fPtr->holelist[j];
					}
				}
			}

			if (withFacetMarker)
			{
				if (tmpObj.facetmarkerlist != 0)
				{
					for (int i = 0; i < tmpObj.numberoffacets; i++)
						obj->facetmarkerlist[i] = tmpObj.facetmarkerlist[i];
				}
				else
				{
					for (int i = 0; i < tmpObj.numberoffacets; i++)
						obj->facetmarkerlist[i] = 0;
				}
			}

			// *** add additional outward cylinder facets

			{
				TetGen::tetgenio::polygon* pPtr;
				TetGen::tetgenio::facet* fPtr;

				{
					// top
					fPtr = &obj->facetlist[tmpObj.numberoffacets];
					fPtr->numberofpolygons = 1;
					fPtr->polygonlist = new TetGen::tetgenio::polygon[fPtr->numberofpolygons];
					fPtr->numberofholes = 0;
					fPtr->holelist = 0;

					pPtr = &fPtr->polygonlist[0];
					pPtr->numberofvertices = _outwardPhiSegments;
					pPtr->vertexlist = new int[pPtr->numberofvertices];

					for (unsigned int i = 0; i < _outwardPhiSegments; i++)
						pPtr->vertexlist[i] = tmpObj.numberofpoints + _outwardPhiSegments + i;

					if (withFacetMarker) obj->facetmarkerlist[tmpObj.numberoffacets] = _outwardTopMarker;
				}

				{
					// sides

					for (unsigned int i = 0; i < _outwardPhiSegments; i++)
					{
						fPtr = &obj->facetlist[tmpObj.numberoffacets+1+i];
						fPtr->numberofpolygons = 1;
						fPtr->polygonlist = new TetGen::tetgenio::polygon[fPtr->numberofpolygons];
						fPtr->numberofholes = 0;
						fPtr->holelist = 0;

						pPtr = &fPtr->polygonlist[0];
						pPtr->numberofvertices = 4;
						pPtr->vertexlist = new int[pPtr->numberofvertices];

						pPtr->vertexlist[0] = tmpObj.numberofpoints + i;
						pPtr->vertexlist[1] = tmpObj.numberofpoints + ((i+1) % _outwardPhiSegments);
						pPtr->vertexlist[2] = tmpObj.numberofpoints + _outwardPhiSegments + ((i+1) % _outwardPhiSegments);
						pPtr->vertexlist[3] = tmpObj.numberofpoints + _outwardPhiSegments + i;

						if (withFacetMarker) obj->facetmarkerlist[tmpObj.numberoffacets+1+i] = _outwardSidesMarker;
					}
				}
			}

			// *** common bottom facet: inner hull + outward cylinder
			
			{
				TetGen::tetgenio::polygon* pPtr;
				TetGen::tetgenio::facet* fPtr;

				fPtr = &obj->facetlist[0];

				fPtr->numberofpolygons = 2;
				fPtr->polygonlist = new TetGen::tetgenio::polygon[fPtr->numberofpolygons];
				fPtr->numberofholes = 1;
				fPtr->holelist = new REAL[3];

				pPtr = &fPtr->polygonlist[0];
				pPtr->numberofvertices = _outwardPhiSegments;
				pPtr->vertexlist = new int[pPtr->numberofvertices];

				for (unsigned int i = 0; i < _outwardPhiSegments; i++)
					pPtr->vertexlist[i] = tmpObj.numberofpoints + i;

				pPtr = &fPtr->polygonlist[1];
				pPtr->numberofvertices = tmpObj.facetlist[0].polygonlist[0].numberofvertices;
				pPtr->vertexlist = new int[pPtr->numberofvertices];

				for (int i = 0; i < pPtr->numberofvertices; i++)
					pPtr->vertexlist[i] = tmpObj.facetlist[0].polygonlist[0].vertexlist[i];

				for (int i = 0; i < 3; i++)
					fPtr->holelist[i] = REAL(0.0);

				if (withFacetMarker) obj->facetmarkerlist[0] = _outwardBottomMarker;
			}

			// *** insert inner hole at the place of the sample

			REAL zCoordinate(0.0);
			for (int i = 0; i < tmpObj.numberofpoints; i++)
				zCoordinate += tmpObj.pointlist[i*3+2];

			zCoordinate /= tmpObj.numberofpoints;
			zCoordinate /= 2.0;

			obj->numberofholes = tmpObj.numberofholes + 1;
			obj->holelist = new REAL[obj->numberofholes*3];

			if (tmpObj.holelist != 0)
			{
				for (int i = 0; i < tmpObj.numberofholes*3; i++)
					obj->holelist[i] = tmpObj.holelist[i];
			}
			else
			{
				for (int i = 0; i < tmpObj.numberofholes*3; i++)
					obj->holelist[i] = REAL(0.0);
			}	
			
			obj->holelist[tmpObj.numberofholes*3] = REAL(0.0);
			obj->holelist[tmpObj.numberofholes*3+1] = REAL(0.0);
			obj->holelist[tmpObj.numberofholes*3+2] = zCoordinate;
		}
		
	private:
		unsigned int trifacePLC(TetGen::tetgenio* obj, const TetGen::tetgenio& input) const
		{
			obj->initialize();
			obj->firstnumber = 0;

			// *** unique list of node indizes of the convex hull

			std::set<int> trifaceNodes;
			for (int i = 0; i < input.numberoftrifaces*3; i++)
				trifaceNodes.insert(input.trifacelist[i]);

			// *** determine convex hull constraints

			REAL convexRadiusXY(0.0);
			REAL convexRadiusZ(0.0);

			{
				for (std::set<int>::const_iterator i = trifaceNodes.begin(); i != trifaceNodes.end(); i++)
				{
					REAL tmpRadius;
					tmpRadius = std::pow(input.pointlist[(*i)*3],2.0);
					tmpRadius += std::pow(input.pointlist[(*i)*3+1],2.0);

					if (tmpRadius > convexRadiusXY)
						convexRadiusXY = tmpRadius;

					if (input.pointlist[(*i)*3+2] > convexRadiusZ)
						convexRadiusZ = input.pointlist[(*i)*3+2];
				}

				convexRadiusXY = std::sqrt(convexRadiusXY);
			}

			// *** look for approximate z-contour of the convex hull

			std::vector<REAL> contourPoints(_zContourResolution*2,0.0); // allocate memory for pairs of z and r(z)

			{
				// *** find contour radius at the bottom

				contourPoints[0] = 0.0;
				contourPoints[1] = 0.0;

				for (std::set<int>::const_iterator i = trifaceNodes.begin(); i != trifaceNodes.end(); i++)
				{
					if (input.pointlist[(*i)*3+2] != 0.0) continue;

					REAL radius = std::pow(input.pointlist[(*i)*3],2.0);
					radius += std::pow(input.pointlist[(*i)*3+1],2.0);

					if (contourPoints[1] < radius) contourPoints[1] = radius;
				}

				contourPoints[1] = std::sqrt(contourPoints[1]);

				// *** find contour radius at the top

				contourPoints[(_zContourResolution-1)*2] = convexRadiusZ;
				contourPoints[(_zContourResolution-1)*2+1] = 0.0;

				for (std::set<int>::const_iterator i = trifaceNodes.begin(); i != trifaceNodes.end(); i++)
				{
					if (input.pointlist[(*i)*3+2] != convexRadiusZ) continue;

					REAL radius = std::pow(input.pointlist[(*i)*3],2.0);
					radius += std::pow(input.pointlist[(*i)*3+1],2.0);

					if (contourPoints[(_zContourResolution-1)*2+1] < radius) contourPoints[(_zContourResolution-1)*2+1] = radius;
				}

				contourPoints[(_zContourResolution-1)*2+1] = std::sqrt(contourPoints[(_zContourResolution-1)*2+1]);

				// ***

				for (unsigned int i = 1; i < _zContourResolution - 1; i++)
				{
					REAL samplePoint[3];
					samplePoint[2] = convexRadiusZ;
					samplePoint[2] /= _zContourResolution;
					samplePoint[2] *= i;

					contourPoints[i*2] = samplePoint[2];
					contourPoints[i*2+1] = REAL(0.0);

					for (unsigned int j = 0; j < _phiContourResolution; j++)
					{
						const REAL deltaPhi = 2.0 * pi / _phiContourResolution;
						const REAL phi = j * deltaPhi;

						samplePoint[0] = convexRadiusXY;
						samplePoint[0] /= _rContourResolution;
						samplePoint[0] *= std::cos(phi);
						samplePoint[1] = convexRadiusXY;
						samplePoint[1] /= _rContourResolution;
						samplePoint[1] *= std::sin(phi);

						REAL testPoint[3];

						testPoint[2] = samplePoint[2];

						int lowerBound(0);
						int upperBound(_rContourResolution);

						while (upperBound - lowerBound > 1)
						{
							int testBound = lowerBound;
							testBound += upperBound;
							testBound /= 2;

							for (int k = 0; k < 2; k++)
							{
								testPoint[k] = samplePoint[k];
								testPoint[k] *= testBound;
							}

							if (testLocation(input,testPoint,-1,1,0) != 1)
								lowerBound = testBound;
							else
								upperBound = testBound;
						}

						REAL radius(upperBound);
						radius *= convexRadiusXY;
						radius /= _rContourResolution;

						contourPoints[i*2+1] = std::max(contourPoints[i*2+1],radius);
					}
				}
			}

			// *** find final z-coordinates relying on the contour data

			std::vector<REAL> zCoords;

			{
				zCoords.push_back(contourPoints[0]); // bottom z-coordinate

				unsigned int index = 1;
				while (index < _zContourResolution)
				{
					REAL theta;
					theta = contourPoints[index*2+1]-contourPoints[(index-1)*2+1];
					theta /= contourPoints[index*2]-contourPoints[(index-1)*2];
					theta = std::atan(theta);

					REAL deltaZ;
					deltaZ = std::cos(theta);
					deltaZ *= _inner_zResolution;

					zCoords.push_back(zCoords.back()+deltaZ);

					while (zCoords.back() > contourPoints[index*2] && index < _zContourResolution)
						index++;
				}

				while (zCoords.back() > contourPoints[(_zContourResolution-1)*2])
					zCoords.pop_back();

				zCoords.push_back(contourPoints[(_zContourResolution-1)*2]);
			}

			// *** compute final coordinates

			obj->numberofpoints = 1 + zCoords.size() * _inner_phiResolution;
			obj->pointlist = new REAL[obj->numberofpoints*3];

			if (_point_marker >= 0) obj->pointmarkerlist = new int[obj->numberofpoints];

			{
				const REAL deltaPhi = 2.0 * pi / _inner_phiResolution;

				// intermediate points in between "]bottom,top]"
				
				for (unsigned int i = 1; i < zCoords.size(); i++)
				{
					REAL point[3];
					point[2] = zCoords.at(i);

					for (unsigned int j = 0; j < _inner_phiResolution; j++)
					{
						const REAL phi = j * deltaPhi;

						point[0] = convexRadiusXY;
						point[0] /= _inner_rResolution;
						point[0] *= std::cos(phi);
						point[1] = convexRadiusXY;
						point[1] /= _inner_rResolution;
						point[1] *= std::sin(phi);

						REAL testPoint[3];
						testPoint[2] = point[2];

						int lowerBound(0);
						int upperBound(_inner_rResolution);

						while (upperBound - lowerBound > 1)
						{
							int testBound = lowerBound;
							testBound += upperBound;
							testBound /= 2;

							for (int k = 0; k < 2; k++)
							{
								testPoint[k] = point[k];
								testPoint[k] *= testBound;
							}

							if (testLocation(input,testPoint,-1,1,0) != 1)
								lowerBound = testBound;
							else
								upperBound = testBound;
						}

						testPoint[0] = point[0];
						testPoint[0] *= upperBound;

						testPoint[1] = point[1];
						testPoint[1] *= upperBound;

						obj->pointlist[(i*_inner_phiResolution+j)*3] = testPoint[0];
						obj->pointlist[(i*_inner_phiResolution+j)*3+1] = testPoint[1];
						obj->pointlist[(i*_inner_phiResolution+j)*3+2] = testPoint[2];
					}
				}

				// construct bottom points => x,y coordinates are copied from next z-layer above the bottom one

				for (unsigned int i = 0; i < _inner_phiResolution; i++)
				{
					obj->pointlist[i*3] = obj->pointlist[(_inner_phiResolution+i)*3];
					obj->pointlist[i*3+1] = obj->pointlist[(_inner_phiResolution+i)*3+1];
					obj->pointlist[i*3+2] = zCoords.front();
				}
				
				// single point at top-center
				
				obj->pointlist[(zCoords.size()*_inner_phiResolution)*3] = 0.0;
				obj->pointlist[(zCoords.size()*_inner_phiResolution)*3+1] = 0.0;
				obj->pointlist[(zCoords.size()*_inner_phiResolution)*3+2] = zCoords.back();
			}

			if (obj->pointmarkerlist != 0)
			{
				for (int i = 0; i < obj->numberofpoints; i++)
					obj->pointmarkerlist[i] = _point_marker;
			}
			
			return zCoords.size();
		}


		// ***
		
		REAL _outwardRadius;
		REAL _outwardHeight;

		unsigned int _outwardPhiSegments;

		int _outwardTopMarker;
		int _outwardSidesMarker;
		int _outwardBottomMarker;

		unsigned int _rContourResolution;
		unsigned int _phiContourResolution;
		unsigned int _zContourResolution;
		
		unsigned int _inner_rResolution;
		unsigned int _inner_phiResolution;
		REAL _inner_zResolution;

		int _innerMarker;

		int _point_marker;
};

class plc_boundedCylinder
{
	public:
		plc_boundedCylinder()
			: _innerRadius(0.0),
			  _innerHeight(0.0),
			  _innerPhiSegments(32),
			  _innerMarker(),
			  _outwardRadius(0.0),
			  _outwardHeight(0.0),
			  _outwardPhiSegments(32),
			  _outwardTopMarker(-1),
			  _outwardSidesMarker(-1),
			  _outwardBottomMarker(-1),
			  _point_marker(-1)
		{}

		void defineInnerCylinder(const REAL radius, const REAL height)
		{
			_innerRadius = radius;
			_innerHeight = height;
		}

		void setInnerCylinderResolution(const unsigned int phi)
		{
			_innerPhiSegments = phi;
		}

		void defineOutwardCylinder(const REAL radius, const REAL height)
		{
			_outwardRadius = radius;
			_outwardHeight = height;
		}

		void setOutwardCylinderResolution(const unsigned int phi)
		{
			_outwardPhiSegments = phi;
		}

		void setMarker_at_top(const int marker)
		{
			_outwardTopMarker = marker;
		}

		void setMarker_at_bottom(const int marker)
		{
			_outwardBottomMarker = marker;
		}

		void setMarker_at_sides(const int marker)
		{
			_outwardSidesMarker = marker;
		}

		void setMarker_at_innerCylinder(const int marker)
		{
			_innerMarker = marker;
		}

		void setPointMarker(const int marker)
		{
			_point_marker = marker;
		}

		void generate(TetGen::tetgenio* obj) const
		{
			if (obj == 0) throw std::runtime_error("plc_boundedCylinder::generate()");

			// *** check geometry

			if (_innerRadius > _outwardRadius || _innerHeight > _outwardHeight)
				throw std::runtime_error("plc_boundedCylinder::generate(): invalid geometry!");

			// ***
			
			obj->initialize();

			obj->firstnumber = 0;

			const int innerPoints = 2 * _innerPhiSegments;
			const int outwardPoints = 2 * _outwardPhiSegments;

			obj->numberofpoints = innerPoints;
			obj->numberofpoints += outwardPoints;
			obj->pointlist = new REAL[obj->numberofpoints*3];
			obj->pointmarkerlist = new int[obj->numberofpoints];

			// *** inner cylinder points

			{
				const REAL deltaPhi = REAL(2.0) * pi / _innerPhiSegments;

				for (unsigned int i = 0; i < _innerPhiSegments; i++)
				{
					// bottom points
					obj->pointlist[i*3] = _innerRadius * std::cos(deltaPhi*i);
					obj->pointlist[i*3+1] = _innerRadius * std::sin(deltaPhi*i);
					obj->pointlist[i*3+2] = 0.0;

					// top points
					obj->pointlist[(_innerPhiSegments+i)*3] = _innerRadius * std::cos(deltaPhi*i);
					obj->pointlist[(_innerPhiSegments+i)*3+1] = _innerRadius * std::sin(deltaPhi*i);
					obj->pointlist[(_innerPhiSegments+i)*3+2] = _innerHeight;
				}
			}

			// *** outward cylinder points

			{
				const int point_offset = innerPoints;

				const REAL deltaPhi = REAL(2.0) * pi / _outwardPhiSegments;

				for (unsigned int i = 0; i < _outwardPhiSegments; i++)
				{
					// bottom points
					obj->pointlist[(point_offset+i)*3] = _outwardRadius * std::cos(deltaPhi*i);
					obj->pointlist[(point_offset+i)*3+1] = _outwardRadius * std::sin(deltaPhi*i);
					obj->pointlist[(point_offset+i)*3+2] = 0.0;

					// top points
					obj->pointlist[(point_offset+_outwardPhiSegments+i)*3] = _outwardRadius * std::cos(deltaPhi*i);
					obj->pointlist[(point_offset+_outwardPhiSegments+i)*3+1] = _outwardRadius * std::sin(deltaPhi*i);
					obj->pointlist[(point_offset+_outwardPhiSegments+i)*3+2] = _outwardHeight;
				}
			}

			if (_point_marker >= 0)
			{
				for (int i = 0; i < obj->numberofpoints; i++)
					obj->pointmarkerlist[i] = _point_marker;
			}

			// ***

			const int innerFaces = 1+ _innerPhiSegments;
			const int outwardFaces = 1 + _outwardPhiSegments;

			const bool withFacetMarker = _innerMarker >= 0 && _outwardBottomMarker >= 0 && _outwardTopMarker >= 0 && _outwardSidesMarker != 0;
			
			obj->numberoffacets = innerFaces;
			obj->numberoffacets += outwardFaces;
			obj->numberoffacets += 1;
			obj->facetlist = new TetGen::tetgenio::facet[obj->numberoffacets];

			if (withFacetMarker) obj->facetmarkerlist = new int[obj->numberoffacets];

			// *** define inner cylinder facets

			{
				TetGen::tetgenio::polygon* pPtr;
				TetGen::tetgenio::facet* fPtr;

				{
					// top
					fPtr = &obj->facetlist[0];
					fPtr->numberofpolygons = 1;
					fPtr->polygonlist = new TetGen::tetgenio::polygon[fPtr->numberofpolygons];
					fPtr->numberofholes = 0;
					fPtr->holelist = 0;

					pPtr = &fPtr->polygonlist[0];
					pPtr->numberofvertices = _innerPhiSegments;
					pPtr->vertexlist = new int[pPtr->numberofvertices];

					for (unsigned int i = 0; i < _innerPhiSegments; i++)
						pPtr->vertexlist[i] = _innerPhiSegments + i;

					if (withFacetMarker) obj->facetmarkerlist[0] = _innerMarker;
				}

				// sides
				for (unsigned int i = 0; i < _innerPhiSegments; i++)
				{
					fPtr = &obj->facetlist[1+i];
					fPtr->numberofpolygons = 1;
					fPtr->polygonlist = new TetGen::tetgenio::polygon[fPtr->numberofpolygons];
					fPtr->numberofholes = 0;
					fPtr->holelist = 0;

					pPtr = &fPtr->polygonlist[0];
					pPtr->numberofvertices = 4;
					pPtr->vertexlist = new int[pPtr->numberofvertices];

					pPtr->vertexlist[0] = i;
					pPtr->vertexlist[1] = (i+1) % _innerPhiSegments;
					pPtr->vertexlist[2] = _innerPhiSegments + (i+1) % _innerPhiSegments;
					pPtr->vertexlist[3] = _innerPhiSegments + i;

					if (withFacetMarker) obj->facetmarkerlist[1+i] = _innerMarker;
				}
			}

			// *** define outward cylinder facets

			{
				const int face_offset = innerFaces;
				const int point_offset = innerPoints;

				TetGen::tetgenio::polygon* pPtr;
				TetGen::tetgenio::facet* fPtr;

				{
					// top
					fPtr = &obj->facetlist[face_offset];
					fPtr->numberofpolygons = 1;
					fPtr->polygonlist = new TetGen::tetgenio::polygon[fPtr->numberofpolygons];
					fPtr->numberofholes = 0;
					fPtr->holelist = 0;

					pPtr = &fPtr->polygonlist[0];
					pPtr->numberofvertices = _outwardPhiSegments;
					pPtr->vertexlist = new int[pPtr->numberofvertices];

					for (unsigned int i = 0; i < _outwardPhiSegments; i++)
						pPtr->vertexlist[i] = point_offset + _outwardPhiSegments + i;

					if (withFacetMarker) obj->facetmarkerlist[face_offset] = _outwardTopMarker;
				}

				// sides
				for (unsigned int i = 0; i < _outwardPhiSegments; i++)
				{
					fPtr = &obj->facetlist[face_offset+1+i];
					fPtr->numberofpolygons = 1;
					fPtr->polygonlist = new TetGen::tetgenio::polygon[fPtr->numberofpolygons];
					fPtr->numberofholes = 0;
					fPtr->holelist = 0;

					pPtr = &fPtr->polygonlist[0];
					pPtr->numberofvertices = 4;
					pPtr->vertexlist = new int[pPtr->numberofvertices];

					pPtr->vertexlist[0] = point_offset + i;
					pPtr->vertexlist[1] = point_offset + (i+1) % _outwardPhiSegments;
					pPtr->vertexlist[2] = point_offset + _outwardPhiSegments + (i+1) % _outwardPhiSegments;
					pPtr->vertexlist[3] = point_offset + _outwardPhiSegments + i;

					if (withFacetMarker) obj->facetmarkerlist[face_offset+1+i] = _outwardSidesMarker;
				}
			}

			// *** common bottom facet: inner cylinder + outward cylinder

			{
				const int face_offset = innerFaces + outwardFaces;

				TetGen::tetgenio::polygon* pPtr;
				TetGen::tetgenio::facet* fPtr;

				fPtr = &obj->facetlist[face_offset];
				fPtr->numberofpolygons = 2;
				fPtr->polygonlist = new TetGen::tetgenio::polygon[fPtr->numberofpolygons];
				fPtr->numberofholes = 1;
				fPtr->holelist = new REAL[3];

				pPtr = &fPtr->polygonlist[0];
				pPtr->numberofvertices = _outwardPhiSegments;
				pPtr->vertexlist = new int[pPtr->numberofvertices];

				for (unsigned int i = 0; i < _outwardPhiSegments; i++)
					pPtr->vertexlist[i] = innerPoints + i;

				pPtr = &fPtr->polygonlist[1];
				pPtr->numberofvertices = _innerPhiSegments;
				pPtr->vertexlist = new int[pPtr->numberofvertices];

				for (unsigned int i = 0; i < _innerPhiSegments; i++)
					pPtr->vertexlist[i] = i;

				for (unsigned int i = 0; i < 3; i++)
					fPtr->holelist[i] = REAL(0.0);

				if (withFacetMarker) obj->facetmarkerlist[face_offset] = _outwardBottomMarker;
			}

			// *** insert inner hole at the place of the sample

			obj->numberofholes = 1;
			obj->holelist = new REAL[3];

			obj->holelist[0] = REAL(0.0);
			obj->holelist[1] = REAL(0.0);
			obj->holelist[2] = _innerHeight / REAL(2.0);
		}

	private:
			REAL _innerRadius;
			REAL _innerHeight;

			unsigned int _innerPhiSegments;

			int _innerMarker;

			REAL _outwardRadius;
			REAL _outwardHeight;

			unsigned int _outwardPhiSegments;

			int _outwardTopMarker;
			int _outwardSidesMarker;
			int _outwardBottomMarker;

			int _point_marker;
};

void trifaceExtents(const TetGen::tetgenio& obj, REAL* xy, REAL* z)
{
	*xy = 0.0;
	*z = 0.0;

	for (int i = 0; i < obj.numberoftrifaces; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			const int nodeIndex = obj.trifacelist[i*3+j];

			REAL tmpRadiusXY(0.0);
			for (int k = 0; k < 2; k++)
				tmpRadiusXY += obj.pointlist[nodeIndex*3+k] * obj.pointlist[nodeIndex*3+k];

			if (tmpRadiusXY > *xy) *xy = tmpRadiusXY;

			if (obj.pointlist[nodeIndex*3+2] > *z) *z = obj.pointlist[nodeIndex*3+2];
		}
	}

	*xy = std::sqrt(*xy);
}

class DistanceMap
{
	public:
		DistanceMap(const unsigned int radius =1000, const unsigned int theta =32, const unsigned int phi =64)
			: _radiusSegments(radius),
			  _thetaSegments(theta),
			  _phiSegments(phi),
			  _index(),
			  _zeroVolumeMeasure(1.0),
			  _decayExponent(1.0),
			  _maxConstraint(std::numeric_limits<REAL>::max())
		{}

		void setRadiusResolution(const unsigned int value)
		{
			_radiusSegments = value;
		}

		void setThetaResolution(const unsigned int value)
		{
			_thetaSegments = value;
		}
		
		void setPhiResolution(const unsigned int value)
		{
			_phiSegments = value;
		}

		void setVolumeMeasure(const REAL value)
		{
			_zeroVolumeMeasure = value;
		}

		void setVolumeDecayExponent(const REAL value)
		{
			_decayExponent = value;
		}

		void setVolumeMaxConstraint(const REAL value)
		{
			_maxConstraint = value;
		}
		
		void setReferenceObject(const TetGen::tetgenio& obj)
		{
			_index.clear();

			REAL maxRadius;

			{
				REAL extentsXY;
				REAL extentsZ;

				trifaceExtents(obj,&extentsXY,&extentsZ);

				maxRadius = std::max(extentsXY,extentsZ);
			}

			const REAL deltaTheta = pi / 2.0 / _thetaSegments;
			const REAL deltaPhi = 2.0 * pi / _phiSegments;
			
			for (unsigned int i = 1; i <= _thetaSegments; i++)
			{
				const REAL theta = i * deltaTheta;
				
				rrMap localIndex;
				for (unsigned int j = 0; j <= _phiSegments; j++)
				{
					const REAL phi = j * deltaPhi - pi;

					REAL srcPoint[3];
					srcPoint[0] = maxRadius / _radiusSegments * std::cos(theta) * std::cos(phi);
					srcPoint[1] = maxRadius / _radiusSegments * std::cos(theta) * std::sin(phi);
					srcPoint[2] = maxRadius / _radiusSegments * std::sin(theta);
					
					int lowerBound(0);
					int upperBound(_radiusSegments);

					REAL point[3];
					while (upperBound - lowerBound > 1)
					{
						int testBound = lowerBound;
						testBound += upperBound;
						testBound /= 2;

						for (int k = 0; k < 3; k++)
						{
							point[k] = srcPoint[k];
							point[k] *= testBound;
						}

						if (testLocation(obj,point,-1,1,0) != 1)
							lowerBound = testBound;
						else
							upperBound = testBound;
					}
					
					REAL distance(0.0);
					for (int j = 0; j < 3; j++)
						distance += std::pow(srcPoint[j]*lowerBound,2.0);

					distance = std::sqrt(distance);

					const std::pair<REAL,REAL> localEntry(phi,distance);
					localIndex.insert(localEntry);
				}

				if (i == 1)
				{
					const std::pair<REAL,rrMap> zeroEntry(0.0,localIndex);
					_index.insert(zeroEntry);
				}

				const std::pair<REAL,rrMap> entry(theta,localIndex);
				_index.insert(entry);
			}
		}

		REAL findDistance(const REAL* point, REAL* pointRadius =0) const
		{
			REAL radius, xyRadius;
			xyRadius = std::pow(point[0],2.0);
			xyRadius += std::pow(point[1],2.0);

			radius = xyRadius;
			radius += std::pow(point[2],2.0);
			
			xyRadius = std::sqrt(xyRadius);
			radius = std::sqrt(radius);

			// ***
			
			const REAL theta = std::atan(point[2]/xyRadius);
			if (theta < 0.0 || theta > pi / 2.0)
				throw std::runtime_error("DistanceMap::findDistance(): theta of of range!");

			const REAL phi = std::atan2(point[1],point[0]);

			std::map<REAL,rrMap>::const_iterator thetaPtr = _index.lower_bound(theta);
			if (_index.end() == thetaPtr)
				throw std::runtime_error("DistanceMap::findDistance(): Bad theta search!");

			std::map<REAL,REAL>::const_iterator phiPtr = thetaPtr->second.lower_bound(phi);
			if (thetaPtr->second.end() == phiPtr)
				throw std::runtime_error("DistanceMap::findDistance(): Bad phi search!");
			
			if (pointRadius != 0) *pointRadius = radius;

			return phiPtr->second;
		}

		int getVolumeConstraints(TetGen::tetgenio* obj)
		{
			if (obj->tetrahedronvolumelist == 0)
				obj->tetrahedronvolumelist = new REAL[obj->numberoftetrahedra];

			// ***

			int nonConforming(0);
			for (int i = 0; i < obj->numberoftetrahedra; i++)
			{
				const REAL actualValue = tetVolume(i,*obj);

				// ***

				REAL center[3];
				tetCenter(i,*obj,center);

				REAL r,r0;
				r0 = findDistance(center,&r);
				r = std::max(r,r0);

				// ***

				REAL nominalValue(_zeroVolumeMeasure);
				nominalValue *= std::pow(r/r0,_decayExponent);

				nominalValue = std::min<REAL>(nominalValue,_maxConstraint);

				if (nominalValue < actualValue)
				{
					obj->tetrahedronvolumelist[i] = nominalValue;
					nonConforming++;
				}
				else
					obj->tetrahedronvolumelist[i] = REAL(-1.0);
			}

			return nonConforming;
		}
		
	private:
		inline const REAL* tetCoords(const int tetIndex, const int tetVertex, const TetGen::tetgenio& obj)
		{
			return &obj.pointlist[obj.tetrahedronlist[tetIndex*4+tetVertex]*3];
		}

		void tetCenter(const int tetIndex, const TetGen::tetgenio& obj, REAL* center)
		{
			for (int i = 0; i < 3; i++)
			{
				center[i] = REAL(0.0);
				for (int j = 0; j < 4; j++)
					center[i] += tetCoords(tetIndex,j,obj)[i];
				center[i] /= REAL(4.0);
			}
		}

		REAL tetVolume(const int tetIndex, const TetGen::tetgenio& obj)
		{
			REAL a[3],b[3],c[3];

			for (int i = 0; i < 3; i++)
			{
				a[i] = tetCoords(tetIndex,0,obj)[i];
				a[i] -= tetCoords(tetIndex,3,obj)[i];

				b[i] = tetCoords(tetIndex,1,obj)[i];
				b[i] -= tetCoords(tetIndex,3,obj)[i];

				c[i] = tetCoords(tetIndex,2,obj)[i];
				c[i] -= tetCoords(tetIndex,3,obj)[i];
			}

			REAL value;
			value = a[0]*(b[1]*c[2]-b[2]*c[1]);
			value -= a[1]*(b[0]*c[2]-b[2]*c[0]);
			value += a[2]*(b[0]*c[1]-b[1]*c[0]);
			value /= REAL(6.0);

			return std::fabs(value);
		}

		typedef std::map<REAL,REAL> rrMap;

		unsigned int _radiusSegments;
		unsigned int _thetaSegments;
		unsigned int _phiSegments;
		
		std::map<REAL,rrMap> _index;

		// ***

		REAL _zeroVolumeMeasure;
		REAL _decayExponent;
		REAL _maxConstraint;
};

struct InputParams
{
	std::string filename;
	
	REAL xyExpansion;
	REAL zExpansion;

	unsigned int plcContourRadius;
	unsigned int plcContourPhi;
	unsigned int plcContourHeight;

	unsigned int plcRecoveryRadius;
	unsigned int plcRecoveryPhi;
	REAL plcRecoveryHeight;

	REAL cylinder_xyExpansion;
	REAL cylinder_zExpansion;
	unsigned int cylinder_resolution;
};

struct InputData
{
	InputData()
		: mesh(),
		  extentsXY(),
		  extentsZ(),
		  inflatedMesh(),
		  inflatedExtentsXY(),
		  inflatedExtentsZ()
	{}

	TetGen::tetgenio mesh;
	REAL extentsXY;
	REAL extentsZ;

	TetGen::tetgenio inflatedMesh;
	REAL inflatedExtentsXY;
	REAL inflatedExtentsZ;

	TetGen::tetgenio boundedPLC;
	TetGen::tetgenio unboundedPLC;
};

struct LevelData
{
	TetGen::tetgenio mesh;
	REAL extentsXY;
	REAL extentsZ;
};

int main(int argc, char** argv)
{
	const bool forceSinglePrecision(true);

	// enable additional file output

	bool meshOutput(false);
	bool edgeFilterOutput(false);

	// ***

	std::cout << "MeshGen 1.0b" << std::endl;
	std::cout << std::endl;
	std::cout << "Copyright (C) 2012 Christian Oberdorfer" << std::endl;
	std::cout << std::endl;
	
	// *** initialization

	unsigned int levelLimit;

	InputParams inputParams;

	std::string firstLevel_tetGenParams;
	REAL firstLevel_volumeConstraints;

	std::string secondLevel_tetGenParams;
	REAL secondLevel_maxConstraints;
	
	REAL thirdLevel_cylinderRadius;
	REAL thirdLevel_cylinderHeight;
	unsigned int thirdLevel_innerCylinderResolution;
	unsigned int thirdLevel_outwardCylinderResolution;
	REAL thirdLevel_maxConstraints;
	
	REAL fourthLevel_cylinderRadius;
	REAL fourthLevel_cylinderHeight;
	unsigned int fourthLevel_innerCylinderResolution;
	unsigned int fourthLevel_outwardCylinderResolution;
	REAL fourthLevel_maxConstraints;
	
	REAL zeroVolumeMeasure;
	REAL volumeDecayExponent;

	unsigned int refineCycles;
	REAL refineThreshold;

	unsigned int edgeFilter_cycles;
	REAL edgeFilter_length;
	
	unsigned int distance_radiusResolution;
	unsigned int distance_thetaResolution;
	unsigned int distance_phiResolution;

	int final_vacuumMarker;
	int final_topMarker;
	int final_bottomMarker;
	int final_sidesMarker;
	int final_outputMode;

	// ***

	bool skipIniFile(false);

	if (argc < 3)
		skipIniFile = true;
	else
	{
		for (int i = 1; i < argc; i++)
		{
			const char* aboutCmd = "--about";
			const char* helpCmd = "--help";

			if (std::strcmp(argv[i],aboutCmd) == 0 || std::strcmp(argv[i],helpCmd) == 0)
			{
				skipIniFile = true;
				break;
			}
		}
	}
	
	if (!skipIniFile)
	{
		const char* iniFilename = "meshgen.ini";
		const char* iniDelimiter = " = ";
		
		if (!std::ifstream(iniFilename,std::ifstream::in).good())
		{
			std::cout << "Cannot read initialization file!" << std::endl;
			std::cout << std::endl;

			while (true)
			{
				std::cout << "Create default file (y/n)? ";

				char value;
				std::cin >> value;

				if (value == 'y')
				{
					writeInitialization(iniFilename,iniDelimiter);

					std::cout << std::endl;
					std::cout << "Default initialization file was successfully written!" << std::endl;
					std::cout << std::endl;
					break;
				}
				else if (value == 'n')
				{
					std::cout << std::endl;
					return 0;
				}
			}

			while (true)
			{
				std::cout << "Quit (y/n)? ";

				char value;
				std::cin >> value;

				if (value == 'y')
				{
					std::cout << std::endl;
					return 0;
				}
				else if (value == 'n')
				{
					std::cout << std::endl;
					break;
				}
			}
		}
		
		std::map<std::string,std::string> iniParams;
		readInitialization(iniFilename,&iniParams,iniDelimiter);

		std::map<std::string,std::string>::iterator entry;

		// ***
		
		entry = iniParams.find("LEVEL_LIMIT");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%u",&levelLimit) != 1)
			throw std::runtime_error("Error in initialization data at entry 'LEVEL_LIMIT'!");

		entry = iniParams.find("INPUT_XY_EXPANSION");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%lf",&inputParams.xyExpansion) != 1)
			throw std::runtime_error("Error in initialization data at entry 'INPUT_XY_EXPANSION'!");

		entry = iniParams.find("INPUT_Z_EXPANSION");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%lf",&inputParams.zExpansion) != 1)
			throw std::runtime_error("Error in initialization data at entry 'INPUT_Z_EXPANSION'!");

		entry = iniParams.find("INPUT_PLC_CONTOUR_RADIUS");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%u",&inputParams.plcContourRadius) != 1)
			throw std::runtime_error("Error in initialization data at entry 'INPUT_PLC_CONTOUR_RADIUS'!");

		entry = iniParams.find("INPUT_PLC_CONTOUR_PHI");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%u",&inputParams.plcContourPhi) != 1)
			throw std::runtime_error("Error in initialization data at entry 'INPUT_PLC_CONTOUR_PHI'!");

		entry = iniParams.find("INPUT_PLC_CONTOUR_HEIGHT");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%u",&inputParams.plcContourHeight) != 1)
			throw std::runtime_error("Error in initialization data at entry 'INPUT_PLC_CONTOUR_HEIGHT'!");

		entry = iniParams.find("INPUT_PLC_RECOVERY_RADIUS");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%u",&inputParams.plcRecoveryRadius) != 1)
			throw std::runtime_error("Error in initialization data at entry 'INPUT_PLC_RECOVERY_RADIUS'!");

		entry = iniParams.find("INPUT_PLC_RECOVERY_PHI");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%u",&inputParams.plcRecoveryPhi) != 1)
			throw std::runtime_error("Error in initialization data at entry 'INPUT_PLC_RECOVERY_PHI'!");

		entry = iniParams.find("INPUT_PLC_RECOVERY_HEIGHT");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%lf",&inputParams.plcRecoveryHeight) != 1)
			throw std::runtime_error("Error in initialization data at entry 'INPUT_PLC_RECOVERY_HEIGHT'!");

		entry = iniParams.find("INPUT_CYLINDER_XY_EXPANSION");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%lf",&inputParams.cylinder_xyExpansion) != 1)
			throw std::runtime_error("Error in initialization data at entry 'INPUT_CYLINDER_XY_EXPANSION'!");

		entry = iniParams.find("INPUT_CYLINDER_Z_EXPANSION");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%lf",&inputParams.cylinder_zExpansion) != 1)
			throw std::runtime_error("Error in initialization data at entry 'INPUT_CYLINDER_Z_EXPANSION'!");

		entry = iniParams.find("INPUT_CYLINDER_RESOLUTION");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%u",&inputParams.cylinder_resolution) != 1)
			throw std::runtime_error("Error in initialization data at entry 'INPUT_CYLINDER_RESOLUTION'!");

		// ***

		entry = iniParams.find("FIRST_LEVEL_TETGEN_PARAMS");
		if (iniParams.end() == entry)
			throw std::runtime_error("Error in initialization data at entry 'FIRST_LEVEL_TETGEN_PARAMS'!");
		else
			firstLevel_tetGenParams = entry->second;

		entry = iniParams.find("FIRST_LEVEL_VOLUME_CONSTRAINTS");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%lf",&firstLevel_volumeConstraints) != 1)
			throw std::runtime_error("Error in initialization data at entry 'FIRST_LEVEL_VOLUME_CONSTRAINTS'!");

		// ***

		entry = iniParams.find("SECOND_LEVEL_TETGEN_PARAMS");
		if (iniParams.end() == entry)
			throw std::runtime_error("Error in initialization data at entry 'SECOND_LEVEL_TETGEN_PARAMS'!");
		else
			secondLevel_tetGenParams = entry->second;

		entry = iniParams.find("SECOND_LEVEL_MAXIMUM_CONSTRAINTS");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%lf",&secondLevel_maxConstraints) != 1)
			throw std::runtime_error("Error in initialization data at entry 'SECOND_LEVEL_MAXIMUM_CONSTRAINTS'!");

		// ***

		entry = iniParams.find("THIRD_LEVEL_CYLINDER_RADIUS");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%lf",&thirdLevel_cylinderRadius) != 1)
			throw std::runtime_error("Error in initialization data at entry 'THIRD_LEVEL_CYLINDER_RADIUS'!");

		entry = iniParams.find("THIRD_LEVEL_CYLINDER_HEIGHT");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%lf",&thirdLevel_cylinderHeight) != 1)
			throw std::runtime_error("Error in initialization data at entry 'THIRD_LEVEL_CYLINDER_HEIGHT'!");

		entry = iniParams.find("THIRD_LEVEL_INNER_CYLINDER_RESOLUTION");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%u",&thirdLevel_innerCylinderResolution) != 1)
			throw std::runtime_error("Error in initialization data at entry 'THIRD_LEVEL_INNER_CYLINDER_RESOLUTION'!");

		entry = iniParams.find("THIRD_LEVEL_OUTWARD_CYLINDER_RESOLUTION");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%u",&thirdLevel_outwardCylinderResolution) != 1)
			throw std::runtime_error("Error in initialization data at entry 'THIRD_LEVEL_OUTWARD_CYLINDER_RESOLUTION'!");

		entry = iniParams.find("THIRD_LEVEL_MAXIMUM_CONSTRAINTS");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%lf",&thirdLevel_maxConstraints) != 1)
			throw std::runtime_error("Error in initialization data at entry 'THIRD_LEVEL_MAXIMUM_CONSTRAINTS'!");

		// ***

		entry = iniParams.find("FOURTH_LEVEL_CYLINDER_RADIUS");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%lf",&fourthLevel_cylinderRadius) != 1)
			throw std::runtime_error("Error in initialization data at entry 'FOURTH_LEVEL_CYLINDER_RADIUS'!");

		entry = iniParams.find("FOURTH_LEVEL_CYLINDER_HEIGHT");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%lf",&fourthLevel_cylinderHeight) != 1)
			throw std::runtime_error("Error in initialization data at entry 'FOURTH_LEVEL_CYLINDER_HEIGHT'!");

		entry = iniParams.find("FOURTH_LEVEL_INNER_CYLINDER_RESOLUTION");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%u",&fourthLevel_innerCylinderResolution) != 1)
			throw std::runtime_error("Error in initialization data at entry 'FOURTH_LEVEL_INNER_CYLINDER_RESOLUTION'!");

		entry = iniParams.find("FOURTH_LEVEL_OUTWARD_CYLINDER_RESOLUTION");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%u",&fourthLevel_outwardCylinderResolution) != 1)
			throw std::runtime_error("Error in initialization data at entry 'FOURTH_LEVEL_OUTWARD_CYLINDER_RESOLUTION'!");

		entry = iniParams.find("FOURTH_LEVEL_MAXIMUM_CONSTRAINTS");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%lf",&fourthLevel_maxConstraints) != 1)
			throw std::runtime_error("Error in initialization data at entry 'FOURTH_LEVEL_MAXIMUM_CONSTRAINTS'!");

		// ***

		entry = iniParams.find("MINIMUM_VOLUME");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%lf",&zeroVolumeMeasure) != 1)
			throw std::runtime_error("Error in initialization data at entry 'DISTANCE_MEASURE_INITIAL_VOLUME'!");

		entry = iniParams.find("VOLUME_DECAY_EXPONENT");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%lf",&volumeDecayExponent) != 1)
			throw std::runtime_error("Error in initialization data at entry 'DISTANCE_MEASURE_VOLUME_DECAY_EXPONENT'!");
	
		entry = iniParams.find("DISTANCE_MEASURE_RADIUS_RESOLUTION");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%u",&distance_radiusResolution) != 1)
			throw std::runtime_error("Error in initialization data at entry 'DISTANCE_MEASURE_RADIUS_RESOLUTION'!");

		entry = iniParams.find("DISTANCE_MEASURE_THETA_RESOLUTION");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%u",&distance_thetaResolution) != 1)
			throw std::runtime_error("Error in initialization data at entry 'DISTANCE_MEASURE_THETA_RESOLUTION'!");

		entry = iniParams.find("DISTANCE_MEASURE_PHI_RESOLUTION");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%u",&distance_phiResolution) != 1)
			throw std::runtime_error("Error in initialization data at entry 'DISTANCE_MEASURE_PHI_RESOLUTION'!");

		// ***

		entry = iniParams.find("MAXIMUM_REFINE_CYCLES");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%u",&refineCycles) != 1)
			throw std::runtime_error("Error in initialization data at entry 'MAXIMUM_REFINE_CYCLES'!");

		entry = iniParams.find("REFINE_THRESHOLD");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%lf",&refineThreshold) != 1)
			throw std::runtime_error("Error in initialization data at entry 'REFINE_THRESHOLD'!");

		// ***


		entry = iniParams.find("MAXIMUM_EDGE_FILTER_CYCLES");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%u",&edgeFilter_cycles) != 1)
			throw std::runtime_error("Error in initialization data at entry 'EDGE_FILTER_CYCLES'!");

		entry = iniParams.find("EDGE_FILTER_LENGTH");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%lf",&edgeFilter_length) != 1)
			throw std::runtime_error("Error in initialization data at entry 'REFINE_THRESHOLD'!");

		// ***

		entry = iniParams.find("VACUUM_MARKER");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%d",&final_vacuumMarker) != 1)
			throw std::runtime_error("Error in initialization data at entry 'VACUUM_MARKER'!");

		entry = iniParams.find("TOP_MARKER");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%d",&final_topMarker) != 1)
			throw std::runtime_error("Error in initialization data at entry 'TOP_MARKER'!");

		entry = iniParams.find("BOTTOM_MARKER");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%d",&final_bottomMarker) != 1)
			throw std::runtime_error("Error in initialization data at entry 'BOTTOM_MARKER'!");

		entry = iniParams.find("SIDES_MARKER");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%d",&final_sidesMarker) != 1)
			throw std::runtime_error("Error in initialization data at entry 'SIDES_MARKER'!");

		// ***

		entry = iniParams.find("ASCII_OUTPUT_MODE");
		if (iniParams.end() == entry || std::sscanf(entry->second.c_str(),"%d",&final_outputMode) != 1)
			throw std::runtime_error("Error in initialization data at entry 'OUTPUT_MODE'!");
	}

	// *** control of mesh generation

	const bool withFirstLevel(true);

	bool withSecondLevel(true);
	bool withThirdLevel(true);
	bool withFourthLevel(true);
	
	if (levelLimit > 0)
	{
		if (levelLimit > 1) withSecondLevel = true; else withSecondLevel = false;
		if (levelLimit > 2) withThirdLevel = true; else withThirdLevel = false;
		if (levelLimit > 3) withFourthLevel = true; else withFourthLevel = false;
	}

	// *** parse command line parameters

	int verboseLevel(1);
	int writeSubgrids(0);
	
	std::string outputFile;
	std::string configFile;
	std::string csvFile;

	{
		int showHelp(0);
		int aboutFlag(0);

		struct option optionList[] =
		{
			{ "help", no_argument, &showHelp, 1},
			{ "about", no_argument, &aboutFlag, 1 },
			{ "quiet", no_argument, &verboseLevel, 0 },
			{ "write-ascii", no_argument, &final_outputMode, 1 },
			{ "write-binary", no_argument, &final_outputMode, 0 },
			{ "write-csv-subgrids", no_argument, &writeSubgrids, 1 },
			{ "write-csv", required_argument, 0, 1 },
			{ "verbose", optional_argument, 0, 2 },
			{ "initial-volume", required_argument, 0, 3 },
			{ "volume-decay-exponent", required_argument, 0, 4 },
			{ "length", required_argument, 0, 5 },
			{ "radius", required_argument, 0, 6 },
			{ "vacuum-marker", required_argument, 0, 7 },
			{ "top-marker", required_argument, 0, 8 },
			{ "bottom-marker", required_argument, 0, 9 },
			{ "sides-marker", required_argument, 0, 10 },
			{ "create-config-template", required_argument, 0, 11 },
			{ "debug-mode", no_argument, 0, 12 },
			{ 0, 0, 0, 0 }
		};

		while (true)
		{
			int key(0);

			do
				key = getopt_long(argc,argv,"",optionList,0);
			while (key == 0 || key == '?');

			if (key == -1) break;

			switch (key)
			{
				case 1: // --write-csv
					if (optarg == 0 || std::strlen(optarg) == 0)
						throw std::runtime_error("Error parsing 'write-csv' argument!");
					csvFile = optarg;
					break;
				case 2: // --verbose
					if (optarg == 0)
						verboseLevel = 2;
					else if (std::sscanf(optarg,"%d",&verboseLevel) != 1)
						throw std::runtime_error("Error parsing 'verbose' argument!");
					break;
				case 3: // --intial-volume
					if (optarg == 0 || std::sscanf(optarg,"%lf",&zeroVolumeMeasure) != 1)
						throw std::runtime_error("Error parsing 'initial-volume' argument!");
					break;
				case 4: // --volume-decay-exponent
					if (optarg == 0 || std::sscanf(optarg,"%lf",&volumeDecayExponent) != 1)
						throw std::runtime_error("Error parsing 'volume-decay-exponent' argument!");
					break;
				case 5: // --length
					if (optarg == 0 || std::sscanf(optarg,"%lf",&fourthLevel_cylinderHeight) != 1)
						throw std::runtime_error("Error parsing 'length' argument!");
					break;
				case 6: // --radius
					if (optarg == 0 || std::sscanf(optarg,"%lf",&fourthLevel_cylinderRadius) != 1)
						throw std::runtime_error("Error parsing 'radius' argument!");
					break;
				case 7: // --vacuum-marker
					if (optarg == 0 || std::sscanf(optarg,"%d",&final_vacuumMarker) != 1)
						throw std::runtime_error("Error parsing 'vacuum-marker' argument!");
					break;
				case 8: // --top-marker
					if (optarg == 0 || std::sscanf(optarg,"%d",&final_topMarker) != 1)
						throw std::runtime_error("Error parsing 'top-marker' argument!");
					break;
				case 9: // --bottom-marker
					if (optarg == 0 || std::sscanf(optarg,"%d",&final_bottomMarker) != 1)
						throw std::runtime_error("Error parsing 'bottom-marker' argument!");
					break;
				case 10: // --sides-marker
					if (optarg == 0 || std::sscanf(optarg,"%d",&final_sidesMarker) != 1)
						throw std::runtime_error("Error parsing 'sides-marker' argument!");
					break;
				case 11: // --create-config-template
					if (optarg == 0 || std::strlen(optarg) == 0) 
						throw std::runtime_error("Error parsing 'create-config-template' argument!");
					configFile = optarg;
					break;
				case 12: // --debug-mode
					meshOutput = true;
					edgeFilterOutput = true;
					verboseLevel = 99;
					break;
				default:
					std::string errorString("Unknown argument in command line: \"");
					errorString += optarg;
					errorString += "\"!";

					throw std::runtime_error(errorString);
			}
		}
		
		if (aboutFlag)
		{
			std::cout << "Compilation time-stamp: " << __DATE__ << ", " << __TIME__ << std::endl;
			std::cout << std::endl;
			std::cout << "MeshGen is a mesh generation utility for TAPSim - an atom probe (3DAP) data" << std::endl;
			std::cout << "simulation  program written by  Christian  Oberdorfer.  Its source  code is" << std::endl;
			std::cout << "licensed under the  terms of the  GNU  General Public License  Version 3 or" << std::endl;
			std::cout << "later ( GNU GPL v3+ ) as  published by  the Free Software Foundation.  This" << std::endl;
			std::cout << "program is distributed  in the hope that it will be useful, but WITHOUT ANY" << std::endl;
			std::cout << "WARRANTY;  without even the implied warranty  of MERCHANTABILITY or FITNESS" << std::endl;
			std::cout << "FOR A PARTICULAR PURPOSE. See the GPL for more details." << std::endl;

			std::cout << "You should have receeived a copy of the GNU GPL along with this program. If" << std::endl;
			std::cout << "not, look at 'http://www.gnu.org/licenses/gpl.txt'." << std::endl;
			std::cout << std::endl;
			std::cout << "More information about the program and  its application is available in the" << std::endl;
			std::cout << "internet, see 'https://imperia.uni-muenster.de/Physik.MP/Schmitz/tapsim'."   << std::endl;
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

			return 0;
		}

		if (optind+2 != argc || showHelp || skipIniFile)
		{
			// ^^^ checking the 'skipIniFile' parameter maybe redundant but it ensures that program execution
			// *** will only continue in the case any file based initialization has been carried out!
			
			std::cout << "Calling convention: ./meshgen [...] input-file output-file" << std::endl;
			std::cout << std::endl;

			std::cout << "Parameters:" << std::endl;
			std::cout << "\t--help (prints this output)" << std::endl;
			std::cout << "\t--about" << std::endl;
			std::cout << "\t--quiet, --verbose, --verbose=?" << std::endl;
			std::cout << "\t--write-ascii, --write-binary, --write-csv-subgrids, --write-csv=?" << std::endl;
			std::cout << "\t--initial-volume=?, --volume-decay-exponent=?" << std::endl;
			std::cout << "\t--length=?, --radius=?" << std::endl;
			std::cout << "\t--vacuum-marker=?, --top-marker=?, --bottom-marker=?, --sides-marker=?" << std::endl;
			std::cout << "\t--create-config-template=?" << std::endl;
			std::cout << "\t--debug-mode (enables additional output)" << std::endl;
			std::cout << std::endl;

			return 0;
		}

		inputParams.filename = argv[optind];
		outputFile = argv[optind+1];
	}
	
	// *** I. Process input data
	
	InputData input;

	{
		// *** read file

		readSample(inputParams.filename.c_str(),&input.mesh);

		{
			const int sizeBefore = input.mesh.numberofpoints;

			TetGen::tetgenio tmp;
			TetGen::tetrahedralize(const_cast<char*>("zQ"),&input.mesh,&tmp);

			if (sizeBefore != tmp.numberofpoints)
				throw std::runtime_error("Input contains duplicated points!");
			
			input.mesh.deinitialize();
			input.mesh = tmp;

			tmp.initialize();
		}

		trifaceExtents(input.mesh,&input.extentsXY,&input.extentsZ);

		// *** compute inflated mesh

		{
			std::set<int> trifaceNodes;
			for (int i = 0; i < input.mesh.numberoftrifaces*3; i++)
				trifaceNodes.insert(input.mesh.trifacelist[i]);

			TetGen::tetgenio tmp;

			tmp.firstnumber = 0;
			tmp.numberofpoints = trifaceNodes.size();
			tmp.pointlist = new REAL[tmp.numberofpoints*3];

			int index(0);
			for (std::set<int>::const_iterator i = trifaceNodes.begin(); i != trifaceNodes.end(); i++)
			{
				tmp.pointlist[index*3] = input.mesh.pointlist[(*i)*3];
				tmp.pointlist[index*3] *= inputParams.xyExpansion;

				tmp.pointlist[index*3+1] = input.mesh.pointlist[(*i)*3+1];
				tmp.pointlist[index*3+1] *= inputParams.xyExpansion;

				tmp.pointlist[index*3+2] = input.mesh.pointlist[(*i)*3+2];
				tmp.pointlist[index*3+2] *= inputParams.zExpansion;

				index++;
			}
			
			// ***

			TetGen::tetrahedralize(const_cast<char*>("zQ"),&tmp,&input.inflatedMesh);
		}

		trifaceExtents(input.inflatedMesh,&input.inflatedExtentsXY,&input.inflatedExtentsZ);

		// *** construct PLCs from the inflated mesh convex hull

		plc_from_triface_hull plc;

		plc.defineInnerContour(inputParams.plcContourRadius,inputParams.plcContourPhi,inputParams.plcContourHeight);
		plc.defineInnerMesh(inputParams.plcRecoveryRadius,inputParams.plcRecoveryPhi,inputParams.plcRecoveryHeight);

		plc.defineOutwardCylinder(input.inflatedExtentsXY*inputParams.cylinder_xyExpansion,input.extentsZ*inputParams.cylinder_zExpansion);
		plc.setOutwardCylinderResolution(inputParams.cylinder_resolution);

		plc.generateUnbounded(&input.unboundedPLC,input.inflatedMesh);
		plc.generateBounded(&input.boundedPLC,input.inflatedMesh);
	}
	
	if (verboseLevel > 0)
	{
		std::cout << "I) Input processed" << std::endl;
		std::cout << "\t-> filename: '" << inputParams.filename << "'" << std::endl;
		std::cout << "\t-> xy-extents: " << input.extentsXY << std::endl;
		std::cout << "\t-> z-extents: " << input.extentsZ << std::endl;

		{
			std::cout << "\t-> point marker statistics: " << std::endl;

			std::map<int,int> pointMarkers;
			for (int i = 0; i < input.mesh.numberofpoints; i++)
				pointMarkers[input.mesh.pointmarkerlist[i]]++;

			for (std::map<int,int>::const_iterator i = pointMarkers.begin(); i != pointMarkers.end(); i++)
				std::cout << "\t\t" << i->first << " => occurrence: " << i->second << std::endl;
		}

		std::cout << "\t-> xy-inflation-factor: " << inputParams.xyExpansion << std::endl;
		std::cout << "\t-> z-inflation-factor: " << inputParams.zExpansion << std::endl;

		std::cout << "\t-> inflated xy-extents: " << input.inflatedExtentsXY << std::endl;
		std::cout << "\t-> inflated z-extents: " << input.inflatedExtentsZ << std::endl;

		std::cout << "\t-> plc contour parameters (r/phi/z): ";
		std::cout << inputParams.plcContourRadius << "/";
		std::cout << inputParams.plcContourPhi << "/";
		std::cout << inputParams.plcContourHeight;
		std::cout << std::endl;

		std::cout << "\t-> plc recovery parameters (r/phi/z): ";
		std::cout << inputParams.plcRecoveryRadius << "/";
		std::cout << inputParams.plcRecoveryPhi << "/";
		std::cout << inputParams.plcRecoveryHeight;
		std::cout << std::endl;

		std::cout << "\t-> plc bounding cylinder xy-expansion: " << inputParams.cylinder_xyExpansion << std::endl;
		std::cout << "\t-> plc bounding cylinder z-expansion: " << inputParams.cylinder_zExpansion << std::endl;
		std::cout << "\t-> plc bounding cylinder resolution: " << inputParams.cylinder_resolution << std::endl;

		std::cout << "\t=> number of mesh points: " << input.mesh.numberofpoints << std::endl;
		std::cout << std::endl;
	}

	if (meshOutput)
	{
		writeCSV("raw_input.csv",input.mesh);
		writeCSV("raw_inflatedInput.csv",input.inflatedMesh);
		writeCSV("raw_inflatedInput_unboundedPLC.csv",input.unboundedPLC);
		writeCSV("raw_inflatedInput_boundedPLC.csv",input.boundedPLC);
		
		if (verboseLevel > 0)
		{
			std::cout << "Wrote input data to file 'raw_input.csv'!" << std::endl;
			std::cout << "Wrote inflated input data to file 'raw_inflatedInput.csv'!" << std::endl;
			std::cout << "Wrote inflated input unbounded plc data to file 'raw_inflatedInput_unboundedPLC.csv'!" << std::endl;
			std::cout << "Wrote inflated input bounded plc data to file 'raw_inflatedInput_boundedPLC.csv'!" << std::endl;
			std::cout << std::endl;
		}
	}
	
	// *** II. Create first level mesh
	// *** => mesh in between the inflated sample hull and the sample

	LevelData firstLevel;

	if (withFirstLevel)
	{
		std::ostringstream params;
		params << firstLevel_tetGenParams << firstLevel_volumeConstraints;
		
		TetGen::tetrahedralize(const_cast<char*>(params.str().c_str()),&input.unboundedPLC,&firstLevel.mesh);
		
		// *** remove points within the source domain
		
		TetGen::tetgenio tmp = firstLevel.mesh;
		firstLevel.mesh.initialize();

		int marker[3] = { 1, 0, 1 };
		markPoints(input.mesh,&tmp,marker,0);
		removePoints(tmp,&tmp,0,1);
		
		TetGen::tetrahedralize(const_cast<char*>("zQ"),&tmp,&firstLevel.mesh);

		// *** set extents

		trifaceExtents(firstLevel.mesh,&firstLevel.extentsXY,&firstLevel.extentsZ);
	}
	else
	{
		firstLevel.extentsXY = input.inflatedExtentsXY;
		firstLevel.extentsZ = input.inflatedExtentsZ;
	}

	if (verboseLevel > 0)
	{
		std::cout << "II) First level mesh constructed" << std::endl;

		if (verboseLevel > 1) std::cout << "\t-> invoked parameters for TetGen: \"" << firstLevel_tetGenParams << "\"" << std::endl;

		std::cout << "\t-> volume constraints: " << firstLevel_volumeConstraints << std::endl;
		std::cout << "\t-> xy-extents: " << firstLevel.extentsXY << std::endl;
		std::cout << "\t-> z-extents: " << firstLevel.extentsZ << std::endl;
		std::cout << "\t=> number of mesh points: " << firstLevel.mesh.numberofpoints << std::endl;
		std::cout << std::endl;
	}

	if (meshOutput)
	{
		writeCSV("raw_firstLevel.csv",firstLevel.mesh);

		if (verboseLevel > 1)
		{
			std::cout << "Wrote first level data to file 'raw_firstLevel.csv'!" << std::endl;
			std::cout << std::endl;
		}
	}

	// *** III. Create second level mesh
	// *** => mesh between the inflated sample hull and a somewhat bigger cylinder
	
	LevelData secondLevel;

	if (withSecondLevel)
	{
		TetGen::tetgenio tmp;
		TetGen::tetrahedralize(const_cast<char*>(secondLevel_tetGenParams.c_str()),&input.boundedPLC,&tmp);
		secondLevel.mesh = tmp;
		tmp.initialize();

		int lastCycleSize(-1);
		for (unsigned int i = 0; i < refineCycles; i++)
		{
			DistanceMap distObj(distance_radiusResolution,distance_phiResolution,distance_thetaResolution);
			distObj.setVolumeMeasure(zeroVolumeMeasure);
			distObj.setVolumeDecayExponent(volumeDecayExponent);
			distObj.setVolumeMaxConstraint(secondLevel_maxConstraints);
			distObj.setReferenceObject(input.inflatedMesh);
			
			REAL threshold = distObj.getVolumeConstraints(&secondLevel.mesh);
			threshold /= secondLevel.mesh.numberoftetrahedra;

			if (threshold <= refineThreshold)
			{
				if (verboseLevel > 2)
				{
					std::cout << "Threshold (";
					std::cout << std::setprecision(2) << std::fixed << (refineThreshold*1e2) << "%) reached!" << std::endl;
					std::cout << std::resetiosflags(std::ios_base::floatfield);
				}
				break;
			}
			tmp = secondLevel.mesh;
			secondLevel.mesh.initialize();

			TetGen::tetrahedralize(const_cast<char*>("zrqaQ"),&tmp,&secondLevel.mesh);

			tmp.deinitialize();
			tmp.initialize();

			if (verboseLevel > 2)
			{
				std::cout << "Mesh refinement: cycle #" << 1+i << ", " << secondLevel.mesh.numberofpoints << "  points ";
				std::cout << "(unconstrained tetrahedra: " << std::setprecision(2) << std::fixed << (threshold*1e2) << "%)" << std::endl;
				std::cout << std::resetiosflags(std::ios_base::floatfield);
			}

			if (lastCycleSize == secondLevel.mesh.numberofpoints)
			{
				if (verboseLevel > 2) std::cout << "Limit reached! Outstanding cycles are skipped." << std::endl;
				break;
			}

			lastCycleSize = secondLevel.mesh.numberofpoints;
		}

		if (verboseLevel > 2) std::cout << std::endl;
		
		// ***

		tmp = secondLevel.mesh;
		secondLevel.mesh.initialize();

		TetGen::tetrahedralize(const_cast<char*>("zQ"),&tmp,&secondLevel.mesh);

		trifaceExtents(secondLevel.mesh,&secondLevel.extentsXY,&secondLevel.extentsZ);
	}
	else
	{
		secondLevel.extentsXY = input.inflatedExtentsXY*inputParams.cylinder_xyExpansion;
		secondLevel.extentsZ = input.extentsZ*inputParams.cylinder_zExpansion;
	}

	if (verboseLevel > 0)
	{
		std::cout << "III) Second level mesh constructed" << std::endl;

		if (verboseLevel > 1)
		{
			std::cout << "\t-> invoked parameters for TetGen: \"" << firstLevel_tetGenParams << "\"" << std::endl;
			std::cout << "\t-> volume measure: " << zeroVolumeMeasure << std::endl;
			std::cout << "\t-> volume decay exponent: " << volumeDecayExponent << std::endl;
			std::cout << "\t-> volume constraint: " << secondLevel_maxConstraints << std::endl;
		}

		std::cout << "\t-> xy-extents: " << secondLevel.extentsXY << std::endl;
		std::cout << "\t-> z-extents: " << secondLevel.extentsZ << std::endl;
		std::cout << "\t=> number of mesh points: " << secondLevel.mesh.numberofpoints << std::endl;
		std::cout << std::endl;
	}

	if (meshOutput)
	{
		writeCSV("raw_secondLevel.csv",secondLevel.mesh);

		if (verboseLevel > 0)
		{
			std::cout << "Wrote second level data to file 'raw_secondLevel.csv'!" << std::endl;
			std::cout << std::endl;
		}
	}

	// *** IV. Create third level mesh
	// *** => mesh between the first inner cylinder and another even bigger cylinder

	LevelData thirdLevel;

	if (withThirdLevel)
	{
		plc_boundedCylinder boundedCylinder;
		boundedCylinder.defineInnerCylinder(secondLevel.extentsXY,secondLevel.extentsZ);
		boundedCylinder.defineOutwardCylinder(thirdLevel_cylinderRadius,thirdLevel_cylinderHeight);
		boundedCylinder.setInnerCylinderResolution(thirdLevel_innerCylinderResolution);
		boundedCylinder.setOutwardCylinderResolution(thirdLevel_outwardCylinderResolution);
		boundedCylinder.generate(&thirdLevel.mesh);
		
		TetGen::tetgenio tmp;
		TetGen::tetrahedralize(const_cast<char*>("zpqMQ"),&thirdLevel.mesh,&tmp);

		thirdLevel.mesh.deinitialize();
		thirdLevel.mesh = tmp;
		tmp.initialize();
		
		int lastCycleSize(-1);
		for (unsigned int i = 0; i < refineCycles; i++)
		{
			DistanceMap distObj(distance_radiusResolution,distance_phiResolution,distance_thetaResolution);
			distObj.setVolumeMeasure(zeroVolumeMeasure);
			distObj.setVolumeDecayExponent(volumeDecayExponent);
			distObj.setVolumeMaxConstraint(thirdLevel_maxConstraints);
			distObj.setReferenceObject(input.inflatedMesh);

			REAL threshold = distObj.getVolumeConstraints(&thirdLevel.mesh);
			threshold /= thirdLevel.mesh.numberoftetrahedra;

			if (threshold <= refineThreshold)
			{
				if (verboseLevel > 2)
				{
					std::cout << "Threshold (";
					std::cout << std::setprecision(2) << std::fixed << (refineThreshold*1e2) << "%) reached!" << std::endl;
					std::cout << std::resetiosflags(std::ios_base::floatfield);
				}
				break;
			}

			tmp = thirdLevel.mesh;
			thirdLevel.mesh.initialize();

			TetGen::tetrahedralize(const_cast<char*>("zrqaQ"),&tmp,&thirdLevel.mesh);

			tmp.deinitialize();
			tmp.initialize();

			if (verboseLevel > 2)
			{
				std::cout << "Mesh refinement: cycle #" << 1+i << ", " << thirdLevel.mesh.numberofpoints << "  points ";
				std::cout << "(unconstrained tetrahedra: " << std::setprecision(2) << std::fixed << (threshold*1e2) << "%)" << std::endl;
				std::cout << std::resetiosflags(std::ios_base::floatfield);
			}

			if (lastCycleSize == thirdLevel.mesh.numberofpoints)
			{
				if (verboseLevel > 2) std::cout << "Limit reached! Outstanding cycles are skipped." << std::endl;
				break;
			}

			lastCycleSize = thirdLevel.mesh.numberofpoints;
		}

		if (verboseLevel > 2) std::cout << std::endl;

		// ***

		tmp = thirdLevel.mesh;
		thirdLevel.mesh.initialize();

		TetGen::tetrahedralize(const_cast<char*>("zQ"),&tmp,&thirdLevel.mesh);

		trifaceExtents(thirdLevel.mesh,&thirdLevel.extentsXY,&thirdLevel.extentsZ);
	}
	else
	{
		thirdLevel.extentsXY = thirdLevel_cylinderRadius;
		thirdLevel.extentsZ = thirdLevel_cylinderHeight;
	}

	if (verboseLevel > 0)
	{
		std::cout << "IV) Third level mesh constructed " << std::endl;

		if (verboseLevel > 1)
		{
			std::cout << "\t-> plc inner cylinder radius: " << secondLevel.extentsXY << std::endl;
			std::cout << "\t-> plc inner cylinder height: " << secondLevel.extentsZ << std::endl;
			std::cout << "\t-> plc inner cylinder resolution: " << thirdLevel_innerCylinderResolution << std::endl;
			std::cout << "\t-> plc outward cylinder radius: " << thirdLevel_cylinderRadius << std::endl;
			std::cout << "\t-> plc outward cylinder height: " << thirdLevel_cylinderHeight << std::endl;
			std::cout << "\t-> plc outward cylinder resolution: " << thirdLevel_outwardCylinderResolution << std::endl;
			std::cout << "\t-> volume measure: " << zeroVolumeMeasure << std::endl;
			std::cout << "\t-> volume decay exponent: " << volumeDecayExponent << std::endl;
			std::cout << "\t-> volume constraint: " << thirdLevel_maxConstraints << std::endl;
		}

		std::cout << "\t-> xy-extents: " << thirdLevel.extentsXY << std::endl;
		std::cout << "\t-> z-extents: " << thirdLevel.extentsZ << std::endl;
		std::cout << "\t=> number of mesh points: " << thirdLevel.mesh.numberofpoints << std::endl;
		std::cout << std::endl;
	}

	if (meshOutput)
	{
		writeCSV("raw_thirdLevel.csv",thirdLevel.mesh);

		if (verboseLevel > 0)
		{
			std::cout << "Wrote third level data to file 'raw_thirdLevel.csv'!" << std::endl;
			std::cout << std::endl;
		}
	}

	// *** V. Create fourth level mesh
	// *** => mesh between the foremost inner cylinder and yet another even more bigger cylinder
	
	LevelData fourthLevel;
	
	if (withFourthLevel)
	{
		plc_boundedCylinder boundedCylinder;
		boundedCylinder.defineInnerCylinder(thirdLevel.extentsXY,thirdLevel.extentsZ);
		boundedCylinder.defineOutwardCylinder(fourthLevel_cylinderRadius,fourthLevel_cylinderHeight);
		boundedCylinder.setInnerCylinderResolution(fourthLevel_innerCylinderResolution);
		boundedCylinder.setOutwardCylinderResolution(fourthLevel_outwardCylinderResolution);
		boundedCylinder.generate(&fourthLevel.mesh);
		
		TetGen::tetgenio tmp;
		TetGen::tetrahedralize(const_cast<char*>("zpqMQ"),&fourthLevel.mesh,&tmp);

		fourthLevel.mesh.deinitialize();
		fourthLevel.mesh = tmp;
		tmp.initialize();

		int lastCycleSize(-1);
		for (unsigned int i = 0; i < refineCycles; i++)
		{
			DistanceMap distObj(distance_radiusResolution,distance_phiResolution,distance_thetaResolution);
			distObj.setVolumeMeasure(zeroVolumeMeasure);
			distObj.setVolumeDecayExponent(volumeDecayExponent);
			distObj.setVolumeMaxConstraint(fourthLevel_maxConstraints);
			distObj.setReferenceObject(input.inflatedMesh);

			REAL threshold = distObj.getVolumeConstraints(&fourthLevel.mesh);
			threshold /= fourthLevel.mesh.numberoftetrahedra;

			if (threshold <= refineThreshold)
			{
				if (verboseLevel > 2)
				{
					std::cout << "Threshold (";
					std::cout << std::setprecision(2) << std::fixed << (refineThreshold*1e2) << "%) reached!" << std::endl;
					std::cout << std::resetiosflags(std::ios_base::floatfield);
				}
				break;
			}

			tmp = fourthLevel.mesh;
			fourthLevel.mesh.initialize();

			TetGen::tetrahedralize(const_cast<char*>("zrqaQ"),&tmp,&fourthLevel.mesh);

			tmp.deinitialize();
			tmp.initialize();

			if (verboseLevel > 2)
			{
				std::cout << "Mesh refinement: cycle #" << 1+i << ", " << fourthLevel.mesh.numberofpoints << "  points ";
				std::cout << "(unconstrained tetrahedra: " << std::setprecision(2) << std::fixed << (threshold*1e2) << "%)" << std::endl;
				std::cout << std::resetiosflags(std::ios_base::floatfield);
			}

			if (lastCycleSize == fourthLevel.mesh.numberofpoints)
			{
				if (verboseLevel > 2) std::cout << "Limit reached! Outstanding cycles are skipped." << std::endl;
				break;
			}

			lastCycleSize = fourthLevel.mesh.numberofpoints;
		}

		if (verboseLevel > 2) std::cout << std::endl;

		// ***

		tmp = fourthLevel.mesh;
		fourthLevel.mesh.initialize();

		TetGen::tetrahedralize(const_cast<char*>("zQ"),&tmp,&fourthLevel.mesh);
		
		trifaceExtents(fourthLevel.mesh,&fourthLevel.extentsXY,&fourthLevel.extentsZ);
	}

	if (verboseLevel > 0)
	{
		std::cout << "V) Fourth level mesh constructed " << std::endl;

		if (verboseLevel > 1)
		{
			std::cout << "\t-> plc inner cylinder radius: " << thirdLevel.extentsXY << std::endl;
			std::cout << "\t-> plc inner cylinder height: " << thirdLevel.extentsZ << std::endl;
			std::cout << "\t-> plc inner cylinder resolution: " << fourthLevel_innerCylinderResolution << std::endl;
			std::cout << "\t-> plc outward cylinder radius: " << fourthLevel_cylinderRadius << std::endl;
			std::cout << "\t-> plc outward cylinder height: " << fourthLevel_cylinderHeight << std::endl;
			std::cout << "\t-> plc outward cylinder resolution: " << fourthLevel_outwardCylinderResolution << std::endl;
			std::cout << "\t-> volume measure: " << zeroVolumeMeasure << std::endl;
			std::cout << "\t-> volume decay exponent: " << volumeDecayExponent << std::endl;
			std::cout << "\t-> volume constraint: " << fourthLevel_maxConstraints << std::endl;
		}

		std::cout << "\t-> xy-extents: " << fourthLevel.extentsXY << std::endl;
		std::cout << "\t-> z-extents: " << fourthLevel.extentsZ << std::endl;
		std::cout << "\t=> number of mesh points: " << fourthLevel.mesh.numberofpoints << std::endl;
		std::cout << std::endl;
	}

	if (meshOutput)
	{
		writeCSV("raw_fourthLevel.csv",fourthLevel.mesh);

		if (verboseLevel > 0)
		{
			std::cout << "Wrote third level data to file 'raw_fourthLevel.csv'!" << std::endl;
			std::cout << std::endl;
		}
	}

	// *** VI. Merges all meshes

	if (verboseLevel > 0)
	{
		std::cout << "VI) Post processing" << std::endl;
		std::cout << std::endl;
	}
	
	LevelData final;

	{
		// remove coinciding mesh points

		{
			int marker[3] = { 1, 0, 1 };

			if (withThirdLevel) 
				markPoints(thirdLevel.mesh,&fourthLevel.mesh,marker,0);
			else
				setPointMarker(&fourthLevel.mesh,0);
			
			if (withSecondLevel)
				markPoints(secondLevel.mesh,&thirdLevel.mesh,marker,0);
			else
				setPointMarker(&thirdLevel.mesh,0);
			
			if (withFirstLevel)
				markPoints(firstLevel.mesh,&secondLevel.mesh,marker,0);
			else
				setPointMarker(&secondLevel.mesh,0);
			
			markPoints(input.mesh,&firstLevel.mesh,marker,0);

			// ***
			
			int numPoints_1;
			int numPoints_2;
			int numPoints_3;
			int numPoints_4;
			
			if (withFirstLevel)
				numPoints_1 = removePoints(firstLevel.mesh,&firstLevel.mesh,0,1);
			else
				numPoints_1 = 0;
			
			if (withSecondLevel)
				numPoints_2 = removePoints(secondLevel.mesh,&secondLevel.mesh,0,1);
			else
				numPoints_2 = 0;
			
			if (withThirdLevel)
				numPoints_3 = removePoints(thirdLevel.mesh,&thirdLevel.mesh,0,1);
			else
				numPoints_3 = 0;
			
			if (withFourthLevel)
				numPoints_4 = removePoints(fourthLevel.mesh,&fourthLevel.mesh,0,1);
			else
				numPoints_4 = 0;

			int totalNum(0);
			totalNum += numPoints_4;
			totalNum += numPoints_3;
			totalNum += numPoints_2;
			totalNum += numPoints_1;
			
			if (verboseLevel > 1)
			{
				std::cout << "Removed coinciding points" << std::endl;
				std::cout << "\t-> first level mesh: " << numPoints_1 << std::endl;
				std::cout << "\t-> second level mesh: " << numPoints_2 << std::endl;
				std::cout << "\t-> third level mesh: " << numPoints_3 << std::endl;
				std::cout << "\t-> fourth level mesh: " << numPoints_4 << std::endl;
				std::cout << std::endl;

				std::cout << "Total number of removed coinciding points: " << totalNum << std::endl;
				std::cout << std::endl;
			}

		}
		
		// set pending points attributes (equivalent to tapsim - invalid - point numbers)

		if (withFirstLevel)
		{
			firstLevel.mesh.numberofpointattributes = 1;
			firstLevel.mesh.pointattributelist = new REAL[firstLevel.mesh.numberofpointattributes*firstLevel.mesh.numberofpoints];
			for (int i = 0; i < firstLevel.mesh.numberofpoints*firstLevel.mesh.numberofpointattributes; i++)
				firstLevel.mesh.pointattributelist[i] = static_cast<REAL>(0);
		}

		if (withSecondLevel)
		{
			secondLevel.mesh.numberofpointattributes = 1;
			secondLevel.mesh.pointattributelist = new REAL[secondLevel.mesh.numberofpointattributes*secondLevel.mesh.numberofpoints];
			for (int i = 0; i < secondLevel.mesh.numberofpoints*secondLevel.mesh.numberofpointattributes; i++)
				secondLevel.mesh.pointattributelist[i] = static_cast<REAL>(0);
		}

		if (withThirdLevel)
		{
			thirdLevel.mesh.numberofpointattributes = 1;
			thirdLevel.mesh.pointattributelist = new REAL[thirdLevel.mesh.numberofpointattributes*thirdLevel.mesh.numberofpoints];
			for (int i = 0; i < thirdLevel.mesh.numberofpoints*thirdLevel.mesh.numberofpointattributes; i++)
				thirdLevel.mesh.pointattributelist[i] = static_cast<REAL>(0);
		}

		if (withFourthLevel)
		{
			fourthLevel.mesh.numberofpointattributes = 1;
			fourthLevel.mesh.pointattributelist = new REAL[fourthLevel.mesh.numberofpointattributes*fourthLevel.mesh.numberofpoints];
			for (int i = 0; i < fourthLevel.mesh.numberofpoints*fourthLevel.mesh.numberofpointattributes; i++)
				fourthLevel.mesh.pointattributelist[i] = static_cast<REAL>(0);
		}

		if (verboseLevel > 3)
		{
			std::cout << "Set point attributes ('point numbers')." << std::endl;
			std::cout << std::endl;
		}

		// ***
		
		if (writeSubgrids)
		{
			setPointMarker(&firstLevel.mesh,1);
			writeCSV("firstLevel.csv",firstLevel.mesh);

			if (verboseLevel > 1) std::cout << "Wrote final first level data to file 'firstLevel.csv'!" << std::endl;

			// ***
			setPointMarker(&secondLevel.mesh,2);
			writeCSV("secondLevel.csv",secondLevel.mesh);

			if (verboseLevel > 1) std::cout << "Wrote final second level data to file 'secondLevel.csv'!" << std::endl;

			// ***

			setPointMarker(&thirdLevel.mesh,3);
			writeCSV("thirdLevel.csv",thirdLevel.mesh);

			if (verboseLevel > 1) std::cout << "Wrote final third level data to file 'thirdLevel.csv'!" << std::endl;

			// ***

			setPointMarker(&fourthLevel.mesh,4);
			writeCSV("fourthLevel.csv",fourthLevel.mesh);

			if (verboseLevel > 1) std::cout << "Wrote final fourth level data to file 'fourthLevel.csv'!" << std::endl;

			if (verboseLevel > 1) std::cout << std::endl;
		}
		
		// merge exterior meshes

		if (withFirstLevel)
		{
			addPoints(&final.mesh,firstLevel.mesh);
			firstLevel.mesh.deinitialize();
			firstLevel.mesh.initialize();
		}

		if (withSecondLevel)
		{
			addPoints(&final.mesh,secondLevel.mesh);
			secondLevel.mesh.deinitialize();
			secondLevel.mesh.initialize();
		}

		if (withThirdLevel)
		{
			addPoints(&final.mesh,thirdLevel.mesh);
			thirdLevel.mesh.deinitialize();
			thirdLevel.mesh.initialize();
		}

		if (withFourthLevel)
		{
			addPoints(&final.mesh,fourthLevel.mesh);
			fourthLevel.mesh.deinitialize();
			fourthLevel.mesh.initialize();
		}

		if (verboseLevel > 2)
		{
			std::cout << "Merged first to fourth level meshes." << std::endl;
			std::cout << std::endl;
		}

		// *** do explicit conversion to single floating point precision

		if (forceSinglePrecision)
		{
			for (int i = 0; i <final.mesh.numberofpoints*3; i++)
				final.mesh.pointlist[i] = static_cast<float>(final.mesh.pointlist[i]);

			if (verboseLevel > 2)
			{
				std::cout << "Force casting of point coordinates to single precision." << std::endl;
				std::cout << std::endl;
			}
		}

		// *** do consistency check

		{
			const int sizeBefore = final.mesh.numberofpoints;

			TetGen::tetgenio tmp;
			TetGen::tetrahedralize(const_cast<char*>("zCQ"),&final.mesh,&tmp); 

			if (verboseLevel > 0 && sizeBefore != tmp.numberofpoints)
			{
				std::cout << "!!! ATTENTION" << std::endl;
				std::cout << "!!! =========" << std::endl;
				std::cout << "!!!" << std::endl;
				std::cout << "!!! Due to the limited precision of single floating point arithmetic" << std::endl;
				std::cout << "!!! redundant points have been removed!" << std::endl;
				std::cout << "!!!" << std::endl;
				std::cout << "!!!\t=> size before removal: " << sizeBefore << std::endl;
				std::cout << "!!!\t=> current size: " << final.mesh.numberofpoints << std::endl;
				std::cout << "!!!\t=> removed points: " << sizeBefore - final.mesh.numberofpoints << std::endl;
				std::cout << std::endl;
			}
			
			final.mesh.deinitialize();
			final.mesh = tmp;
			tmp.initialize();
		}

		// *** apply edge length filter and do mesh coarsening

		if (edgeFilter_length > 0.0)
		{
			int removeCntSum(0);
			for (unsigned int i = 0; i < edgeFilter_cycles; i++)
			{
				TetGen::tetgenio tmp;
				TetGen::tetrahedralize(const_cast<char*>("zQ"),&final.mesh,&tmp);

				final.mesh.deinitialize();
				final.mesh.initialize();

				final.mesh = tmp;
				tmp.initialize();

				// ***

				TetGen::tetgenio removedPoints;
				const int removeCnt = filterEdges(&final.mesh,edgeFilter_length,&removedPoints);

				if (verboseLevel > 2)
				{
					std::cout << "Removing short edges: cycle #" << i+1;
					std::cout << ", " << removeCnt << " points removed" << std::endl;
				}

				removeCntSum += removeCnt;

				if (removeCnt == 0)
				{
					if (verboseLevel > 2) std::cout << "No unconstrained edges. Outstanding cycles are skipped." << std::endl;
					break;
				}

				if (edgeFilterOutput)
				{
					std::ostringstream filename;
					filename << "removedPoints_";
					filename << i;
					filename << ".csv";
					
					writeCSV(filename.str().c_str(),removedPoints);

					if (verboseLevel > 2) std::cout << "Wrote removed points to file '" << filename.str() << "'." << std::endl;
				}
			}

			if (verboseLevel > 2) std::cout << std::endl;

			if (verboseLevel > 1)
			{
				std::cout << "Edge length filter " << std::endl;
				std::cout << "\t-> length constraint: " << edgeFilter_length << std::endl;
				std::cout << "\t-> total number of removed points: " << removeCntSum << std::endl;
				std::cout << std::endl;
			}
		}

		// *** set final extents

		{
			TetGen::tetgenio tmp;
			TetGen::tetrahedralize(const_cast<char*>("zQ"),&final.mesh,&tmp);

			final.mesh.deinitialize();
			final.mesh = tmp;
			tmp.initialize();

			trifaceExtents(final.mesh,&final.extentsXY,&final.extentsZ);
		}

		if (verboseLevel > 1)
		{
			std::cout << "Final extents: xy = " << final.extentsXY << ", z = " << final.extentsZ << std::endl;
			std::cout << std::endl;
		}

		// set final point markers

		setPointMarker(&final.mesh,final_vacuumMarker);
		
		for (int i = 0; i < final.mesh.numberoftrifaces*3; i++)
		{
			const int nodeIndex = final.mesh.trifacelist[i];
			
			if (final.mesh.pointlist[nodeIndex*3+2] == 0.0)
				final.mesh.pointmarkerlist[nodeIndex] = final_bottomMarker;
			else if (final.mesh.pointlist[nodeIndex*3+2] == final.extentsZ)
				final.mesh.pointmarkerlist[nodeIndex] = final_topMarker;
			else
				final.mesh.pointmarkerlist[nodeIndex] = final_sidesMarker;
		}

		for (int i = 0; i < final.mesh.numberoftrifaces; i++)
		{
			int bottomCnt(0);
			int topCnt(0);
			
			for (int j = 0; j < 3; j++)
			{
				const int nodeIndex = final.mesh.trifacelist[i*3+j];

				if (final.mesh.pointlist[nodeIndex*3+2] == 0.0)
					bottomCnt++;
				else if (final.mesh.pointlist[nodeIndex*3+2] == final.extentsZ)
					topCnt++;
			}

			
			if ((bottomCnt > 0 && bottomCnt < 3) || (topCnt > 0 && topCnt < 3))
			{
				for (int j = 0; j < 3; j++)
					final.mesh.pointmarkerlist[final.mesh.trifacelist[i*3+j]] = final_sidesMarker;
			}
		}

		if (verboseLevel > 0)
		{
			std::cout << "Final point markers are set." << std::endl;
			std::cout << std::endl;
			std::cout << "Used marker values: " << std::endl;
			std::cout << "\t-> vacuum: " << final_vacuumMarker << std::endl;
			std::cout << "\t-> bottom: " << final_bottomMarker << std::endl;
			std::cout << "\t-> top: " << final_topMarker << std::endl;
			std::cout << "\t-> sides: " << final_sidesMarker << std::endl;
			std::cout << std::endl;
		}
		
		// *** finally merge with user provided input mesh

		addPoints(&final.mesh,input.mesh);

		if (verboseLevel > 1)
		{
			std::cout << "Merged generated mesh with input." << std::endl;
			std::cout << std::endl;
		}
	}

	// *** VII. write final mesh to file

	if (verboseLevel > 0)
	{
		std::cout << "VII) Finishing" << std::endl;
		std::cout << std::endl;

		std::map<int,int> pointMarkers;
		for (int i = 0; i < final.mesh.numberofpoints; i++)
			pointMarkers[final.mesh.pointmarkerlist[i]]++;

		std::cout << "Point marker statistics: " << std::endl;
		for (std::map<int,int>::const_iterator i = pointMarkers.begin(); i != pointMarkers.end(); i++)
			std::cout << '\t' << i->first << " => occurrence: " << i->second << std::endl;

		std::cout << std::endl;

		std::cout << "Total number of points: " << final.mesh.numberofpoints << " points." << std::endl;
		std::cout << std::endl;
	}

	writeTapSim_Nodefile(outputFile.c_str(),final.mesh,final_outputMode);

	if (verboseLevel > 0)
	{
		std::cout << "Wrote TapSim nodefile!" << std::endl;
		std::cout << "Filename: '" << outputFile.c_str() << "'." << std::endl;
		std::cout << std::endl;
	}

	if (!csvFile.empty())
	{
		writeCSV(csvFile.c_str(),final.mesh);

		if (verboseLevel > 0)
		{
			std::cout << "Wrote final data to csv-file '" << csvFile << "'!" << std::endl;
			std::cout << std::endl;
		}
	}
	
	if (!configFile.empty())
	{
		ConfigTemplate_Params params;
		{
			params.bottomId = final_bottomMarker;
			params.sidesId = final_sidesMarker;
			params.topId = final_topMarker;
			params.vacuumId = final_vacuumMarker;

			for (int i = 0; i < final.mesh.numberofpoints; i++)
				params.userIds.insert(final.mesh.pointmarkerlist[i]);
			
			params.userIds.erase(params.bottomId);
			params.userIds.erase(params.sidesId);
			params.userIds.erase(params.topId);
			params.userIds.erase(params.vacuumId);
		}
		
		writeConfig(configFile.c_str(),params);
		
		if (verboseLevel > 0)
		{
			std::cout << "Wrote configuration template to file '" << configFile << "'!" << std::endl;
			std::cout << std::endl;
		}
	}

	return 0;
}
