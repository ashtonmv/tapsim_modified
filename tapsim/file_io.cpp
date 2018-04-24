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

#include "file_io.h"
#include "file_util.h"

#include <list>
#include <string>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <cstring>

#include <iostream>

#include "debug.h"

namespace
{	
	const char* configKeys[][2] =
	{
		{ "TEMPERATURE", "[K]" },
		{ "ID", "[1]" },
		{ "NAME", "[1]", },
		{ "CHARGE_DENSITY", "[C/m³]" },
		{ "DIELECTRICITY", "[1]" },
		{ "REMOVABLE", "[1]" },
		{ "NEUMANN_BOUNDARY", "[1]" },
		{ "DIRICHLET_BOUNDARY", "[1]" },
		{ "POTENTIAL", "[V]" },
		{ "MASS", "[u]" },
		{ "EVAPORATION_CHARGE_STATE", "[1]" },
		{ "EVAPORATION_FIELD_STRENGTH", "[V/m]" },
		{ "EVAPORATION_ACTIVATION_ENERGY", "[eV]" }
	};
}

// ***** ===== FILE_IO::READ_INITIALIZATION() ===== *****

void File_Io::readInitialization(const char* filename, std::map<std::string,std::string>* values, const char* delimiter)
{	
	if (values == 0) throw std::runtime_error("File_Io::readInitialization()");

	values->clear();

	// ***

	std::ifstream file(filename,std::ifstream::in);
	if (!file.good()) throw std::runtime_error("File_Io::readInitialization(): cannot open file!");


	while (file.good())
	{
		std::string buffer;

		while (file.good() && (buffer.empty() || buffer.at(0) == '#'))
			std::getline(file,buffer);

		if (!file.good()) break;

		std::pair<std::string,std::string> dataEntry;

		const size_t pos = buffer.find(delimiter);
		if (pos == std::string::npos) throw std::runtime_error("File_Io::readInitialization()");

		dataEntry.first = File_Utils::trim(buffer.substr(0,pos));
		dataEntry.second = File_Utils::trim(buffer.substr(pos+std::strlen(delimiter)));
		
		const std::pair<std::map<std::string,std::string>::iterator,bool> status = values->insert(dataEntry);
		if (!status.second) throw std::runtime_error("File_Io::readInitialization()");
	}

	file.close();
}

// ***** ===== FILE_IO::WRITE_INITIALIZATION() ===== *****

void File_Io::writeInitialization(const char* filename, const std::list< File_Io::KeyValue >& values, const char* delimiter)
{
	std::ofstream file(filename,std::ifstream::out);
	if (!file.good()) throw std::runtime_error("File_Io::writeInitialization()");

	for (std::list<KeyValue>::const_iterator i = values.begin(); i != values.end(); i++)
	{
		if (i->first.empty())
			file << '\n';
		else if (i->first.at(0) == '#')
			file << File_Utils::trim(i->first) << '\n';
		else
			file << File_Utils::trim(i->first) << delimiter << File_Utils::trim(i->second) << '\n';
	}

	file.close();
}

// ***** ===== FILE_IO::WRITE_CONFIG() ===== *****

void File_Io::writeConfig(const char* filename, const Configuration::Table& configTable, const int ioMode)
{
	if (ioMode == ASCII)
	{
		const char* spacer = " = ";
		const char endline = '\n';
	
		std::ofstream file(filename,std::ofstream::out|std::ofstream::trunc);
		if (!file.good()) throw std::runtime_error("File_Io::writeConfig()");

		file.precision(8);

		// *** header

		file << "# global parameters" << endline;
		file << "#" << endline;
		
		file << "# " << configKeys[0][0] << " ";
		file << configKeys[0][1] << endline;
		file << "#" << endline;

		file << "# species dependant parameters" << endline;
		file << "#" << endline;
		
		for (int i = 1; i < 13; i++)
		{
			file << "# ";
			file << configKeys[i][0] << " ";
			file << configKeys[i][1] << endline;
		}

		file << "ASCII" << endline;
		file << endline;

		// *** data section

		file << configKeys[0][0] << spacer << configTable.temperature() << endline;
		file << endline;
		
		std::list<Configuration::NodeId> ids = configTable.ids();
		for (std::list<Configuration::NodeId>::const_iterator i = ids.begin(); i != ids.end(); i++)
		{
			const Configuration::NodeData& data = configTable[*i];

			file << configKeys[1][0] << spacer << data.id().toValue() << endline;
			file << configKeys[2][0] << spacer << data.name() << endline;
			file << configKeys[3][0] << spacer << data.chargeDensity() << endline;
			file << configKeys[4][0] << spacer << data.epsilon() << endline;
			file << configKeys[5][0] << spacer << data.removable() << endline;
			file << configKeys[6][0] << spacer << data.neumannBoundary() << endline;
			file << configKeys[7][0] << spacer << data.dirichletBoundary() << endline;
			file << configKeys[8][0] << spacer << data.phi() << endline;
			file << configKeys[9][0] << spacer << data.mass() << endline;
			file << configKeys[10][0] << spacer << data.evapCharge() << endline;
			file << configKeys[11][0] << spacer << data.evapField() << endline;
			file << configKeys[12][0] << spacer << data.evapEnergy() << endline;
			file << endline;
		}

		file.close();
	}
	else if (ioMode == BINARY)
	{
		std::ofstream file(filename,std::ofstream::out|std::ofstream::trunc);
		if (!file.good()) throw std::runtime_error("File_Io::writeConfig()");

		file << "BINARY" << '\n';

		File_Utils::binWrite<float>(file,configTable.temperature());
		
		File_Utils::binWrite<int>(file,configTable.data().size());

		for (std::map<Configuration::NodeId,Configuration::NodeData>::const_iterator i = configTable.data().begin(); i != configTable.data().end(); i++)
		{
			File_Utils::binWrite<Configuration::NodeId>(file,i->second.id());
			File_Utils::binWrite<std::string>(file,i->second.name());
			File_Utils::binWrite<float>(file,i->second.chargeDensity());
			File_Utils::binWrite<float>(file,i->second.epsilon());
			File_Utils::binWrite<bool>(file,i->second.removable());
			File_Utils::binWrite<bool>(file,i->second.neumannBoundary());
			File_Utils::binWrite<bool>(file,i->second.dirichletBoundary());
			File_Utils::binWrite<float>(file,i->second.phi());
			File_Utils::binWrite<float>(file,i->second.mass());
			File_Utils::binWrite<short>(file,i->second.evapCharge());
			File_Utils::binWrite<float>(file,i->second.evapField());
			File_Utils::binWrite<float>(file,i->second.evapEnergy());
		}

		file.close();
	}
	else
		throw std::runtime_error("File_Io::writeConfig()");
}

// ***** ===== FILE_IO::READ_CONFIG() ===== *****

void File_Io::readConfig(const char* filename, Configuration::Table* configTable)
{
	if (configTable == 0) throw std::runtime_error("File_Io::readConfig()");

	configTable->clear();

	// ***

	std::ifstream file(filename,std::ifstream::in);
	if (!file.good()) throw std::runtime_error("File_Io::readConfig()");

	std::string firstLine;
	std::getline(file,firstLine);

	while (firstLine.empty() || firstLine.at(0) == '#')
		std::getline(file,firstLine);

	bool asciiMode;

	if (firstLine == "ASCII")
		asciiMode = true;
	else
	{
		if (firstLine == "BINARY")
 			asciiMode = false;
		else
			throw std::runtime_error("File_Io::readConfig()");
	}

	// ***

	if (asciiMode)
	{
		configTable->setTemperature(0.0); // initialize default value
		
		// ***
		
		const char* delimiter = " = ";
		
		while (file.good())
		{
			std::list<KeyValue> dataset;

			{
				std::map<std::string,std::string> uniqueDataset;

				std::string buffer;
				while (file.good() && (buffer.empty() || buffer.at(0) == '#'))
					std::getline(file,buffer);

				while (file.good() && !buffer.empty())
				{
					unsigned int pos = buffer.find(delimiter);
					if (pos == std::string::npos) throw std::runtime_error("File_Io::readConfig()");

					KeyValue dataObj;
					dataObj.first = File_Utils::trim(buffer.substr(0,pos));
					dataObj.second = File_Utils::trim(buffer.substr(pos+3));

					if (!uniqueDataset.insert(dataObj).second) throw std::runtime_error("File_Io::readConfig()");

					do
						std::getline(file,buffer);
					while (!buffer.empty() && buffer.at(0) == '#');
				}

				for (std::map<std::string,std::string>::const_iterator i = uniqueDataset.begin(); i != uniqueDataset.end(); i++)
					dataset.push_back(*i);
			}
			
			
			// ***
			
			if (dataset.size() == 0)
				continue;
			else if (dataset.size() == 1)
			{
				// temperature
				if (dataset.front().first != configKeys[0][0])
					throw std::runtime_error("File_Io::readConfig()");

				const float temperature = File_Utils::read_element<float>(dataset.front().second);
				configTable->setTemperature(temperature);
			}
			else if (dataset.size() == 12)
			{
				Configuration::NodeData data;
				
				while (!dataset.empty())
				{
					const std::string& key = dataset.front().first;
					const std::string& value = dataset.front().second;
					
					int index(1);
					while (index < 13 && key != configKeys[index][0])
						index++;

					if (index == 0 || index >= 13)  throw std::runtime_error("File_Io::readConfig()");

					switch (index)
					{
						case 1: // id
						{
							const short id = File_Utils::read_element<short>(value);
							data.setId(Configuration::NodeId::fromValue(id));
							break;
						}
						case 2: // name
						{
							data.setName(value);
							break;
						}
						case 3: // charge density
						{
							const float chargeDensity = File_Utils::read_element<float>(value);
							data.setChargeDensity(chargeDensity);
							break;
						}
						case 4: // dielectricity
						{
							const float epsilon = File_Utils::read_element<float>(value);
							data.setEpsilon(epsilon);
							break;
						}
						case 5: // removable
						{
							const bool removable = File_Utils::read_element<bool>(value);
							data.setRemovable(removable);
							break;
						}
						case 6: // neumann boundary
						{
							const bool neumannBoundary = File_Utils::read_element<bool>(value);
							data.setNeumannBoundary(neumannBoundary);
							break;
						}
						case 7: // dirichlet boundary
						{
							const bool dirichletBoundary = File_Utils::read_element<bool>(value);
							data.setDirichletBoundary(dirichletBoundary);
							break;
						}
						case 8: // potential
						{
							const float phi = File_Utils::read_element<float>(value);
							data.setPhi(phi);
							break;
						}
						case 9: // mass
						{
							const float mass = File_Utils::read_element<float>(value);
							data.setMass(mass);
							break;
						}
						case 10: // evaporation charge state
						{
							const short evapCharge = File_Utils::read_element<short>(value);
							data.setEvapCharge(evapCharge);
							break;
						}
						case 11: // evaporation field strength
						{
							const float evapField = File_Utils::read_element<float>(value);
							data.setEvapField(evapField);
							break;
						}
						case 12: // evaporation activiation energy
						{
							const float evapEnergy = File_Utils::read_element<float>(value);
							data.setEvapEnergy(evapEnergy);
							break;
						}
					}

					dataset.pop_front();
				}

				configTable->insert(data.id(),data);
			}
			else
				throw std::runtime_error("File_Io::readConfig()");
		}
	}
	else
	{
		float temperature;
		File_Utils::binRead<float>(file,&temperature);
		
		int size;
		File_Utils::binRead<int>(file,&size);

		for (int i = 0; i < size; i++)
		{
			Configuration::NodeData data;
			
			Configuration::NodeId id;
			File_Utils::binRead<Configuration::NodeId>(file,&id);
			data.setId(id);
			
			std::string name;
			File_Utils::binRead<std::string>(file,&name);
			data.setName(name);
			
			float chargeDensity;
			File_Utils::binRead<float>(file,&chargeDensity);
			data.setChargeDensity(chargeDensity);
			
			float epsilon;
			File_Utils::binRead<float>(file,&epsilon);
			data.setEpsilon(epsilon);
			
			bool removable;
			File_Utils::binRead<bool>(file,&removable);
			data.setRemovable(removable);
			
			bool neumannBoundary;
			File_Utils::binRead<bool>(file,&neumannBoundary);
			data.setNeumannBoundary(neumannBoundary);
			
			bool dirichletBoundary;
			File_Utils::binRead<bool>(file,&dirichletBoundary);
			data.setDirichletBoundary(dirichletBoundary);
			
			float phi;
			File_Utils::binRead<float>(file,&phi);
			data.setPhi(phi);
			
			float mass;
			File_Utils::binRead<float>(file,&mass);
			data.setMass(mass);
			
			short evapCharge;
			File_Utils::binRead<short>(file,&evapCharge);
			data.setEvapCharge(evapCharge);
			
			float evapField;
			File_Utils::binRead<float>(file,&evapField);
			data.setEvapField(evapField);
			
			float evapEnergy;
			File_Utils::binRead<float>(file,&evapEnergy);
			data.setEvapEnergy(evapEnergy);
			
			configTable->insert(data.id(),data);
		}
	}

	file.close();
}

// ***** ===== FILE_IO::WRITE_GEOMETRY() ===== *****

void File_Io::writeGeometry(const char* filename, const Geometry_3d::Table& geomTable, const unsigned char contentFlags, const int ioMode)
{
	std::ofstream file(filename,std::ofstream::out|std::ofstream::trunc);
	if (!file.good()) throw std::runtime_error("File_Io::writeGeometry()");

	file << "TAPSIM GEOMETRY DATA" << '\n';
	file << "VERSION " << "1.0" << '\n';

	if (ioMode == ASCII)
	{
		file.setf(std::ios_base::scientific|std::ios_base::showpoint);
		file.precision(8);
	
		// ***

		const char spacer = '\t';
		const char endline = '\n';

		
		file << "ASCII" << endline;
		file << endline;
		
		// ***

		if (contentFlags & File_Io::NODES)
		{
			file << "NODES " << geomTable.numNodes() << endline;
			for (int i = 0; i < geomTable.numNodes(); i++)
			{
				const Geometry_3d::Point& coords = geomTable.nodeCoords(i);
				
				file << coords.x() << spacer;
				file << coords.y() << spacer;
				file << coords.z() << spacer;

				file << (geomTable.isBoundaryNode(i) ? 1 : 0) << spacer;
			}
			file << endline;
		}
		
		// ***

		if (contentFlags & File_Io::TETRAHEDRA)
		{
			file << "TETRAHEDRA " << geomTable.numTetrahedra() << endline;
			for(int i = 0; i < geomTable.numTetrahedra(); i++)
			{
				const Geometry_3d::Table::TetrahedronData& tet = geomTable.tetrahedron(i);
				
				for (int j = 0; j < 4; j++)
					file << tet.vertices[j] << spacer;
				
				for (int j = 0; j < 3; j++)
					file << tet.neighbours[j] << spacer;

				file << tet.neighbours[3] << endline;
			}
			file << endline;
		}
		
		// ***

		{
			std::multimap<int,int> voronoiFacets;

			if (contentFlags & File_Io::VORONOI_CELLS)
			{
				file << "VORONOI_CELLS " << geomTable.numNodes() << endline;
				for (int i = 0; i < geomTable.numNodes(); i++)
				{
					std::set<int> delaunayNeighbours;
					geomTable.adjacentNodes(i,&delaunayNeighbours);

					std::multimap<int,int>::iterator hint = voronoiFacets.find(i);

					// ***

					file << delaunayNeighbours.size();
					for (std::set<int>::const_iterator j = delaunayNeighbours.begin(); j != delaunayNeighbours.end(); j++)
					{
						file << spacer << *j;
						
						// ***
						
						if (i < *j && contentFlags & File_Io::VORONOI_FACETS)
						{
							const std::pair<int,int> obj(i,*j);
							hint = voronoiFacets.insert(hint,obj);
						}
					}

					file << endline;
				}
				
				file << endline;
			}
			else if (contentFlags & File_Io::VORONOI_FACETS)
			{
				for (int i = 0; i < geomTable.numNodes(); i++)
				{
					std::set<int> delaunayNeighbours;
					geomTable.adjacentNodes(i,&delaunayNeighbours);

					std::multimap<int,int>::iterator hint = voronoiFacets.find(i);

					for (std::set<int>::const_iterator j = delaunayNeighbours.begin(); j != delaunayNeighbours.end(); j++)
					{
						if (i < *j)
						{
							const std::pair<int,int> obj(i,*j);
							hint = voronoiFacets.insert(hint,obj);
						}
					}
				}
			}
			
			// ***

			if (contentFlags & VORONOI_FACETS)
			{
				file << "VORONOI_FACETS " << voronoiFacets.size() << endline;
				
				for (std::multimap<int,int>::iterator i = voronoiFacets.begin(); i != voronoiFacets.end(); i++)
				{
					Geometry_3d::VoronoiFace face;
					geomTable.get_voronoiFace(i->first,i->second,&face);
					
					// ***
					
					file << face.n1 << spacer;
					file << face.n2 << spacer;

					for (int j = 0; j < 3; j++)
						file << face.d1[j] << spacer;

					for (int j = 0; j < 3; j++)
						file << face.d2[j] << spacer;
					
					file << face.tets.size(); 
					for (std::list<int>::const_iterator j = face.tets.begin(); j != face.tets.end();  j++)
						file << spacer << *j;
					
					file << std::endl;
				}
				file << endline;
				
				voronoiFacets.clear();
			}
			
			// ***
			
			if (contentFlags & File_Io::VORONOI_VERTICES)
			{
				file << "VORONOI_VERTICES " << geomTable.numTetrahedra() << endline;
				for (int i = 0; i < geomTable.numTetrahedra(); i++)
				{
					const Geometry_3d::Point center = geomTable.tetCenter(i);
					
					file << center.x() << spacer;
					file << center.y() << spacer;
					file << center.z() << endline;
				}
				file << endline;
			}
		}

		file.close();
	}
	else if (ioMode == BINARY)
	{
		const char endline = '\n';
		
		file << "BINARY" << endline;

		// ***

		if (contentFlags & File_Io::NODES)
		{
			file << "NODES " << geomTable.numNodes() << endline;
			for (int i = 0; i < geomTable.numNodes(); i++)
			{
				File_Utils::binWrite(file,geomTable.nodeCoords(i));

				const char boundary = (geomTable.isBoundaryNode(i) ? 1 : 0);
				File_Utils::binWrite<char>(file,boundary);
			}
			file << endline;
		}

		// ***

		if (contentFlags & File_Io::TETRAHEDRA)
		{
			file << "TETRAHEDRA " << geomTable.numTetrahedra() << endline;
			for(int i = 0; i < geomTable.numTetrahedra(); i++)
				File_Utils::binWrite(file,geomTable.tetrahedron(i));
			file << endline;
		}
		
		// ***

		{
			std::multimap<int,int> voronoiFacets;

			if (contentFlags & File_Io::VORONOI_CELLS)
			{
				file << "VORONOI_CELLS " << geomTable.numNodes() << endline;
				for (int i = 0; i < geomTable.numNodes(); i++)
				{
					std::set<int> delaunayNeighbours;
					geomTable.adjacentNodes(i,&delaunayNeighbours);

					std::multimap<int,int>::iterator hint = voronoiFacets.find(i);
				
					// ***
					
					File_Utils::binWrite<unsigned int>(file,delaunayNeighbours.size());
					for (std::set<int>::const_iterator j = delaunayNeighbours.begin(); j != delaunayNeighbours.end(); j++)
					{
						File_Utils::binWrite<int>(file,*j);

						// ***
						
						if (i < *j && contentFlags & File_Io::VORONOI_FACETS)
						{
							const std::pair<int,int> obj(i,*j);
							hint = voronoiFacets.insert(hint,obj);
						}
					}
				}
				
				file << endline;
			}
			else if (contentFlags & File_Io::VORONOI_FACETS)
			{
				for (int i = 0; i < geomTable.numNodes(); i++)
				{
					std::set<int> delaunayNeighbours;
					geomTable.adjacentNodes(i,&delaunayNeighbours);

					std::multimap<int,int>::iterator hint = voronoiFacets.find(i);

					for (std::set<int>::const_iterator j = delaunayNeighbours.begin(); j != delaunayNeighbours.end(); j++)
					{
						if (i < *j)
						{
							const std::pair<int,int> obj(i,*j);
							hint = voronoiFacets.insert(hint,obj);
						}
					}
				}
			}
			
			// ***

			if (contentFlags & File_Io::VORONOI_FACETS)
			{
				file << "VORONOI_FACETS " << voronoiFacets.size() << endline;
				
				for (std::multimap<int,int>::iterator i = voronoiFacets.begin(); i != voronoiFacets.end(); i++)
				{
					Geometry_3d::VoronoiFace face;
					geomTable.get_voronoiFace(i->first,i->second,&face);
					
					// ***

					File_Utils::binWrite<int>(file,face.n1);
					File_Utils::binWrite<int>(file,face.n2);
					File_Utils::binWrite(file,face.d1);
					File_Utils::binWrite(file,face.d2);
					
					File_Utils::binWrite<unsigned int>(file,face.tets.size());
					for (std::list<int>::const_iterator j = face.tets.begin(); j != face.tets.end();  j++)
						File_Utils::binWrite<int>(file,*j);
				}
				file << endline;
				
				voronoiFacets.clear();
			}
			
			// ***
			
			if (contentFlags & File_Io::VORONOI_VERTICES)
			{
				file << "VORONOI_VERTICES " << geomTable.numTetrahedra() << endline;
				for (int i = 0; i < geomTable.numTetrahedra(); i++)
					File_Utils::binWrite(file,geomTable.tetCenter(i));
				file << endline;
			}
		}

		file.close();
	}
	else
		throw std::runtime_error("File_Io::writeGeometry()");

	file.close();
}

// ***** ===== FILE_IO::WRITE_GEOM_GRID() ===== *****

void File_Io::writeGeomGrid(const char* filename, const Geometry_3d::Table& geomTable, const Grid_3d::Table& gridTable, const int ioMode, const int withNumbers, const int withPotentials)
{
	if (geomTable.numNodes() != gridTable.numNodes()) throw std::runtime_error("File_Io::writeGeomGrid()");

	// ***

	std::ofstream file(filename,std::ofstream::out|std::ofstream::trunc);
	if (!file.good()) throw std::runtime_error("File_Io::writeGeomGrid()");

	if (ioMode == ASCII)
	{
		file.setf(std::ios_base::scientific|std::ios_base::showpoint);
		file.precision(8);
	
		// ***

		const char spacer = '\t';
		const char endline = '\n';

		file << "ASCII " << geomTable.numNodes() << " ";
		file << withNumbers << " " << withPotentials << endline;

		for (int i = 0; i < geomTable.numNodes(); i++)
		{
			for (short j = 0; j < 3; j++)
				file << geomTable.node(i).coords[j] << spacer;
	
			file << gridTable.node(i).id().toValue();

			if (withNumbers) file << spacer << gridTable.node(i).number().toValue();
			if (withPotentials) file << spacer << gridTable.potential(i);

			file << endline;
		}
	}
	else if (ioMode == BINARY)
	{
		file << "BINARY " << geomTable.numNodes() << " ";
		file << withNumbers << " " << withPotentials << '\n';

		for (int index = 0; index < geomTable.numNodes();index++)
		{
			File_Utils::binWrite(file,geomTable.node(index).coords);
			File_Utils::binWrite<short>(file,gridTable.node(index).id().toValue());

			if (withNumbers) File_Utils::binWrite<unsigned int>(file,gridTable.node(index).number().toValue());
			if (withPotentials) File_Utils::binWrite<float>(file,gridTable.potential(index));
		}
	}
	else
		throw std::runtime_error("File_Io::writeGeomGrid()");
	
	file.close();
}

// ***** ===== FILE_IO::READ_GEOM_GRID() ===== *****

void File_Io::readGeomGrid(const char* filename, Geometry_3d::Table* geomTable, Grid_3d::Table* gridTable, int* withNumbers, int* withPotentials)
{
	std::ifstream file(filename,std::ifstream::in);
	if (!file.good()) throw std::runtime_error("File_Io::readGeomGrid()");

	// ***

	std::string headerLine;
	
	do
		std::getline(file,headerLine);
	while (headerLine.empty() || headerLine.at(0) == '#');

	bool asciiMode;

	if (headerLine.substr(0,5) == "ASCII")
	{
		asciiMode = true;
		headerLine = headerLine.substr(5,std::string::npos);
	}
	else if (headerLine.substr(0,6) == "BINARY")
	{
		asciiMode = false;
		headerLine = headerLine.substr(6,std::string::npos);
	}
	else
		throw std::runtime_error("File_Io::readGeomGrid()");

	int size(0);
	
	bool readNumbers(false);
	bool readPotentials(false);
	
	{
		std::istringstream line(headerLine);

		line >> size;
		if (line.fail()) throw std::runtime_error("File_Io::readGeomGrid()");
		
		line >> std::noboolalpha >> readNumbers;
		if (line.fail()) throw std::runtime_error("File_Io::readGeomGrid()");
		
		line >> std::noboolalpha >> readPotentials;
		if  (line.fail()) throw std::runtime_error("File_Io::readGeomGrid()");
	}
	
	// ***

	geomTable->allocate(size);
	gridTable->allocate(size);

	for (int index = 0; index < size; index++)
	{
		Geometry_3d::Point coords(0.0f);
		short id = Configuration::NodeId::invalid().toValue();

		unsigned int number = Configuration::NodeNumber::invalid().toValue();
		float phi = 0.0f;
		
		// ***
		
		if (asciiMode)
		{
			if (!file.good()) throw std::runtime_error("File_Io::readGeomGrid()");
	
			std::string buffer;
			
			do
				std::getline(file,buffer);
			while (buffer.at(0) == '#');
	
			std::list<std::string> values;
			std::istringstream line(buffer);
		
			do
			{
				std::string tmp;
				std::getline(line,tmp,'\t');
				values.push_back(tmp);
			}
			while (line.good());
	
			std::list<std::string>::const_iterator i = values.begin();
	
			// *** mandatory content (coordinate information, identifier)
	
			if (values.end() != i)
				coords.x() = File_Utils::read_element<float>(*i++);
			else
				throw std::runtime_error("File_Io::readGeomGrid(): coords.x");
	
			if (values.end() != i)
				coords.y() = File_Utils::read_element<float>(*i++);
			else
				throw std::runtime_error("File_Io::readGeomGrid(): coords.y");
	
			if (values.end() != i)
				coords.z() = File_Utils::read_element<float>(*i++);
			else
				throw std::runtime_error("File_Io::readGeomGrid(): coords.z");

			if (values.end() != i)
				id = File_Utils::read_element<short>(*i++);
			else
				throw std::runtime_error("File_Io::readGeomGrid(): id");
	
			// *** optional content (number, potential)
				
			if (readNumbers)
			{
				if (values.end() != i) 
					number = File_Utils::read_element<unsigned int>(*i++);
				else
					throw std::runtime_error("File_Io::readGeomGrid(): number");
			}

			if (readPotentials)
			{
				if (values.end() != i)
					phi = File_Utils::read_element<float>(*i++);
				else
					throw std::runtime_error("File_Io::readGeomGrid(): potential");
			}
		}
		else
		{
			File_Utils::binRead(file,&coords);
			File_Utils::binRead(file,&id);
			
			if (readNumbers) File_Utils::binRead(file,&number);
			if (readPotentials) File_Utils::binRead(file,&phi);
		}

		// ***

		geomTable->node(index).coords = coords;
		geomTable->node(index).boundary = false;
		
		gridTable->node(index).setId(Configuration::NodeId::fromValue(id));
		gridTable->node(index).setNumber(Configuration::NodeNumber::fromValue(number));
		gridTable->node(index).setPhi(0,phi);
		gridTable->node(index).setPhi(1,phi);
		
		// ***
		
		if (!file.good()) throw std::runtime_error("File_Io::readGeomGrid()");
	}

	// *** do renumbering of invalid unique number entries
	// *** (and also checks for unique numbering of the valid number input)

	// => initially invalid numbers are assigned a valid value
	// => initially valid numbers remain unchanged

	{
		std::set<Configuration::NodeNumber> numberTable;
		for (int i = 0; i < gridTable->numNodes(); i++)
		{
			if (!gridTable->number(i).isValid()) continue;
			std::pair<std::set<Configuration::NodeNumber>::iterator,bool> tmp = numberTable.insert(gridTable->number(i));
			if (!tmp.second) throw std::runtime_error("File_Io::readGeomGrid(): non unique numbering detected");
		}
	}

	{
		Configuration::NodeNumber number = Configuration::NodeNumber::first();

		for (int i = 0; i < gridTable->numNodes(); i++)
		{
			if (!gridTable->number(i).isValid()) continue;
			if (gridTable->number(i) >= number) number = gridTable->number(i).next();
		}

		for (int i = 0; i < gridTable->numNodes(); i++)
		{
			if (!gridTable->number(i).isValid())
			{
				gridTable->node(i).setNumber(number);
				number = number.next();
			}
		}
	}

	// ***

	if (withNumbers != 0) *withNumbers = readNumbers;
	if (withPotentials != 0) *withPotentials = readNumbers;

	file.close();
}

void File_Io::writeSystem(const char* filename, const System_3d& system, const unsigned int index)
{
	std::ofstream file(filename,std::ofstream::out);
	if (!file.good()) throw std::runtime_error("File_Io::write()");

	file << "TAPSIM DUMP INFORMATION" << '\n';
	file << "LAST_INDEX " << index << '\n';

	system.configTable >> file;
	system.geomTable >> file;
	system.gridTable >> file;

	file.close();
}
void File_Io::readSystem(const char* filename, System_3d* system, unsigned int* index)
{
	std::ifstream file(filename,std::ifstream::in);
	if (!file.good()) throw std::runtime_error("File_Io::read()");

	std::string buffer;

	std::getline(file,buffer);
	if (buffer.substr(0,23) != "TAPSIM DUMP INFORMATION")
		throw std::runtime_error("File_Io::readSystem(): wrong file type");

	std::getline(file,buffer);
	if (buffer.substr(0,10) != "LAST_INDEX") 
		throw std::runtime_error("File_Io::readSystem()");
	
	try
	{
		buffer = buffer.substr(10,std::string::npos);
		buffer = File_Utils::trim(buffer);
	}
	catch (std::out_of_range&)
	{
		throw std::runtime_error("File_Io::readSystem()");
	}
	
	if (index != 0) 
		*index = File_Utils::read_element<unsigned int>(buffer);
	else
		File_Utils::read_element<unsigned int>(buffer); // test for conforming file format

	system->configTable << file;
	system->geomTable << file;
	system->gridTable << file;

	file.close();
}
