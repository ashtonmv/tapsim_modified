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

#include "debug.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <vector>
#include <string>
#include <map>
#include <set>

#include <cstdlib>

#include "fieldline_3d.h"
#include "trajectory_3d.h"
#include "geometry_3d.h"

// ***

void Debug::printNode(std::ostream* stream, const int index, const System_3d& system)
{
	const Geometry_3d::Table::NodeData& geomObj = system.geomTable.node(index);
	const Grid_3d::Node& gridObj = system.gridTable.node(index);

	*stream << "*** Node: #" << index << std::endl;

	if (index < 0)
	{
		*stream << "\tINVALID" << std::endl;
		return;
	}

	*stream << "\tcoords: (";
	*stream << geomObj.coords[0] << ", ";
	*stream << geomObj.coords[1] << ", ";
	*stream << geomObj.coords[2] << ")" << std::endl;

	*stream << "\tboundary: " << geomObj.boundary;
	*stream << "/" << gridObj.boundary() << std::endl;

	*stream <<"\tassociated tetrahedron: " << geomObj.associatedTetrahedron << std::endl;

	// ***

	*stream << "\tid: " << gridObj.id().toValue() << std::endl;
	*stream << "\tnumber: " << gridObj.number().toValue() << std::endl;
	*stream << "\tphi0/phi1: " << gridObj.phi(0) << "/" << gridObj.phi(1) << std::endl;
	*stream << "\tcharge: " << gridObj.charge() << std::endl;
	*stream << "\tboundary: " << gridObj.boundary() << std::endl;
	*stream << "\tdirichletBoundary: " << gridObj.dirichletBoundary() << std::endl;
	*stream << "\tneumannBoundary: " << gridObj.neumannBoundary() << std::endl;
	*stream << "\tneighbours (size = " << gridObj.numNeighbours() << "):" << std::endl;
	
	for (int i = 0; i < gridObj.numNeighbours(); i++)
	{
		const int neighIndex = gridObj.neighbour(i);

		*stream << "\t\t=> #" << i << ": ";
		*stream << "index = " << gridObj.neighbour(i) << ", ";
		*stream << "coords = (" << system.geomTable.node(neighIndex).coords[0] << "/";
		*stream << system.geomTable.node(neighIndex).coords[1] << "/";
		*stream << system.geomTable.node(neighIndex).coords[2] << "), ";
		*stream << "phi = " << system.gridTable.node(neighIndex).phi(0) << "/";
		*stream << system.gridTable.node(neighIndex).phi(1) << ", ";
		*stream << "coupling = " << gridObj.coupling(i) << std::endl;
	}

	// ***

	*stream << "=> potential: " << system.gridTable.potential(index) << std::endl;

	const MathVector3d<float> fieldValue = system.gridTable.field_o1(index,system.geomTable);
	*stream << "=> field: (" << fieldValue.x() << ", " <<fieldValue.y() << ", " << fieldValue.z() << ")" << std::endl;
}

void Debug::printTetrahedron(std::ostream* stream, const int index, const System_3d& system)
{
	*stream << "*** Tetrahedron: #" << index << std::endl;

	if (index < 0)
	{
		*stream << "\tINVALID" << std::endl;
		return;
	}

	*stream << "\tvertices: ";
	*stream << system.geomTable.tetrahedron(index).vertices[0] << ", ";
	*stream << system.geomTable.tetrahedron(index).vertices[1] << ", ";
	*stream << system.geomTable.tetrahedron(index).vertices[2] << ", ";
	*stream << system.geomTable.tetrahedron(index).vertices[3] << std::endl;

	*stream << "\tneighbours: ";
	*stream << system.geomTable.tetrahedron(index).neighbours[0] << ", ";
	*stream << system.geomTable.tetrahedron(index).neighbours[1] << ", ";
	*stream << system.geomTable.tetrahedron(index).neighbours[2] << ", ";
	*stream << system.geomTable.tetrahedron(index).neighbours[3] << std::endl;

	*stream << "\tcoords: ";
	for (int i = 0; i < 4; i++)
	{
		const int nodeIndex = system.geomTable.tetrahedron(index).vertices[i];

		*stream << "(" << system.geomTable.node(nodeIndex).coords.x() << "/";
		*stream << system.geomTable.node(nodeIndex).coords.y() << "/";
		*stream << system.geomTable.node(nodeIndex).coords.z() << ")";

		if (i < 3) *stream << ", "; else *stream << std::endl;
	}

	*stream << "\tcenter: (" << system.geomTable.tetCenter(index).x() << "/";
	*stream << system.geomTable.tetCenter(index).y() << "/";
	*stream << system.geomTable.tetCenter(index).z() << ")";
	*stream << std::endl;

	*stream << "\torient: " << system.geomTable.tetOrient(index) << std::endl;
}

void Debug::printVoronoiFace(std::ostream* stream, const Geometry_3d::VoronoiFace& face)
{
	*stream << "*** Voronoi-face: #" << face.n1 << ", #" << face.n2 << std::endl;
	*stream << "\tn1-coords: " << face.p1.x() << "/" << face.p1.y() << "/" << face.p1.z() << std::endl;
	*stream << "\tn2-coords: " << face.p2.x() << "/" << face.p2.y() << "/" << face.p2.z() << std::endl;

	*stream << "\tfinite: " << face.isFinite() << std::endl;
	*stream << "\td1-coords: " << face.d1.x() << "/" << face.d1.y() << "/" << face.d1.z() << std::endl;
	*stream << "\td2-coords: " << face.d2.x() << "/" << face.d2.y() << "/" << face.d2.z() << std::endl;

	*stream << "\tcenter: " << face.center().x() << "/" << face.center().y() << "/" << face.center().z() << std::endl;
	*stream << "\tarea: " << face.area() << std::endl;
	*stream << "\ttetrahedra => vertices: " << std::endl;

	int index = 1;
	std::list<int>::const_iterator i = face.tets.begin();
	std::list<Geometry_3d::Point>::const_iterator j = face.vertices.begin();

	while (i != face.tets.end())
	{
		*stream << "\t\t" << index << " -> tet: " << *i;
		*stream << ", center: (" << j->x() << ", " << j->y() << ", " << j->z() << ")" << std::endl;

		i++;
		j++;
		index++;
	}
}

void Debug::printVoronoiCell(std::ostream* stream, const Geometry_3d::VoronoiCell& cell)
{
	*stream << "*** Voronoi-cell: #" << cell.n << std::endl;
	*stream << "\tcenter-coords: " << cell.p.x() << "/" << cell.p.y() << "/" << cell.p.z() << std::endl;
	*stream << "\tfinite: " << cell.isFinite() << std::endl;
	*stream << "\tvolume: " << cell.volume() << std::endl;

	*stream << "\tfaces (n1 -> n2: finite, area):" << std::endl;

	for (std::list<Geometry_3d::VoronoiFace>::const_iterator i = cell.faces.begin(); i != cell.faces.end(); i++)
		*stream << "\t\t" << i->n1 << " -> " << i->n2 << ": finite = " << i->isFinite() << ", " << i->area() << std::endl;
}

// ***

void Debug::printEvapOptions(std::ostream* stream, const Process::EvaporationOptions& options)
{
	*stream << "results data file: \"" << options.resultsFile << "\"" << std::endl;
	*stream << "results data file mode: \"" << options.resultsMode << "\"" << std::endl;
	*stream << "results data chunk size: \"" << options.resultsChunkSize << "\"" << std::endl;

	*stream << "trajectory data file: \"" << options.trajectoryFile << "\"" << std::endl;
	*stream << "trajectory data file mode: \"" << options.trajectoryMode << "\"" << std::endl;
	*stream << "trajectory data chunk size: \"" << options.trajectoryChunkSize << "\"" << std::endl;

	*stream << "grid data file: \"" << options.gridFile << "\"" << std::endl;
	*stream << "grid data file mode: \"" << options.gridMode << "\"" << std::endl;
	*stream << "grid data snapshot interval: \"" << options.gridInterval << "\"" << std::endl;

	*stream << "surface data file: \"" << options.surfaceFile << "\"" << std::endl;
	*stream << "surface data file mode: \"" << options.surfaceMode << "\"" << std::endl;
	*stream << "surface data snapshot interval: \"" << options.surfaceInterval << "\"" << std::endl;

	*stream << "geometry data file: \"" << options.geometryFile << "\"" << std::endl;
	*stream << "geometry data file mode: \"" << options.geometryMode << "\"" << std::endl;
	*stream << "geometry data contents: \"" << options.geometryParams << "\"" << std::endl;
	
	*stream << "dump file: \"" << options.dumpFile << "\"" << std::endl;
	*stream << "dump interval: \"" << options.dumpInterval << "\"" << std::endl;

	*stream << "output delay time: \"" << options.delayTime << "\"" << std::endl;

	*stream << "method used for computing evaporation probabilities: \"" << options.probMode << "\"" << std::endl;
	*stream << "method used for choosing atoms next in evaporation: \"" << options.evapMode << "\"" << std::endl;
	
	*stream << "vacuum identifier name: \"" << options.vacuumName << "\"" << std::endl;
	
	*stream << "trajectory integrator-type: \"" << options.trajectoryIntegrator << "\"" << std::endl;
	*stream << "trajectory stepper-type: \"" << options.trajectoryStepper << "\"" << std::endl;
	*stream << "trajectory initial time step: \"" << options.trajectoryTimeStep << "\"" << std::endl;
	*stream << "trajectory time-step-limit: \"" << options.trajectoryTimeStepLimit << "\"" << std::endl;
	*stream << "trajectory error threshdol: \"" << options.trajectoryErrorThreshold << "\"" << std::endl;
	*stream << "evaporation fixed initial position: \"" << options.fixedInitialPosition << "\"" << std::endl;
	*stream << "no initial after desorption: \"" << options.noInitialVelocity << "\"" << std::endl;

	*stream << "initial shrinkage: \"" << options.initShrinkage << "\"" << std::endl;
	*stream << "shrink limit: \"" << options.shrinkLimit << "\"" << std::endl;

	*stream << "initial event count: \"" << options.initEventCnt << "\"" << std::endl;
	*stream << "event count limit: \"" << options.eventCntLimit << "\"" << std::endl;

	*stream << "voltage queue limit: \"" << options.voltageQueueSize << "\"" << std::endl;

	*stream << "local relaxation threshold: \"" << options.relaxThreshold << "\"" << std::endl;
	*stream << "local relaxation shell size: \"" << options.relaxShellSize << "\"" << std::endl;
	*stream << "local relaxation cycle size: \"" << options.relaxCycleSize << "\"" << std::endl;
	*stream << "local relaxation queue size: \"" << options.relaxQueueSize << "\"" << std::endl;
	*stream << "global relaxation iteration number: \"" << options.relaxGlobalCycles << "\"" << std::endl;

	*stream << "refresh relaxation interval: \"" << options.refreshInterval << "\"" << std::endl;
	*stream << "refresh relaxation threshold: \"" << options.refreshThreshold << "\"" << std::endl;
	*stream << "refresh relaxation cycle size: \"" << options.refreshCycleSize << "\"" << std::endl;
	*stream << "refresh relaxation queue size: \"" << options.refreshQueueSize << "\"" << std::endl;
}

// ***

namespace
{
	class pv_pointSet
	{
		public:
			pv_pointSet(const int dataSize =0)
				: _nextIndex(0),
				  _points(&less),
				  _data(dataSize)
			{}

			// ***

			int insert(const Geometry_3d::Point& point)
			{
				if (_points.find(point) == _points.end()) _points[point] = _nextIndex++;
				return _points[point];
			}

			void setData(const int pointIndex, const int dataIndex, const float value, bool overwrite =false)
			{
				if (!overwrite && _data.at(dataIndex).find(pointIndex) != _data.at(dataIndex).end()) return;
				_data.at(dataIndex)[pointIndex] = value;
			}

			// ***

			std::map<int,Geometry_3d::Point> points()
			{
				std::map<int,Geometry_3d::Point> results;
				for (std::map<Geometry_3d::Point,int>::const_iterator i = _points.begin(); i != _points.end(); i++)
					results[i->second] = i->first;

				return results;
			}

			const std::map<int,float>& data(const int index)
			{
				if (_data.at(index).size() != _points.size())
				{
					for (std::map<Geometry_3d::Point,int>::const_iterator i = _points.begin(); i != _points.end(); i++)
						if (_data.at(index).find(i->second) == _data.at(index).end()) _data.at(index)[i->second] = 0.0f;
				}

				return _data.at(index);
			}
	
			int size() const { return _points.size(); }

			int dataSize() const { return _data.size(); }

		private:
			static bool less(const Geometry_3d::Point& a, const Geometry_3d::Point& b)
			{
				for (int i = 0; i < 3; i++)
				{
					if (a[i] < b[i])
						return true;
					else
						if (a[i] > b[i]) return false;
				}
		
				return false;
			}

			int _nextIndex;

			std::map<Geometry_3d::Point,int,bool(*)(const Geometry_3d::Point&,const Geometry_3d::Point&)> _points;
			std::vector<std::map<int,float> > _data;
	};

	enum { POINT = 1, LINE = 2, POLYLINE =4, POLYGON =7, TETRAHEDRON = 10 };

	class pv_cellData
	{
		public:
			pv_cellData(const int type =1, const int dataSize =0)
				: _type(type),
				  _nodes(),
				  _data(dataSize)
			{}

			pv_cellData& operator<<(const int index)
			{
				_nodes.push_back(index);
				return *this;
			}

			int type() const { return _type; }

			const std::list<int>& nodes() const { return _nodes; }

			void setData(const int index, const float value) { _data.at(index) = value; }
			float data(const int index) const { return _data.at(index); }

			void remove() { _nodes.pop_back(); }

			int size() const { return _nodes.size(); }
			int dataSize() const { return _data.size(); }

		private:
			int _type;
			std::list<int> _nodes;
			std::vector<float> _data;
	};

	class pv_cellSet
	{
		public:
			pv_cellSet(const int dataSize =0)
				: _dataSize(dataSize),
				  _cells()
			{}
	
			void add(const pv_cellData& obj)
			{
				if (obj.dataSize() != _dataSize) throw std::runtime_error("pv_cellSet(): datasize mismatch!");
				_cells.push_back(obj);
			}
	
			const std::list<pv_cellData>& cells() const { return _cells; }

			int numCells() { return _cells.size(); }
	
			int numNodes()
			{
				int totalCnt(0);
				for (std::list<pv_cellData>::const_iterator i = _cells.begin(); i != _cells.end(); i++)
					totalCnt += i->nodes().size();
	
				return totalCnt;
			}

		private:
			int _dataSize;
			std::list<pv_cellData> _cells;
	};

	int swapEndian(const int value)
	{
		int result;

		const char* vPtr = reinterpret_cast<const char*>(&value);
		char* rPtr = reinterpret_cast<char*>(&result);

		for (unsigned int i = 0; i < sizeof(int); i++)
			rPtr[i] = vPtr[sizeof(int)-1-i];

		return result;
	}

	float swapEndian(const float value)
	{
		float result;

		const char* vPtr = reinterpret_cast<const char*>(&value);
		char* rPtr = reinterpret_cast<char*>(&result);

		for (unsigned int i = 0; i < sizeof(int); i++)
			rPtr[i] = vPtr[sizeof(int)-1-i];

		return result;
	}
}

void Debug::writeData_PARAVIEW(const char* filename, const System_3d& system)
{
	if (system.geomTable.numNodes() != system.gridTable.numNodes())
		throw std::runtime_error("Debug::writeData_PARAVIEW() / DELAUNAY NODES");

	// ***

	std::cout << "*** Writing delaunay network in PARAVIEW format:" << std::endl;

	// *** GENERATE ASCII OUTPUT

	const char spacer = ' ';
	const char endline = '\n';

	std::ofstream file(filename,std::ofstream::out);
	file.setf(std::ios_base::scientific|std::ios_base::showpoint);

	std::cout << "\t-> filename: \"" << filename << "\"" << std::endl;

	// *** HEADER

	file << "# vtk DataFile Version 1.0" << endline;
	file << "Vorfinite-Tap => Delaunay Network" << endline;
	file << "ASCII" << endline;
	file << endline;

	file << "DATASET" << spacer << "UNSTRUCTURED_GRID" << endline;
	file << endline;

	// *** POINT DEFINITION

	file << "POINTS" << spacer << system.geomTable.numNodes() << spacer << "float" << endline;
	for (int i = 0; i < system.geomTable.numNodes(); i++)
	{
		file << system.geomTable.nodeCoords(i).x() << spacer;
		file << system.geomTable.nodeCoords(i).y() << spacer;
		file << system.geomTable.nodeCoords(i).z() << endline;
	}
	file << endline;

	std::cout << "\t-> number of points: " << system.geomTable.numNodes() << std::endl;

	// *** CELL DEFINITION

	file << "CELLS" << spacer << system.geomTable.numTetrahedra() << spacer << system.geomTable.numTetrahedra() * 5 << endline;
	for (int i = 0; i < system.geomTable.numTetrahedra(); i++)
	{
		file << 4;
		for (int j = 0; j < 4; j++)
			file << spacer << system.geomTable.tetrahedron(i).vertices[j];
		file << endline;
	}
	file << endline;

	file << "CELL_TYPES" << spacer << system.geomTable.numTetrahedra() << endline;
	for (int i = 0; i < system.geomTable.numTetrahedra(); i++)
		file << TETRAHEDRON << endline;
	file << endline;

	std::cout << "\t-> number of cells: " << system.geomTable.numNodes() << " + " << system.geomTable.numTetrahedra() << std::endl;

	// ***  POINT DATA

	file << "POINT_DATA" << spacer << system.geomTable.numNodes() << endline;
	file << endline;

// 	file << "SCALARS" << spacer << "nodeId" << spacer << "int" << spacer << 1 << endline;
// 	file << "LOOKUP_TABLE" << spacer << "default" << endline;
// 	for (int i = 0; i < system.geomTable.numNodes(); i++)
// 		file << i << endline;
// 	file << endline;

	file << "SCALARS" << spacer << "phi" << spacer << "float" << spacer << 1 << endline;
	file << "LOOKUP_TABLE" << spacer << "default" << endline;
	for (int i = 0; i < system.gridTable.numNodes(); i++)
		file << system.gridTable.potential(i) << endline;
	file << endline;

// 	file << "SCALARS" << spacer << "flux" << spacer << "float" << spacer << 1 << endline;
// 	file << "LOOKUP_TABLE" << spacer << "default" << endline;
// 	for (int i = 0; i < system.gridTable.numNodes(); i++)
// 		file << system.gridTable.flux(i,system.geomTable) << endline;
// 	file << endline;

// 	file << "VECTORS" << spacer << "field" << spacer << "float" << endline;
// 	for (int i = 0; i < system.gridTable.numNodes(); i++)
// 	{
// 		MathVector3d<float> field = system.gridTable.field_o1(i,system.geomTable);
// 
// 		file << field.x() << spacer;
// 		file << field.y() << spacer;
// 		file << field.z() << endline;
// 	}
// 	file << endline;

	std::cout << "\t-> point-data: node-ids, potentials, fluxes, field-vectors" << std::endl;
	std::cout << std::endl;

	file.close();


	// ***

	// 	int index;
	// 
	// 	pv_pointSet pointSet(6);
	// 
	// 	// *** generators
	// 
	// 	for (int i = 0; i < system.geomTable.numNodes(); i++)
	// 	{
	// 		if (i % (1+system.geomTable.numNodes()/100) == 0)
	// 		{
	// 			std::cout << "\t-> processing generators: " << i << " / " << system.geomTable.numNodes() << '\r';
	// 			std::cout.flush();
	// 		}
	// 
	// 		// ***
	// 
	// 		index = pointSet.insert(system.geomTable.nodeCoords(i));
	// 		pointSet.setData(index,0,static_cast<float>(system.gridTable.node(i).id().toValue()));
	// 		pointSet.setData(index,1,system.gridTable.potential(i));
	// 		pointSet.setData(index,2,system.gridTable.flux(i,system.geomTable));
	// 
	// 		const MathVector3d<float> fieldValue = system.gridTable.field_o1(i,system.geomTable);
	// 		pointSet.setData(index,3,fieldValue.x());
	// 		pointSet.setData(index,4,fieldValue.y());
	// 		pointSet.setData(index,5,fieldValue.z());
	// 	}
	// 
	// 	std::cout << "\t-> processing generators: finished" << std::string(25,' ') << std::endl;
	// 
	// 	// *** tetrahedra
	// 
	// 	pv_cellSet tetrahedra(1);
	// 
	// 	for (int i = 0; i < system.geomTable.numTetrahedra(); i++)
	// 	{
	// 		if (i % (1+system.geomTable.numTetrahedra()/100) == 0)
	// 		{
	// 			std::cout << "\t-> processing tetrahedra: " << i << " / " << system.geomTable.numTetrahedra() << '\r';
	// 			std::cout.flush();
	// 		}
	// 
	// 		// ***
	// 
	// 		Geometry_3d::Tetrahedron obj;
	// 
	// 		for (int j = 0; j < 4; j++)
	// 			obj[j] = system.geomTable.nodeCoords(system.geomTable.tetrahedron(i).vertices[j]);
	// 
	// 		pv_cellData tetrahedron(TETRAHEDRON,1);
	// 		tetrahedron.setData(0,obj.volume());
	// 
	// 		for (int j = 0; j < 4; j++)
	// 		{
	// 			index = pointSet.insert(obj[j]);
	// 			tetrahedron << index;
	// 		}
	// 
	// 		tetrahedra.add(tetrahedron);
	// 	}
	// 
	// 	std::cout << "\t-> processing tetrahedra: finished" << std::string(25,' ') << std::endl;
	// 
	// 	// *** GENERATE ASCII OUTPUT
	// 
	// 	const char spacer = ' ';
	// 	const char endline = '\n';
	// 
	// 	std::cout << "\t-> writing data to file \"" << filename << "\" ..." << '\r';
	// 	std::cout.flush();
	// 
	// 	std::ofstream file(filename,std::ofstream::out);
	// 	file.setf(std::ios_base::scientific|std::ios_base::showpoint);
	// 
	// 	// *** HEADER
	// 
	// 	file << "# vtk DataFile Version 1.0" << endline;
	// 	file << "Vorfinite-Tap => Delaunay Network" << endline;
	// 	file << "ASCII" << endline;
	// 	file << endline;
	// 
	// 	file << "DATASET" << spacer << "UNSTRUCTURED_GRID" << endline;
	// 	file << endline;
	// 
	// 	// *** POINT DEFINITION
	// 
	// 	file << "POINTS" << spacer << pointSet.size() << spacer << "float" << endline;
	// 	const std::map<int,Geometry_3d::Point> sortedPoints = pointSet.points();
	// 	for (std::map<int,Geometry_3d::Point>::const_iterator i = sortedPoints.begin(); i != sortedPoints.end(); i++)
	// 	{
	// 		file << i->second.x() << spacer;
	// 		file << i->second.y() << spacer;
	// 		file << i->second.z() << endline;
	// 	}
	// 	file << endline;
	// 
	// 	// *** CELL DEFINITION (1/2)
	// 
	// 	file << "CELLS" << spacer << tetrahedra.numCells() << spacer << tetrahedra.numCells() + tetrahedra.numNodes() << endline;
	// 	for (std::list<pv_cellData>::const_iterator j = tetrahedra.cells().begin(); j != tetrahedra.cells().end(); j++)
	// 	{
	// 		if (j->nodes().size() == 0) continue;
	// 
	// 		file << j->nodes().size();
	// 		for (std::list<int>::const_iterator k = j->nodes().begin(); k != j->nodes().end(); k++)
	// 			file << spacer << *k;
	// 		file << endline;
	// 	}
	// 	file << endline;
	// 
	// 	// *** CELL DEFINITION (2/2)
	// 
	// 	file << "CELL_TYPES" << spacer << tetrahedra.numCells() << endline;
	// 	for (std::list<pv_cellData>::const_iterator j = tetrahedra.cells().begin(); j != tetrahedra.cells().end(); j++)
	// 		if (j->nodes().size() != 0) file << j->type() << endline;
	// 	file << endline;
	// 
	// 	// ***  POINT DATA
	// 
	// 	file << "POINT_DATA" << spacer << pointSet.size() << endline;
	// 	file << endline;
	// 
	// 	file << "SCALARS" << spacer << "nodeId" << spacer << "int" << spacer << 1 << endline;
	// 	file << "LOOKUP_TABLE" << spacer << "default" << endline;
	// 	for (std::map<int,float>::const_iterator i = pointSet.data(0).begin(); i != pointSet.data(0).end(); i++)
	// 		file << static_cast<int>(i->second) << endline;
	// 	file << endline;
	// 
	// 	file << "SCALARS" << spacer << "phi" << spacer << "float" << spacer << 1 << endline;
	// 	file << "LOOKUP_TABLE" << spacer << "default" << endline;
	// 	for (std::map<int,float>::const_iterator i = pointSet.data(1).begin(); i != pointSet.data(1).end(); i++)
	// 		file << i->second << endline;
	// 	file << endline;
	// 
	// 	file << "SCALARS" << spacer << "flux" << spacer << "float" << spacer << 1 << endline;
	// 	file << "LOOKUP_TABLE" << spacer << "default" << endline;
	// 	for (std::map<int,float>::const_iterator i = pointSet.data(2).begin(); i != pointSet.data(2).end(); i++)
	// 		file << i->second << endline;
	// 	file << endline;
	// 
	// 	file << "VECTORS" << spacer << "field" << spacer << "float" << endline;
	// 	{
	// 		std::map<int,float>::const_iterator i = pointSet.data(3).begin();
	// 		std::map<int,float>::const_iterator j = pointSet.data(4).begin();
	// 		std::map<int,float>::const_iterator k = pointSet.data(5).begin();
	// 
	// 		while (pointSet.data(3).end() != i)
	// 		{
	// 			file << i->second << spacer;
	// 			file << j->second << spacer;
	// 			file << k->second << endline;
	// 
	// 			i++; j++; k++;
	// 		}
	// 	}
	// 	file << endline;
	// 
	// 	// *** CELL DATA
	// 
	// 	file << "CELL_DATA" << spacer << tetrahedra.numCells() << endline;
	// 	file << endline;
	// 	
	// 	file << "SCALARS" << spacer << "delaunayVolume" << spacer << "float" << spacer << 1 << endline;
	// 	file << "LOOKUP_TABLE" << spacer << "default" << spacer << endline;
	// 	for (std::list<pv_cellData>::const_iterator j = tetrahedra.cells().begin(); j != tetrahedra.cells().end(); j++)
	// 		file << j->data(0) << endline;
	// 	file << endline;
	// 
	// 	// ***
	// 
	// 	std::cout << "\t-> writing data to file \"" << filename << "\" (finished)" << std::endl;
	// 	std::cout << std::endl;
	// 
	// 	file.close();
}

void Debug::writeData_PARAVIEW(const char* filename, const System_3d& system, const Geometry_3d::VoronoiFace& face)
{
	if (system.geomTable.numNodes() != system.gridTable.numNodes())
		throw std::runtime_error("Debug::writeData_PARAVIEW() / SINGLE FACE");

	// ***

	const float infiniteLength = 0.01f * (system.geomTable.max() - system.geomTable.min()).length();

	int index;

	pv_pointSet pointSet(1);
	pv_cellSet cellSet(1);

	// ***

	{
		pv_cellData delaunayCell(POLYLINE,1);
		delaunayCell.setData(0,0.0f);

		index = pointSet.insert(face.p1);
		pointSet.setData(index,0,0.0f);
		delaunayCell << index;

		index = pointSet.insert(face.p2);
		pointSet.setData(index,0,0.0f);
		delaunayCell << index;

		cellSet.add(delaunayCell);
	}

	if (!face.isFinite())
	{
		pv_cellData faceCell(POLYGON,1);
		faceCell.setData(0,2.0f);

		index = pointSet.insert(face.vertices.front() + infiniteLength * face.d1);
		pointSet.setData(index,0,2);
		faceCell << index;

		for (std::list<Geometry_3d::Point>::const_iterator i = face.vertices.begin(); i != face.vertices.end(); i++)
		{
			index = pointSet.insert(*i);
			pointSet.setData(index,0,1);
			faceCell << index;
		}

		index = pointSet.insert(face.vertices.back() + infiniteLength * face.d2);
		pointSet.setData(index,0,2);
		faceCell << index;

		cellSet.add(faceCell);
	}
	else
	{
		pv_cellData faceCell(POLYGON,1);
		faceCell.setData(0,1.0f);

		for (std::list<Geometry_3d::Point>::const_iterator i = face.vertices.begin(); i != face.vertices.end(); i++)
		{
			index = pointSet.insert(*i);
			pointSet.setData(index,0,1);
			faceCell << index;
		}

		cellSet.add(faceCell);
	}

	// *** GENERATE ASCII OUTPUT

	const char spacer = ' ';
	const char endline = '\n';

	std::cout << "*** Writing voronoi face in PARAVIEW format to file \"" << filename << "\" ... ";

	std::ofstream file(filename,std::ofstream::out|std::ofstream::trunc);
	file.setf(std::ios_base::scientific|std::ios_base::showpoint);

	// *** HEADER

	file << "# vtk DataFile Version 1.0" << endline;
	file << "Vorfinite-Tap => Voronoi Face: #" << face.n1 << " -> #" << face.n2 << endline;
	file << "ASCII" << endline;
	file << endline;

	file << "DATASET" << spacer << "UNSTRUCTURED_GRID" << endline;
	file << endline;

	// *** POINT DEFINITION

	file << "POINTS" << spacer << pointSet.size() << spacer << "float" << endline;
	const std::map<int,Geometry_3d::Point> sortedPoints = pointSet.points();
	for (std::map<int,Geometry_3d::Point>::const_iterator i = sortedPoints.begin(); i != sortedPoints.end(); i++)
	{
		file << i->second.x() << spacer;
		file << i->second.y() << spacer;
		file << i->second.z() << endline;
	}
	file << endline;

	// *** CELL DEFINITION (1/2)

	file << "CELLS" << spacer << cellSet.numCells() << spacer << cellSet.numCells() + cellSet.numNodes() << endline;
	for (std::list<pv_cellData>::const_iterator i = cellSet.cells().begin(); i != cellSet.cells().end(); i++)
	{
		if (i->nodes().size() == 0) continue;

		file << i->nodes().size();
		for (std::list<int>::const_iterator j = i->nodes().begin(); j != i->nodes().end(); j++)
			file << spacer << *j;
		file << endline;
	}
	file << endline;

	// *** CELL DEFINITION (2/2)

	file << "CELL_TYPES" << spacer << cellSet.numCells() << endline;
	for (std::list<pv_cellData>::const_iterator i = cellSet.cells().begin(); i != cellSet.cells().end(); i++)
		if (i->nodes().size() != 0) file << i->type() << endline;
	file << endline;

	// ***  POINT DATA

	file << "POINT_DATA" << spacer << pointSet.size() << endline;
	file << endline;

	file << "SCALARS" << spacer << "pointData" << spacer << "float" << spacer << 1 << endline;
	file << "LOOKUP_TABLE" << spacer << "default" << endline;
	for (std::map<int,float>::const_iterator i = pointSet.data(0).begin(); i != pointSet.data(0).end(); i++)
		file << i->second << endline;
	file << endline;

	// *** CELL DATA

	file << "CELL_DATA" << spacer << cellSet.numCells() << endline;
	file << endline;
	
	file << "SCALARS" << spacer << "cellData" << spacer << "float" << spacer << 1 << endline;
	file << "LOOKUP_TABLE" << spacer << "default" << spacer << endline;
	for (std::list<pv_cellData>::const_iterator i = cellSet.cells().begin(); i != cellSet.cells().end(); i++)
		file << i->data(0) << endline;
	file << endline;

	// ***

	std::cout << "finished." << std::endl;

	file.close();
}

void Debug::writeData_PARAVIEW(const char* filename, const System_3d& system, const Geometry_3d::VoronoiCell& obj)
{
	if (system.geomTable.numNodes() != system.gridTable.numNodes())
		throw std::runtime_error("Debug::writeData_PARAVIEW() / SINGLE CELL");

	const float infiniteLength = 0.01f * (system.geomTable.max() - system.geomTable.min()).length();

	const bool withGenerator(true);
	const bool withNeighbours(true);
	const bool withTetrahedra(false);

	// ***

	const float GENERATOR = 0.0f;
	const float DELAUNAY = 1.0f;
	const float FACE = 2.0f;

	int index;

	pv_pointSet pointSet(1);

	// *** generator

	pv_cellSet generator(2);

	if (withGenerator)
	{
		float volume = obj.volume();
		if (volume != volume || std::numeric_limits<float>::infinity() == volume) volume = -1.0f;

		// ***

		pv_cellData cell(POINT,2);
		cell.setData(0,GENERATOR);
		cell.setData(1,volume);

		index = pointSet.insert(obj.p);
		pointSet.setData(index,0,GENERATOR);
		cell << index;

		generator.add(cell);
	}

	// *** delaunay neighbours

	pv_cellSet neighbours(2);

	if (withNeighbours)
	{
		float volume = obj.volume();
		if (volume != volume || std::numeric_limits<float>::infinity() == volume) volume = -1.0f;

		// ***

		pv_cellData cell(POLYLINE,2);
		cell.setData(0,DELAUNAY);
		cell.setData(1,volume);

		for (std::list<Geometry_3d::VoronoiFace>::const_iterator i = obj.faces.begin(); i != obj.faces.end(); i++)
		{
			index = pointSet.insert(system.geomTable.nodeCoords(i->n1));
			pointSet.setData(index,0,GENERATOR);
			cell << index;

			index = pointSet.insert(system.geomTable.nodeCoords(i->n2));
			pointSet.setData(index,0,DELAUNAY);
			cell << index;
		}

		neighbours.add(cell);
	}

	// *** voronoi faces

	pv_cellSet voronoiFaces(2);

	for (std::list<Geometry_3d::VoronoiFace>::const_iterator i = obj.faces.begin(); i != obj.faces.end(); i++)
	{
		float area = i->area();
		if (area != area || std::numeric_limits<float>::infinity() == area) area = -1.0f;

		if (!i->isFinite())
		{
			pv_cellData face(POLYGON,2);
			face.setData(0,FACE);
			face.setData(1,area);
	
			index = pointSet.insert(i->vertices.front() + infiniteLength * i->d1);
			pointSet.setData(index,0,FACE);
			face << index;
	
			for (std::list<Geometry_3d::Point>::const_iterator j = i->vertices.begin(); j != i->vertices.end(); j++)
			{
				index = pointSet.insert(*j);
				pointSet.setData(index,0,FACE);
				face << index;
			}
	
			index = pointSet.insert(i->vertices.back() + infiniteLength * i->d2);
			pointSet.setData(index,0,FACE);
			face << index;
	
			voronoiFaces.add(face);
		}
		else
		{
			pv_cellData face(POLYGON,2);
			face.setData(0,FACE);
			face.setData(1,area);
	
			for (std::list<Geometry_3d::Point>::const_iterator j = i->vertices.begin(); j != i->vertices.end(); j++)
			{
				index = pointSet.insert(*j);
				pointSet.setData(index,0,FACE);
				face << index;
			}
	
			voronoiFaces.add(face);
		}
	}

	// *** delaunay tetrahedra

	pv_cellSet tetrahedra(2);

	if (withTetrahedra)
	{
		std::set<int> appearedTets;

		for (std::list<Geometry_3d::VoronoiFace>::const_iterator i = obj.faces.begin(); i != obj.faces.end(); i++)
		{
			std::list<int>::const_iterator j = i->tets.begin();

			int tetCnt(0);
			while (i->tets.end() != j)
			{
				if (appearedTets.count(*j) == 0)
				{
					appearedTets.insert(*j);

					pv_cellData tetCell(TETRAHEDRON,2);
					tetCell.setData(0,DELAUNAY);
					tetCell.setData(1,tetCnt++);
		
					for (int l = 0; l < 4; l++)
					{
						index = pointSet.insert(system.geomTable.nodeCoords(system.geomTable.tetrahedron(*j).vertices[l]));
						pointSet.setData(index,0,DELAUNAY);
						tetCell << index;
					}

					tetrahedra.add(tetCell);
				}
	
				j++;
			}
		}
	}

	// *** GENERATE ASCII OUTPUT

	const char spacer = ' ';
	const char endline = '\n';

	std::cout << "*** Writing voronoi cell in PARAVIEW format to file \"" << filename << "\" ... ";

	std::ofstream file(filename,std::ofstream::out|std::ofstream::trunc);
	file.setf(std::ios_base::scientific|std::ios_base::showpoint);

	// *** HEADER

	file << "# vtk DataFile Version 1.0" << endline;
	file << "Vorfinite-Tap => Voronoi Cell: #" << obj.n << endline;
	file << "ASCII" << endline;
	file << endline;

	file << "DATASET" << spacer << "UNSTRUCTURED_GRID" << endline;
	file << endline;

	// *** POINT DEFINITION

	file << "POINTS" << spacer << pointSet.size() << spacer << "float" << endline;
	const std::map<int,Geometry_3d::Point> sortedPoints = pointSet.points();
	for (std::map<int,Geometry_3d::Point>::const_iterator i = sortedPoints.begin(); i != sortedPoints.end(); i++)
	{
		file << i->second.x() << spacer;
		file << i->second.y() << spacer;
		file << i->second.z() << endline;
	}
	file << endline;

	// *** CELL DEFINITION (1/2)

	int numCells(0);
	numCells += generator.numCells();
	numCells += neighbours.numCells();
	numCells += voronoiFaces.numCells();
	numCells += tetrahedra.numCells();

	int numNodes(0);
	numNodes += generator.numNodes();
	numNodes += neighbours.numNodes();
	numNodes += voronoiFaces.numNodes();
	numNodes += tetrahedra.numNodes();

	const pv_cellSet* cellsets[4] = { &generator, &neighbours, &voronoiFaces, &tetrahedra };

	file << "CELLS" << spacer << numCells << spacer << numCells + numNodes << endline;
	for (int i = 0; i < 4; i++)
	{
		for (std::list<pv_cellData>::const_iterator j = cellsets[i]->cells().begin(); j != cellsets[i]->cells().end(); j++)
		{
			if (j->nodes().size() == 0) continue;
	
			file << j->nodes().size();
			for (std::list<int>::const_iterator k = j->nodes().begin(); k != j->nodes().end(); k++)
				file << spacer << *k;
			file << endline;
		}
	}
	file << endline;

	// *** CELL DEFINITION (2/2)

	file << "CELL_TYPES" << spacer << numCells << endline;
	for (int i = 0; i < 4; i++)
	{
		for (std::list<pv_cellData>::const_iterator j = cellsets[i]->cells().begin(); j != cellsets[i]->cells().end(); j++)
			if (j->nodes().size() != 0) file << j->type() << endline;
	}
	file << endline;

	// ***  POINT DATA

	file << "POINT_DATA" << spacer << pointSet.size() << endline;
	file << endline;

	file << "SCALARS" << spacer << "pointData" << spacer << "float" << spacer << 1 << endline;
	file << "LOOKUP_TABLE" << spacer << "default" << endline;
	for (std::map<int,float>::const_iterator i = pointSet.data(0).begin(); i != pointSet.data(0).end(); i++)
		file << i->second << endline;
	file << endline;

	// *** CELL DATA

	file << "CELL_DATA" << spacer << numCells << endline;
	file << endline;
	
	file << "SCALARS" << spacer << "cellData_1" << spacer << "float" << spacer << 1 << endline;
	file << "LOOKUP_TABLE" << spacer << "default" << spacer << endline;
	for (int i = 0; i < 4; i++)
	{
		for (std::list<pv_cellData>::const_iterator j = cellsets[i]->cells().begin(); j != cellsets[i]->cells().end(); j++)
			file << j->data(0) << endline;
	}
	file << endline;

	file << "SCALARS" << spacer << "cellData_2" << spacer << "float" << spacer << 1 << endline;
	file << "LOOKUP_TABLE" << spacer << "default" << spacer << endline;
	for (int i = 0; i < 4; i++)
	{
		for (std::list<pv_cellData>::const_iterator j = cellsets[i]->cells().begin(); j != cellsets[i]->cells().end(); j++)
			file << j->data(1) << endline;
	}
	file << endline;

	// ***

	std::cout << "finished." << std::endl;
	std::cout << std::endl;

	file.close();
}

// ***

void Debug::writeData_PARAVIEW_potentials(const char* filename, const System_3d& system)
{
	if (system.geomTable.numNodes() != system.gridTable.numNodes())
		throw std::runtime_error("Debug::writeData_PARAVIEW_potentials()");

	// ***

	std::cout << "*** Writing potential data in PARAVIEW format:" << std::endl;

	int index;

	pv_pointSet pointSet(1);

	// *** generators

	for (int i = 0; i < system.geomTable.numNodes(); i++)
	{
		if (i % (1+system.geomTable.numNodes()/100) == 0)
		{
			std::cout << "\t-> processing generators: " << i << " / " << system.geomTable.numNodes() << '\r';
			std::cout.flush();
		}

		// ***

		index = pointSet.insert(system.geomTable.nodeCoords(i));
		pointSet.setData(index,0,system.gridTable.potential(i));
	}

	std::cout << "\t-> processing generators: finished" << std::string(25,' ') << std::endl;

	// *** tetrahedra

	pv_cellSet tetrahedra(0);

	for (int i = 0; i < system.geomTable.numTetrahedra(); i++)
	{
		if (i % (1+system.geomTable.numTetrahedra()/100) == 0)
		{
			std::cout << "\t-> processing tetrahedra: " << i << " / " << system.geomTable.numTetrahedra() << '\r';
			std::cout.flush();
		}

		// ***

		Geometry_3d::Tetrahedron obj;

		for (int j = 0; j < 4; j++)
			obj[j] = system.geomTable.nodeCoords(system.geomTable.tetrahedron(i).vertices[j]);

		pv_cellData tetrahedron(TETRAHEDRON,0);

		for (int j = 0; j < 4; j++)
		{
			index = pointSet.insert(obj[j]);
			tetrahedron << index;
		}

		tetrahedra.add(tetrahedron);
	}

	std::cout << "\t-> processing tetrahedra: finished" << std::string(25,' ') << std::endl;

	// *** GENERATE ASCII OUTPUT

	const char spacer = ' ';
	const char endline = '\n';

	std::cout << "\t-> writing data to file \"" << filename << "\" ..." << '\r';
	std::cout.flush();

	std::ofstream file(filename,std::ofstream::out|std::ofstream::trunc);
	file.setf(std::ios_base::scientific|std::ios_base::showpoint);

	// *** HEADER

	file << "# vtk DataFile Version 1.0" << endline;
	file << "Vorfinite-Tap => potential distribution" << endline;
	file << "BINARY" << endline;

	file << "DATASET" << spacer << "UNSTRUCTURED_GRID" << endline;

	// *** POINT DEFINITION

	file << "POINTS" << spacer << pointSet.size() << spacer << "float" << endline;
	const std::map<int,Geometry_3d::Point> sortedPoints = pointSet.points();
	for (std::map<int,Geometry_3d::Point>::const_iterator i = sortedPoints.begin(); i != sortedPoints.end(); i++)
	{
		const float x = swapEndian(i->second.x());
		file.write(reinterpret_cast<const char*>(&x),sizeof(float));

		const float y = swapEndian(i->second.y());
		file.write(reinterpret_cast<const char*>(&y),sizeof(float));

		const float z = swapEndian(i->second.z());
		file.write(reinterpret_cast<const char*>(&z),sizeof(float));
	}
	file << endline;

	// *** CELL DEFINITION (1/2)

	file << "CELLS" << spacer << tetrahedra.numCells() << spacer << tetrahedra.numCells() + tetrahedra.numNodes() << endline;
	for (std::list<pv_cellData>::const_iterator j = tetrahedra.cells().begin(); j != tetrahedra.cells().end(); j++)
	{
		if (j->nodes().size() == 0) continue;

		const int size = swapEndian(static_cast<int>(j->nodes().size()));
		file.write(reinterpret_cast<const char*>(&size),sizeof(int));

		for (std::list<int>::const_iterator k = j->nodes().begin(); k != j->nodes().end(); k++)
		{
			const int node = swapEndian(*k);
			file.write(reinterpret_cast<const char*>(&node),sizeof(int));
		}
	}
	file << endline;

	// *** CELL DEFINITION (2/2)

	file << "CELL_TYPES" << spacer << tetrahedra.numCells() << endline;
	for (std::list<pv_cellData>::const_iterator j = tetrahedra.cells().begin(); j != tetrahedra.cells().end(); j++)
	{
		if (j->nodes().size() != 0)
		{
			const int type = swapEndian(j->type());
			file.write(reinterpret_cast<const char*>(&type),sizeof(int));
		}
	}
	file << endline;

	// ***  POINT DATA

	file << "POINT_DATA" << spacer << pointSet.size() << endline;

	file << "SCALARS" << spacer << "phi" << spacer << "float" << spacer << 1 << endline;
	file << "LOOKUP_TABLE" << spacer << "default" << endline;
	for (std::map<int,float>::const_iterator i = pointSet.data(0).begin(); i != pointSet.data(0).end(); i++)
	{
		const float phi = swapEndian(i->second);
		file.write(reinterpret_cast<const char*>(&phi),sizeof(float));
	}
	file << endline;

	// ***

	file.close();

	std::cout << "\t-> writing data to file \"" << filename << "\" (finished)" << std::endl;
	std::cout << std::endl;
}

void Debug::writeData_PARAVIEW_voronoiCells(const char* filename, const System_3d& system)
{
	if (system.geomTable.numNodes() != system.gridTable.numNodes())
		throw std::runtime_error("Debug::writeData_PARAVIEW() / VORONOI CELLS");

	const float infiniteLength = 0.01f * (system.geomTable.max() - system.geomTable.min()).length();

	// ***

	std::cout << "*** Writing voronoi cells in PARAVIEW format:" << std::endl;
	std::cout << "\t-> filename: \"" << filename << "\"" << std::endl;
	std::cout << "\t-> coordinate constraints: " << system.geomTable.min() << " - " << system.geomTable.max() << std::endl;

	int index;

	pv_pointSet pointSet(0);
	pv_cellSet voronoiCells(4);

	for (int nodeIndex = 0; nodeIndex < system.geomTable.numNodes(); nodeIndex++)
	{
		if (nodeIndex % (1+system.geomTable.numNodes()/100) == 0)
		{
			std::cout << "\t-> computing cell geometries: " << nodeIndex << "/" << system.geomTable.numNodes() << '\r';
			std::cout.flush();
		}

		Geometry_3d::VoronoiCell obj;

		system.geomTable.get_voronoiCell(nodeIndex,&obj);

		for (std::list<Geometry_3d::VoronoiFace>::const_iterator i = obj.faces.begin(); i != obj.faces.end(); i++)
		{
			// *** process each face only once

			if (nodeIndex > i->n2) continue;

			// ***

			float area = i->finiteArea(infiniteLength);

			float id;
			{
				short id1 = system.gridTable.node(nodeIndex).id().toValue();
				short id2 = system.gridTable.node(i->n2).id().toValue();

				if (id1 > id2)
					id = static_cast<float>(id1);
				else
					id = static_cast<float>(id2);
			}

			if (!i->isFinite())
			{
				pv_cellData face(POLYGON,4);

				face.setData(0,area);
				face.setData(1,id);

				float phi = system.gridTable.potential(nodeIndex);
				face.setData(2,phi);

				MathVector3d<float> field = system.gridTable.field_o1(nodeIndex,system.geomTable);
				face.setData(3,field.length());
		
				index = pointSet.insert(i->vertices.front() + infiniteLength * i->d1);
				face << index;
		
				for (std::list<Geometry_3d::Point>::const_iterator j = i->vertices.begin(); j != i->vertices.end(); j++)
				{
					index = pointSet.insert(*j);
					face << index;
				}
		
				index = pointSet.insert(i->vertices.back() + infiniteLength * i->d2);
				face << index;
		
				voronoiCells.add(face);
			}
			else
			{
				pv_cellData face(POLYGON,4);

				face.setData(0,area);
				face.setData(1,id);

				float phi = system.gridTable.potential(nodeIndex);
				phi += system.gridTable.potential(i->n2);
				phi /= 2.0f;

				face.setData(2,phi);

				MathVector3d<float> field = system.gridTable.field_o1(nodeIndex,system.geomTable);
				field += system.gridTable.field_o1(i->n2,system.geomTable);
				field /= 2.0f;

				face.setData(3,field.length());
		
				for (std::list<Geometry_3d::Point>::const_iterator j = i->vertices.begin(); j != i->vertices.end(); j++)
				{
					index = pointSet.insert(*j);
					face << index;
				}
		
				voronoiCells.add(face);
			}
		}
	}

	std::cout << "\t-> computing cell geometries: finished" << std::string(50,' ') << std::endl;

	// *** GENERATE ASCII OUTPUT

	const char spacer = ' ';
	const char endline = '\n';

	std::ofstream file(filename,std::ofstream::out|std::ofstream::trunc);
	file.setf(std::ios_base::scientific|std::ios_base::showpoint);

	std::cout << "\t-> writing data to file: processing\r";
	std::cout.flush();

	// *** HEADER

	file << "# vtk DataFile Version 1.0" << endline;
	file << "Vorfinite-Tap => Voronoi Cells" << endline;
	file << "ASCII" << endline;
	file << endline;

	file << "DATASET" << spacer << "UNSTRUCTURED_GRID" << endline;
	file << endline;

	// *** POINT DEFINITION

	file << "POINTS" << spacer << pointSet.size() << spacer << "float" << endline;
	const std::map<int,Geometry_3d::Point> sortedPoints = pointSet.points();
	for (std::map<int,Geometry_3d::Point>::const_iterator i = sortedPoints.begin(); i != sortedPoints.end(); i++)
	{
		file << i->second.x() << spacer;
		file << i->second.y() << spacer;
		file << i->second.z() << endline;
	}
	file << endline;

	// *** CELL DEFINITION (1/2)

	const int numCells = voronoiCells.numCells();
	const int numNodes = voronoiCells.numNodes();

	file << "CELLS" << spacer << numCells << spacer << numCells + numNodes << endline;
	for (std::list<pv_cellData>::const_iterator j = voronoiCells.cells().begin(); j != voronoiCells.cells().end(); j++)
	{
		if (j->nodes().size() == 0) continue;

		file << j->nodes().size();
		for (std::list<int>::const_iterator k = j->nodes().begin(); k != j->nodes().end(); k++)
			file << spacer << *k;
		file << endline;
	}
	file << endline;

	// *** CELL DEFINITION (2/2)

	file << "CELL_TYPES" << spacer << numCells << endline;
	for (std::list<pv_cellData>::const_iterator j = voronoiCells.cells().begin(); j != voronoiCells.cells().end(); j++)
		if (j->nodes().size() != 0) file << j->type() << endline;
	file << endline;

	// *** CELL DATA

	file << "CELL_DATA" << spacer << numCells << endline;
	file << endline;
	
	file << "SCALARS" << spacer << "area" << spacer << "float" << spacer << 1 << endline;
	file << "LOOKUP_TABLE" << spacer << "default" << spacer << endline;
	for (std::list<pv_cellData>::const_iterator j = voronoiCells.cells().begin(); j != voronoiCells.cells().end(); j++)
		file << j->data(0) << endline;
	file << endline;

	file << "SCALARS" << spacer << "id" << spacer << "int" << spacer << 1 << endline;
	file << "LOOKUP_TABLE" << spacer << "default" << spacer << endline;
	for (std::list<pv_cellData>::const_iterator j = voronoiCells.cells().begin(); j != voronoiCells.cells().end(); j++)
		file << static_cast<int>(j->data(1)) << endline;
	file << endline;

	file << "SCALARS" << spacer << "phi" << spacer << "float" << spacer << 1 << endline;
	file << "LOOKUP_TABLE" << spacer << "default" << spacer << endline;
	for (std::list<pv_cellData>::const_iterator j = voronoiCells.cells().begin(); j != voronoiCells.cells().end(); j++)
		file << j->data(2) << endline;
	file << endline;

	file << "SCALARS" << spacer << "field" << spacer << "float" << spacer << 1 << endline;
	file << "LOOKUP_TABLE" << spacer << "default" << spacer << endline;
	for (std::list<pv_cellData>::const_iterator j =voronoiCells.cells().begin(); j != voronoiCells.cells().end(); j++)
		file << j->data(3) << endline;
	file << endline;

	// ***

	file.close();

	std::cout << "\t-> writing data to file: finished" << std::string(25,' ') << std::endl;
	std::cout << std::endl;
}

/*void Debug::writeData_PARAVIEW_constrainedVoronoiCells(const char* filename, const System_3d& system)
{
	if (system.geomTable.numNodes() != system.gridTable.numNodes())
		throw std::runtime_error("Debug::writeData_PARAVIEW() / CONSTRAINED VORONOI CELLS");

	// ***

	std::cout << "*** Writing constrained voronoi cells in PARAVIEW format:" << std::endl;
	std::cout << "\t-> filename: \"" << filename << "\"" << std::endl;
	std::cout << "\t-> coordinate extents: " << system.geomTable.min() << " - " << system.geomTable.max() << std::endl;
	std::cout << "\t-> coordinate constraints: " << system.geomTable.constrainedMin() << " - " << system.geomTable.constrainedMax() << std::endl;

	int index;

	pv_pointSet pointSet(0);
	pv_cellSet voronoiCells(4);

	for (int nodeIndex = 0; nodeIndex < system.geomTable.numNodes(); nodeIndex++)
	{
		if (nodeIndex % (1+system.geomTable.numNodes()/100) == 0)
		{
			std::cout << "\t-> computing cell geometries: " << nodeIndex << "/" << system.geomTable.numNodes() << '\r';
			std::cout.flush();
		}

		Geometry_3d::VoronoiCell obj;

		system.geomTable.get_voronoiCell(nodeIndex,&obj);

		for (std::list<Geometry_3d::VoronoiFace>::const_iterator i = obj.faces.begin(); i != obj.faces.end(); i++)
		{
			// *** process each face only once

			if (nodeIndex > i->n2) continue;

			// ***

			float area = i->constrainedArea(system.geomTable.constrainedMax(),system.geomTable.constrainedMax());

			float id;
			{
				short id1 = system.gridTable.node(nodeIndex).id().toValue();
				short id2 = system.gridTable.node(i->n2).id().toValue();

				if (id1 > id2)
					id = static_cast<float>(id1);
				else
					id = static_cast<float>(id2);
			}

			if (!i->isFinite())
			{
				pv_cellData face(POLYGON,4);

				face.setData(0,area);
				face.setData(1,id);

				float phi = system.gridTable.potential(nodeIndex);
				face.setData(2,phi);

				//MathVector3d<float> field = system.gridTable.field_o1(nodeIndex,system.geomTable);

				MathVector3d<float> field(0.0f);
				face.setData(3,field.length());
		
				std::list<Geometry_3d::Point> points;
				i->constrainedVertices(system.geomTable.constrainedMax(),system.geomTable.constrainedMax(),&points);
		
				for (std::list<Geometry_3d::Point>::const_iterator j = points.begin(); j != points.end(); j++)
				{
					index = pointSet.insert(*j);
					face << index;
				}
		
				voronoiCells.add(face);
			}
			else
			{
				pv_cellData face(POLYGON,4);

				face.setData(0,area);
				face.setData(1,id);

				float phi = system.gridTable.potential(nodeIndex);
				phi += system.gridTable.potential(i->n2);
				phi /= 2.0f;

				face.setData(2,phi);

				//MathVector3d<float> field = system.gridTable.field_o1(nodeIndex,system.geomTable);
				//field += system.gridTable.field_o1(i->n2,system.geomTable);
				//field /= 2.0f;
				
				MathVector3d<float> field(0.0f);

				face.setData(3,field.length());
		
				for (std::list<Geometry_3d::Point>::const_iterator j = i->vertices.begin(); j != i->vertices.end(); j++)
				{
					index = pointSet.insert(*j);
					face << index;
				}
		
				voronoiCells.add(face);
			}
		}
	}

	std::cout << "\t-> computing cell geometries: finished" << std::string(50,' ') << std::endl;

	// *** GENERATE ASCII OUTPUT

	const char spacer = ' ';
	const char endline = '\n';

	std::ofstream file(filename,std::ofstream::out|std::ofstream::trunc);
	file.setf(std::ios_base::scientific|std::ios_base::showpoint);

	std::cout << "\t-> writing data to file: processing\r";
	std::cout.flush();

	// *** HEADER

	file << "# vtk DataFile Version 1.0" << endline;
	file << "Vorfinite-Tap => Voronoi Cells" << endline;
	file << "ASCII" << endline;
	file << endline;

	file << "DATASET" << spacer << "UNSTRUCTURED_GRID" << endline;
	file << endline;

	// *** POINT DEFINITION

	file << "POINTS" << spacer << pointSet.size() << spacer << "float" << endline;
	const std::map<int,Geometry_3d::Point> sortedPoints = pointSet.points();
	for (std::map<int,Geometry_3d::Point>::const_iterator i = sortedPoints.begin(); i != sortedPoints.end(); i++)
	{
		file << i->second.x() << spacer;
		file << i->second.y() << spacer;
		file << i->second.z() << endline;
	}
	file << endline;

	// *** CELL DEFINITION (1/2)

	const int numCells = voronoiCells.numCells();
	const int numNodes = voronoiCells.numNodes();

	file << "CELLS" << spacer << numCells << spacer << numCells + numNodes << endline;
	for (std::list<pv_cellData>::const_iterator j = voronoiCells.cells().begin(); j != voronoiCells.cells().end(); j++)
	{
		if (j->nodes().size() == 0) continue;

		file << j->nodes().size();
		for (std::list<int>::const_iterator k = j->nodes().begin(); k != j->nodes().end(); k++)
			file << spacer << *k;
		file << endline;
	}
	file << endline;

	// *** CELL DEFINITION (2/2)

	file << "CELL_TYPES" << spacer << numCells << endline;
	for (std::list<pv_cellData>::const_iterator j = voronoiCells.cells().begin(); j != voronoiCells.cells().end(); j++)
		if (j->nodes().size() != 0) file << j->type() << endline;
	file << endline;

	// *** CELL DATA

	file << "CELL_DATA" << spacer << numCells << endline;
	file << endline;
	
	file << "SCALARS" << spacer << "area" << spacer << "float" << spacer << 1 << endline;
	file << "LOOKUP_TABLE" << spacer << "default" << spacer << endline;
	for (std::list<pv_cellData>::const_iterator j = voronoiCells.cells().begin(); j != voronoiCells.cells().end(); j++)
		file << j->data(0) << endline;
	file << endline;

	file << "SCALARS" << spacer << "id" << spacer << "int" << spacer << 1 << endline;
	file << "LOOKUP_TABLE" << spacer << "default" << spacer << endline;
	for (std::list<pv_cellData>::const_iterator j = voronoiCells.cells().begin(); j != voronoiCells.cells().end(); j++)
		file << static_cast<int>(j->data(1)) << endline;
	file << endline;

	file << "SCALARS" << spacer << "phi" << spacer << "float" << spacer << 1 << endline;
	file << "LOOKUP_TABLE" << spacer << "default" << spacer << endline;
	for (std::list<pv_cellData>::const_iterator j = voronoiCells.cells().begin(); j != voronoiCells.cells().end(); j++)
		file << j->data(2) << endline;
	file << endline;

	file << "SCALARS" << spacer << "field" << spacer << "float" << spacer << 1 << endline;
	file << "LOOKUP_TABLE" << spacer << "default" << spacer << endline;
	for (std::list<pv_cellData>::const_iterator j =voronoiCells.cells().begin(); j != voronoiCells.cells().end(); j++)
		file << j->data(3) << endline;
	file << endline;

	// ***

	file.close();

	std::cout << "\t-> writing data to file: finished" << std::string(25,' ') << std::endl;
	std::cout << std::endl;
}*/

// ***

void Debug::writeData_PARAVIEW_surfaceNodes(const char* filename, const Surface_3d::Table& surfaceTable, const System_3d& system)
{
	if (system.geomTable.numNodes() != system.gridTable.numNodes())
		throw std::runtime_error("Debug::writeData_PARAVIEW_surfaceNodes()");

	std::cout << "*** Writing surface node-data in PARAVIEW format to file \"" << filename << "\":" << std::endl;

	const bool withNodes(true);
	const bool withTets(false);

	int progressCnt;
	int index;

	// ***

	pv_pointSet points(7);

	pv_cellSet nodes(0);
	pv_cellSet tetrahedra(0);

	progressCnt = 0;

	std::set<int> delaunayTets;
	for (Surface_3d::Nodeset::const_iterator i = surfaceTable.nodes().begin(); i != surfaceTable.nodes().end(); i++)
	{
		if (++progressCnt % (1+surfaceTable.nodes().size()/100) == 0)
		{
			std::cout << "\t-> collecting tetrahedra: " << progressCnt << " / " << surfaceTable.nodes().size() << "\r";
			std::cout.flush();
		}

		std::set<int> commonTets;
		system.geomTable.tetsWithCommonNode(i->index(),&commonTets);

		delaunayTets.insert(commonTets.begin(),commonTets.end());
	}

	std::cout << "\t-> collecting tetrahedra: finished" << std::string(25,' ') << std::endl;

	progressCnt = 0;
	for (std::set<int>::const_iterator i = delaunayTets.begin(); i != delaunayTets.end(); i++)
	{
		if (++progressCnt % (1+delaunayTets.size()/100) == 0)
		{
			std::cout << "\t-> processing tetrahedra: " << progressCnt << " / " << delaunayTets.size() << "\r";
			std::cout.flush();
		}

		pv_cellData tetCell(TETRAHEDRON,0);
	
		for (int j = 0; j < 4; j++)
		{
			const int tetVertex = system.geomTable.tetrahedron(*i).vertices[j];

			index = points.insert(system.geomTable.nodeCoords(tetVertex));

			points.setData(index,0,static_cast<float>(system.gridTable.node(tetVertex).id().toValue()));
			points.setData(index,1,system.gridTable.potential(tetVertex));
			points.setData(index,2,system.gridTable.flux(tetVertex,system.geomTable));

			float probability(0.0f);
			for (Surface_3d::Nodeset::const_iterator i = surfaceTable.nodes().begin(); i != surfaceTable.nodes().end(); i++)
			{
				if (i->index() == tetVertex)
				{
					probability = i->probability();
					break;
				}
			}

			points.setData(index,3,probability);

			const MathVector3d<float> fieldValue = system.gridTable.field_o1(tetVertex,system.geomTable);

			points.setData(index,4,fieldValue.x());
			points.setData(index,5,fieldValue.y());
			points.setData(index,6,fieldValue.z());

			// ***

			if (withNodes)
			{
				pv_cellData nodeCell(POINT,0);
				nodeCell <<  index;

				nodes.add(nodeCell);
			}

			if (withTets) tetCell << index;
		}

		if (withTets) tetrahedra.add(tetCell);
	}

	std::cout << "\t-> processing tetrahedra: finished" << std::string(25,' ') << std::endl;

	// *** GENERATE ASCII OUTPUT

	const char spacer = ' ';
	const char endline = '\n';

	std::ofstream file(filename,std::ofstream::out|std::ofstream::trunc);
	file.setf(std::ios_base::scientific|std::ios_base::showpoint);

	std::cout << "\t-> writing data: processing" << '\r';
	std::cout.flush();

	// *** HEADER

	file << "# vtk DataFile Version 1.0" << endline;
	file << "Vorfinite-Tap => Surface Nodes" << endline;
	file << "ASCII" << endline;
	file << endline;

	file << "DATASET" << spacer << "UNSTRUCTURED_GRID" << endline;
	file << endline;

	// *** POINT DEFINITION

	file << "POINTS" << spacer << points.size() << spacer << "float" << endline;
	const std::map<int,Geometry_3d::Point> sortedPoints = points.points();
	for (std::map<int,Geometry_3d::Point>::const_iterator i = sortedPoints.begin(); i != sortedPoints.end(); i++)
	{
		file << i->second.x() << spacer;
		file << i->second.y() << spacer;
		file << i->second.z() << endline;
	}
	file << endline;

	// *** CELL DEFINITION (1/2)

	const int numCells = tetrahedra.numCells() + nodes.numCells();
	const int numNodes = tetrahedra.numNodes() + nodes.numNodes() + numCells;

	file << "CELLS" << spacer <<  numCells << spacer << numCells + numNodes << endline;

	for (std::list<pv_cellData>::const_iterator i = tetrahedra.cells().begin(); i != tetrahedra.cells().end(); i++)
	{
		if (i->nodes().size() == 0) continue;

		file << i->nodes().size();
		for (std::list<int>::const_iterator j = i->nodes().begin(); j != i->nodes().end(); j++)
			file << spacer << *j;
		file << endline;
	}
	for (std::list<pv_cellData>::const_iterator i = nodes.cells().begin(); i != nodes.cells().end(); i++)
	{
		if (i->nodes().size() == 0) continue;

		file << i->nodes().size();
		for (std::list<int>::const_iterator j = i->nodes().begin(); j != i->nodes().end(); j++)
			file << spacer << *j;
		file << endline;
	}
	file << endline;

	// *** CELL DEFINITION (2/2)

	file << "CELL_TYPES" << spacer << numCells << endline;
	for (std::list<pv_cellData>::const_iterator j = tetrahedra.cells().begin(); j != tetrahedra.cells().end(); j++)
		if (j->nodes().size() != 0) file << j->type() << endline;
	for (std::list<pv_cellData>::const_iterator j = nodes.cells().begin(); j != nodes.cells().end(); j++)
		if (j->nodes().size() != 0) file << j->type() << endline;
	file << endline;

	// ***  POINT DATA

	file << "POINT_DATA" << spacer << points.size() << endline;
	file << endline;

	file << "SCALARS" << spacer << "id" << spacer << "int" << spacer << 1 << endline;
	file << "LOOKUP_TABLE" << spacer << "default" << endline;
	for (std::map<int,float>::const_iterator i = points.data(0).begin(); i != points.data(0).end(); i++)
		file << static_cast<int>(i->second) << endline;
	file << endline;

	file << "SCALARS" << spacer << "potential" << spacer << "float" << spacer << 1 << endline;
	file << "LOOKUP_TABLE" << spacer << "default" << endline;
	for (std::map<int,float>::const_iterator i = points.data(1).begin(); i != points.data(1).end(); i++)
		file << i->second << endline;
	file << endline;

	file << "SCALARS" << spacer << "flux" << spacer << "float" << spacer << 1 << endline;
	file << "LOOKUP_TABLE" << spacer << "default" << endline;
	for (std::map<int,float>::const_iterator i = points.data(2).begin(); i != points.data(2).end(); i++)
		file << i->second << endline;
	file << endline;

	file << "SCALARS" << spacer << "probability" << spacer << "float" << spacer << 1 << endline;
	file << "LOOKUP_TABLE" << spacer << "default" << endline;
	for (std::map<int,float>::const_iterator i = points.data(3).begin(); i != points.data(3).end(); i++)
		file << i->second << endline;
	file << endline;

	file << "VECTORS" << spacer << "field" << spacer << "float" << endline;
	{
		std::map<int,float>::const_iterator i = points.data(4).begin();
		std::map<int,float>::const_iterator j = points.data(5).begin();
		std::map<int,float>::const_iterator k = points.data(6).begin();

		while (points.data(4).end() != i)
		{
			file << i->second << spacer;
			file << j->second << spacer;
			file << k->second << endline;

			i++, j++, k++;
		}
	}
	file << endline;

	// ***

	std::cout << "\t-> writing data: finished" << std::string(25,' ') << std::endl;
	std::cout << std::endl;
}

void Debug::writeData_PARAVIEW_surfaceCells(const char* filename, const Surface_3d::Table& surfaceTable, const System_3d& system)
{
	float infiniteLength;

	{
		if (surfaceTable.nodes().empty()) throw std::runtime_error("Debug::writeData_PARAVIEW_surfaceCells()");

		Geometry_3d::Point min = system.geomTable.nodeCoords(surfaceTable.nodes().begin()->index());
		Geometry_3d::Point max = system.geomTable.nodeCoords(surfaceTable.nodes().begin()->index());

		for (Surface_3d::Nodeset::const_iterator i = surfaceTable.nodes().begin(); i != surfaceTable.nodes().end(); i++)
		{
			const Geometry_3d::Point& coords = system.geomTable.nodeCoords(i->index());

			for (int j = 0; j < 3; j++)
			{
				if (coords[j] < min[j])
					min[j] = coords[j];
				else
					if (coords[j] > max[j]) max[j] = coords[j];
			}
		}

		infiniteLength = 0.01f * (max - min).length();
	}

	// ***

	std::cout << "*** Writing surface cell-data in PARAVIEW format to file \"" << filename << "\":" << std::endl;

	int progressCnt(0);

	int index;

	pv_pointSet pointSet(0);
	pv_cellSet surfaceCells(7);

	for (Surface_3d::Nodeset::const_iterator i = surfaceTable.nodes().begin(); i != surfaceTable.nodes().end(); i++)
	{
		if (++progressCnt % (1+surfaceTable.nodes().size()/100) == 0)
		{
			std::cout << "\t-> processing cells: " << progressCnt << " / " << surfaceTable.nodes().size() << '\r';
			std::cout.flush();
		}
	
		Geometry_3d::VoronoiCell obj;
		system.geomTable.get_voronoiCell(i->index(),&obj);

		for (std::list<Geometry_3d::VoronoiFace>::const_iterator j = obj.faces.begin(); j != obj.faces.end(); j++)
		{
			pv_cellData face(POLYGON,7);

			const int nodeId = system.gridTable.node(i->index()).id().toValue();

			face.setData(0,static_cast<float>(nodeId));
			face.setData(1,i->probability());

			float phi = system.gridTable.potential(j->n1);
			phi += system.gridTable.potential(j->n2);
			phi /= 2.0f;
			face.setData(2,phi);

			float flux = system.gridTable.flux(j->n1,system.geomTable);
			flux += system.gridTable.flux(j->n2,system.geomTable);
			flux /= 2.0f;
			face.setData(3,flux);
	
			MathVector3d<float> field = system.gridTable.field_o1(j->n1,system.geomTable);
			field += system.gridTable.field_o1(j->n2,system.geomTable);
			field /= 2.0f;

			// ***

			for (int k = 0; k < 3; k++)
				face.setData(4+k,field[k]);

			if (!j->isFinite())
			{
				index = pointSet.insert(j->vertices.front() + infiniteLength * j->d1);
				face << index;
		
				for (std::list<Geometry_3d::Point>::const_iterator k = j->vertices.begin(); k != j->vertices.end(); k++)
				{
					index = pointSet.insert(*k);
					face << index;
				}
		
				index = pointSet.insert(j->vertices.back() + infiniteLength * j->d2);
				face << index;
			}
			else
			{
				for (std::list<Geometry_3d::Point>::const_iterator k = j->vertices.begin(); k != j->vertices.end(); k++)
				{
					index = pointSet.insert(*k);
					face << index;
				}
			}

			surfaceCells.add(face);
		}
	}

	std::cout << "\t-> processing cells: finished" << std::string(25,' ') << std::endl;

	// *** GENERATE ASCII OUTPUT

	const char spacer = ' ';
	const char endline = '\n';

	std::cout << "\t-> writing data: processing" << '\r';

	std::ofstream file(filename,std::ofstream::out|std::ofstream::trunc);
	file.setf(std::ios_base::scientific|std::ios_base::showpoint);

	// *** HEADER

	file << "# vtk DataFile Version 1.0" << endline;
	file << "Vorfinite-Tap => Voronoi Surface Cells" << endline;
	file << "ASCII" << endline;
	file << endline;

	file << "DATASET" << spacer << "UNSTRUCTURED_GRID" << endline;
	file << endline;

	// *** POINT DEFINITION

	file << "POINTS" << spacer << pointSet.size() << spacer << "float" << endline;
	const std::map<int,Geometry_3d::Point> sortedPoints = pointSet.points();
	for (std::map<int,Geometry_3d::Point>::const_iterator i = sortedPoints.begin(); i != sortedPoints.end(); i++)
	{
		file << i->second.x() << spacer;
		file << i->second.y() << spacer;
		file << i->second.z() << endline;
	}
	file << endline;

	// *** CELL DEFINITION (1/2)

	const int numCells = surfaceCells.numCells();
	const int numNodes = surfaceCells.numNodes();

	file << "CELLS" << spacer << numCells << spacer << numCells + numNodes << endline;
	for (std::list<pv_cellData>::const_iterator j = surfaceCells.cells().begin(); j != surfaceCells.cells().end(); j++)
	{
		if (j->nodes().size() == 0) continue;

		file << j->nodes().size();
		for (std::list<int>::const_iterator k = j->nodes().begin(); k != j->nodes().end(); k++)
			file << spacer << *k;
		file << endline;
	}
	file << endline;

	// *** CELL DEFINITION (2/2)

	file << "CELL_TYPES" << spacer << numCells << endline;
	for (std::list<pv_cellData>::const_iterator j = surfaceCells.cells().begin(); j != surfaceCells.cells().end(); j++)
		if (j->nodes().size() != 0) file << j->type() << endline;
	file << endline;

	// *** CELL DATA

	file << "CELL_DATA" << spacer << numCells << endline;
	file << endline;
	
	file << "SCALARS" << spacer << "id" << spacer << "int" << spacer << 1 << endline;
	file << "LOOKUP_TABLE" << spacer << "default" << spacer << endline;
	for (std::list<pv_cellData>::const_iterator j = surfaceCells.cells().begin(); j != surfaceCells.cells().end(); j++)
		file << static_cast<int>(j->data(0)) << endline;
	file << endline;

	file << "SCALARS" << spacer << "probability" << spacer << "float" << spacer << 1 << endline;
	file << "LOOKUP_TABLE" << spacer << "default" << spacer << endline;
	for (std::list<pv_cellData>::const_iterator j = surfaceCells.cells().begin(); j != surfaceCells.cells().end(); j++)
		file << j->data(1) << endline;
	file << endline;

	file << "SCALARS" << spacer << "phi" << spacer << "float" << spacer << 1 << endline;
	file << "LOOKUP_TABLE" << spacer << "default" << spacer << endline;
	for (std::list<pv_cellData>::const_iterator j = surfaceCells.cells().begin(); j != surfaceCells.cells().end(); j++)
		file << j->data(2) << endline;
	file << endline;

	file << "SCALARS" << spacer << "flux" << spacer << "float" << spacer << 1 << endline;
	file << "LOOKUP_TABLE" << spacer << "default" << spacer << endline;
	for (std::list<pv_cellData>::const_iterator j = surfaceCells.cells().begin(); j != surfaceCells.cells().end(); j++)
		file << j->data(3) << endline;
	file << endline;

	file << "VECTORS" << spacer << "field" << spacer << "float" << std::endl;
	for (std::list<pv_cellData>::const_iterator j = surfaceCells.cells().begin(); j != surfaceCells.cells().end(); j++)
	{
		file << j->data(4) << spacer;
		file << j->data(5) << spacer;
		file << j->data(6) << endline;
	}
	file << endline;

	// ***

	std::cout << "\t-> writing data: finished" << std::string(25,' ') << std::endl;
	std::cout << std::endl;

	file.close();
}

void Debug::writeData_PARAVIEW_surfaceFieldlines(const char* filename, const Surface_3d::Table& surfaceTable, const System_3d& system, const float length, const float mask)
{
	std::cout << "*** Writing field lines for given surface in PARAVIEW format to file \"" << filename << "\":" << std::endl;;

	// *** compute field lines

	std::list<Fieldline_3d> fieldlines;
	std::list<float> probabilities;

	{
		std::cout << "\t-> mask-factor = " << mask << std::endl;

		// ***

		int progressCnt(0);

		std::cout << "\t-> computing field lines (this may take a very long time): ";
		std::cout << progressCnt << " / " << surfaceTable.nodes().size() << std::string(10, ' ') << '\r';
		std::cout.flush();

		std::map<int,int> statistics;

		std::srand(1); // ***
		
		for (Surface_3d::Nodeset::const_iterator i = surfaceTable.nodes().begin(); i != surfaceTable.nodes().end(); i++)
		{
			if (++progressCnt % 1 == 0)
			{
				std::cout << "\t-> computing fieldlines (this may take a very long time): ";
				std::cout << progressCnt << " / " << surfaceTable.nodes().size() << std::string(10,' ') << '\r';
				std::cout.flush();
			}

			//if (std::rand()/static_cast<float>(RAND_MAX) > mask) continue;
			
			Geometry_3d::Point position = system.geomTable.nodeCoords(i->index());
			if (position.x() < -0.2e-9f || position.x() > 0.2e-9f) continue;
	
			// ***

			Fieldline_3d obj(&system.geomTable,&system.gridTable);
	
			obj.init(system.geomTable.nodeCoords(i->index()));
			obj.compute(length);

			statistics[obj.status()]++;
	
			fieldlines.push_back(obj);
			
			probabilities.push_back(i->probability());
		}

		std::cout << "\t-> computing fieldlines: finished";
		std::cout << std::string(100,' ') << std::endl;
	
		std::cout << "\t-> status statistics: " << std::endl;
		for (std::map<int,int>::const_iterator i = statistics.begin(); i != statistics.end(); i++)
			std::cout << "\t\t=> field lines with status #(" << i->first << "): " << i->second << std::endl;
	}

	// *** GENERATE ASCII OUTPUT

	std::cout << " \t-> writing output:";

	enum { VTK_VERTEX =1, VTK_POLY_LINE =4};

	const char spacer = ' ';
	const char endline = '\n';

	const int cellType = VTK_POLY_LINE;

	std::ofstream file(filename,std::ofstream::out|std::ofstream::trunc);
	file.setf(std::ios_base::scientific|std::ios_base::showpoint);

	// *** HEADER

	file << "# vtk DataFile Version 1.0" << endline;
	file << "TAPSIM => Surface fieldlines" << endline;
	file << "ASCII" << endline;
	file << endline;

	file << "DATASET" << spacer << "UNSTRUCTURED_GRID" << endline;
	file << endline;

	// *** POINT DEFINITION

	int pointNum(0);
	for (std::list<Fieldline_3d>::const_iterator i = fieldlines.begin(); i != fieldlines.end(); i++)
		pointNum += i->data().size();

	file << "POINTS" << spacer << pointNum << spacer << "float" << endline;
	for (std::list<Fieldline_3d>::const_iterator i = fieldlines.begin(); i != fieldlines.end(); i++)
	{
		for (std::vector<Fieldline_3d::phaseVector>::const_iterator j = i->data().begin(); j != i->data().end(); j++)
		{
			file << j->position().x() << spacer;
			file << j->position().y() << spacer;
			file << j->position().z() << endline;
		}
	}
	file << endline;

	// *** CELL DEFINITION (1/2)

	int index(0);

	file << "CELLS" << spacer << fieldlines.size() << spacer << fieldlines.size() + pointNum << endline;
	for (std::list<Fieldline_3d>::const_iterator i = fieldlines.begin(); i != fieldlines.end(); i++)
	{
		file << i->data().size();
		for (std::vector<Fieldline_3d::phaseVector>::const_iterator j = i->data().begin(); j != i->data().end(); j++)
			file << spacer << index++;
		file << endline;
	}
	file << endline;

	// *** CELL DEFINITION (2/2)

	file << "CELL_TYPES" << spacer << fieldlines.size() << endline;
	for (std::list<Fieldline_3d>::const_iterator i = fieldlines.begin(); i != fieldlines.end(); i++)
		file << cellType << endline;
	file << endline;

	// ***  POINT DATA

	file << "POINT_DATA" << spacer << pointNum << endline;
	file << endline;

	file << "VECTORS" << spacer << "field" << spacer << "float" << std::endl;
	for (std::list<Fieldline_3d>::const_iterator i = fieldlines.begin(); i != fieldlines.end(); i++)
	{
		for (std::vector<Fieldline_3d::phaseVector>::const_iterator j = i->data().begin(); j != i->data().end(); j++)
		{
			file << j->field().x() << spacer;
			file << j->field().y() << spacer;
			file << j->field().z() << endline;
		}
	}
	file << endline;

	// *** CELL DATA

	file << "CELL_DATA" << spacer << fieldlines.size() << endline;
	file << endline;
	
	file << "SCALARS" << spacer << "probability" << spacer << "float" << spacer << 1 << endline;
	file << "LOOKUP_TABLE" << spacer << "default" << spacer << endline;
	for (std::list<float>::const_iterator i = probabilities.begin(); i != probabilities.end(); i++)
		file << *i << endline;
	file << endline;

	file << "SCALARS" << spacer << "status" << spacer << "int" << spacer << 1 << endline;
	file << "LOOKUP_TABLE" << spacer << "default" << spacer << endline;
	for (std::list<Fieldline_3d>::const_iterator i = fieldlines.begin(); i != fieldlines.end(); i++)
		file << i->status() << endline;
	file << endline;

	// ***

	std::cout << " finished." << std::endl;

	file.close();
}

void Debug::writeData_PARAVIEW_surfaceTrajectories(const char* filename, const Surface_3d::Table& surfaceTable, const System_3d& system, const float distance, const float accuracy, const float mask)
{
	std::cout << "*** Writing trajectories for given surface in PARAVIEW format to file \"" << filename << "\":" << std::endl;;

	// *** compute trajectories

	std::list<Trajectory_3d> trajectories;
	std::list<float> probabilities;

	{
		std::cout << "\t-> accuracy = " << accuracy << std::endl;
		std::cout << "\t-> mask-factor = " << mask << std::endl;

		// ***

		int progressCnt(0);

		int noFieldCnt(0);

		std::cout << "\t-> computing trajectories (this may take a very long time): ";
		std::cout << progressCnt << " / " << surfaceTable.nodes().size() << std::string(10, ' ') << '\r';
		std::cout.flush();

		std::map<int,int> statistics;

		std::srand(1); // ***
		
		for (Surface_3d::Nodeset::const_iterator i = surfaceTable.nodes().begin(); i != surfaceTable.nodes().end(); i++)
		{
			if (++progressCnt % 1 == 0)
			{
				std::cout << "\t-> computing trajectories (this may take a very long time): ";
				std::cout << progressCnt << " / " << surfaceTable.nodes().size() << std::string(10,' ') << '\r';
				std::cout.flush();
			}

			//if (std::rand()/static_cast<float>(RAND_MAX) > mask) continue;
			
			Geometry_3d::Point position = system.geomTable.nodeCoords(i->index());
			if (position.x() < -0.2e-9f || position.x() > 0.2e-9f) continue;
	
			// ***

			probabilities.push_back(i->probability());

			// ***

			Trajectory_3d obj(&system.geomTable,&system.gridTable);
	
			obj.setIntegratorType(Trajectory_3d::O5_INTEGRATOR);
			obj.setStepperType(Trajectory_3d::ERROR_RESTRICTED);
			obj.setTimeStepLimit(1e-16);

			obj.init(system.geomTable.nodeCoords(i->index()),Trajectory_3d::Vector3d(0.0f),Trajectory_3d::eCharge,63.4f * Trajectory_3d::amu2kg);
			obj.integrate(1e-14,accuracy,distance);

			statistics[obj.status()]++;

			//statistics[obj.status()-Trajectory_3d::INVALID]++;
	
			trajectories.push_back(obj);
		}

		std::cout << "\t-> computing trajectories: finished";
		std::cout << std::string(100,' ') << std::endl;
	
		std::cout << "\t-> status statistics: " << std::endl;
		for (std::map<int,int>::const_iterator i = statistics.begin(); i != statistics.end(); i++)
			std::cout << "\t\t=> trajectories with status #(" << i->first << "): " << i->second << std::endl;

		std::cout << "NO FIELD COUNT => " << noFieldCnt << std::endl;
	}

	// *** GENERATE ASCII OUTPUT

	std::cout << " \t-> writing output:";

	enum { VTK_POLY_LINE =4, VTK_QUADRATIC_EDGE =21 };

	const char spacer = ' ';
	const char endline = '\n';

	const int cellType = VTK_POLY_LINE;

	std::ofstream file(filename,std::ofstream::out|std::ofstream::trunc);
	file.setf(std::ios_base::scientific|std::ios_base::showpoint);

	// *** HEADER

	file << "# vtk DataFile Version 1.0" << endline;
	file << "Vorfinite-Tap => Trajectories computed for each node of a given surface" << endline;
	file << "ASCII" << endline;
	file << endline;

	file << "DATASET" << spacer << "UNSTRUCTURED_GRID" << endline;
	file << endline;

	// *** POINT DEFINITION

	int pointNum(0);
	for (std::list<Trajectory_3d>::const_iterator i = trajectories.begin(); i != trajectories.end(); i++)
		pointNum += i->data().size();

	file << "POINTS" << spacer << pointNum << spacer << "float" << endline;
	for (std::list<Trajectory_3d>::const_iterator i = trajectories.begin(); i != trajectories.end(); i++)
	{
		for (std::vector<Trajectory_3d::phaseVector>::const_iterator j = i->data().begin(); j != i->data().end(); j++)
		{
			file << j->position(0) << spacer;
			file << j->position(1) << spacer;
			file << j->position(2) << endline;
		}
	}
	file << endline;

	// *** CELL DEFINITION (1/2)

	int index(0);

	file << "CELLS" << spacer << trajectories.size() << spacer << trajectories.size() + pointNum << endline;
	for (std::list<Trajectory_3d>::const_iterator i = trajectories.begin(); i != trajectories.end(); i++)
	{
		file << i->data().size();
		for (std::vector<Trajectory_3d::phaseVector>::const_iterator j = i->data().begin(); j != i->data().end(); j++)
			file << spacer << index++;
		file << endline;
	}
	file << endline;

	// *** CELL DEFINITION (2/2)

	file << "CELL_TYPES" << spacer << trajectories.size() << endline;
	for (std::list<Trajectory_3d>::const_iterator i = trajectories.begin(); i != trajectories.end(); i++)
		file << cellType << endline;
	file << endline;

	// ***  POINT DATA

	file << "POINT_DATA" << spacer << pointNum << endline;
	file << endline;

	file << "VECTORS" << spacer << "velocity" << spacer << "float" << std::endl;
	for (std::list<Trajectory_3d>::const_iterator i = trajectories.begin(); i != trajectories.end(); i++)
	{
		for (std::vector<Trajectory_3d::phaseVector>::const_iterator j = i->data().begin(); j != i->data().end(); j++)
		{
			file << j->velocity(0) << spacer;
			file << j->velocity(1) << spacer;
			file << j->velocity(2) << endline;
		}
	}
	file << endline;

	// *** CELL DATA

	file << "CELL_DATA" << spacer << trajectories.size() << endline;
	file << endline;
	
	file << "SCALARS" << spacer << "probability" << spacer << "float" << spacer << 1 << endline;
	file << "LOOKUP_TABLE" << spacer << "default" << spacer << endline;
	for (std::list<float>::const_iterator i = probabilities.begin(); i != probabilities.end(); i++)
		file << *i << endline;
	file << endline;

	file << "SCALARS" << spacer << "status" << spacer << "int" << spacer << 1 << endline;
	file << "LOOKUP_TABLE" << spacer << "default" << spacer << endline;
	for (std::list<Trajectory_3d>::const_iterator i = trajectories.begin(); i != trajectories.end(); i++)
		file << i->status() << endline;
	file << endline;

	// ***

	std::cout << " finished." << std::endl;

	file.close();
}

void Debug::writeData(const char* filename, const System_3d& system, const char spacer, const char endline)
{
	std::cout << "*** Writing node-data to file \"" << filename << "\":" << std::endl;

	std::ofstream file(filename,std::ofstream::out|std::ofstream::trunc);
	file.setf(std::ios_base::scientific|std::ios_base::showpoint);

	// *** header

	file << "# x" << spacer;
	file << "y" << spacer;
	file << "z" << spacer;
	file << "id" << spacer;
	file << "phi" << spacer;
	file << "fieldX" << spacer;
	file << "fieldY" << spacer;
	file << "fieldZ" << spacer;
	file << "fieldStrength" << endline;
	file << endline;

	// *** contents

	for (int i = 0; i < system.geomTable.numNodes(); i++)
	{
		if (i % (1 + system.geomTable.numNodes() / 100) == 0)
		{
			std::cout << "\t -> processing nodes: " << i << " / " << system.geomTable.numNodes();
			std::cout << std::string(15, ' ') << '\r';
			std::cout.flush();
		}

		file << system.geomTable.nodeCoords(i).x() << spacer;
		file << system.geomTable.nodeCoords(i).y() << spacer;
		file << system.geomTable.nodeCoords(i).z() << spacer;
		file << system.gridTable.node(i).id().toValue() << spacer;
		file << system.gridTable.potential(i) << spacer;
		
		const MathVector3d<float> field = system.gridTable.field_o1(i,system.geomTable);

		file << field.x() << spacer;
		file << field.y() << spacer;
		file << field.z() << spacer;
		file << field.length() << endline;
	}

	std::cout << "\t -> processing nodes: finished" << std::string(25,' ') << std::endl;
	std::cout << std::endl;
}

void Debug::writeData_surfaceTrajectories(const char* filename, const Surface_3d::Table& surfaceTable, const System_3d& system, const float accuracy, const float mask)
{
// 	std::cout << "*** Writing trajectories to file \"" << filename << "\":" << std::endl;
// 	std::cout << "\t-> accuracy = " << accuracy << std::endl;
// 	std::cout << "\t-> mask-factor = " << mask << std::endl;
// 
// 	// ***
// 
// 	const int ioMode = 0;
// 
// 	const char spacer = '\t';
// 	const char endline = '\n';
// 	
// 	std::ofstream file(filename,std::ofstream::out|std::ofstream::trunc);
// 	file.setf(std::ios_base::scientific|std::ios_base::showpoint);
// 
// 	// *** header
// 
// 	file << "# time" << spacer;
// 	file << "position x" << spacer;
// 	file << "position y" << spacer;
// 	file << "position z" << spacer;
// 	file << "velocity x" << spacer;
// 	file << "velocity y" << spacer;
// 	file << "velocity z" << spacer;
// 	file << endline;
// 	file << endline;
// 
// 	// *** contents
// 
// 	int statistics[5];
// 	for (int i = 0; i < 5; i++)
// 		statistics[i] = 0;
// 
// 	int progressCnt(0);
// 	for (SurfaceNodes::const_iterator i = surface.begin(); i != surface.end(); i++)
// 	{
// 		std::cout << "\t-> computing trajectories (this may take a very long time): ";
// 		std::cout << ++progressCnt << " / " << surface.size() << std::string(10, ' ') << '\r';
// 		std::cout.flush();
// 
// 		// skip trajectories starting at the bottom
// 		if (system.geomTable.nodeCoords(i->index()).z() == 0.0f) continue;
// 
// 		if (std::rand()/static_cast<float>(RAND_MAX) > mask) continue;
// 
// 		// *** compute trajectory
// 
// 		const Trajectory_3d::Vector3d velocity = Trajectory_3d::Vector3d(0.0f);
// 
// 		const float charge = Trajectory_3d::eCharge;
// 		const float mass = 63.4f * Trajectory_3d::amu2kg;
// 
// 		Trajectory_3d trajectory(&system.geomTable,&system.gridTable,Trajectory_3d::O5_INTEGRATOR,Trajectory_3d::ERROR_RESTRICTED);
// 		trajectory.setDeltaLimit(1e-16);
// 		trajectory.init(system.geomTable.nodeCoords(i->index()),velocity,charge,mass);
// 		trajectory.integrate(1e-14,accuracy);
// 
// 		statistics[trajectory.status()]++;
// 
// 		// *** write to file
// 
// 		file << "# index = " << progressCnt << endline;
// 		file << "# integrator status = " << trajectory.status() << " => " << Trajectory_3d::status_str(trajectory.status()) << endline;
// 		file << "# charge [C] / mass [kg] = " << trajectory.charge() << " / " << trajectory.mass() << endline;
// 	
// 		file << "# node index = " << i->index() << endline;
// 
// 		file << "# evaporation probability = " << i->probability() << endline;
// 	
// 		if (ioMode == File_Io::BINARY)
// 		{
// 			stream << "BINARY " << trajectory.data().size() << '\n';
// 			for (std::vector<Trajectory_3d::phaseVector>::const_iterator i = trajectory.data().begin(); i != trajectory.data().end(); i++)
// 				stream.write(reinterpret_cast<const char*>(&*i),sizeof(Trajectory_3d::phaseVector));
// 		}
// 		else
// 		{
// 			stream << "ASCII " << trajectory.data().size() << '\n';
// 			for (std::vector<Trajectory_3d::phaseVector>::const_iterator i = trajectory.data().begin(); i != trajectory.data().end(); i++)
// 			{
// 				stream << i->time() << '\t';
// 				for (int j = 0; j < 3; j++)
// 					stream << i->position(j) << '\t';
// 				for (int j = 0; j < 3; j++)
// 					stream << i->velocity(j) << '\t';
// 				stream << i->tetIndex() << '\n';
// 			}
// 		} 
// 	
// 		file << "# final error estimate (position) = (";
// 		file << trajectory.error_estimate().position(Trajectory_3d::Vector3d::X) << "/";
// 		file << trajectory.error_estimate().position(Trajectory_3d::Vector3d::Y) << "/";
// 		file << trajectory.error_estimate().position(Trajectory_3d::Vector3d::Z) << ")";
// 		file << endline;
// 	
// 		file << "# final error estimate (velocity) = (";
// 		file << trajectory.error_estimate().velocity(Trajectory_3d::Vector3d::X) << "/";
// 		file << trajectory.error_estimate().velocity(Trajectory_3d::Vector3d::Y) << "/";
// 		file << trajectory.error_estimate().velocity(Trajectory_3d::Vector3d::Z) << ")";
// 		file << endline;
// 
// 		file << endline;
// 	}
// 
// 	file.close();
// 
// 	// ***
// 
// 	std::cout << "\t-> computing trajectories: finished, " << progressCnt << " records(s) processed";
// 	std::cout << std::string(75,' ') << std::endl;
// 
// 	std::cout << "\t-> statistics: " << std::endl;
// 	for (int i = 0; i < 5; i++)
// 		std::cout << "\t\t-> trajectories with status #" << i << ": " << statistics[i] << std::endl;
// 
// 	std::cout << std::endl;
}

// ***

void Debug::testCoupling(std::ostream& stream, const Grid_3d::Table& gridTable)
{
	int errorCnt(0);
	std::ostringstream errorStream;

	stream << "*** Testing coupling factors:" << std::endl;

	for (int i = 0; i < gridTable.numNodes(); i++)
	{
		if (i % (1+gridTable.numNodes()/100) == 0)
		{
			stream << "\t-> processing: " << i << " / " << gridTable.numNodes() << '\r';
			stream.flush();
		}

		for (int j = 0; j < gridTable.node(i).numNeighbours(); j++)
		{
			const int neighIndex = gridTable.node(i).neighbour(j);

			for (int k = 0; k < gridTable.node(neighIndex).numNeighbours(); k++)
			{
				if (gridTable.node(neighIndex).neighbour(k) == i)
				{
					if (gridTable.node(i).coupling(j) != gridTable.node(neighIndex).coupling(k))
					{
						errorStream << "\t\tnode: " << i << ", neighbour: " << neighIndex;
						errorStream << " => coupling-factors: " << gridTable.node(i).coupling(j) << "/";
						errorStream << gridTable.node(neighIndex).coupling(k) << std::endl;

						errorCnt++;
					}

					break;
				}
			}
		}
	}

	stream << "\t-> processing: finished" << std::string(25,' ') << std::endl;

	if (errorCnt > 0)
	{
		stream << "\t-> suspicious nodes detected: " << std::endl;
		stream << errorStream.str();
	}

	stream << "\t-> " << errorCnt << " problems(s) found" << std::endl;
	stream << std::endl;
}

void Debug::testNeighbourDistance(const char* filename, const Grid_3d::Table& gridTable)
{
	std::ofstream file(filename,std::ofstream::out|std::ofstream::trunc);
	file.setf(std::ios_base::scientific|std::ios_base::showpoint);

	float meanDiff(0.0);

	for (int i = 0; i < gridTable.numNodes(); i++)
	{
		file << "*** node: #" << i << std::endl;

		float diff(0.0);

		file << "\tneighbours (size = " << gridTable.node(i).numNeighbours() << ": ";

		for (int j = 0; j < gridTable.node(i).numNeighbours(); j++)
		{
			if (j < gridTable.node(i).numNeighbours() - 1)
				file << "#" << gridTable.node(i).neighbour(j) << ", ";
			else
				file << "#" << gridTable.node(i).neighbour(j) << ", ";

			diff += (gridTable.node(i).neighbour(j) - i) * (gridTable.node(i).neighbour(j) - i);
		}

		file << std::endl;

		file << "\t diff-Value: " << diff << std::endl;

		meanDiff += diff;
	}

	file << "### MEAN DIFF-VALUE = " << meanDiff << std::endl;
}

void Debug::surfaceSize_statistics(const System_3d& system)
{
	std::map<short,float> idArea;
	std::map<short,unsigned int> idCount;

	for (int i = 0; i < system.geomTable.numNodes(); i++)
	{
		const float epsilon_1 = 1.0f / system.configTable[system.gridTable.node(i).id()].epsilon();

		float surfaceArea(0.0f);
		for (int j = 0; j < system.gridTable.node(i).numNeighbours(); j++)
		{
			const int neighIndex = system.gridTable.node(i).neighbour(j);

			const float epsilon_2 = 1.0f / system.configTable[system.gridTable.node(neighIndex).id()].epsilon();

			float tmp = system.gridTable.node(i).coupling(j);
			tmp *= system.geomTable.distance(i,neighIndex);
			tmp /= 2.0f;
			tmp *= epsilon_1 + epsilon_2;

			surfaceArea += tmp;
		}

		idArea[system.gridTable.node(i).id().toValue()] += surfaceArea;
		idCount[system.gridTable.node(i).id().toValue()] += 1;
	}

	std::cout << "Surface size statistics:" << std::endl;
	for (std::map<short,float>::const_iterator iter = idArea.begin(); iter != idArea.end(); iter++)
	{
		float meanValue = iter->second;
		meanValue /= idCount[iter->first];

		float devValue(0.0f);
		int count(0);

		for (int i = 0; i < system.geomTable.numNodes(); i++)
		{
			if (system.gridTable.node(i).id().toValue() != iter->first) continue;

			const float epsilon_1 = 1.0f / system.configTable[system.gridTable.node(i).id()].epsilon();

			float surfaceArea(0.0f);
			for (int j = 0; j < system.gridTable.node(i).numNeighbours(); j++)
			{
				const int neighIndex = system.gridTable.node(i).neighbour(j);

				const float epsilon_2 = 1.0f / system.configTable[system.gridTable.node(neighIndex).id()].epsilon();

				float tmp = system.gridTable.node(i).coupling(j);
				tmp *= system.geomTable.distance(i,neighIndex);
				tmp /= 2.0f;
				tmp *= epsilon_1 + epsilon_2;

				surfaceArea += tmp;
			}

			devValue += std::pow(surfaceArea-meanValue,2.0f);
			count++;
		}

		if (count > 1) devValue /= count;
		devValue = std::sqrt(devValue);

		std::cout << iter->first << "\tarea: " << iter->second;
	       	std::cout << "\tcounts: " << idCount[iter->first];
		std::cout << "\tmean area: " << meanValue;
		std::cout << "\tdeviation: " << devValue;
		std::cout << std::endl;
	}
}
