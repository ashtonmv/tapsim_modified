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

#include "geometry_3d.h"

#include <algorithm>
#include <stdexcept>
#include <fstream>
#include <limits>
#include <queue>
#include <set>

#include <cmath>
#include <cstdlib>

#include "../tetgen/dist_tetgen.h"
#include "../predicates/predicates.h"
#include "../vector/vector.h"

#include "utils.h"

#include <iostream>

#include "debug.h"

// ***** ===== GEOMETRY_3D::TRIANGLE ===== *****

Geometry_3d::Triangle::Triangle()
{
	for (int index = 0; index < 3; index++)
		_vertices[index] = Point(0.0f);
}

Geometry_3d::Triangle::Triangle(const Point& a, const Point& b, const Point& c)
{
	_vertices[0] = a;
	_vertices[1] = b;
	_vertices[2] = c;
}

Geometry_3d::Point Geometry_3d::Triangle::normal() const
{
	Point p[2];
	for (int index = 0; index < 2; index++)
	{
		p[index] = _vertices[index+1];
		p[index] -= _vertices[0];
	}

	return vecProduct(p[0],p[1]).normalize();
}

float Geometry_3d::Triangle::area() const
{
	float lumat[2*2];
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
			lumat[i*2+j] = _vertices[i][j] - _vertices[2][j];
	}

	float d;
	if (!lu_decmp<float,2>(lumat,0,&d)) return 0.0f;

	float value = lu_determinant<float,2>(lumat,d);
	value /= 2.0f;

	if (value != value) throw std::runtime_error("Geometry_3d::Triangle::area()");

	return std::fabs(value);
	
// 	Point p[2];
// 	for (int index = 0; index < 2; index++)
// 	{
// 		p[index] = _vertices[index+1];
// 		p[index] -= _vertices[0];
// 	}
// 
// 	float value = vecProduct(p[0],p[1]).length();
// 	value /= 2.0f;
// 
// 	if (value != value) throw std::runtime_error("Geometry_3d::Triangle::area()");
// 	
// 	return std::fabs(value);
}

float Geometry_3d::Triangle::orient(const Point &d) const
{
	const float* a = _vertices[0].raw();
	const float* b = _vertices[1].raw();
	const float* c = _vertices[2].raw();

	return Predicates::orient3d(b,a,c,d.raw());
}

Geometry_3d::Triangle& Geometry_3d::Triangle::flip()
{
	Point tmp(_vertices[1]);
	_vertices[1] = _vertices[0];
	_vertices[0] = tmp;

	return *this;
}

Geometry_3d::Point Geometry_3d::Triangle::barycenter() const
{
	Point tmp(0.0f);
	for (int i = 0; i < 3; i++)
		tmp += _vertices[i];

	return tmp / 3.0f;
}

// ***** ===== GEOMETRY_3D::TETRAHEDRON ===== *****

Geometry_3d::Tetrahedron::Tetrahedron()
{
	for (int index = 0; index < 4; index++)
		_vertices[index] = Point(0.0f);
}

Geometry_3d::Tetrahedron::Tetrahedron(const Point& a, const Point& b, const Point& c, const Point& d)
{
	_vertices[0] = a;
	_vertices[1] = b;
	_vertices[2] = c;
	_vertices[3] = d;
}

float Geometry_3d::Tetrahedron::volume() const
{
	float lumat[3*3];
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
			lumat[i*3+j] = _vertices[i][j] - _vertices[3][j];
	}

	float d;
	if (!lu_decmp<float,3>(lumat,0,&d)) return 0.0f;

	float value = lu_determinant<float,3>(lumat,d);
	value /= 6.0f;

	if (value != value) throw std::runtime_error("Geometry_3d::Tetrahedron::volume()");

	return std::fabs(value);

// 	Point p[3];
// 	for (int index = 0; index < 3; index++)
// 	{
// 		p[index] = _vertices[index+1];
// 		p[index] -=  _vertices[0];
// 	}
// 
// 	float value = tripProduct(p[0],p[1],p[2]);
// 	value /= 6.0f;
// 
// 	if (value != value) throw std::runtime_error("Geometry3d::Tetrahedron::volume()");
// 
// 	return std::fabs(value);
}

float Geometry_3d::Tetrahedron::orient() const
{
	const float* p1 = _vertices[0].raw();
	const float* p2 = _vertices[1].raw();
	const float* p3 = _vertices[2].raw();
	const float* p4 = _vertices[3].raw();

	return Predicates::orient3d(p2,p1,p3,p4);
}

Geometry_3d::Tetrahedron& Geometry_3d::Tetrahedron::flip()
{
	Point tmp(_vertices[3]);
	_vertices[3] = _vertices[2];
	_vertices[2] = tmp;

	return *this;
}

Geometry_3d::Point Geometry_3d::Tetrahedron::barycenter() const
{
	Point tmp(0.0f);
	for (int i = 0; i < 4; i++)
		tmp += _vertices[i];

	return tmp / 4.0f;
}

// ***** ===== GEOMETRY_3D::VORONOIFACE ===== *****

bool Geometry_3d::VoronoiFace::isFinite() const
{
	if (tets.empty()) throw std::runtime_error("Geometry_3d::VoronoiFace::isFinite()");

	if (tets.size() > 1)
	{
		const int p1 = tets.front();
		const int p2 = tets.back();

		return (p1 == p2  ? true : false);
	}
	else
		return false;
}

Geometry_3d::Point Geometry_3d::VoronoiFace::center() const
{
	if (vertices.size() != tets.size()) throw std::runtime_error("Geometry_3d::VoronoiFace::center()");

	// ***
	
	Point center(0.0f);
	for (std::list<Point>::const_iterator i = vertices.begin(); i != vertices.end(); i++)
		center += *i;
	
	center /= vertices.size();
	
	return center;
}

Geometry_3d::Point Geometry_3d::VoronoiFace::finiteCenter(const float infinityFactor) const
{
	if (vertices.size() != tets.size()) throw std::runtime_error("Geometry_3d::VoronoiFace::finiteCenter()");

	// ***

	Point value(0.0f);

	for (std::list<Point>::const_iterator i = vertices.begin(); i != vertices.end(); i++)
		value += *i;

	int cnt = vertices.size();

	if (!isFinite())
	{
		value += vertices.front();
		value += infinityFactor * d1;

		value += vertices.back();
		value += infinityFactor * d2;

		cnt += 2;
	}

	value /= cnt;

	return value;
}

// Geometry_3d::Point Geometry_3d::VoronoiFace::constrainedCenter(const Point& min, const Point& max) const
// {
// 	if (vertices.size() != tets.size()) throw std::runtime_error("Geometry_3d::VoronoiFace::constrainedCenter()");
// 
// 	// ***
// 	
// 	std::list<Point> polygon;
// 	constrainedVertices(min,max,&polygon);
// 	
// 	Point center(0.0f);
// 	for (std::list<Point>::const_iterator i = polygon.begin(); i != polygon.end(); i++)
// 		center += *i;
// 	
// 	center /= polygon.size();
// 	
// 	return center;
// }

float Geometry_3d::VoronoiFace::area() const
{
	if (vertices.size() != tets.size()) throw std::runtime_error("Geometry_3d::VoronoiFace::area()");

	// ***

	if (!isFinite())
		return std::numeric_limits<float>::infinity();
	else
		return finiteArea(0.0f);
}

float Geometry_3d::VoronoiFace::finiteArea(const float infinityFactor) const
{
	if (vertices.size() != tets.size()) throw std::runtime_error("Geometry_3d::VoronoiFace::finiteArea()");

	// ***

	const Point areaCenter = finiteCenter(infinityFactor);

	float value(0.0f);

	{
		std::list<Point>::const_iterator currentVertex = vertices.begin();
		std::list<Point>::const_iterator nextVertex = vertices.begin();

		advance(nextVertex,1);
		
		while (vertices.end() != nextVertex)
		{
			Point e1 = *currentVertex;
			e1 -= areaCenter;
	
			Point e2 = *nextVertex;
			e2 -= areaCenter;
	
			value += vecProduct(e1,e2).length() / 2.0f;
	
			currentVertex++;
			nextVertex++;
		}
	}

	if (!isFinite())
	{
		Point e1,e2;

		e1 = vertices.front();
		e1 += infinityFactor * d1;
		e1 -= areaCenter;
	
		e2 = vertices.front();
		e2 -= areaCenter;
	
		value += vecProduct(e1,e2).length() / 2.0f;

		// ***

		e1 = vertices.back();
		e1 -= areaCenter;

		e2 = vertices.back();
		e2 += infinityFactor * d2;
		e2 -= areaCenter;

		value += vecProduct(e1,e2).length() / 2.0f;
	}

	if (value != value) throw std::runtime_error("Geometry_3d::VoronoiFace::finiteArea()");
	
	return value;
}

// float Geometry_3d::VoronoiFace::constrainedArea(const Point& min, const Point& max) const
// {
// 	if (vertices.size() != tets.size()) throw std::runtime_error("Geometry_3d::VoronoiFace::constrainedArea()");
// 
// 	// ***
// 	
// 	std::list<Point> polygon;
// 	constrainedVertices(min,max,&polygon);
// 	
// 	Point center(0.0f);
// 	for (std::list<Point>::const_iterator i = polygon.begin(); i != polygon.end(); i++)
// 		center += *i;
// 	
// 	center /= polygon.size();
// 	
// 	// ***
// 
// 	float value(0.0f);
// 
// 	{
// 		std::list<Point>::const_iterator currentVertex = vertices.begin();
// 		std::list<Point>::const_iterator nextVertex = vertices.begin();
// 
// 		advance(nextVertex,1);
// 		
// 		while (vertices.end() != nextVertex)
// 		{
// 			Point e1 = *currentVertex;
// 			e1 -= center;
// 	
// 			Point e2 = *nextVertex;
// 			e2 -= center;
// 	
// 			value += vecProduct(e1,e2).length() / 2.0f;
// 	
// 			currentVertex++;
// 			nextVertex++;
// 		}
// 	}
// 	
// 	if (value != value) throw std::runtime_error("Geometry_3d::VoronoiFace::constrainedArea()");
// 
// 	return value;
// }

Geometry_3d::Point Geometry_3d::VoronoiFace::normal() const
{
	if (vertices.size() != tets.size()) throw std::runtime_error("Geometry_3d::VoronoiFace::normal()");

	return (p2-p1).norm();
}


// void Geometry_3d::VoronoiFace::constrainedVertices(const Point& min, const Point& max, std::list<Point>* points) const
// {
// 	if (vertices.size() != tets.size()) throw std::runtime_error("Geometry_3d::VoronoiFace::constrainedVertices()");
// 
// 	// ***
// 
// 	points->clear();
// 	
// 	const Point directions = (max - min).norm();
// 	const Point* planeOffsets[2] = { &min, & max };
// 
// 	if (!isFinite())
// 	{
// 		float x(0.0f);
// 		
// 		for (int i = 0; i < 2; i++)
// 		{
// 			for (int j = 0; j < 3; j++)
// 			{
// 				// *** add a test in order to check if vertices are within the box constraints...
// 
// 				Point n1(0.0);
// 				n1[j] = directions[j];
// 
// 				Point n2(0.0);
// 				n2[(j+1)%3] = directions[(j+1)%3];
// 
// 				Point offset(0.0f);
// 				offset[(j+2)%3] = (*planeOffsets[i])[(j+2)%3];
// 				offset = vertices.front() - offset;
// 
// 				Point solVector;
// 
// 				const bool result = intersect_plane_with_line(n1,n2,d1,offset,&solVector,0);
// 
// 				if (result && solVector.z() < x) x = solVector.z();
// 			}
// 		}
// 
// 		if (x > 0.0f) points->push_back(vertices.front()+x*d1);
// 	}
// 	
// 	for (std::list<Point>::const_iterator i = vertices.begin(); i != vertices.end(); i++)
// 		points->push_back(*i);
// 
// 	if (!isFinite())
// 	{
// 		float x(0.0f);
// 		
// 		for (int i = 0; i < 2; i++)
// 		{
// 			for (int j = 0; j < 3; j++)
// 			{
// 				// *** add a test in order to check if vertices are within the box constraints...
// 
// 				Point n1(0.0);
// 				n1[j] = directions[j];
// 
// 				Point n2(0.0);
// 				n2[(j+1)%3] = directions[(j+1)%3];
// 
// 				Point offset(0.0f);
// 				offset[(j+2)%3] = (*planeOffsets[i])[(j+2)%3];
// 				offset = vertices.front() - offset;
// 
// 				Point solVector;
// 
// 				const bool result = intersect_plane_with_line(n1,n2,d2,offset,&solVector,0);
// 				
// 				if (result && solVector.z() < x) x = solVector.z();
// 			}
// 		}
// 		
// 		if (x > 0.0f) points->push_back(vertices.back()+x*d2);
// 	}
// }

// ***** ===== GEOMETRY_3D::VORONOICELL ===== *****

bool Geometry_3d::VoronoiCell::isFinite() const
{
	for (std::list<VoronoiFace>::const_iterator i = faces.begin(); i != faces.end(); i++)
		if (!i->isFinite()) return false;

	return true;
}

float Geometry_3d::VoronoiCell::volume() const
{
	float volume(0.0);
	for (std::list<VoronoiFace>::const_iterator face = faces.begin(); face != faces.end(); face++)
	{
		if (!face->isFinite()) return std::numeric_limits<float>::infinity();

		float tmp = face->area();
		tmp *= (face->p1 - face->p2).length();
		tmp /= 3.0f * 2.0f;

		volume += tmp;
	}

	if (volume != volume) throw std::runtime_error("Geometry3d::VoronoiCell::volume()");
	
	return volume;
}

float Geometry_3d::VoronoiCell::finiteVolume(const float infinityFactor) const
{
	float volume(0.0);
	for (std::list<VoronoiFace>::const_iterator face = faces.begin(); face != faces.end(); face++)
	{
		float tmp = face->finiteArea(infinityFactor);
		tmp *= (face->p1 - face->p2).length();
		tmp /= 3.0f * 2.0f;
		
		volume += tmp;
	}
	
	if (volume != volume) throw std::runtime_error("Geometry_3d::VoronoiCell::finiteVolume()");

	return volume;
}

// float Geometry_3d::VoronoiCell::constrainedVolume(const Point& min, const Point& max) const
// {
// 	float volume(0.0);
// 	for (std::list<VoronoiFace>::const_iterator face = faces.begin(); face != faces.end(); face++)
// 	{
// 		float tmp = face->constrainedArea(min,max);
// 		tmp *= (face->p1 - face->p2).length();
// 		tmp /= 3.0f * 2.0f;
// 
// 		volume += tmp;
// 	}
// 	
// 	if (volume != volume) throw std::runtime_error("Geometry_3d::VoronoiCell::constrainedVolume()");
// 
// 	return volume;
// }


// ***** ===== GEOMETRY_3D::TABLE ===== *****

const std::string Geometry_3d::Table::_binaryVersion("1.0");

Geometry_3d::Table::Table()
	: _minCoords(),
	  _maxCoords(),
	  _constrainedMinCoords(),
	  _constrainedMaxCoords(),
	  _nodes(),
	  _tetrahedra()
{
	Predicates::exactinit();
}

void Geometry_3d::Table::operator<<(std::ifstream& stream)
{
	std::string version;
	{
		char c;
		stream.read(reinterpret_cast<char*>(&c),sizeof(char));


		while (c != '\0')
		{
			version.push_back(c);
			stream.read(reinterpret_cast<char*>(&c),sizeof(char));
		}
	}

	if (version != _binaryVersion) throw std::runtime_error("Geometry_3d::Table<<(): wrong version!");

	// ***

	for (int i = 0; i < 3; i++)
		stream.read(reinterpret_cast<char*>(&_minCoords[i]),sizeof(float));

	for (int i = 0; i < 3; i++)
		stream.read(reinterpret_cast<char*>(&_maxCoords[i]),sizeof(float));

	// ***

	for (int i = 0; i < 3; i++)
		stream.read(reinterpret_cast<char*>(&_constrainedMinCoords[i]),sizeof(float));

	for (int i = 0; i < 3; i++)
		stream.read(reinterpret_cast<char*>(&_constrainedMaxCoords[i]),sizeof(float));

	// ***

	{
		int nodeSize;
		stream.read(reinterpret_cast<char*>(&nodeSize),sizeof(int));
		_nodes.allocate(nodeSize);
	}

	for (unsigned int i = 0; i < _nodes.size(); i++)
	{
		for (int j = 0; j < 3; j++)
			stream.read(reinterpret_cast<char*>(&_nodes[i].coords[j]),sizeof(float));

		stream.read(reinterpret_cast<char*>(&_nodes[i].boundary),sizeof(bool));
		stream.read(reinterpret_cast<char*>(&_nodes[i].associatedTetrahedron),sizeof(int));
	}

	// ***

	{
		int tetSize;
		stream.read(reinterpret_cast<char*>(&tetSize),sizeof(tetSize));
		_tetrahedra.allocate(tetSize);
	}

	for (unsigned int i = 0; i < _tetrahedra.size(); i++)
	{
		for (int j = 0; j < 4; j++)
			stream.read(reinterpret_cast<char*>(&_tetrahedra[i].vertices[j]),sizeof(int));

		for (int j = 0; j < 4; j++)
			stream.read(reinterpret_cast<char*>(&_tetrahedra[i].neighbours[j]),sizeof(int));
	}
}

void Geometry_3d::Table::operator>>(std::ofstream& stream) const
{
	const char* version = _binaryVersion.c_str();
	stream.write(version,_binaryVersion.size()+1);

	// ***

	for (int i = 0; i < 3; i++)
		stream.write(reinterpret_cast<const char*>(&_minCoords[i]),sizeof(float));

	for (int i = 0; i < 3; i++)
		stream.write(reinterpret_cast<const char*>(&_maxCoords[i]),sizeof(float));

	// ***

	for (int i = 0; i < 3; i++)
		stream.write(reinterpret_cast<const char*>(&_constrainedMinCoords[i]),sizeof(float));

	for (int i = 0; i < 3; i++)
		stream.write(reinterpret_cast<const char*>(&_constrainedMaxCoords[i]),sizeof(float));

	// ***

	const unsigned int nodeSize = _nodes.size();
	stream.write(reinterpret_cast<const char*>(&nodeSize),sizeof(unsigned int));

	for (unsigned int i = 0; i < _nodes.size(); i++)
	{
		for (int j = 0; j < 3; j++)
			stream.write(reinterpret_cast<const char*>(&_nodes[i].coords[j]),sizeof(float));

		stream.write(reinterpret_cast<const char*>(&_nodes[i].boundary),sizeof(bool));
		stream.write(reinterpret_cast<const char*>(&_nodes[i].associatedTetrahedron),sizeof(int));
	}

	// ***

	const unsigned int tetSize = _tetrahedra.size();
	stream.write(reinterpret_cast<const char*>(&tetSize),sizeof(unsigned int));

	for (unsigned int i = 0; i < _tetrahedra.size(); i++)
	{
		for (int j = 0; j < 4; j++)
			stream.write(reinterpret_cast<const char*>(&_tetrahedra[i].vertices[j]),sizeof(int));

		for (int j = 0; j < 4; j++)
			stream.write(reinterpret_cast<const char*>(&_tetrahedra[i].neighbours[j]),sizeof(int));
	}
}

void Geometry_3d::Table::allocate(const int nodeSize, const int tetSize)
{
	allocateNodes(nodeSize);
	allocateTets(tetSize);
}

void Geometry_3d::Table::set_metaInformation()
{
	// *** update adjacency and additional information
	// *** => keep this sequence of function calls!

	setExtents();
	setNeighbours();
	setBoundary();
	setAssociatedTetrahedra();
	setConstrainedExtents();
}

void Geometry_3d::Table::reset_metaInformation()
{
	for (unsigned int i = 0; i < _nodes.size(); i++)
	{
		_nodes[i].boundary = false;
		_nodes[i].associatedTetrahedron = -1;
	}

	for (unsigned int i = 0; i < _tetrahedra.size(); i++)
	{
		for (int j = 0; j < 4; j++)
			_tetrahedra[i].neighbours[j] = -1;
	}

	_minCoords = Geometry_3d::Point(0.0f);
	_maxCoords = Geometry_3d::Point(0.0f);

	_constrainedMinCoords = Geometry_3d::Point(0.0f);
	_constrainedMaxCoords = Geometry_3d::Point(0.0f);
}

void Geometry_3d::Table::allocateNodes(const int nodeSize)
{
	if (nodeSize < 0)
		_nodes.release();
	else
	{
		_nodes.allocate(nodeSize);
		for (unsigned int index = 0; index < _nodes.size(); index++)
			_nodes[index] = NodeData();
	}
	
	_minCoords = Geometry_3d::Point(0.0f);
	_maxCoords = Geometry_3d::Point(0.0f);
}

void Geometry_3d::Table::allocateTets(const int tetSize)
{
	if (tetSize < 0)
		_tetrahedra.release();
	else
	{
		_tetrahedra.allocate(tetSize);
		for (unsigned int index = 0; index < _tetrahedra.size(); index++)
			_tetrahedra[index] = TetrahedronData();
	}
	
	_constrainedMinCoords = Geometry_3d::Point(0.0f);
	_constrainedMaxCoords = Geometry_3d::Point(0.0f);	
}


void Geometry_3d::Table::tetrahedralize()
{
	// *** compute tetrahedralization
	
	::TetGen::tetgenio input;

	input.firstnumber = 0;
	input.numberofpoints = _nodes.size();
	input.pointlist = new REAL[3 * input.numberofpoints];

	for (int i = 0; i < input.numberofpoints; i++)
	{
		for (int j = 0; j < 3; j++)
			input.pointlist[i*3+j] = static_cast<REAL>(_nodes[i].coords[j]);
	}

	::TetGen::tetgenio output;

	const char* params = "znNFCQ";
	::TetGen::tetrahedralize(const_cast<char*>(params),&input,&output);

	if (static_cast<int>(_nodes.size()) != input.numberofpoints)
		throw std::runtime_error("Geometry_3d::Table::tetrahedralize()");

	input.deinitialize();
	input.initialize();
	
	// *** copy results

	_tetrahedra.allocate(output.numberoftetrahedra);

	int tetIndex(0);
	for (int i = 0; i < output.numberoftetrahedra; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			_tetrahedra[tetIndex].vertices[j] = output.tetrahedronlist[i*4+j];
			_tetrahedra[tetIndex].neighbours[j] = output.neighborlist[i*4+j];
		}

		tetIndex++;
	}

	// *** update adjacency and additional information
	// *** => keep this sequence of function calls!

	setExtents();
	setOrientation();
	setBoundary();
	setAssociatedTetrahedra();
	setConstrainedExtents();
}

const Geometry_3d::Point& Geometry_3d::Table::tetCoords(const int tetIndex, const int index) const
{
	const int nodeIndex = _tetrahedra[tetIndex].vertices[index];

	return _nodes[nodeIndex].coords;
}

bool Geometry_3d::Table::isBoundaryTetrahedron(const int index) const
{
	for (int i = 0; i < 4; i++)
	{
		const int nodeIndex = _tetrahedra[index].vertices[i];
		if (_nodes[nodeIndex].boundary) return true;
	}

	return false;
}

Geometry_3d::Point Geometry_3d::Table::tetCenter(const int tetIndex) const
{
	double lumat[3*3];
	double rhs[3];
	
	for (int i = 0; i < 3; i++)
	{
		double tmp(0.0);
		for (int j = 0; j < 3; j++)
		{
			double diff = tetCoords(tetIndex,i)[j];
			diff -= tetCoords(tetIndex,3)[j];
			
			tmp += diff * diff;
			lumat[i*3+j] = diff;
		}

		rhs[i] = 0.5 * tmp;
	}

	// ***

	double d;
	int index[3];
	
	if (!lu_decmp<double,3>(lumat,index,&d)) throw std::runtime_error("Geometry_3d::Table::tetCenter()");
	lu_solve<double,3>(lumat,index,rhs);

	// ***
	
	Geometry_3d::Point center;
	for (int i = 0; i < 3; i++)
	{
		center[i] = tetCoords(tetIndex,3)[i];
		center[i] += rhs[i];
		
		if (center[i] != center[i]) throw std::runtime_error("Geometry_3d::Table::tetCenter()");
	}
	
	return center;
}

float Geometry_3d::Table::tetOrient(const int index) const
{
	const int n1 = _tetrahedra[index].vertices[0];
	const int n2 = _tetrahedra[index].vertices[1];
	const int n3 = _tetrahedra[index].vertices[2];
	const int n4 = _tetrahedra[index].vertices[3];

	return tetOrient(n1,n2,n3,n4);
}

float Geometry_3d::Table::tetVolume(const int tetIndex) const
{
	float lumat[3*3];
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
			lumat[i*3+j] = tetCoords(tetIndex,i)[j] - tetCoords(tetIndex,3)[j];
	}
	
	float d;
	if (!lu_decmp<float,3>(lumat,0,&d)) return 0.0f;

	float value = lu_determinant<float,3>(lumat,d);
	value /= 6.0f;

	if (value != value) throw std::runtime_error("Geometry_3d::Table::tetVolume()");

	return std::fabs(value);

// 	Point p[3];
// 	for (int i = 0; i < 3; i++)
// 	{
// 		p[i] = tetCoords(tetIndex,i+1);
// 		p[i] -= tetCoords(tetIndex,0);
// 	}
// 
// 	float value = tripProduct(p[0],p[1],p[2]);
// 	value /= 6.0f;
// 	
// 	if (value != value) throw std::runtime_error("Geometry_3d::Table::tetVolume()");
// 
// 	return value;
}

void Geometry_3d::Table::tetWeights(const int tetIndex, const Point& p, Quadruple<float>* weights) const
{
	const float volume = tetVolume(tetIndex);

	if (volume == 0.0f)
	{
		std::cerr << "*** Geometry_3d::Table::tetWeights(): zero volume" << std::endl;
		std::cerr << "\t=> tetrahedron:" << tetIndex << std::endl;
		for (int i = 0; i < 4; i++)
		{
			std::cerr << "\t\t-> " << i << " = " << tetCoords(tetIndex,i) << " ";
			std::cerr << "(vertex = " << _tetrahedra[tetIndex].vertices[i] << ")" << std::endl;
		}
		std::cerr << "\t=> point: " << p << std::endl;

		//throw std::runtime_error("Geometry_3d::Table::tetWeights(): zero volume");

		// *** preliminary fix !!!

		for (int i = 0; i < 4; i++)
			(*weights)[i] = 1.0f/4.0f;
		
		return;
	}

	for (int i = 0; i < 4; i++)
	{
		float lumat[3*3];
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
				lumat[j*3+k] = tetCoords(tetIndex,(i+j)%4)[k] - p[k];
		}

		float d;
		if (!lu_decmp<float,3>(lumat,0,&d))
		{
			(*weights)[(i+3)%4] = 0.0f;
			continue;
		};

		(*weights)[(i+3)%4] = lu_determinant<float,3>(lumat,d);
		(*weights)[(i+3)%4] /= 6.0f;
		(*weights)[(i+3)%4] = std::fabs((*weights)[(i+3)%4]);

		(*weights)[(i+3)%4] /= volume;

		if ((*weights)[(i+3)%4] != (*weights)[(i+3)%4])
		{
			std::cerr << "*** Geometry_3d::Table::tetWeights(): nan" << std::endl;
			std::cerr << "\t=> tetrahedron:" << tetIndex << std::endl;
			for (int j = 0; j < 4; j++)
			{
				std::cerr << "\t\t-> " << j << " = " << tetCoords(tetIndex,j) << " ";
				std::cerr << "(vertex = " << _tetrahedra[tetIndex].vertices[j] << ")" << std::endl;
			}
			std::cerr << "\t=> point: " << p << std::endl;
			std::cerr << "\t=> location: i = " << i << std::endl;

			throw std::runtime_error("Geometry_3d::Table::tetWeights(): nan");
		}
	}

	return;
}

int Geometry_3d::Table::findTetrahedron(const Point& pos) const
{
	const int guessIndex = tetGuess(pos);
	return findTetrahedron(pos,guessIndex);
}

int Geometry_3d::Table::findTetrahedron(const Point& pos, int tetIndex) const
{
	/*for (int loopCnt = 0; loopCnt < _tetrahedra.size(); loopCnt++)
	{
		if (tetIndex < 0) return -1;
		
		// ***

		const TetrahedronData& obj = _tetrahedra[tetIndex];

		int index;
		for (index = 0; index < 4; index++)
		{
			if (tetOrient(obj.vertices[index],obj.vertices[3-(2+index)%4],obj.vertices[(2+index)%4],pos) < 0.0f)
			{
				tetIndex = obj.neighbours[3-index];
				break;
			}
		}

		if (index == 4) return tetIndex;
	}*/

	// *** fall back to stochastic walk which should prevent from running in cycles

	//std::cerr << "Geometry_3d::Table::findTetrahedron(): fall back to stochastic walk" << std::endl;
	
	for (unsigned int loopCnt = 0; loopCnt < _tetrahedra.size(); loopCnt++)
	{
		if (tetIndex < 0) return -1;

		// ***

		const int minIndex = std::rand() % 4;
		const int maxIndex = minIndex + 4;

		int index;
		for (index = minIndex; index < maxIndex; index++)
		{
			const TetrahedronData& obj = _tetrahedra[tetIndex];
			if (tetOrient(obj.vertices[index%4],obj.vertices[3-(2+index)%4],obj.vertices[(2+index)%4],pos) < 0.0f)
			{
				tetIndex = obj.neighbours[3-index%4];
				break;
			}
		}

		if (index == maxIndex) return tetIndex;
	}
	
	// *** error handling

	for (int i = 0; i < 3; i++)
		if (pos[i] != pos[i]) throw std::runtime_error("Geometry_3d::Table::findTetrahedron(): bad position");

	throw std::runtime_error("Geometry_3d::Table::findTetrahedron()");
}

void Geometry_3d::Table::adjacentTetrahedra(const int primaryTetrahedron, std::set<int>* secondaryTetrahedra) const
{
	if (primaryTetrahedron < 0) throw std::runtime_error("Geometry_3d::Table::adjacentTetrahedra()");

	// ***

	secondaryTetrahedra->insert(primaryTetrahedron);

	std::queue<int> untestedTets;
	for (int i = 0; i < 4; i++)
	{
		if (_tetrahedra[primaryTetrahedron].neighbours[i] < 0)
			continue;
		else
			untestedTets.push(_tetrahedra[primaryTetrahedron].neighbours[i]);
	}

	// ***

	while (!untestedTets.empty())
	{
		if (secondaryTetrahedra->find(untestedTets.front()) == secondaryTetrahedra->end())
		{
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					if (_tetrahedra[primaryTetrahedron].vertices[i] == _tetrahedra[untestedTets.front()].vertices[j])
					{
						for (int k = 0; k < 4; k++)
						{
							if (_tetrahedra[untestedTets.front()].neighbours[k] < 0)
								continue;
							else
								untestedTets.push(_tetrahedra[untestedTets.front()].neighbours[k]);
						}
	
						secondaryTetrahedra->insert(untestedTets.front());
	
						i = 4;
						j = 4;
					}
				}
			}
		}

		untestedTets.pop();
	}

	secondaryTetrahedra->erase(primaryTetrahedron);
}

void Geometry_3d::Table::adjacentNodes(const int primaryNode, std::set<int>* secondaryNodes) const
{
	std::set<int> associatedTets;
	tetsWithCommonNode(primaryNode,&associatedTets);
	
	for (std::set<int>::const_iterator i = associatedTets.begin(); i != associatedTets.end(); i++)
	{
		for (int j = 0; j < 4; j++)
			secondaryNodes->insert(_tetrahedra[*i].vertices[j]);
	}

	secondaryNodes->erase(primaryNode);
}

void Geometry_3d::Table::adjacentNodes(const int primaryNode, std::set<int>* secondaryNodes, const int limit) const
{
	if (primaryNode < 0) throw std::runtime_error("Geometry_3d::Table::find_localNodes()");

	// ***

	std::list<int> untestedNodes;
	untestedNodes.push_back(primaryNode);

	for (int i = 0; i < limit; i++)
	{
		std::list<int> nextCandidates;

		while (!untestedNodes.empty())
		{
			std::set<int> localNeighbours;
			adjacentNodes(untestedNodes.front(),&localNeighbours);

			for (std::set<int>::const_iterator j = localNeighbours.begin(); j != localNeighbours.end(); j++)
			{
				if (secondaryNodes->find(*j) == secondaryNodes->end())
				{
					secondaryNodes->insert(*j);
					nextCandidates.push_back(*j);
				}
			}

			untestedNodes.pop_front();
		}

		untestedNodes.splice(untestedNodes.end(),nextCandidates);
	}
	
	secondaryNodes->erase(primaryNode);
}

void Geometry_3d::Table::tetsWithCommonNode(const int nodeIndex, std::set<int>* secondaryTetrahedra) const
{
	const int primaryTetrahedron = _nodes[nodeIndex].associatedTetrahedron;

	int localIndex(0);
	while (_tetrahedra[primaryTetrahedron].vertices[localIndex] != nodeIndex)
		localIndex++;

	tetsWithCommonNode(localIndex,primaryTetrahedron,secondaryTetrahedra);
}


void Geometry_3d::Table::tetsWithCommonNode(const int localIndex, const int primaryTetrahedron, std::set<int>* secondaryTetrahedra) const
{
	secondaryTetrahedra->insert(primaryTetrahedron);

	std::queue<int> untestedTets;
	for (int i = 0; i < 4; i++)
	{
		if (i == localIndex || _tetrahedra[primaryTetrahedron].neighbours[i] < 0)
			continue;
		else
			untestedTets.push(_tetrahedra[primaryTetrahedron].neighbours[i]);
	}

	// ***

	while (!untestedTets.empty())
	{
		const int localNode = _tetrahedra[primaryTetrahedron].vertices[localIndex];

		if (secondaryTetrahedra->find(untestedTets.front()) == secondaryTetrahedra->end())
		{
			for (int i = 0; i < 4; i++)
			{
				if (localNode == _tetrahedra[untestedTets.front()].vertices[i])
				{
					for (int k = 0; k < 4; k++)
					{
						if (_tetrahedra[untestedTets.front()].neighbours[k] < 0)
							continue;
						else
							untestedTets.push(_tetrahedra[untestedTets.front()].neighbours[k]);
					}

					secondaryTetrahedra->insert(untestedTets.front());
					break;
				}
			}
		}

		untestedTets.pop();
	}
}

float Geometry_3d::Table::distance(const int n1, const int n2) const
{
	Point value = _nodes[n1].coords;
	value -= _nodes[n2].coords;

	return value.length();
}

void Geometry_3d::Table::get_voronoiFace(const int n1, const int n2, Geometry_3d::VoronoiFace* obj) const
{
	if (obj == 0) throw std::runtime_error("Geometry_3d::Table::get_voronoiFacet()");

	obj->n1 = n1;
	obj->n2 = n2;

	obj->p1 = _nodes[obj->n1].coords;
	obj->p2 = _nodes[obj->n2].coords;

	bool finiteFace(true);

	// *** I) find indices of any two tetrahedra with common nodes n1 and n2

	{
		bool extendedSearch(false);

		std::list<int> untestedTets;
		untestedTets.push_back(_nodes[obj->n1].associatedTetrahedron);

		obj->tets.clear();
		while (obj->tets.empty())
		{
			Quadruple<int> localIndex;
	
			localIndex[0] = 0;
			while (_tetrahedra[untestedTets.front()].vertices[localIndex[0]] != obj->n1)
				localIndex[0]++;
	
			localIndex[1] = 0;
			while (_tetrahedra[untestedTets.front()].vertices[localIndex[1]] != obj->n2 && localIndex[1] < 4)
				localIndex[1]++;
	
			// ***
	
			if (localIndex[1] == 4)
			{
				if (untestedTets.size() == 1)
				{
					if (!extendedSearch)
					{
						std::set<int> tmpTets;
						tetsWithCommonNode(localIndex[0],untestedTets.front(),&tmpTets);

						untestedTets.clear();
						for (std::set<int>::const_iterator i = tmpTets.begin(); i != tmpTets.end(); i++)
							untestedTets.push_back(*i);

						extendedSearch = true;
					}
					else
						throw std::runtime_error("Geometry_3d::Table::get_voronoiFace()");
				}
				else
					untestedTets.pop_front();
			}
			else
			{
				int j = 2;
				for (int i = 0; i < 4; i++)
				{
					if (i == localIndex[0] || i == localIndex[1])
						continue;
					else
						localIndex.at(j++) = i;
				}
		
				if (sign(localIndex) > 0.0f) std::swap<int>(localIndex[2],localIndex[3]);

				obj->tets.push_back(untestedTets.front());

				// ***

				if (_tetrahedra[untestedTets.front()].neighbours[localIndex[2]] != -1)
					obj->tets.push_back(_tetrahedra[untestedTets.front()].neighbours[localIndex[2]]);
				else
				{
					obj->d2 = vecProduct(nodeCoords(_tetrahedra[obj->tets.back()].vertices[localIndex[3]]) - obj->p1,obj->p2 - obj->p1);
					if (obj->d2.length() > 0.0f) obj->d2 /= obj->d2.length();
		
					finiteFace = false;
				}
			}
		}
	}

	// *** II) traverse clockwise all adjacent tetrahedra with common nodes n1 and n2

	if (finiteFace)
	{
		while (obj->tets.front() != obj->tets.back())
		{
			Quadruple<int> localIndex;

			localIndex[0] = 0;
			while (_tetrahedra[obj->tets.back()].vertices[localIndex[0]] != obj->n1)
				localIndex[0]++;
	
			localIndex[1] = 0;
			while (_tetrahedra[obj->tets.back()].vertices[localIndex[1]] != obj->n2)
				localIndex[1]++;
	
			int j = 2;
			for (int i = 0; i < 4; i++)
			{
				if (i == localIndex[0] || i == localIndex[1])
					continue;
				else
					localIndex.at(j++) = i;
			}
	
			if (sign(localIndex) > 0.0f) std::swap<int>(localIndex[2],localIndex[3]);
	
			if (_tetrahedra[obj->tets.back()].neighbours[localIndex[2]] != -1)
				obj->tets.push_back(_tetrahedra[obj->tets.back()].neighbours[localIndex[2]]);
			else
			{
				obj->d2 = vecProduct(nodeCoords(_tetrahedra[obj->tets.back()].vertices[localIndex[3]]) - obj->p1, obj->p2 - obj->p1);
				if (obj->d2.length() > 0.0f) obj->d2 /= obj->d2.length();
	
				finiteFace = false;

				break;
			}
		}
	}

	// *** III) additionally reverse direction in the case of an unclosed voronoi facet (neighbour index == -1)

	if (!finiteFace)
	{
		while (true)
		{
			Quadruple<int> localIndex;

			localIndex[0] = 0;
			while (_tetrahedra[obj->tets.front()].vertices[localIndex[0]] != obj->n1)
				localIndex[0]++;
	
			localIndex[1] = 0;
			while (_tetrahedra[obj->tets.front()].vertices[localIndex[1]] != obj->n2)
				localIndex[1]++;
	
			int j = 2;
			for (int i = 0; i < 4; i++)
			{
				if (i == localIndex[0] || i == localIndex[1])
					continue;
				else
					localIndex.at(j++) = i;
			}
	
			if (sign(localIndex) > 0.0f) std::swap<int>(localIndex[2],localIndex[3]);
	
			if (_tetrahedra[obj->tets.front()].neighbours[localIndex[3]] != -1)
				obj->tets.push_front(_tetrahedra[obj->tets.front()].neighbours[localIndex[3]]);
			else
			{
				obj->d1 = vecProduct(obj->p2 - obj->p1, nodeCoords(_tetrahedra[obj->tets.front()].vertices[localIndex[2]]) - obj->p1);
				if (obj->d1.length() > 0.0f) obj->d1 /= obj->d1.length();

				break;
			}
		}
	}

	// *** IV) compute voronoi vertices for each entry in the list of tetrahedra

	obj->vertices.clear();
	for (std::list<int>::const_iterator i = obj->tets.begin(); i != obj->tets.end(); i++)
		obj->vertices.push_back(tetCenter(*i));
}

float Geometry_3d::Table::voronoiArea(const int n1, const int n2, const float infinityFactor) const
{
	VoronoiFace face;
	get_voronoiFace(n1,n2,&face);

	if (infinityFactor > 0.0f)
		return face.finiteArea(infinityFactor);
	else
		return face.area();
}

// float Geometry_3d::Table::constrainedVoronoiArea(const int n1, const int n2) const
// {
// 	VoronoiFace face;
// 	get_voronoiFace(n1,n2,&face);
// 	
// 	return face.constrainedArea(_constrainedMinCoords,_constrainedMaxCoords);
// }


void Geometry_3d::Table::get_voronoiCell(const int node, VoronoiCell* obj) const
{
	obj->n = node;
	obj->p = _nodes[obj->n].coords;

	// ***

	std::set<int> delaunayNeighbours;
	adjacentNodes(node,&delaunayNeighbours);

	for (std::set<int>::const_iterator i = delaunayNeighbours.begin(); i != delaunayNeighbours.end(); i++)
	{
		VoronoiFace face;
		get_voronoiFace(node,*i,&face);

		obj->faces.push_back(face);
	}
}

float Geometry_3d::Table::voronoiVolume(const int primaryNode, const float infinityFactor) const
{
	VoronoiCell cell;
	get_voronoiCell(primaryNode,&cell);

	if (infinityFactor > 0.0f)
		return cell.finiteVolume(infinityFactor);
	else
		return cell.volume();

	// 	std::set<int> secondaryNodes;
	// 	adjacentNodes(primaryNode,&secondaryNodes);
	// 
	// 	// ***
	// 
	// 	float volume(0.0);
	// 	for (std::set<int>::const_iterator secondaryNode = secondaryNodes.begin(); secondaryNode != secondaryNodes.end(); secondaryNode++)
	// 	{
	// 		float tmp = voronoiArea(primaryNode,*secondaryNode,infinityFactor);
	// 		tmp *= distance(primaryNode,*secondaryNode);
	// 		tmp /= 3.0f * 2.0f;
	// 
	// 		volume += tmp;
	// 	}
	// 
	// 	return volume;
}

// float Geometry_3d::Table::constrainedVoronoiVolume(const int primaryNode) const
// {
// 	VoronoiCell cell;
// 	get_voronoiCell(primaryNode,&cell);
// 	
// 	return cell.constrainedVolume(_constrainedMinCoords,_constrainedMaxCoords);
// }


// void Geometry_3d::Table::test_pointsWithin(const std::vector<Point>& points, std::vector<bool>* flags) const
// {
// 	std::vector<int> boundaryFacets;
// 	get_boundaryFacets(&boundaryFacets);
// 
// 	flags->resize(points.size());
// 
// 	for (unsigned int i = 0; i < points.size(); i++)
// 	{
// 		for (unsigned int j = 0; j < boundaryFacets.size(); j += 2)
// 		{
// 			int facet[3];
// 
// 			int index(0);
// 			for (int k = 0; k < 4; k++)
// 			{
// 				if (k == boundaryFacets[j+1])
// 					continue;
// 				else
// 					facet[index++] = _tetrahedra[boundaryFacets[j]].vertices[k];
// 			}
// 
// 			const float orient = tetOrient(facet[0],facet[1],facet[2],points[i]);
// 
// 			if (orient == 0.0f)
// 				(*flags)[i] = true;
// 			else
// 			{
// 				if (orient < 0.0f)
// 					(*flags)[i] = true;
// 				else
// 					(*flags)[i] = false;
// 
// 				if (boundaryFacets[j+1] % 2 != 0) (*flags)[i] = !(*flags)[i];
// 			}
// 
// 			if (!(*flags)[i]) break;
// 		}
// 	}
// }

void Geometry_3d::Table::get_boundaryFacets(std::vector<int>* boundaryFacets) const
{
	boundaryFacets->clear();

	for (unsigned int i = 0; i < _tetrahedra.size(); i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if (_tetrahedra[i].neighbours[j] < 0)
			{
				boundaryFacets->push_back(i);
				boundaryFacets->push_back(j);
			}
		}
	}
}

void Geometry_3d::Table::setBoundary()
{
	for (unsigned int i = 0; i < _nodes.size(); i++)
		_nodes[i].boundary = false;

	// ***

	for (unsigned int i = 0; i < _tetrahedra.size(); i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if (_tetrahedra[i].neighbours[j] < 0)
			{
				for (int k = 0; k < 4; k++)
				{
					if (j == k) continue;
	
					const int nodeIndex = _tetrahedra[i].vertices[k];
					_nodes[nodeIndex].boundary = true;
				}
			}
		}
	}
}

void Geometry_3d::Table::setAssociatedTetrahedra()
{
	for (unsigned int i = 0; i < _tetrahedra.size(); i++)
	{
		for (int j = 0; j < 4; j++)
		{
			const int nodeIndex = _tetrahedra[i].vertices[j];
			_nodes[nodeIndex].associatedTetrahedron = i;
		}
	}
}

void Geometry_3d::Table::setExtents()
{
	_minCoords = _nodes[0].coords;
	_maxCoords = _nodes[0].coords;

	for (unsigned int i = 1; i < _nodes.size(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (_nodes[i].coords[j] < _minCoords[j])
				_minCoords[j] = _nodes[i].coords[j];
			else
			{
				if (_nodes[i].coords[j] > _maxCoords[j])
					_maxCoords[j] = _nodes[i].coords[j];
			}
		}
	}
}

void Geometry_3d::Table::setConstrainedExtents()
{
	_constrainedMinCoords = _minCoords;
	_constrainedMaxCoords = _maxCoords;

	for (unsigned int i = 0; i < _nodes.size(); i++)
	{
		if (!_nodes[i].boundary) continue;

		std::set<int> neighNodes;
		adjacentNodes(i,&neighNodes);
		
		for (std::set<int>::const_iterator j = neighNodes.begin(); j != neighNodes.end(); j++)
		{
			VoronoiFace face;
			get_voronoiFace(i,*j,&face);

			if (face.isFinite()) continue;

			for (int k = 0; k < 3; k++)
			{
				if (_constrainedMinCoords[k] > face.vertices.front()[k])
					_constrainedMinCoords[k] = face.vertices.front()[k];
				else if (_constrainedMaxCoords[k] < face.vertices.front()[k])
					_constrainedMaxCoords[k] = face.vertices.front()[k];
			}

			for (int k = 0; k < 3; k++)
			{
				if (_constrainedMinCoords[k] > face.vertices.back()[k])
					_constrainedMinCoords[k] = face.vertices.back()[k];
				else if (_constrainedMaxCoords[k] < face.vertices.back()[k])
					_constrainedMaxCoords[k] = face.vertices.back()[k];
			}
		}
	}
}

void Geometry_3d::Table::setOrientation()
{
	for (unsigned int tetIndex = 0; tetIndex < _tetrahedra.size(); tetIndex++)
	{
		if (tetOrient(tetIndex) < 0.0f)
		{
			std::swap<int>(_tetrahedra[tetIndex].vertices[0],_tetrahedra[tetIndex].vertices[1]);
			std::swap<int>(_tetrahedra[tetIndex].neighbours[0],_tetrahedra[tetIndex].neighbours[1]);
		}
	}
}

void Geometry_3d::Table::setNeighbours()
{
	typedef std::multimap<int,int> Multimap;

	Multimap indexTable;
	for (unsigned int i = 0; i < _tetrahedra.size(); i++)
	{
		for (int j = 0; j < 4; j++)
			indexTable.insert(std::pair<int,int>(_tetrahedra[i].vertices[j],i));

		for (int j = 0; j < 4; j++)
			_tetrahedra[i].neighbours[j] = -1;
	}

	// ***

	for (unsigned int i = 0; i < _nodes.size(); i++)
	{
		const std::pair<Multimap::const_iterator,Multimap::const_iterator> range = indexTable.equal_range(i);

		if (range.first == range.second) throw std::runtime_error("Geometry_3d::Table::setNeighbours()");

		// ***

		int jCnt(0);

		for (Multimap::const_iterator j = range.first; j != range.second; j++)
		{
			++jCnt;

			for (Multimap::const_iterator k = range.first; k != range.second; k++)
			{
				if (j == k) continue;

				int jNeighIndex(6);
				int kNeighIndex(6);

				int equalCnt(0);

				for  (int vj = 0; vj < 4; vj++)
				{
					for (int vk = 0; vk < 4; vk++)
					{
						if (_tetrahedra[j->second].vertices[vj] == _tetrahedra[k->second].vertices[vk])
						{
							equalCnt++;
							jNeighIndex -= vj;
							kNeighIndex -= vk;
							break;
						}
					}
				}

				if (equalCnt == 3)
				{
					_tetrahedra[j->second].neighbours[jNeighIndex] = k->second;
					_tetrahedra[k->second].neighbours[kNeighIndex] = j->second;
				}
			}
		}

		indexTable.erase(i);
	}
}

float Geometry_3d::Table::sign(const Quadruple<int>& n) const
{
	int value(0);
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if (n[i] == j)
			{
				for (int k = 0; k < j; k++)
					if (n[k] > n[j]) value++;

				break;
			}
		}
	}

	value %= 2;

	if (value) return -1.0f; else return 1.0f;
}

int Geometry_3d::Table::tetGuess(const Point& pos) const
{
	int refIndex = -1; 
	float refDistance = std::numeric_limits<float>::max();

	int sampleCnt = std::pow(_tetrahedra.size(),0.25);

	while (sampleCnt-- > 0)
	{
		int tetIndex = std::rand();
		tetIndex %= _tetrahedra.size();
		
		for (int vertex = 0; vertex < 4; vertex++)
		{
			const Geometry_3d::Point& coords = tetCoords(tetIndex,vertex);
			
			float distance(0.0f);
			for (int i = 0; i < 3; i++)
			{
				float tmp = coords[i] - pos[i];
				tmp *= tmp;
				
				distance += tmp;
			}
			
			if (distance < refDistance)
			{
				refDistance = distance;
				refIndex = tetIndex;
			}
		}
	}
	
	if (refIndex < 0) throw std::runtime_error("Geometry_3d::Table::tetGuess()");
	
	return refIndex;
}

float Geometry_3d::Table::tetOrient(const int n1, const int n2, const int n3, const Point& v) const
{
	const float* p1 = _nodes[n1].coords.raw();
	const float* p2 = _nodes[n2].coords.raw();
	const float* p3 = _nodes[n3].coords.raw();
	const float* p4 = v.raw();

	return Predicates::orient3d(p2,p1,p3,p4); // returned value is >0 if p1,p2,p3,p4 are in clockwise order
}

float Geometry_3d::Table::tetOrient(const int n1, const int n2, const int n3, const int n4) const
{
	const float* p1 = _nodes[n1].coords.raw();
	const float* p2 = _nodes[n2].coords.raw();
	const float* p3 = _nodes[n3].coords.raw();
	const float* p4 = _nodes[n4].coords.raw();

	return Predicates::orient3d(p2,p1,p3,p4); // returned value is >0 if p1,p2,p3,p4 are in clockwise order
}

// ***** ===== GEOMETRY_3D::TOBARYCENTRIC() ===== *****

// void  tetWeights(const int, const Point&, Quadruple<float>* weights) const;
// Quadruple<float> Geometry_3d::toBarycentric(const Tetrahedron& obj, Point p)
// {
// 	const float volume = obj.volume();
// 
// 	if (volume == 0.0f) 
// 	{
// 		// write a message to stderr
// 		cerr << "*** Geometry_3d::toBarycentric(): volume(a,b,c,d) == 0.0f" << std::endl;
// 		cerr << "\t=> tetrahedron:" << std::endl;
// 		cerr << "\t\ta = " << obj.a() << std::endl;
// 		cerr << "\t\tb = " << obj.b() << std::endl;
// 		cerr << "\t\tc = " << obj.c() << std::endl;
// 		cerr << "\t\td = " << obj.d() << std::endl;
// 		cerr << "\t=> point: " << p < <std::endl;
// 
// 		//throw std::runtime_error("Geometry_3d::toBarycentric(): volume(a,b,c,d) == 0.0f");
// 		return Quadruple<float>(0.25f,0.25f,0.25f,0.25f); // !!! preliminary fix !!!
// 	}
// 
// 	Quadruple<float> value;
// 	for (int index = 0; index < 4; index++)
// 	{
// 		value[(index+3)%4] = Tetrahedron(obj[index],obj[(index+1)%4],obj[(index+2)%4],p).volume();
// 		value[(index+3)%4] /= volume;
// 
// 		if (value[(index+3)%4] != value[(index+3)%4]) throw std::runtime_error("Geometry3d::toBarycentric()");
// 	}
// 
// 	return value;
// }

// ***** ===== GEOMETRY_3D::FIND_BEST_PLANE() ===== *****

Geometry_3d::Point Geometry_3d::find_best_plane(const std::vector<Geometry_3d::Point>& coords)
{
	float a[3*3];

	a[0*3+0] = 0.0;
	for (unsigned int i = 0; i < coords.size(); i++)
		a[0*3+0] += coords[i].x() * coords[i].x();

	a[0*3+1] = 0.0;
	for (unsigned int i = 0; i < coords.size(); i++)
		a[0*3+1] += coords[i].x() * coords[i].y();

	a[0*3+2] = 0.0;
	for (unsigned int i = 0; i < coords.size(); i++)
		a[0*3+2] += coords[i].x() * coords[i].z();

	a[1*3+1] = 0.0;
	for (unsigned int i = 0; i < coords.size(); i++)
		a[1*3+1] += coords[i].y() * coords[i].y();

	a[1*3+2] = 0.0;
	for (unsigned int i = 0; i < coords.size(); i++)
		a[1*3+2] += coords[i].y() * coords[i].z();

	a[2*3+2] = 0.0;
	for (unsigned int i = 0; i < coords.size(); i++)
		a[2*3+2] += coords[i].z() * coords[i].z();

	// ***

	a[1*3+0] = a[0*3+1];
	a[2*3+0] = a[0*3+2];
	a[2*3+1] = a[1*3+2];

	// ***

	float b[3],c[3];
	for (int i = 0; i < 3; i++)
	{
		b[i] = 0.0f;
		c[i] = 0.0f;
	}

	householder<float,3>(a,b,c);
	tri_eigen_qli<float,3>(b,c,a);

	Geometry_3d::Point value;

	{
		int index(0);
		for (int i = 1; i < 3; i++)
			if (b[i] < b[index]) index = i;
			    
		for (int i = 0; i < 3; i++)
		{
			if (a[i*3+index] != a[i*3+index]) throw std::runtime_error("Geometry_3d::find_best_plane()");
			
			value[i] = a[i*3+index];
		}
	}

	if (Geometry_3d::Point(0.0f,0.0f,0.0f) != value) value.normalize();

	return value;
}

Geometry_3d::Point Geometry_3d::find_best_plane(const Geometry_3d::Point& reference, const std::vector<Geometry_3d::Point>& source)
{
	std::vector<Geometry_3d::Point> coords;
	for (std::vector<Geometry_3d::Point>::const_iterator i = source.begin(); i != source.end(); i++)
		coords.push_back(reference - *i);

	return find_best_plane(coords);
}

/*bool Geometry_3d::intersect_plane_with_line(const Point& e1, const Point& e2, const Point& n, const Point& offset, Point* solVector, Point* solCoords)
{
	LOCAL_REAL a[4][4];

	for (int i = 0; i < 3; i++)
	{
		a[0][i] = e1[i];
		a[1][i] = e2[i];
		a[2][i] = n[i];
	}

	LOCAL_REAL b[4];

	for (int i = 0; i < 3; i++)
		b[i] = offset[i];

	int index[4];
	LOCAL_REAL d;

	const bool hasResult = lu_decmp<LOCAL_REAL,3>(a,index,&d,0);

	// ***

	if (hasResult)
	{
		lu_solve<LOCAL_REAL,3>(a,index,b,0);

		if (solVector != 0)
		{
			for (int i = 0; i < 3; i++)
				(*solVector)[i] = b[i];
		}

		if (solCoords != 0)
		{
			for (int i = 0; i < 3; i++)
			{
				(*solCoords)[i] = 0.0;
				for (int j = 0; j< 3; j++)
					(*solCoords)[i] += a[i][j] * b[i];
			}
		}
	}
	else
	{
		if (solVector != 0)
		{
			for (int i = 0; i < 3; i++)
				(*solVector)[i] = 0.0;
		}
		
		if (solCoords != 0)
		{
			for (int i = 0; i < 3; i++)
				(*solCoords)[i] = 0.0;
		}
	}
	
	return hasResult;
}*/
