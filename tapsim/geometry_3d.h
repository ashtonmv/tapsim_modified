#ifndef TAPSIM_GEOMETRY_3D_H
#define TAPSIM_GEOMETRY_3D_H

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
#include <vector>
#include <string>
#include <set>

#include "../vector/vector.h"

namespace Geometry_3d
{
	typedef MathVector3d<float> Point;

	class Triangle
	{
		public:
			Triangle();
			Triangle(const Point&, const Point&, const Point&);
	
			float area() const;
			float orient(const Point&) const;

			Triangle& flip();

			Point normal() const;
	
			const Point& operator[](const int index) const { return _vertices[index]; }
			Point& operator[](const int index) { return _vertices[index]; }
	
			const Point& a() const { return _vertices[0]; }
			Point& a() { return _vertices[0]; }
		
			const Point& b() const { return _vertices[1]; }
			Point& b() { return _vertices[1]; }
		
			const Point& c() const { return _vertices[2]; }
			Point& c() { return _vertices[2]; }

			Point barycenter() const;

		private:
			Point _vertices[3];
	};

	class Tetrahedron
	{
		public:
			Tetrahedron();
			Tetrahedron(const Point&, const Point&, const Point&, const Point&);

			float volume() const;
			float orient() const;

			Tetrahedron& flip();
	
			const Point& operator[](const int index) const { return _vertices[index]; }
			Point& operator[](const int index) { return _vertices[index]; }
	
			const Point& a() const { return _vertices[0]; }
			Point& a() { return _vertices[0]; }
	
			const Point& b() const { return _vertices[1]; }
			Point& b() { return _vertices[1]; }
	
			const Point& c() const { return _vertices[2]; }
			Point& c() { return _vertices[2]; }

			const Point& d() const { return _vertices[3]; }
			Point& d() { return _vertices[3]; }

			Point barycenter() const;

		private:
			Point _vertices[4];
	};

	// ***

	struct VoronoiFace
	{
		/* !!! no constructor !!! */

		bool isFinite() const;

		Point center() const;
		Point finiteCenter(const float) const;
		//Point constrainedCenter(const Point&, const Point&) const;

		float area() const;
		float finiteArea(const float) const;
		//float constrainedArea(const Point&, const Point&) const;

		Point normal() const;

		int n1;
		int n2;

		Point p1;
		Point p2;

		Point d1;
		Point d2;

		std::list<int> tets;
		std::list<Point> vertices;
		
		//void constrainedVertices(const Point&, const Point&, std::list<Point>*) const;
	};

	struct VoronoiCell
	{
		/* !!! no constructor !!! */

		bool isFinite() const;

		float volume() const;
		float finiteVolume(const float) const;
		//float constrainedVolume(const Point& min, const Point& max) const;

		int n;
		Point p;

		std::list<VoronoiFace> faces;
	};

	// ***

	class Table
	{
		public:
			struct NodeData
			{
				NodeData() : coords(0.0f), boundary(false), associatedTetrahedron(-1) {}
		
				Point coords;
				bool boundary;

				int index;
				int associatedTetrahedron;
			};
	
			struct TetrahedronData
			{
				TetrahedronData() : vertices(-1), neighbours(-1) {}
				TetrahedronData(const int a, const int b, const int c, const int d) : vertices(a,b,c,d), neighbours(-1) {}

				bool hasNeighbour(const int index) const { return neighbours[index] != -1; }
	
				Quadruple<int> vertices;
				Quadruple<int> neighbours;
			};

			// ***
	
			Table();
			~Table() {}

			void operator<<(std::ifstream&);
			void operator>>(std::ofstream&) const;

			void allocate(const int =-1, const int =-1);
			void reset() { allocate(); }

			void set_metaInformation();
			void reset_metaInformation();

			// ***

			void allocateNodes(const int =-1);
			void resetNodes() { allocateNodes(); }
			
			int numNodes() const { return _nodes.size(); }

			const NodeData& node(const int index) const { return _nodes[index]; }
			NodeData& node(const int index) { return _nodes[index]; }

			const Point& nodeCoords(const int index) const { return _nodes[index].coords; }

			bool isBoundaryNode(const int index) const { return _nodes[index].boundary; }

			const Point& min() const { return _minCoords; }
			const Point& max() const { return _maxCoords; }

			const Point& constrainedMin() const { return _constrainedMinCoords; }
			const Point& constrainedMax() const { return _constrainedMaxCoords; }

			// ***

			void allocateTets(const int =-1);
			void resetTets() { allocateTets(); }
			
			void tetrahedralize();

			int numTetrahedra() const { return _tetrahedra.size(); }
			const TetrahedronData& tetrahedron(const int index) const { return _tetrahedra[index]; }
			TetrahedronData& tetrahedron(const int index) { return _tetrahedra[index]; }

			const Point& tetCoords(const int, const int) const;

			bool isBoundaryTetrahedron(const int index) const;

			Point tetCenter(const int) const;
			float tetOrient(const int) const;
			float tetVolume(const int) const;

			void  tetWeights(const int, const Point&, Quadruple<float>*) const;

			// ***

			int findTetrahedron(const Point&) const;
			int findTetrahedron(const Point&, int) const;

			// ***

			void adjacentTetrahedra(const int, std::set<int>*) const;

			void adjacentNodes(const int, std::set<int>*) const;
			void adjacentNodes(const int, std::set<int>*, const int) const;

			void tetsWithCommonNode(const int, std::set<int>*) const;
			void tetsWithCommonNode(const int, const int, std::set<int>*) const;

			// ***

			float distance(const int, const int) const;

			void get_voronoiFace(const int, const int, VoronoiFace*) const;
			float voronoiArea(const int, const int, const float =0.0f) const;
			//float constrainedVoronoiArea(const int, const int) const;

			void get_voronoiCell(const int, VoronoiCell*) const;
			float voronoiVolume(const int, const float =0.0f) const;
			//float constrainedVoronoiVolume(const int) const;

			// ***

		private:
			static const std::string _binaryVersion;

			// ***

			void setExtents();
			void setBoundary();
			void setConstrainedExtents();
			void setAssociatedTetrahedra();

			void setOrientation();
			void setNeighbours();

			// ***

			void get_boundaryFacets(std::vector<int>*) const;

			// ***

			float sign(const Quadruple<int>&) const;
			
			int tetGuess(const Point&) const;

			float tetOrient(const int, const int, const int, const Point&) const;
			float tetOrient(const int, const int, const int, const int) const;

			// ***

			Point _minCoords;
			Point _maxCoords;
			
			Point _constrainedMinCoords;
			Point _constrainedMaxCoords;

			DynamicVector<NodeData> _nodes;
			DynamicVector<TetrahedronData> _tetrahedra;
	};

	//Quadruple<float> toBarycentric(const Tetrahedron&, Point);
	
	Point find_best_plane(const std::vector<Point>&);
	Point find_best_plane(const Point&, const std::vector<Point>&);
	
	bool intersect_plane_with_line(const Point&, const Point&, const Point&, const Point&, Point* =0, Point* =0);
}

#endif
