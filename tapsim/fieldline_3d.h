#ifndef TAPSIM_FIELDLINE_3D_H
#define TAPSIM_FIELDLINE_3D_H

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

#include <vector>
#include <limits>

#include "mathVector_3d.h"
#include "geometry_3d.h"
#include "grid_3d.h"

class Fieldline_3d
{
	public:
		typedef MathVector3d<float> Vector3d;

		// ***

		class phaseVector
		{
			public:
				phaseVector()
					: _position(),
					  _field(),
					  _tetIndex(-1)
				{}

				phaseVector(const Vector3d& p)
					: _position(p),
					  _field(),
					  _tetIndex(-1)
				{}

				const Vector3d& position() const { return _position; }
				Vector3d& position() { return _position; }
				
				const Vector3d& field() const { return _field; }
				Vector3d& field() { return _field; }
				
				int tetIndex() const { return _tetIndex; }
				int& tetIndex() { return _tetIndex; }

			private:
				Vector3d _position;	// [_position] = m
				Vector3d _field;	// [_field] = V/m
				int _tetIndex;		// [_tetIndex] = 1
		};

		// ***

		static const char* status_str(const int);
		
		enum { INVALID =-1, NO_COMPUTATION, NO_FIELD, STEPPER_LIMIT_PASSED, ITERATION_LIMIT_PASSED, GEOMETRY_LIMIT_PASSED, SYSTEM_LIMIT_PASSED };

		// ***

		Fieldline_3d(const Geometry_3d::Table* =0, const Grid_3d::Table* =0);

		void setGeometry(const Geometry_3d::Table* geomTable) {_fastField.setGeometry(geomTable); }
		const Geometry_3d::Table* geometry() const { return _fastField.geometry(); }

		void setGrid(const Grid_3d::Table* gridTable) { _fastField.setGrid(gridTable); }
		const Grid_3d::Table* grid() const { return _fastField.grid(); }

		// ***

		void init(const Vector3d&);

		int compute(const float =-1.0f);
		void reset();

		// ***

		int status() const { return _status; }
		
		const std::vector<phaseVector>& data() const { return _data; }

	private:
		bool ready() const;

		phaseVector stepper(float* delta, const std::set<int>& allowedTetrahedra);

		// ***

		mutable Grid_3d::FastField _fastField;

		unsigned long _stepperLimit;
		unsigned long _iterationLimit;

		// ***

		int _status;

		std::vector<phaseVector> _data;
};

#endif
