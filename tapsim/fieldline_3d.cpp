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

#include "fieldline_3d.h"

#include <iostream>

#include <limits>
#include <stdexcept>
#include <cmath>
#include <cstring>

const char* Fieldline_3d::status_str(const int value)
{
	switch (value)
	{
		case NO_COMPUTATION:
			return "NO_COMPUTATION";
		case NO_FIELD:
			return "NO_FIELD";
		case STEPPER_LIMIT_PASSED:
			return "STEPPER_LIMIT_PASSED";
		case ITERATION_LIMIT_PASSED:
			return "ITERATION_LIMIT_PASSED";
		case GEOMETRY_LIMIT_PASSED:
			return "GEOMETRY_LIMIT_PASSED";
		case SYSTEM_LIMIT_PASSED:
			return "SYSTEM_LIMIT_PASSED";
		default:
			return "INVALID";
	}
}

Fieldline_3d::Fieldline_3d(const Geometry_3d::Table* geomTable, const Grid_3d::Table* gridTable)
	: _fastField(gridTable,geomTable)
{}

void Fieldline_3d::init(const Vector3d& position)
{
	phaseVector obj(position);
	obj.field() = _fastField.compute(position);
	obj.tetIndex() = _fastField.geometry()->findTetrahedron(obj.position());
	if (obj.tetIndex() < 0) throw std::runtime_error("Fieldline_3d::init()");
	
	_status = NO_COMPUTATION;
	
	_data.clear();
	_data.push_back(obj);
}

// ***

int Fieldline_3d::compute(const float distance)
{
	if (!ready()) throw std::runtime_error("Fieldline_3d::compute()");

	// ***

	_status = INVALID;

	float delta = _fastField.geometry()->tetVolume(_data.back().tetIndex());
	delta = std::pow(delta,1.0f/3.0f);
	delta /= 2.0f;
	
	unsigned long iteration(_iterationLimit);
	
	while (iteration--)
	{
		if (_data.back().field() == Vector3d(0.0f))
		{
			_status = NO_FIELD;
			return _status;
		}
		
		try
		{
			phaseVector next;

			std::set<int> allowedTetrahedra;
			_fastField.geometry()->adjacentTetrahedra(_data.back().tetIndex(),&allowedTetrahedra);
			next = stepper(&delta,allowedTetrahedra);

			_data.push_back(next);
		}
		catch (int status)
		{
			_status = status;
			return status;
		}

		// check global limits

		if (_data.back().tetIndex() < 0)
		{
			_data.pop_back();
			_status = SYSTEM_LIMIT_PASSED;
			return _status;
		}

		// check user limits

		if (distance >= 0.0f)
		{
			Vector3d position = _data.front().position();
			position -= _data.back().position();

			if (position.length() > distance)
			{
				_data.pop_back();
				_status = GEOMETRY_LIMIT_PASSED;
				return _status;
			}
		}
	}

	_status = ITERATION_LIMIT_PASSED;
	return _status;
}

void Fieldline_3d::reset()
{
	const phaseVector tmp(_data.front());

	_status = NO_COMPUTATION;

	_data.clear();
	_data.push_back(tmp);
}

bool Fieldline_3d::ready() const
{
	if (_fastField.geometry() == 0 || _fastField.grid() == 0) return false;
	if (_status != NO_COMPUTATION || _data.back().tetIndex() < 0) return false;

	return true;
}

// ***

Fieldline_3d::phaseVector Fieldline_3d::stepper(float *delta, const std::set<int>& allowedTetrahedra)
{
	unsigned long iteration(_stepperLimit);

	while (iteration--)
	{
		Fieldline_3d::phaseVector value;
		value.position() =  _data.back().position();
		value.position() += (*delta) * _data.back().field().norm();
		value.tetIndex() = _fastField.geometry()->findTetrahedron(value.position(),_data.back().tetIndex());

		// ***

		if (value.tetIndex() < 0 && _fastField.geometry()->isBoundaryTetrahedron(_data.back().tetIndex()))
		{
			value.field() = _fastField.compute(value.position());
			return value;
		}

		if (allowedTetrahedra.end() != allowedTetrahedra.find(value.tetIndex()))
		{
			value.field() = _fastField.compute(value.position());
			return value;
		}

		if (value.tetIndex() == _data.back().tetIndex())
			*delta *= 1.1f;
		else
			*delta *= 0.9f;
	}

	throw int(STEPPER_LIMIT_PASSED);
}
