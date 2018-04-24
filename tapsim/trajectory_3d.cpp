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

#include "trajectory_3d.h"

#include <iostream>
#include "debug.h"

#include <limits>
#include <stdexcept>
#include <cmath>
#include <cstring>

#include "../vector/vector.h"

const char* Trajectory_3d::integrator_str(const int value)
{
	switch (value)
	{
		case O1_INTEGRATOR:
			return "O1_INTEGRATOR";
		case O5_INTEGRATOR:
			return "O5_INTEGRATOR";
		default:
			return "INVALID";
	}
}

const char* Trajectory_3d::stepper_str(const int value)
{
	switch (value)
	{
		case ERROR_RESTRICTED:
			return "ERROR_RESTRICTED";
		case GEOMETRY_RESTRICTED:
			return "GEOMETRY_RESTRICTED";
		case ADAPTIVE:
			return "ADAPTIVE";

		default:
		case NONE:
			return "INVALID / NONDE";
	}
}

const char* Trajectory_3d::status_str(const int value)
{
	switch (value)
	{
		case NO_INTEGRATION:
			return "NO_INTEGRATION";
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

Trajectory_3d::Trajectory_3d(const Geometry_3d::Table* geomTable, const Grid_3d::Table* gridTable,const int intType, const int stepType)
	: _fastField(gridTable,geomTable),
	  _integratorType(intType),
	  _stepperType(stepType),
	  _stepperLimit(std::numeric_limits<unsigned long>::max()),
	  _iterationLimit(std::numeric_limits<unsigned long>::max()),
	  _timeStepLimit(std::numeric_limits<float>::min()),
	  _status(INVALID),
	  _charge(),
	  _mass(),
	  _data(),
	  _error_estimate()
{}

void Trajectory_3d::init(phaseVector tvp, const float charge, const float mass)
{
	if (_fastField.geometry() == 0)
		throw std::runtime_error("Trajectory_3d::init()");
	else
	{
		// set and check initial index for the tetrahedron

		if (tvp.tetIndex() < 0)
			tvp.tetIndex() = _fastField.geometry()->findTetrahedron(tvp.position());
		else
			tvp.tetIndex() = _fastField.geometry()->findTetrahedron(tvp.position(),tvp.tetIndex());

		if (tvp.tetIndex() < 0) throw std::runtime_error("Trajectory_3d::init(): initial point is outside");
	}

	// ***

	_status = NO_INTEGRATION;

	_data.clear();
	_data.push_back(tvp);

	_charge = charge;
	_mass = mass;

	_error_estimate = phaseVector();
}

// ***

int Trajectory_3d::integrate(float delta, const float errorLimit)
{
	if (!ready()) throw std::runtime_error("Trajectory_3d::integrate()");

	// ***

	if (_fastField.compute(_data.front().position(),_data.front().tetIndex()) == MathVector3d<float>(0.0f))
	{
		_status = NO_FIELD;
		return _status;
	}

	// ***

	_status = INVALID;

	unsigned long iteration(_iterationLimit);

	while (iteration--)
	{
		phaseVector error;

		try
		{
			phaseVector next;

			switch (_stepperType)
			{
				case NONE:
				{
					next = doIntegration(delta,&error);
					break;
				}
				case ERROR_RESTRICTED:
				{
					next = errorRestricted_stepper(&delta,errorLimit,&error);
					break;
				}
				case GEOMETRY_RESTRICTED:
				{
					std::set<int> allowedTetrahedra;
					_fastField.geometry()->adjacentTetrahedra(_data.back().tetIndex(),&allowedTetrahedra);
					next = geometryRestricted_stepper(&delta,allowedTetrahedra,&error);
					break;
				}
				case ADAPTIVE:
				{
					std::set<int> allowedTetrahedra;
					_fastField.geometry()->adjacentTetrahedra(_data.back().tetIndex(),&allowedTetrahedra);
					next = adaptive_stepper(&delta,errorLimit,allowedTetrahedra,&error);
					break;
				}
				default:
				{
					throw std::logic_error("Trajectory_3d::integrate() +++ invalid stepper-type!");
				}	
			}

			_data.push_back(next);
		}
		catch (int status)
		{
			_status = status;
			return status;
		}

		_error_estimate = cumulative_error(_error_estimate, error);

		// check global limits

		if (_data.back().tetIndex() < 0)
		{
			_status = SYSTEM_LIMIT_PASSED;
			return _status;
		}
	}

	_status = ITERATION_LIMIT_PASSED;
	return _status;
}

int Trajectory_3d::integrate(float delta, const float errorLimit, const float distance)
{
	if (!ready()) throw std::runtime_error("Trajectory_3d::integrate(const float)");

	// ***

	if (_fastField.compute(_data.front().position(),_data.front().tetIndex()) == MathVector3d<float>(0.0f))
	{
		_status = NO_FIELD;
		return _status;
	}

	// ***

	_status = INVALID;

	unsigned long iteration(_iterationLimit);

	while (iteration--)
	{
		phaseVector error;

		try
		{
			phaseVector next;

			switch (_stepperType)
			{
				case NONE:
				{
					next = doIntegration(delta,&error);
					break;
				}
				case ERROR_RESTRICTED:
				{
					next = errorRestricted_stepper(&delta,errorLimit,&error);
					break;
				}
				case GEOMETRY_RESTRICTED:
				{
					std::set<int> allowedTetrahedra;
					_fastField.geometry()->adjacentTetrahedra(_data.back().tetIndex(),&allowedTetrahedra);
					next = geometryRestricted_stepper(&delta,allowedTetrahedra,&error);
					break;
				}
				case ADAPTIVE:
				{
					std::set<int> allowedTetrahedra;
					_fastField.geometry()->adjacentTetrahedra(_data.back().tetIndex(),&allowedTetrahedra);
					next = adaptive_stepper(&delta,errorLimit,allowedTetrahedra,&error);
					break;
				}
				default:
				{
					throw std::logic_error("Trajectory_3d::integrate(const float) +++ invalid stepper-type!");
				}	
			}

			_data.push_back(next);
		}
		catch (int status)
		{
			_status = status;
			return status;
		}

		_error_estimate = cumulative_error(_error_estimate, error);

		// check global limits

		if (_data.back().tetIndex() < 0)
		{
			_status = SYSTEM_LIMIT_PASSED;
			return _status;
		}

		// check user limits

		Vector3d position = _data.front().position();
		position -= _data.back().position();

		if (position.length() > distance)
		{
			_status = GEOMETRY_LIMIT_PASSED;
			return _status;
		}
	}

	_status = ITERATION_LIMIT_PASSED;
	return _status;
}

bool Trajectory_3d::extrapolate(const int vecIndex, const float value)
{
	if (_data.empty()) return false;

	// ***

	float deltaTime(value);
	deltaTime -= _data.back().position(vecIndex);
	deltaTime /= _data.back().velocity(vecIndex);

	phaseVector item;

	item.time() = _data.back().time();
	item.time() += deltaTime;

	item.position() = _data.back().position();
	item.position() += deltaTime * _data.back().velocity();

	item.velocity() = _data.back().velocity();

	item.tetIndex() = _fastField.geometry()->findTetrahedron(item.position(),_data.back().tetIndex());

	_data.push_back(item);

	// ***

	if (deltaTime > 0.0f) _error_estimate.position() += deltaTime * _error_estimate.velocity();

	return true;
}

void Trajectory_3d::reset()
{
	const phaseVector tmp(_data.front());

	_status = NO_INTEGRATION;

	_data.clear();
	_data.push_back(tmp);

	_error_estimate = phaseVector();

	_fastField.reset();
}

bool Trajectory_3d::ready() const
{
	if (_fastField.geometry() == 0 || _fastField.grid() == 0) return false;
	if (_status != NO_INTEGRATION || _data.back().tetIndex() < 0) return false;

	return true;
}

// ***

Trajectory_3d::phaseVector Trajectory_3d::errorRestricted_stepper(float *delta, const float errorLimit, phaseVector* error)
{
	unsigned long iteration(_stepperLimit);

	while (iteration--)
	{
		if (*delta < _timeStepLimit) *delta = _timeStepLimit;

		Trajectory_3d::phaseVector value = doIntegration(*delta,error);

		// find maximum relative position error

		float relative_error(0.0f);

		for (short i = Vector3d::X; i <= Vector3d::Z; i++)
		{
			float tmp;

			if (value.position(i) != 0.0f)
			{
				tmp = error->position(i);
				tmp /= value.position(i);
				tmp = std::fabs(tmp);

				if (relative_error < tmp) relative_error = tmp;
			}

			if (value.velocity(i) != 0.0f)
			{
				tmp = error->velocity(i);
				tmp /= value.velocity(i);
				tmp = std::fabs(tmp);

				if (relative_error < tmp) relative_error = tmp;
			}
		}

		// compare error with preset limit

		if (relative_error < errorLimit || *delta == _timeStepLimit)
		{
			// readjust _delta and accept the current result

			if (relative_error != 0.0f)
			{
				float tmp(*delta);
				tmp *= std::pow(errorLimit/relative_error, 0.25f);
				tmp *= _safetyFactor;

				(tmp < *delta) ? *delta *= 1.1f : *delta = tmp;
			}

			return value;
		}
		else
		{
			// recompute/decrease _delta and try again

			float tmp(*delta);
			tmp *= std::pow(errorLimit/relative_error, 0.25f);
			tmp *= _safetyFactor;

			(tmp > *delta) ? *delta *= 0.9f : *delta = tmp;
		}
	}

	throw int(STEPPER_LIMIT_PASSED);
}

Trajectory_3d::phaseVector Trajectory_3d::geometryRestricted_stepper(float *delta, const std::set<int>& allowedTetrahedra, phaseVector* error)
{
	unsigned long iteration(_stepperLimit);

	while (iteration--)
	{
		if (*delta < _timeStepLimit) *delta = _timeStepLimit;

		Trajectory_3d::phaseVector value = doIntegration(*delta,error);
		if (value.tetIndex() < 0 && _fastField.geometry()->isBoundaryTetrahedron(_data.back().tetIndex()))  return value;

		if (allowedTetrahedra.end() != allowedTetrahedra.find(value.tetIndex())) return value;

		if (value.tetIndex() == _data.back().tetIndex())
			*delta *= 1.1f;
		else
		{
			if (*delta == _timeStepLimit) return value;

			*delta *= 0.9f;
		}
	}

	throw int(STEPPER_LIMIT_PASSED);
}

Trajectory_3d::phaseVector Trajectory_3d::adaptive_stepper(float* delta, const float errorLimit, const std::set<int>& allowedTetrahedra, Trajectory_3d::phaseVector* error)
{
	unsigned long iteration(_stepperLimit);

	while (iteration--)
	{
		if (*delta < _timeStepLimit) *delta = _timeStepLimit;

		Trajectory_3d::phaseVector value = doIntegration(*delta,error);

		// I) check spatial limits

		if (value.tetIndex() >= 0)
		{
			if (value.tetIndex() != _data.back().tetIndex())
			{
				if (allowedTetrahedra.end() == allowedTetrahedra.find(value.tetIndex()))
				{
					*delta *= 0.1f;
					continue;
				}
			}
		}

		// II) compare error with preset limit

		float relative_error(0.0f);

		for (short i = Vector3d::X; i <= Vector3d::Z; i++)
		{
			float tmp;

			if (value.position(i) != 0.0f)
			{
				tmp = error->position(i);
				tmp /= value.position(i);
				tmp = std::fabs(tmp);

				if (relative_error < tmp) relative_error = tmp;
			}

			if (value.velocity(i) != 0.0f)
			{
				tmp = error->velocity(i);
				tmp /= value.velocity(i);
				tmp = std::fabs(tmp);

				if (relative_error < tmp) relative_error = tmp;
			}
		}

		if (relative_error < errorLimit  || *delta == _timeStepLimit)
		{
			//  readjust _delta and accept the current result
	
			if (relative_error != 0.0f)
			{
				float tmp(*delta);
				tmp *= std::pow(errorLimit/relative_error, 0.25f);
				tmp *= _safetyFactor;

				(tmp < *delta) ? *delta *= 1.1f : *delta = tmp;
			}

			return value;
		}
		else
		{
			// recompute/decrease _delta and try again

			float tmp(*delta);
			tmp *= std::pow(errorLimit/relative_error, 0.25f);
			tmp *= _safetyFactor;

			(tmp > *delta) ? *delta *= 0.9f : *delta = tmp;
		}
	}

	throw int(STEPPER_LIMIT_PASSED);
}

// ***

Trajectory_3d::phaseVector Trajectory_3d::doIntegration(const float delta, phaseVector* error) const
{
	switch (_integratorType)
	{
		case O1_INTEGRATOR:
			return o1_integrator(delta,error);
		case O5_INTEGRATOR:
			return o5_integrator(delta,error);
		default:
			throw std::logic_error("Trajectory_3d::doIntegration() +++ invalid integrator-type!");
	}
}

Trajectory_3d::phaseVector Trajectory_3d::o1_integrator(const float delta, phaseVector* error) const
{
	Vector3d field = _fastField.compute(_data.back().position(),_data.back().tetIndex());
	field *= _charge;
	field /= _mass;

	phaseVector next = _data.back();

	for (short i = Vector3d::X; i <= Vector3d::Z; i++)
	{
		next.velocity(i) += delta * field[i];
		next.position(i) += delta * _data.back().velocity(i);
	}
	
	next.time() += delta;

	if (error != 0)
	{
		error->time() = 0.0f;

		error->position() = next.position();
		error->position() *= std::pow(delta,2.0f);

		error->velocity() = next.velocity();
		error->velocity() *= std::pow(delta,2.0f);
	}

	// ***

	next.tetIndex() = _fastField.geometry()->findTetrahedron(next.position(),_data.back().tetIndex());

	return next;
}

Trajectory_3d::phaseVector Trajectory_3d::o5_integrator(const float delta, phaseVector* error) const
{
	// !!! implementation of the (embedded) fifth-order runge-kutta method !!!
	// !!! => for details see: "Numerical Recipes, Integration of Ordinary Differential Equations" !!!

	// static const float cashKarp_a[6] =
	// 	{ 0.0, 1.0 / 5.0, 3.0 / 10.0, 3.0 / 5.0, 1.0, 7.0 / 8.0 };

	static const float cashKarp_b[6][5] =
		{
			{               0.0f,            0.0f,              0.0f,                 0.0f,             0.0f },
			{    1.0f /     5.0f,            0.0f,              0.0f,                 0.0f,             0.0f },
			{    3.0f /    40.0f,   9.0f /  40.0f,              0.0f,                 0.0f,             0.0f },
			{    3.0f /    10.0f,  -9.0f /  10.0f,   6.0f /     5.0f,                 0.0f,             0.0f },
			{  -11.0f /    54.0f,   5.0f /   2.0f, -70.0f /    27.0f,    35.0f /     27.0f,             0.0f },
			{ 1631.0f / 55296.0f, 175.0f / 512.0f, 575.0f / 13824.0f, 44275.0f / 110592.0f, 253.0f / 4096.0f }
		};

	static const float cashKarp_c1[6] =
		{ 37.0f / 378.0f, 0.0f, 250.0f / 621.0f, 125.0f / 594.0f, 0.0f, 512.0f / 1771.0f };

	static const float cashKarp_c2[6] =
		{ 2825.0f / 27648.0f, 0.0f, 18575.0f / 48384.0f, 13525.0f / 55296.0f, 277.0f / 14336.0f, 1.0f / 4.0f };

	const float charge_mass_ratio = _charge / _mass;

	struct
	{
		Vector3d velocity;
		Vector3d position;
	} dydx[6];

	for (short i = 0; i < 6; i++)
	{
		Vector3d velocity = _data.back().velocity();
		Vector3d position = _data.back().position();

		for (short j = 0; j < i; j++)
		{
			for (short l = Vector3d::X; l <= Vector3d::Z; l++)
			{
				velocity[l] += cashKarp_b[i][j] * dydx[j].velocity[l];
				position[l] += cashKarp_b[i][j] * dydx[j].position[l];
			}
		}

		dydx[i].velocity = _fastField.compute(position,_data.back().tetIndex());
		dydx[i].velocity *= charge_mass_ratio;
		dydx[i].velocity *= delta;

		dydx[i].position = velocity;
		dydx[i].position *= delta;
	}

	// *** final solution

	phaseVector next = _data.back();

	for (short i = Vector3d::X; i <= Vector3d::Z; i++)
	{
		for (short j = 0; j < 6; j++)
		{
			next.position(i) += cashKarp_c1[j] * dydx[j].position[i];
			next.velocity(i) += cashKarp_c1[j] * dydx[j].velocity[i];
		}
	}

	next.time() += delta;

	// *** error estimation

	if (error != 0)
	{
		error->time() = 0.0f;

		for (short i = Vector3d::X; i <= Vector3d::Z; i++)
		{
			error->position(i) = 0.0f;
			error->velocity(i) = 0.0f;

			for (short j = 0; j < 6; j++)
			{
				error->position(i) += (cashKarp_c1[j] - cashKarp_c2[j]) * dydx[j].position[i];
				error->velocity(i) += (cashKarp_c1[j] - cashKarp_c2[j]) * dydx[j].velocity[i];
			}
		}
	}

	// ***

	next.tetIndex() = _fastField.geometry()->findTetrahedron(next.position(),_data.back().tetIndex());

	return next;
}

// ***

Trajectory_3d::phaseVector Trajectory_3d::cumulative_error(const Trajectory_3d::phaseVector& current, const Trajectory_3d::phaseVector& next) const
{
	phaseVector value;

	for (short i = Vector3d::X; i <= Vector3d::Z; i++)
	{
		value.position(i) = std::pow(current.position(i),2.0f);
		value.position(i) += std::pow(next.position(i),2.0f);
		value.position(i) = std::sqrt(value.position(i));

		value.velocity(i) = std::pow(current.velocity(i),2.0f);
		value.velocity(i) += std::pow(next.velocity(i),2.0f);
		value.velocity(i) = std::sqrt(value.velocity(i));
	}

	return value;
}
