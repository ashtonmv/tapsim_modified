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

#include "surface_3d.h"

#include <stdexcept>
#include <cstdlib>
#include <cmath>

#include "info.h" /* added for the status informatin in evap_takeAway() */
#include <list>

#include "utils.h"


/* ********** ----------- ********** */

Surface_3d::Table::Table()
	: _nodes(),
	  _vacuumId(-1),
	  _scalingReference(1.0f)
{}

Surface_3d::Table::Table(const System_3d& system, const Configuration::NodeId vacuumId)
	: _nodes(),
	  _vacuumId(vacuumId),
	  _scalingReference(1.0f)
{
	init(system);
}

void Surface_3d::Table::init(const System_3d& system)
{
	// *** local switch, not included in the executable's general parameter list
	// *** => determines distinction between removable and non removable surface nodes

	const bool acceptRemovable(true);

	// ***

	_nodes.clear();

	for (int i = 0; i < system.gridTable.numNodes(); i++)
	{
		const Configuration::NodeId id = system.gridTable.id(i);

		if (id == _vacuumId) continue;

		if (!system.configTable[id].removable() && acceptRemovable) continue;

		for (int j = 0; j < system.gridTable.node(i).numNeighbours(); j++)
		{
			const int neighIndex = system.gridTable.node(i).neighbour(j);

			if (system.gridTable.id(neighIndex) == _vacuumId && system.geomTable.voronoiArea(i,neighIndex,1.0f) > 0.0f)
			{
				const Node obj(i,0.0);
				_nodes.insert(obj);
				break;
			}
		}
	}
}

void Surface_3d::Table::update(const int nodeIndex, const System_3d& system)
{
	// *** local switch, not included in the executable's general parameter list
	// *** => determines distinction between removable and non removable surface nodes

	const bool acceptRemovable(true);

	// ***

	if (_nodes.erase(nodeIndex) != 1) throw std::runtime_error("Surface_3d::Table::update()");

	// ***

	for (int i = 0; i < system.gridTable.node(nodeIndex).numNeighbours(); i++)
	{
		const int neighIndex = system.gridTable.node(nodeIndex).neighbour(i);

		if (system.gridTable.id(neighIndex) == _vacuumId) continue;

		if (!system.configTable[system.gridTable.id(neighIndex)].removable() && acceptRemovable) continue;

		// *** (redundant) test if the neighbour node has at least one vacuum-node in
		// *** its neighbourhood

		for (int j = 0; j < system.gridTable.node(neighIndex).numNeighbours(); j++)
		{
			const int nnIndex = system.gridTable.node(neighIndex).neighbour(j);

			if (system.gridTable.id(nnIndex) == _vacuumId && system.geomTable.voronoiArea(neighIndex,nnIndex,1.0f) > 0.0f)
			{
				const Node obj(neighIndex,0.0);
				_nodes.insert(obj);
				break;
			}
		}
	}
}

/* ********** ----------- ********** */

Geometry_3d::Point Surface_3d::Table::normal(const int nodeIndex, const System_3d& system) const
{
	if (_nodes.find(nodeIndex) == _nodes.end()) return Geometry_3d::Point(0.0f);

	std::set<int> localNodes;
	system.geomTable.adjacentNodes(nodeIndex,&localNodes,3); // limit search to 3rd neighbours

	std::vector<int> vacNodes;
	vacNodes.reserve(localNodes.size());

	{
		std::set<int>::const_iterator i = localNodes.begin();

		while (localNodes.end() != i)
		{
			if (system.gridTable.node(*i).id() == _vacuumId) vacNodes.push_back(*i);

			if (_nodes.find(*i) == _nodes.end())
			{
				const std::set<int>::const_iterator j(i);

				advance(i,1);
				localNodes.erase(j);
			}
			else
				advance(i,1);
		}
	}

	// ***

	std::vector<Geometry_3d::Point> coords;
	for (std::set<int>::const_iterator i = localNodes.begin(); i != localNodes.end(); i++)
		coords.push_back(system.geomTable.nodeCoords(nodeIndex) - system.geomTable.nodeCoords(*i));

	Geometry_3d::Point value = Geometry_3d::find_best_plane(coords);

	// *** adjust orientation

	Geometry_3d::Point tmp;
	for (std::vector<int>::const_iterator i = vacNodes.begin(); i != vacNodes.end(); i++)
		tmp += system.geomTable.nodeCoords(*i) - system.geomTable.nodeCoords(nodeIndex);

	if (value * tmp < 0.0) value *= -1.0f;

	return value;
}

Surface_3d::Nodeset::const_iterator Surface_3d::Table::apex(const System_3d& system) const
{
	Surface_3d::Nodeset::const_iterator topNode = _nodes.begin();
	float topCoordinate = system.geomTable.nodeCoords(topNode->index()).z();

	// ***

	for (Nodeset::const_iterator i = _nodes.begin(); i != _nodes.end(); i++)
	{
		const float coordinate = system.geomTable.nodeCoords(i->index()).z();

		if (coordinate > topCoordinate)
		{
			topCoordinate = coordinate;
			topNode = i;
		}
	}

	return topNode;
}

/* ********** ----------- ********** */

void Surface_3d::evap_compute_specificFields(Surface_3d::Table* surfaceTable, const System_3d& system)
{
	/* probably this is just a bad hack ... */

	if (surfaceTable == 0 || surfaceTable->nodes().empty()) return;

	Surface_3d::Nodeset::iterator i = surfaceTable->nodes().begin();

	float value = system.gridTable.field_o1(i->index(),system.geomTable).length();
	value /= system.configTable[system.gridTable.id(i->index())].evapField();
	advance(i,1);

	float maxValue(value);

	while (surfaceTable->nodes().end() != i)
	{
		value = system.gridTable.field_o1(i->index(),system.geomTable).length();
		value /= system.configTable[system.gridTable.id(i->index())].evapField();

		if (maxValue < value) maxValue = value;

		advance(i,1);
	}

	surfaceTable->setScalingReference(maxValue);
}

namespace
{
	inline void probBoltzmann(Surface_3d::Table* surfaceTable, const System_3d& system)
	{
		/*
		*** The problem of picking the correct distribution for the field dependant evaporation
		*** probability using the Boltzmann Equation is solved by assuming that the surface atom
		*** exposed to the highest field with respect to its specific evaporation field strength
		*** has a evaporation probability of 1.0. This means the activation barrier vanishes.
		*/

		if (system.configTable.temperature() <= 0.0)
			throw std::runtime_error("probBoltzmann(): Cannot compute probability with chosen temperature value!");


		const float kBoltzmann = 8.6173431e-5; // [kBoltzmann] =  eV/K

		float sum(0.0);

		for (Surface_3d::Nodeset::iterator i = surfaceTable->nodes().begin(); i != surfaceTable->nodes().end(); i++)
		{
			float probability(-1.0);
			probability *= system.gridTable.field_o1(i->index(),system.geomTable).length();
			probability /= system.configTable[system.gridTable.id(i->index())].evapField();
			probability += 1.0;
			probability *= system.configTable[system.gridTable.id(i->index())].evapEnergy();
			probability /= system.configTable.temperature();
			probability /= kBoltzmann;
			probability /= -2.0;
			probability = std::exp(probability);
			const_cast<Surface_3d::Node&>(*i).setProbability(probability);

			sum += probability;
		}

		// *** normalize

		for (Surface_3d::Nodeset::iterator i = surfaceTable->nodes().begin(); i != surfaceTable->nodes().end(); i++)
		{
			float value = i->probability();
			value /= sum;

			const_cast<Surface_3d::Node&>(*i).setProbability(value);
		}
	}

	inline void probLinear_field(Surface_3d::Table* surfaceTable, const System_3d& system)
	{
		/*
		*** This method regards different evaporation probabilities for species with
		*** heterogeneous evaporation field strengths. This is achieved by rescaling
		*** of the computed field strength with the inverse of the specific evaporation
		*** field strength. The evaporation probability is computed as the ratio of this.
		*/

		float sum(0.0);

		for (Surface_3d::Nodeset::iterator i = surfaceTable->nodes().begin(); i != surfaceTable->nodes().end(); i++)
		{

			float value = system.gridTable.field_o1(i->index(),system.geomTable).length();
			value /= system.configTable[system.gridTable.id(i->index())].evapField();
			const_cast<Surface_3d::Node&>(*i).setProbability(value);

			sum += i->probability();
		}

		// *** normalize

		for (Surface_3d::Nodeset::iterator i = surfaceTable->nodes().begin(); i != surfaceTable->nodes().end(); i++)
		{
			float value = i->probability();
			value /= sum;

			const_cast<Surface_3d::Node&>(*i).setProbability(value);
		}
	}

	inline void probLinear_force(Surface_3d::Table* surfaceTable, const System_3d& system)
	{
		// computes the probabiliy for field evapoation of each atom by taking into account
		// the acting force on the atom, calculation is started from scratch ignoring any
		// formerly conducted calculations (e.g. from evap_compute_specificFields() above)

		float sum(0.0);

		for (Surface_3d::Nodeset::iterator i = surfaceTable->nodes().begin(); i != surfaceTable->nodes().end(); i++)
		{
			float fieldstrength, charge;
			fieldstrength = system.gridTable.field_o2(i->index(),system.geomTable,&charge).length();

			float probValue = charge;
			probValue *= fieldstrength;
			probValue /= std::pow(system.configTable[system.gridTable.node(i->index()).id()].evapField(),2.0f);

			const_cast<Surface_3d::Node&>(*i).setProbability(probValue);

			sum += probValue;
		}

		// *** normalize

		for (Surface_3d::Nodeset::iterator i = surfaceTable->nodes().begin(); i != surfaceTable->nodes().end(); i++)
		{
			float value = i->probability();
			value /= sum;

			const_cast<Surface_3d::Node&>(*i).setProbability(value);
		}
	}

	inline void probVoronoiFluxForce(Surface_3d::Table* surfaceTable, const System_3d& system)
	{
		const float epsilon0 = 8.85418781762e-12; // [epsilon0] = As/(Vm)

		float sum(0.0);

		for (Surface_3d::Nodeset::iterator i = surfaceTable->nodes().begin(); i != surfaceTable->nodes().end(); i++)
		{
			const Grid_3d::Node& node = system.gridTable.node(i->index());

			// *** BEGIN BUGFIX +++ PROPER SCALING
			float surfaceArea(0.0f);
			// *** END BUGFIX +++ PROPER SCALING

			MathVector3d<float> force(0.0f);
			for (int j = 0; j < node.numNeighbours(); j++)
			{
				const int neighIndex = node.neighbour(j);

				float value = system.gridTable.potential(neighIndex);
				value -= system.gridTable.potential(i->index());
				value *= value;
				value *= -1.0 * epsilon0;
				value *= node.coupling(j);

				const float epsilon_1 = 1.0 / system.configTable[system.gridTable.node(i->index()).id()].epsilon();
				const float epsilon_2 = 1.0 / system.configTable[system.gridTable.node(neighIndex).id()].epsilon();

				value *= epsilon_2 - epsilon_1;

				// *** BEGIN BUGFIX +++ PROPER SCALING
				float corrArea = node.coupling(j);
				corrArea *= system.geomTable.distance(i->index(),neighIndex);
				corrArea /= 2.0;
				corrArea *= epsilon_1 + epsilon_2;
				surfaceArea += corrArea;
				// *** END BUGFIX +++ PROPER SCALING

				MathVector3d<float> localForce = system.geomTable.nodeCoords(neighIndex);
				localForce -= system.geomTable.nodeCoords(i->index());
				localForce /= localForce*localForce;
				localForce *= value;

				force += localForce;
			}

			float probValue = force.length();

			// *** BEGIN BUGFIX +++ PROPER SCALING
			probValue /= surfaceArea;
			// *** END BUGFIX +++ PROPER SCALING

			probValue /= std::pow(system.configTable[system.gridTable.node(i->index()).id()].evapField(),2.0f);
			const_cast<Surface_3d::Node&>(*i).setProbability(probValue);

			sum += probValue;
		}

		// *** normalize

		for (Surface_3d::Nodeset::iterator i = surfaceTable->nodes().begin(); i != surfaceTable->nodes().end(); i++)
		{
			float value = i->probability();
			value /= sum;

			const_cast<Surface_3d::Node&>(*i).setProbability(value);
		}
	}
}

void Surface_3d::evap_compute_probabilities(const int mode, Table* surfaceTable, const System_3d& system)
{
	switch (mode)
	{
		case PROB_BOLTZMANN:
			return probBoltzmann(surfaceTable,system);
		case PROB_LINEAR_FIELD:
			return probLinear_field(surfaceTable,system);
		case PROB_LINEAR_FORCE:
			return probLinear_force(surfaceTable,system);
		case PROB_VORONOI_FLUX_FORCE:
			return probVoronoiFluxForce(surfaceTable,system);
		default:
			throw std::runtime_error("Surface_3d::evap_compute_probabilities()");
	}
}

/* ********** ----------- ********** */

namespace
{
	inline Surface_3d::Nodeset::const_iterator evapMaximum(const Surface_3d::Table& surfaceTable)
	{
		Surface_3d::Nodeset::const_iterator pickNode = surfaceTable.nodes().begin();

		for (Surface_3d::Nodeset::const_iterator i = surfaceTable.nodes().begin(); i != surfaceTable.nodes().end(); i++)
			if (pickNode->probability() < i->probability()) pickNode = i;

		return pickNode;
	}

	inline Surface_3d::Nodeset::const_iterator evapMonteCarlo(const Surface_3d::Table& surfaceTable)
	{
		double pickThreshold = std::rand();
		pickThreshold /= RAND_MAX;

		Surface_3d::Nodeset::const_iterator pickCell = surfaceTable.nodes().begin();

		double threshold = pickCell->probability();

		while (pickThreshold > threshold)
		{
			pickCell++;
			threshold += pickCell->probability();
		}

		return pickCell;
	}
}

Surface_3d::Nodeset::const_iterator Surface_3d::evap_findCandidate(const int mode, const Table& surfaceTable)
{
	switch (mode)
	{
		case EVAP_MAXIMUM:
			return evapMaximum(surfaceTable);
		case EVAP_MONTE_CARLO:
			return evapMonteCarlo(surfaceTable);
		default:
			throw std::runtime_error("Surface_3d::evap_findCandidate()");
	}
}

/* ********** ----------- ********** */

MathVector3d<float> evap_initialPosition(const int nodeIndex, const Geometry_3d::Table& geomTable)
{
	return geomTable.node(nodeIndex).coords;
}

MathVector3d<float> evap_initialVelocity(const Configuration::NodeId& id, const Configuration::Table& configTable)
{
	const double kBoltzmann = 1.3806504e-23;		// [kBoltzmann] = J/K
	const double amu2kg = 1.660538782e-27;			// [amu2kg] = kg/amu

	double max(2.0);
	max *= kBoltzmann;					// [max] = J/K
	max *= configTable.temperature();			// [max] = J
	max /= configTable[id].mass();				// [max] = kg * m² / s² / amu
	max /= amu2kg;						// [max] = m² / s²

	MathVector3d<float> velocity;				// [velocity] = 1

	velocity.x() = std::rand();
	velocity.x() /= RAND_MAX;
	velocity.x() *= max;					// [velocity.x()] = m²/s²

	velocity.y() = std::rand();
	velocity.y() /= RAND_MAX;
	velocity.y() *= max - velocity.x();			// [velocity.y()] = m²/s²

	velocity.z() = max;
	velocity.z() -= velocity.x();
	velocity.z() -= velocity.y();				// [velocity.z()] = m²/s²

	for (short i = MathVector3d<float>::X; i <= MathVector3d<float>::Z; i++)
	{
		velocity[i] = std::sqrt(velocity[i]);
		if (std::rand() % 2) velocity[i] *= -1.0f;
	}
	return velocity;
}

/* ********** ----------- ********** */

void evap_takeAway(System_3d* system, Surface_3d::Table* surfaceTable, const int nodeIndex)
{
	static int lonelyCnt(0);

	// remove atom from surface

	system->gridTable.node(nodeIndex).setId(surfaceTable->vacuumId());
	system->gridTable.fastResync(nodeIndex,system->geomTable,system->configTable);

	surfaceTable->update(nodeIndex,*system);

	// check for detached surface sites in the vicinity of the former atom location

	{
		std::list<int> lonelySites;

		for (int i = 0; i < system->gridTable.node(nodeIndex).numNeighbours(); i++)
		{
			const int neighIndex = system->gridTable.node(nodeIndex).neighbour(i);
			if (system->gridTable.node(neighIndex).id() == surfaceTable->vacuumId()) continue;

			bool flag(true);
			for (int j = 0; j < system->gridTable.node(neighIndex).numNeighbours(); j++)
			{
				const int index = system->gridTable.node(neighIndex).neighbour(j);
				if (system->gridTable.node(index).id() != surfaceTable->vacuumId())
				{
					flag = false;
					break;
				}
			}

			if (flag) lonelySites.push_back(neighIndex);
		}

		if (lonelySites.size() != 0)
		{
			lonelyCnt += lonelySites.size();

			info::out() << "*** Detected detached surface atoms:" << std::endl;
			info::out() << "\t-> node = #" << nodeIndex << std::endl;
			info::out() << "\t-> coordinates = " << system->geomTable.nodeCoords(nodeIndex) << std::endl;
			info::out() << "\t-> type = '" << system->configTable[system->gridTable.id(nodeIndex)].name() << "'" << std::endl;
			info::out() << "\t-> actual number of detached atoms = " << lonelySites.size() << std::endl;
			info::out() << "\t-> total number of detached atoms so far: " << lonelyCnt << std::endl;

			// remove detached surface atoms

			for (std::list<int>::const_iterator i = lonelySites.begin(); i != lonelySites.end(); i++)
			{
				system->gridTable.node(*i).setId(surfaceTable->vacuumId());
				system->gridTable.fastResync(*i,system->geomTable,system->configTable);

				surfaceTable->update(*i,*system);
			}
		}
	}
}

/* ********** ----------- ********** */

void update_coordination(System_3d* system)
{
	int nNodes = system->gridTable.numNodes();
	for (int i=0; i < nNodes; i++)
	{
	    if (system->gridTable.node(i).id().toValue() >= 10)
	    {
	        int numNeighbours = 0;
	        for (int j=0; j < system->gridTable.node(i).numNeighbours(); j++)
	        {
	            int neighborIndex = system->gridTable.node(i).neighbour(j);
	            if (system->gridTable.node(neighborIndex).id().toValue() >= 10)
	            {
	                numNeighbours++;
	            }
	        }
	        int newId = system->gridTable.node(i).id().toValue() + numNeighbours;
	        system->gridTable.node(i).setId(newId);
	        //if (newId > 19)
					//{
					info::out() << newId << std::endl;
					system->gridTable.node(i).setId(19);
					//}
	    }
	}
}
