#ifndef TAPSIM_CONFIGURATION_H
#define TAPSIM_CONFIGURATION_H

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

#include <map>
#include <list>
#include <string>

#include "../vector/vector.h"

namespace Configuration
{
	class NodeId
	{
		public:
			NodeId() : _value(-1) {}
			NodeId(char value) : _value(value) {};

			~NodeId() {}
	
			bool isValid() const { return _value >= 0; }
	
			short toValue() const { return static_cast<short>(_value); }
	
			bool operator<(const NodeId& rhs) const { return _value < rhs._value; }
			bool operator==(const NodeId& rhs) const { return _value == rhs._value; }
			bool operator!=(const NodeId& rhs) const { return _value != rhs._value; }
	
			NodeId& operator++() { _value++; return *this; }
			NodeId& operator++(int) { _value++; return *this; }
	
			NodeId& operator--() { _value--; return *this; }
			NodeId& operator--(int) { _value--; return *this; }
	
			static NodeId invalid() { return NodeId(-1); }
			static NodeId fromValue(const short);

			static NodeId max() { return NodeId(0xff >> 3); } /* maximum id is constrained to 31 - why? */
			static NodeId min() { return NodeId(0); }
	
		private:
			char _value;
	};

	class NodeNumber
	{
		public:
			NodeNumber() : _value(0) {}
			NodeNumber(unsigned int value) { _value = value; }

			~NodeNumber() {}
	
			unsigned int toValue() const { return _value; }

			bool operator<(const NodeNumber& obj) const { return _value < obj._value; }
			bool operator<=(const NodeNumber& obj) const { return _value <= obj._value; }
			bool operator>(const NodeNumber& obj) const { return _value > obj._value; }
			bool operator>=(const NodeNumber& obj) const { return _value >= obj._value; }
			bool operator==(const NodeNumber& obj) const { return _value == obj._value; }
			bool operator!=(const NodeNumber& obj) const { return _value != obj._value; }
			
			bool isValid() const { return _value != 0; }
	
			NodeNumber next() const { if (_value != 0) return NodeNumber(_value + 1); else return invalid(); }
	
			static NodeNumber first() { return NodeNumber(1); }
			static NodeNumber invalid() { return NodeNumber(0); }
			static NodeNumber fromValue(const unsigned int value) { return NodeNumber(value); }
	
		private:
			unsigned int _value;
	};
	
	class NodeData
	{
		public:
			NodeData();
	
			// ***
	
			void setId(const NodeId value) { _id = value; }
			NodeId id() const { return _id; }
	
			void setName(const std::string& value) { _name = value; }
			const std::string& name() const { return _name; }
	
			// ***
	
			void setChargeDensity(const float value) { _chargeDensity = value; }
			float chargeDensity() const { return _chargeDensity; }
	
			void setEpsilon(const float value) { _epsilon = value; }
			float epsilon() const { return _epsilon; }

			void setRemovable(const bool value) { _removable = value; }
			bool removable() const { return _removable; }
	
			void setNeumannBoundary(const bool value) { _neumannBoundary = value; }
			bool neumannBoundary() const { return _neumannBoundary; }

			void setDirichletBoundary(const bool value) { _dirichletBoundary = value; }
			bool dirichletBoundary() const { return _dirichletBoundary; }
	
			void setPhi(const float value) { _phi = value; }
			float phi() const { return _phi; }
	
			void setMass(const float value) { _mass = value; }
			float mass() const { return _mass; }
	
			void setEvapCharge(const short value) { _evapCharge = value; }
			short evapCharge() const { return _evapCharge; }
	
			void setEvapField(const float value) { _evapField = value; }
			float evapField() const { return _evapField; }
	
			void setEvapEnergy(const float value) { _evapEnergy = value; }
			float evapEnergy() const { return _evapEnergy; }
	
			void reset() { *this = NodeData(); }
	
		private:
	
			NodeId _id;
	
			std::string _name;	// short description

			float _chargeDensity;	// [_chargeDensity] = C/m², intrinsic charge-density
			float _epsilon;		// [_epsilon] = 1, relative dielectric constant

			bool _removable;

			bool _neumannBoundary;	// Neumann boundary condition
			bool _dirichletBoundary; // Dirichlet boundary condition

			float _phi;		// [_phi] = V, initial potential value

			float _mass;		// [_mass] = amu, intrinsic atomic mass
		
			short _evapCharge;	// [_evapCharge] = 1, evaporation charge-state
			float _evapField;	// [_evapField] = V/m, evaporation field-strength
			float _evapEnergy;	// [_evapEnergy]] = eV, activation energy for field-evaporation
	};

	class Table
	{
		public:
			Table();

			// ***

			void operator<<(std::ifstream&);
			void operator>>(std::ofstream&) const;

			// ***

			
			void setTemperature (const float value) { _temperature = value; }
			float temperature() const { return _temperature; }
			
			// ***
			
			const NodeData& operator[](const NodeId&) const;
			NodeData& operator[](const NodeId&);
	
			const NodeData& operator[](const std::string&) const;
			NodeData& operator[](const std::string& name);
	
			NodeId insert(const NodeData&);
			NodeId insert(const NodeId, const NodeData&);
	
			void remove(const NodeId&);
			void remove(const std::string&);
	
			bool exists(const NodeId& id) const { return _data.find(id) != _data.end(); }
			bool exists(const std::string& name) const { return _nameIndex.find(name) != _nameIndex.end(); }
	
			std::list<NodeId> ids() const;
			std::list<std::string> names() const;

			const std::map<NodeId,NodeData>& data() const { return _data; }
			
			void clear();

		private:
			NodeId getId() const;

			static const std::string _binaryVersion;

			float _temperature;

			std::map<NodeId,NodeData> _data;
			std::map<std::string,NodeId> _nameIndex;
	};
}

#endif
