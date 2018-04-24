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

#include "configuration.h"

#include <limits>
#include <fstream>
#include <stdexcept>
#include <set>

/* ************* NodeId ************* */

Configuration::NodeId Configuration::NodeId::fromValue(const short value)
{
	if (value < 0 || std::numeric_limits<char>::max() < value)
		return NodeId::invalid();
	else
		return NodeId(value);
}


/* ************* NodeData ************* */

Configuration::NodeData::NodeData()
	:
	  _id(),
	  _name(),
	  _chargeDensity(0.0f),
	  _epsilon(1.0f),
	  _removable(true),
	  _neumannBoundary(false),
	  _dirichletBoundary(false),
	  _phi(0.0f),
	  _mass(0.0f),
	  _evapCharge(0),
	  _evapField(0.0f),
	  _evapEnergy(0.0f)
{}

/* ************* Table ************* */

const std::string Configuration::Table::_binaryVersion("1.0");

Configuration::Table::Table()
	: _temperature(0.0f),
	  _data(),
	  _nameIndex()
{}

void Configuration::Table::operator<<(std::ifstream& stream)
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

	if (version != _binaryVersion) throw std::runtime_error("Configuration::Table<<(): wrong version!");

	// ***
	
	stream.read(reinterpret_cast<char*>(&_temperature),sizeof(float));

	clear();

	int size;
	stream.read(reinterpret_cast<char*>(&size),sizeof(int));

	for (int i = 0; i < size; i++)
	{
		NodeData item;

		// ***

		NodeId id;
		stream.read(reinterpret_cast<char*>(&id),sizeof(NodeId));
		item.setId(id);

		std::string name;
		{
			char c;
			stream.read(reinterpret_cast<char*>(&c),sizeof(char));

			while (c != '\0')
			{
				name.push_back(c);
				stream.read(reinterpret_cast<char*>(&c),sizeof(char));
			}
		}
		item.setName(name);

		float chargeDensity;
		stream.read(reinterpret_cast<char*>(&chargeDensity),sizeof(float));
		item.setChargeDensity(chargeDensity);

		float epsilon;
		stream.read(reinterpret_cast<char*>(&epsilon),sizeof(float));
		item.setEpsilon(epsilon);

		bool removable;
		stream.read(reinterpret_cast<char*>(&removable),sizeof(bool));
		item.setRemovable(removable);

		bool neumannBoundary;
		stream.read(reinterpret_cast<char*>(&neumannBoundary),sizeof(bool));
		item.setNeumannBoundary(neumannBoundary);

		bool dirichletBoundary;
		stream.read(reinterpret_cast<char*>(&dirichletBoundary),sizeof(bool));

		float phi;
		stream.read(reinterpret_cast<char*>(&phi),sizeof(float));
		item.setPhi(phi);

		float mass;
		stream.read(reinterpret_cast<char*>(&mass),sizeof(float));
		item.setMass(mass);

		short evapCharge;
		stream.read(reinterpret_cast<char*>(&evapCharge),sizeof(float));
		item.setEvapCharge(evapCharge);

		float evapField;
		stream.read(reinterpret_cast<char*>(&evapField),sizeof(float));
		item.setEvapField(evapField);

		float evapEnergy;
		stream.read(reinterpret_cast<char*>(&evapEnergy),sizeof(float));
		item.setEvapEnergy(evapEnergy);

		// ***

		if (!insert(item.id(),item).isValid()) throw std::runtime_error("Configuration::Table::operator<<()");
	}
}

void Configuration::Table::operator>>(std::ofstream& stream) const
{
	const char* version = _binaryVersion.c_str();
	stream.write(version,_binaryVersion.size()+1);

	// ***

	stream.write(reinterpret_cast<const char*>(&_temperature),sizeof(float));

	const int size = _data.size();
	stream.write(reinterpret_cast<const char*>(&size),sizeof(int));

	for (std::map<NodeId,NodeData>::const_iterator i = _data.begin(); i != _data.end(); i++)
	{
		NodeId id = i->second.id();
		stream.write(reinterpret_cast<const char*>(&id),sizeof(NodeId));

		const char* name = i->second.name().c_str();
		stream.write(reinterpret_cast<const char*>(name),i->second.name().size()+1);

		const float chargeDensity = i->second.chargeDensity();
		stream.write(reinterpret_cast<const char*>(&chargeDensity),sizeof(float));

		const float epsilon = i->second.epsilon();
		stream.write(reinterpret_cast<const char*>(&epsilon),sizeof(float));

		const bool removable = i->second.removable();
		stream.write(reinterpret_cast<const char*>(&removable),sizeof(bool));

		const bool neumannBoundary = i->second.neumannBoundary();
		stream.write(reinterpret_cast<const char*>(&neumannBoundary),sizeof(bool));

		const bool dirichletBoundary = i->second.dirichletBoundary();
		stream.write(reinterpret_cast<const char*>(&dirichletBoundary),sizeof(bool));

		const float phi = i->second.phi();
		stream.write(reinterpret_cast<const char*>(&phi),sizeof(float));

		const float mass = i->second.mass();
		stream.write(reinterpret_cast<const char*>(&mass),sizeof(float));

		const short evapCharge = i->second.evapCharge();
		stream.write(reinterpret_cast<const char*>(&evapCharge),sizeof(float));

		const float evapField = i->second.evapField();
		stream.write(reinterpret_cast<const char*>(&evapField),sizeof(float));

		const float evapEnergy = i->second.evapEnergy();
		stream.write(reinterpret_cast<const char*>(&evapEnergy),sizeof(float));
	}
}

const Configuration::NodeData& Configuration::Table::operator[](const Configuration::NodeId& id) const
{
	const std::map<NodeId,NodeData>::const_iterator i = _data.find(id);

	if (_data.end() != i)
		return i->second;
	else
		throw std::runtime_error("Configuration::Table::operator[]()");
}

Configuration::NodeData& Configuration::Table::operator[](const Configuration::NodeId& id)
{
	const std::map<NodeId,NodeData>::iterator i = _data.find(id);

	if (_data.end() != i)
		return i->second;
	else
		throw std::runtime_error("Configuration::Table::operator[]()");
}
	
const Configuration::NodeData& Configuration::Table::operator[](const std::string& name) const
{
	const std::map<std::string,NodeId>::const_iterator i = _nameIndex.find(name);

	if (_nameIndex.end() != i)
		return this->operator[](i->second);
	else
		throw std::runtime_error("Configuration::Table::operator[]()");
}

Configuration::NodeData& Configuration::Table::operator[](const std::string& name)
{
	const std::map<std::string,NodeId>::const_iterator i = _nameIndex.find(name);

	if (_nameIndex.end() != i)
		return _data[i->second];
	else
		throw std::runtime_error("Configuration::Table::operator[]()");
}

Configuration::NodeId Configuration::Table::insert(const NodeData& data)
{
	if (exists(data.name())) return NodeId::invalid();

	const NodeId id = getId();

	_data[id] = data;
	_data[id].setId(id);

	_nameIndex[data.name()] = id;

	return id;
}

Configuration::NodeId Configuration::Table::insert(const NodeId id, const NodeData& data)
{
	if (!id.isValid() || exists(id) || exists(data.name())) return NodeId::invalid();

	_data[id] = data;
	_data[id].setId(id);

	_nameIndex[data.name()] = id;

	return id;
}

void Configuration::Table::remove(const NodeId& id)
{
	_nameIndex.erase(_data[id].name());
	_data.erase(id);
}

void Configuration::Table::remove(const std::string& name)
{
	_data.erase(this->operator[](name).id());
	_nameIndex.erase(name);
}

std::list<Configuration::NodeId> Configuration::Table::ids() const
{
	std::list<NodeId> keys;
	for (std::map<NodeId,NodeData>::const_iterator i = _data.begin(); i != _data.end(); i++)
		keys.push_back(i->first);

	return keys;
}

std::list<std::string> Configuration::Table::names() const
{
	std::list<std::string> keys;
	for (std::map<std::string,NodeId>::const_iterator i = _nameIndex.begin(); i != _nameIndex.end(); i++)
		keys.push_back(i->first);

	return keys;
}

void Configuration::Table::clear()
{
	_data.clear();
	_nameIndex.clear();
}

Configuration::NodeId Configuration::Table::getId() const
{
	NodeId id = NodeId::min();

	while (exists(id))
	{
		if (NodeId::max() == id) throw std::runtime_error("Configuration::Table::getId()");
		id++;
	}

	return id;
}
