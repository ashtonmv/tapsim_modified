#ifndef TAPSIM_FILE_UTIL_H
#define TAPSIM_FILE_UTIL_H

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

#include <string>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <typeinfo>

namespace File_Utils
{
	template<class T>
	T read_element(const std::string& str)
	{
		T value;
		std::istringstream obj(str);
		obj >> value;

		if (obj.fail())
		{
			std::string errorMessage;
			errorMessage = "File_Utils::read_element(): value '";
			errorMessage += str;
			errorMessage += "' into '";
			errorMessage += typeid(T).name();
			errorMessage += "' - bad conversion!";

			throw std::runtime_error(errorMessage.c_str());
		}
	
		return value;
	}
	
	template<class T>
	inline void binRead(std::istream& stream, T* value)
	{
		stream.read(reinterpret_cast<char*>(value),sizeof(T));
	}

	template<>
	inline void binRead<std::string>(std::istream& stream, std::string* value)
	{
		value->clear();

		char c;
		stream.get(c);

		while (c != '\0')
		{
			value->push_back(c);
			stream.get(c);
		}
	}
	
	template<class T>
	inline void binWrite(std::ostream& stream, const T& value)
	{
		stream.write(reinterpret_cast<const char*>(&value),sizeof(T));
	}
	
	template<>
	inline void binWrite<std::string>(std::ostream& stream, const std::string& value)
	{
		const char* strValue = value.c_str();
		stream.write(reinterpret_cast<const char*>(strValue),value.size()+1);
	}
	
	std::string trim(const std::string&);
}

#endif
