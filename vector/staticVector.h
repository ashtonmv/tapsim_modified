/************************************************************************
*                                                                       *
* Vector class template library                                         *
*                                                                       *
* Copyright (C) 2011 Christian Oberdorfer                               *
*                                                                       *
* This library is free software: you can redistribute it and/or modify  *
* it under the terms of the GNU General Public License as published by  *
* the Free Software Foundation, either version 3 of the License, or any *
* any later version.                                                    *
*                                                                       *
* This library is distributed in the hope that it will be useful,       *
* but WITHOUT ANY WARRANTY; without even the implied warranty of        *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
* GNU General Public License for more details.                          *
*                                                                       *
* You should have received a copy of the GNU General Public License     *
* along with this library.  If not, see 'http://www.gnu.org/licenses'   *
*                                                                       *
************************************************************************/ 

#ifndef STATIC_VECTOR_H
#define STATIC_VECTOR_H

#include <stdexcept>

template <class elementType, unsigned short _size>
class StaticVector
{
	public:
		StaticVector(const elementType& value =elementType())
		{
			for (unsigned short i = 0; i < _size; i++)
				_data[i] = value;
		}

		StaticVector(const elementType* values)
		{
			for (unsigned short i = 0; i < _size; i++)
				_data[i] = values[i];
		}

		// ***

		const elementType& operator[](unsigned short index) const { return _data[index]; }
		elementType& operator[](unsigned short index) { return _data[index]; }

		const elementType& at(unsigned short index) const
		{
			if (index < _size) 
				return _data[index];
			else
				throw std::out_of_range("StaticVector::at()");
		}

		elementType& at(unsigned short index)
		{
			if (index < _size) 
				return _data[index];
			else
				throw std::out_of_range("StaticVector::at()");
		}

		// ***

		static unsigned short size() { return _size; }

		const elementType* raw() const { return _data; }
		elementType* raw() { return _data; }

	private:
		elementType _data[_size];
};

#endif

