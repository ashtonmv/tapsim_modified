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

#ifndef LIBVECTOR_STATIC_VECTOR_TUPLE_H
#define LIBVECTOR_STATIC_VECTOR_TUPLE_H

#include "staticVector.h"

template<class elementType, unsigned short _size>
class Tuple : public StaticVector<elementType,_size>
{
	public:
		Tuple() : StaticVector<elementType,_size>() {}
		Tuple(const elementType& value) : StaticVector<elementType,_size>(value) {}
		Tuple(const elementType* values) : StaticVector<elementType,_size>(values) {}
	
		// ***
	
		Tuple& operator>>=(const unsigned short value)
		{
			for (unsigned short i = 0; i < _size; i++)
			{
				const elementType tmp = this->operator[](i);
				this->operator[](i) = this->operator[]((i+value) % _size);
				this->operator[]((i+value) % _size) = tmp;
			}
		}
	
		Tuple& operator<<=(const unsigned short value)
		{
			const unsigned short range = value % _size;
			return this->operator>>=(_size - range);
		}
};

// ***

template<class elementType, unsigned short _size>
Tuple<elementType,_size> operator>>(const Tuple<elementType,_size>& value, const unsigned short skip)
{
	Tuple<elementType,_size> result;
	for (unsigned short i = 0; i < _size; i++)
		result[i] = value[(i+skip) % _size];

	return result;
}

template<class elementType, unsigned short _size>
Tuple<elementType,_size> operator<<(const Tuple<elementType,_size>& value, const unsigned short skip)
{
	return value >> (_size - skip % _size);
}

// *** 2-TUPLE => PAIR

template<class elementType>
class Pair : public StaticVector<elementType,2>
{
	public:
		enum tupleIDs { A = 0, B = 1 };
	
		Pair() : StaticVector<elementType,2>() {}
		Pair(const elementType& value) : StaticVector<elementType,2>(value) {}
		Pair(const elementType* values) : StaticVector<elementType,2>(values) {}
	
		Pair(const elementType& a, const elementType& b)
		{
			this->operator[](0) = a;
			this->operator[](1) = b;
		}
	
		const elementType& a() const { return this->operator[](0); }
		elementType& a() { return  this->operator[](0); }
	
		const elementType& b() const { return this->operator[](1); }
		elementType& b() { return this->operator[](1); }
};

// *** 3-TUPLE => TRIPLET

template<class elementType>
class Triplet : public StaticVector<elementType,3>
{
	public:
		enum tupleIDs { A = 0, B = 1, C = 2 };
	
		Triplet() : StaticVector<elementType,3>() {}
		Triplet(const elementType& value) : StaticVector<elementType,3>(value) {}
		Triplet(const elementType* values) : StaticVector<elementType,3>(values) {}
	
		Triplet(const elementType& a, const elementType& b, const elementType& c)
		{
			this->operator[](0) = a;
			this->operator[](1) = b;
			this->operator[](2) = c;
		}
	
		const elementType& a() const { return this->operator[](0); }
		elementType& a() { return  this->operator[](0); }
	
		const elementType& b() const { return this->operator[](1); }
		elementType& b() { return this->operator[](1); }
	
		const elementType& c() const { return this->operator[](2); }
		elementType& c() { return this->operator[](2); }
};

// *** 4-TUPLE => QUADRUPLE

template<class elementType>
class Quadruple : public StaticVector<elementType,4>
{
	public:
		enum tupleIDs { A = 0, B = 1, C = 2, D = 3 };
	
		Quadruple() : StaticVector<elementType,4>() {}
		Quadruple(const elementType& value) : StaticVector<elementType,4>(value) {}
		Quadruple(const elementType* values) : StaticVector<elementType,4>(values) {}
	
		Quadruple(const elementType& a, const elementType& b, const elementType& c, const elementType& d)
		{
			this->operator[](0) = a;
			this->operator[](1) = b;
			this->operator[](2) = c;
			this->operator[](3) = d;
		}
	
		const elementType& a() const { return this->operator[](0); }
		elementType& a() { return  this->operator[](0); }
	
		const elementType& b() const { return this->operator[](1); }
		elementType& b() { return this->operator[](1); }
	
		const elementType& c() const { return this->operator[](2); }
		elementType& c() { return this->operator[](2); }
	
		const elementType& d() const { return this->operator[](3); }
		elementType& d() { return this->operator[](3); }
};

#endif
