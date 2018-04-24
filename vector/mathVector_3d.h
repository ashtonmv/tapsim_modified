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

#ifndef LIBVECTOR_MATH_VECTOR_3D_H
#define LIBVECTOR_MATH_VECTOR_3D_H

#include "mathVector.h"

template<class elementType>
class MathVector3d : public MathVector<elementType,3>
{
	public:
		enum dimIDs { X = 0, Y = 1, Z = 2 };

		MathVector3d()
			: MathVector<elementType,3>()
		{}

		MathVector3d(const MathVector<elementType,3>& value)
			: MathVector<elementType,3>(value)
		{}

		MathVector3d(const elementType* values)
			: MathVector<elementType,3>(values)
		{}

		MathVector3d(const elementType& x, const elementType& y, const elementType& z)
		{
			this->operator[](0) = x;
			this->operator[](1) = y;
			this->operator[](2) = z;
		}

		const elementType& x() const { return this->operator[](0); }
		elementType& x() { return  this->operator[](0); }

		const elementType& y() const { return  this->operator[](1); }
		elementType& y() { return this->operator[](1); }

		const elementType& z() const { return  this->operator[](2); }
		elementType& z() { return this->operator[](2); }
};

template<class elementType>
MathVector3d<elementType> vecProduct(const MathVector<elementType,3>& a, const MathVector<elementType,3>& b)
{
	MathVector3d<elementType> result;

	for (unsigned short i = 0; i < 3; i++)
	{
		const short k = (1+i) % 3;
		const short j = (2+i) % 3;

		result[i] = a[k] * b[j];
		result[i] -= a[j] * b[k];
	}

	return result;
}

template<class elementType>
elementType tripProduct(const MathVector<elementType,3>& a, const MathVector<elementType,3>& b, const MathVector<elementType,3>& c)
{
	return a * vecProduct(b,c);
}

#endif
