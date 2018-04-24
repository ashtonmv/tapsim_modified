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

#ifndef LIBVECTOR_MATH_VECTOR_2D_H
#define LIBVECTOR_MATH_VECTOR_2D_H

#include "mathVector.h"

template<class elementType>
class MathVector2d : public MathVector<elementType,2>
{
	public:
		enum dimIDs { X = 0, Y = 1 };

		MathVector2d()
			: MathVector<elementType,2>()
		{}

		MathVector2d(const MathVector<elementType,2>& value)
			: MathVector<elementType,2>(value)
		{}

		MathVector2d(const elementType* values)
			: MathVector<elementType,2>(values)
		{}

		MathVector2d(const elementType& x, const elementType& y)
		{
			this->operator[](0) = x;
			this->operator[](1) = y;
		}

		const elementType& x() const { return this->operator[](0); }
		elementType& x() { return  this->operator[](0); }

		const elementType& y() const { return  this->operator[](1); }
		elementType& y() { return this->operator[](1); }
};

#endif
