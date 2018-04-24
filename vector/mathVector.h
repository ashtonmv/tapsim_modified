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

#ifndef LIBVECTOR_MATH_VECTOR_H
#define LIBVECTOR_MATH_VECTOR_H

#include <ostream>

#include <cmath>

#include "staticVector.h"

template <class elementType, unsigned short _size>
class MathVector : public StaticVector<elementType,_size>
{
	public:
		MathVector(const elementType& value =elementType())
			: StaticVector<elementType,_size>(value)
		{}

		MathVector(const elementType* values)
			: StaticVector<elementType,_size>(values)
		{}

		// ***

		MathVector& operator*=(const elementType& scalar)
		{
			for (unsigned short i = 0; i < _size; i++)
				this->operator[](i) *= scalar;

			return *this;
		}

		MathVector& operator/=(const elementType& scalar)
		{
			for (unsigned short i = 0; i < _size; i++)
				this->operator[](i) /= scalar;

			return *this;
		}

		MathVector& operator+=(const MathVector& value)
		{
			for (unsigned short i = 0; i < _size; i++)
				this->operator[](i) += value[i];

			return *this;
		}

		MathVector& operator-=(const MathVector& value)
		{
			for (unsigned short i = 0; i < _size; i++)
				this->operator[](i) -= value[i];

			return *this;
		}

		// ***

		MathVector norm() const
		{
			return MathVector(*this).normalize();
		}

		MathVector& normalize()
		{
			*this /= length();

			return *this;
		}

		// ***

		elementType length() const
		{
			elementType length = (*this) * (*this);
			length = std::sqrt(length);

			return length;
		}
};

template<class elementType, unsigned short _size>
elementType angle(const MathVector<elementType, _size>& a, const MathVector<elementType, _size>& b)
{
	elementType result(a * b);

	result /= a.length();
	result /= b.length();

	return result;
}

template <class elementType, unsigned short _size>
elementType operator*(const MathVector<elementType, _size>& a, const MathVector<elementType, _size>& b)
{
	elementType value(0);

	for (unsigned short i = 0; i < _size; i++)
		value += a[i] * b[i];

	return value;
}

template <class elementType, unsigned short _size>
MathVector<elementType,_size> operator*(const MathVector<elementType,_size>& vec, const elementType factor)
{
	MathVector<elementType,_size> value(vec);

	for (unsigned short i = 0; i < _size; i++)
		value[i] *= factor;

	return value;
}

template <class elementType, unsigned short _size>
MathVector<elementType,_size> operator*(const elementType factor, const MathVector<elementType,_size>& vec)
{
	return operator*(vec,factor);
}

template <class elementType, unsigned short _size>
MathVector<elementType,_size> operator/(const MathVector<elementType,_size>& vec, const elementType factor)
{
	MathVector<elementType,_size> value(vec);

	for (unsigned short i = 0; i < _size; i++)
		value[i] /= factor;

	return value;
}

template <class elementType, unsigned short _size>
MathVector<elementType, _size> operator+(const MathVector<elementType, _size>& a, const MathVector<elementType, _size>& b)
{
	MathVector<elementType, _size> value(a);

	for (unsigned short i = 0; i < _size; i++)
		value[i] += b[i];

	return value;
}

template <class elementType, unsigned short _size>
MathVector<elementType, _size> operator-(const MathVector<elementType, _size>& a, const MathVector<elementType, _size>& b)
{
	MathVector<elementType, _size> value(a);

	for (unsigned short i = 0; i < _size; i++)
		value[i] -= b[i];

	return value;
}

template <class elementType, unsigned short _size>
bool operator==(const MathVector<elementType, _size>& a, const MathVector<elementType, _size>& b)
{
	for (unsigned short i = 0; i < _size; i++)
	{
		if (a[i] != b[i]) return false;
	}

	return true;
}

template <class elementType, unsigned short _size>
bool operator!=(const MathVector<elementType, _size>& a, const MathVector<elementType, _size>& b)
{
	return !(a == b);
}

// ***

template <class elementType, unsigned short _size>
inline std::ostream& operator<<(std::ostream& stream, const MathVector<elementType, _size>& obj)
{
	stream << "(";
	if (_size > 0) stream << obj[0];
	for (unsigned short i = 1; i < _size; stream << ", " << obj[i++]);
	stream << ")";

	return stream;
}

#endif
