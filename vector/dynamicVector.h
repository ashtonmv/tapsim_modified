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

#ifndef LIBVECTOR_DYNAMIC_VECTOR_H
#define LIBVECTOR_DYNAMIC_VECTOR_H

#include <stdexcept>

template <class elementType>
class DynamicVector
{
	public:
		DynamicVector(const unsigned int size =0)
		{
			if (size > 0)
			{
				_size = size;
				_data = new elementType[_size];
			}
			else
			{
				_size = 0;
				_data = 0;
			}
		}

		DynamicVector(const unsigned int size, const elementType& value)
		{
			if (size > 0)
			{
				_size = size;
				_data = new elementType[_size];

				assign(value);
			}
			else
			{
				_size = 0;
				_data = 0;
			}
		}

		DynamicVector(const elementType* values, const unsigned int size)
		{
			if (_size > 0)
			{
				_size = size;
				_data = new elementType[_size];

				for (unsigned int i = 0; i < _size; i++)
					_data[i] = values[i];
			}
			else
			{
				_size = 0;
				_data = 0;
			}
		}

		DynamicVector(const DynamicVector& obj)
		{
			_size = obj.size();
			_data = new elementType[_size];

			for (unsigned int i = 0; i < _size; i++)
				_data[i] = obj[i];
		}

		~DynamicVector()
		{
			if (_data != 0) delete[] _data;
		}

		// ***

		DynamicVector& allocate(const unsigned int size)
		{
			if (_data != 0) delete[] _data;

			if (size > 0)
			{
				_size = size;
				_data = new elementType[_size];
			}
			else
			{
				_size = 0;
				_data = 0;
			}

			return *this;
		}

		DynamicVector& allocate(const unsigned int size, const elementType& value)
		{
			allocate(size);
			assign(value);

			return *this;
		}

		DynamicVector& reallocate(const unsigned int size, const elementType& value =elementType())
		{
			unsigned int oldSize = _size;
			elementType* oldData = _data;

			if (size > 0)
			{
				_size = size;
				_data = new elementType[_size];

				const unsigned int copyLimit = (oldSize < size ? oldSize : size);
				
				for (unsigned int i = 0; i < copyLimit; i++)
					_data[i] = oldData[i];

				for (unsigned int i = oldSize; i  <_size; i++)
					_data[i] = value;
			}
			else
			{
				_size = 0;
				_data = 0;
			}

			if (oldData != 0) delete[] oldData;

			return *this;
		}

		void release()
		{
			if (_data != 0) delete[] _data;

			_size = 0;
			_data = 0;
		}

		// ***

		DynamicVector& assign(const elementType& value =elementType())
		{
			for (unsigned int i = 0; i < _size; i++)
				_data[i] = value;

			return *this;
		}

		// ***

		unsigned int size() const { return _size; }

		const elementType& operator[](const unsigned int index) const { return _data[index]; }
		elementType& operator[](const unsigned int index) { return _data[index]; }

		const elementType& at(const unsigned int index) const
		{
			if (index < _size)
				return _data[index];
			else
				throw std::out_of_range("DynamicVector::at()");
		}

		elementType& at(const unsigned int index)
		{
			if (index < _size)
				return _data[index];
			else
				throw std::out_of_range("DynamicVector::at()");
		}

		const elementType* raw() const { return _data; }
		elementType* raw() { return _data; }

		bool isEmpty() const { return _size < 0; }

	private:
		unsigned int _size;
		elementType* _data;
};

#endif
