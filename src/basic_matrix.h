//============================================================================
// Name        : cph
// Author      : Andrew Tausz <atausz@stanford.edu>
// Version     : 1.0
// Copyright   : Copyright © 2011 Andrew Tausz
// Description : A basic package for persistent homology in C++
//============================================================================

#ifndef BASIC_MATRIX_HPP_
#define BASIC_MATRIX_HPP_

#include <math.h>
#include <complex>

#include "metrics.h"

namespace cph
{

	template<class T>
	class basic_matrix
	{
	private:
		T ** const _data;
		const std::size_t _rows, _columns;
	public:
		basic_matrix(const std::size_t rows, const std::size_t columns) :
			_data(new T*[rows]), _rows(rows), _columns(columns)
		{
			for (std::size_t i = 0; i < _rows; ++i)
			{
				this->_data[i] = new T[_columns];
			}
		}

		basic_matrix(const std::size_t rows, const std::size_t columns, const T initial_value) :
			_data(new T*[rows]), _rows(rows), _columns(columns)
		{
			for (std::size_t i = 0; i < _rows; ++i)
			{
				this->_data[i] = new T[_columns];
				for (std::size_t j = 0; j < _columns; j++)
				{
					this->_data[i][j] = initial_value;
				}
			}
		}

		basic_matrix(const basic_matrix<T> & other) :
			_data(new T*[other._rows]), _rows(other._rows), _columns(other._columns)
		{
			for (std::size_t i = 0; i < _rows; ++i)
			{
				this->_data[i] = new T[_columns];
				for (std::size_t j = 0; j < _columns; j++)
				{
					this->_data[i][j] = other._data[i][j];
				}
			}
		}

		virtual ~basic_matrix()
		{
			for (std::size_t i = 0; i < this->_rows; ++i)
			{
				delete[] (this->_data[i]);
			}

			delete[] (this->_data);
		}

		const std::size_t rows() const
		{
			return this->_rows;
		}

		const std::size_t columns() const
		{
			return this->_columns;
		}

		T & operator()(const std::size_t row, const std::size_t column)
		{
			return this->_data[row][column];
		}

		const T operator()(const std::size_t row, const std::size_t column) const
		{
			return this->_data[row][column];
		}

		const T row_norm(const std::size_t i) const
		{
			T sum(0);
			for (std::size_t k = 0; k < this->_columns; k++)
			{
				sum += _data[i][k] * _data[i][k];
			}

			return std::sqrt(sum);
		}

		const T row_distance(const std::size_t i, const std::size_t j, const metric metric_type, const T p = 2) const
		{
			switch (metric_type)
			{
				case euclidean:
					return row_euclidean_distance(i, j);
				case maximum:
					return row_maximum_distance(i, j);
				case manhattan:
					return row_manhattan_distance(i, j);
				case canberra:
					return row_canberra_distance(i, j);
				case binary:
					return row_binary_distance(i, j);
				case minkowski:
					return row_minkowski_distance(i, j, p);
			}

			return 0;
		}

		const T row_euclidean_distance(const std::size_t i, const std::size_t j) const
		{
			T sum(0), distance(0);
			for (std::size_t k = 0; k < this->_columns; k++)
			{
				distance = _data[i][k] - _data[j][k];
				sum += distance * distance;
			}

			return std::sqrt(sum);
		}

		const T row_maximum_distance(const std::size_t i, const std::size_t j) const
		{
			T max(0), difference(0);
			for (std::size_t k = 0; k < this->_columns; k++)
			{
				difference = _data[i][k] - _data[j][k];

				if (difference < 0)
				{
					difference *= -1;
				}

				if (difference > max)
				{
					max = difference;
				}
			}

			return max;
		}

		const T row_manhattan_distance(const std::size_t i, const std::size_t j) const
		{
			T sum(0), difference(0);
			for (std::size_t k = 0; k < this->_columns; k++)
			{
				difference = _data[i][k] - _data[j][k];

				if (difference < 0)
				{
					difference *= -1;
				}

				sum += difference;
			}

			return sum;
		}

		const T row_minkowski_distance(const std::size_t i, const std::size_t j, const T p) const
		{
			T sum(0), difference(0);
			for (std::size_t k = 0; k < this->_columns; k++)
			{
				difference = _data[i][k] - _data[j][k];

				if (difference < 0)
				{
					difference *= -1;
				}

				sum += pow(difference, p);
			}

			return pow(difference, 1.0 / p);
		}

		const T row_canberra_distance(const std::size_t i, const std::size_t j) const
		{
			T sum(0), difference(0), row_sum(0);
			for (std::size_t k = 0; k < this->_columns; k++)
			{
				difference = _data[i][k] - _data[j][k];
				sum = _data[i][k] + _data[j][k];

				if (difference < 0)
				{
					difference *= -1;
				}

				if (sum < 0)
				{
					sum *= -1;
				}

				if (sum != 0)
				{
					row_sum += difference / sum;
				}
			}

			return row_sum;
		}

		const T row_binary_distance(const std::size_t i, const std::size_t j) const
		{
			T numerator(0), denominator(0);
			for (std::size_t k = 0; k < this->_columns; k++)
			{
				bool x = (_data[i][k] != 0);
				bool y = (_data[j][k] != 0);

				if (x || y)
				{
					denominator++;

					if (x ^ y)
					{
						numerator++;
					}
				}
			}

			if (denominator > 0)
			{
				return numerator / denominator;
			}
			else
			{
				return 0;
			}
		}

		template<class S>
		friend S & operator <<(S & s, const basic_matrix<T> & m)
		{

			s << '[';
			for (std::size_t i = 0; i < m._rows; i++)
			{
				if (i > 0)
				{
					s << ';';
					s << std::endl;
				}
				for (std::size_t j = 0; j < m._columns; j++)
				{
					if (j > 0)
					{
						s << ',';
					}
					s << m(i, j);
				}

			}
			s << ']';

			return s;
		}
	};

}

#endif /* BASIC_MATRIX_HPP_ */
