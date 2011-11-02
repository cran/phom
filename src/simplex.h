//============================================================================
// Name        : cph
// Author      : Andrew Tausz <atausz@stanford.edu>
// Version     : 1.0
// Copyright   : Copyright © 2011 Andrew Tausz
// Description : A basic package for persistent homology in C++
//============================================================================

#ifndef SIMPLEX_H_
#define SIMPLEX_H_

#include <vector>
#include <algorithm>
#include <iostream>
#include <set>
#include <list>

namespace cph
{

	template<class V = std::size_t>
	class simplex
	{

	private:
		typedef simplex<V> self;
		std::vector<V> _vertices;

	public:

		static const simplex<V> make_simplex(const V & vertex)
		{
			simplex<V> s(vertex);
			return s;
		}

		static const simplex<V> make_simplex(const V & v1, const V & v2)
		{
			simplex<V> s(v1, v2);
			return s;
		}

		virtual ~simplex()
		{
		}

		simplex(const V & vertex)
		{
			_vertices.push_back(vertex);
		}

		simplex(const V & v1, const V & v2)
		{
			if (v1 < v2)
			{
				_vertices.push_back(v1);
				_vertices.push_back(v2);
			}
			else
			{
				_vertices.push_back(v2);
				_vertices.push_back(v1);
			}
		}

		simplex(const int vertex_array[], int length)
		{
			for (int i = 0; i < length; i++)
			{
				this->_vertices.push_back(vertex_array[i]);
			}
			std::sort(_vertices.begin(), _vertices.end());
		}

		simplex(std::vector<V> & vertices) :
			_vertices(vertices)
		{
			std::sort(vertices.begin(), vertices.end());
		}

		simplex(const self & other) :
			_vertices(other._vertices)
		{
		}

		simplex & operator =(const self & other)
		{
			this->_vertices = other._vertices;
			return *this;
		}

		const V & operator [](const typename std::vector<V>::size_type & index) const
		{
			return _vertices[index];
		}

		typename std::size_t dimension() const
		{
			return _vertices.size() - 1;
		}

		simplex face(const typename std::size_t k) const
		{
			std::vector<V> result;

			for (typename std::size_t i = 0; i < _vertices.size(); i++)
			{
				if ((i < k) || (i > k))
				{
					result.push_back(_vertices[i]);
				}
			}
			return simplex(result);
		}

		std::vector<simplex> boundary() const
		{
			std::vector<simplex> result;
			typename std::size_t n = _vertices.size();

			for (typename std::size_t i = 0; i < n; i++)
			{
				result.push_back(this->face(i));
			}

			return result;
		}

		std::set<simplex> boundary_set() const
		{
			std::set<simplex> result;
			typename std::size_t n = _vertices.size();

			for (typename std::size_t i = 0; i < n; i++)
			{
				result.insert(this->face(i));
			}

			return result;
		}

		std::list<simplex> boundary_list() const
		{
			std::list<simplex> result;
			typename std::size_t n = _vertices.size();

			for (typename std::size_t i = 0; i < n; i++)
			{
				result.push_back(this->face(i));
			}

			return result;
		}

		simplex append_to(const V & v) const
		{
			simplex<V> s(*this);
			s._vertices.push_back(v);
			std::sort(s._vertices.begin(), s._vertices.end());
			return s;
		}

		int compare(const simplex & other) const
		{
			if (this->dimension() > other.dimension())
			{
				return 1;
			}
			else if (this->dimension() < other.dimension())
			{
				return -1;
			}

			typename std::size_t n = _vertices.size();

			for (typename std::size_t i = 0; i < n; i++)
			{
				if (_vertices[i] > other._vertices[i])
				{
					return 1;
				}
				else if (_vertices[i] < other._vertices[i])
				{
					return -1;
				}
			}
			return 0;
		}

		bool operator <(const simplex & other) const
		{
			return (this->compare(other) < 0);
		}

		bool operator >(const simplex & other) const
		{
			return (this->compare(other) > 0);
		}

		bool operator <=(const simplex & other) const
		{
			return (this->compare(other) <= 0);
		}

		bool operator >=(const simplex & other) const
		{
			return (this->compare(other) >= 0);
		}

		bool operator ==(const simplex & other) const
		{
			return (this->compare(other) == 0);
		}

		friend std::ostream & operator <<(std::ostream & s, const simplex<V> & value)
		{
			s << '[';

			typename std::size_t n = value._vertices.size();

			for (typename std::size_t i = 0; i < n; i++)
			{
				if (i > 0)
				{
					s << ',';
				}
				s << value._vertices[i];
			}

			s << ']';

			return s;
		}

	};

}

#endif /* SIMPLEX_H_ */
