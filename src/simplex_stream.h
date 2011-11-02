//============================================================================
// Name        : cph
// Author      : Andrew Tausz <atausz@stanford.edu>
// Version     : 1.0
// Copyright   : Copyright © 2011 Andrew Tausz
// Description : A basic package for persistent homology in C++
//============================================================================

#ifndef SIMPLEX_STREAM_H_
#define SIMPLEX_STREAM_H_

#include <vector>
#include <utility>
#include <functional>
#include <iostream>
#include <algorithm>
#include <map>

#include "utility.h"

namespace cph
{

	template<class X, class Y>
	struct ordered_comparison: public std::binary_function<std::pair<X, Y>, std::pair<X, Y>, bool>
	{
		inline bool operator()(const std::pair<X, Y> & p1, const std::pair<X, Y> & p2)
		{
			if (p1.first < p2.first)
			{
				return true;
			}

			if (p1.first > p2.first)
			{
				return false;
			}

			return (p1.second < p2.second);
		}
	};

	template<class B, class T>
	class auxilary_comparison
	{
		const std::map<B, T> & _map;

	public:
		auxilary_comparison(const std::map<B, T> & map) :
			_map(map)
		{
		}

		inline bool operator()(const B & p1, const B & p2)
		{
			T t1 = this->_map.at(p1);
			T t2 = this->_map.at(p2);

			if (t1 < t2)
			{
				return true;
			}

			if (t1 > t2)
			{
				return false;
			}

			return (p1 < p2);
		}
	};

	template<class X, class Y, class S>
	S & operator <<(S & s, const std::pair<X, Y> & p)
	{
		s << '[' << p.first << ',' << p.second << ']';
		return s;
	}

	template<class B, class T>
	class simplex_stream
	{

	protected:
		std::map<B, T> _filtration_values;
		std::vector<B> _simplices;

	public:
		simplex_stream()
		{
		}
		virtual ~simplex_stream()
		{
		}

		const std::size_t size() const
		{
			return this->_simplices.size();
		}

		const T get_filtration_value(const B & simplex) const
		{
			if (this->_filtration_values.find(simplex) == this->_filtration_values.end())
			{
				return T(0);
			}

			return this->_filtration_values.at(simplex);
		}

		typename std::vector<B>::const_iterator begin() const
		{
			return this->_simplices.begin();
		}

		typename std::vector<B>::const_iterator end() const
		{
			return this->_simplices.end();
		}

		virtual void add_simplex(const B & simplex, const T & filtration_value)
		{
			this->_filtration_values[simplex] = filtration_value;
			this->_simplices.push_back(simplex);
		}

		void print_contents() const
		{

		}

		void ensure_sorted()
		{
			std::sort(this->_simplices.begin(), this->_simplices.end(), this->get_filtered_comparator());
		}

		auxilary_comparison<B, T> get_filtered_comparator() const
		{
			return auxilary_comparison<B, T> (this->_filtration_values);
		}
	};

}

#endif /* SIMPLEX_STREAM_H_ */
