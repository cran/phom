//============================================================================
// Name        : cph
// Author      : Andrew Tausz <atausz@stanford.edu>
// Version     : 1.0
// Copyright   : Copyright © 2011 Andrew Tausz
// Description : A basic package for persistent homology in C++
//============================================================================

#ifndef BARCODE_COLLECTION_HPP_
#define BARCODE_COLLECTION_HPP_

#include <map>
#include <list>

#include "right_open_interval.h"
#include "basic_matrix.h"

namespace cph
{

	template<class T>
	class barcode_collection
	{
	private:
		std::map<std::size_t, std::list<right_open_interval<T> > *> _intervals;
		static const std::size_t ZERO_REPLACEMENT = 9999;
		std::size_t _num_intervals;
	public:
		barcode_collection() :
			_num_intervals(0)
		{
		}

		virtual ~barcode_collection()
		{
			for (typename std::map<std::size_t, std::list<right_open_interval<T> > *>::iterator iter = this->_intervals.begin(); iter != this->_intervals.end(); iter++)
			{
				delete (*iter).second;
			}
		}

		const std::size_t num_intervals() const
		{
			return this->_num_intervals;
		}

		void add_interval(const std::size_t dimension, const T start)
		{
			std::size_t filtered_dimension = barcode_collection::shift(dimension);
			if (_intervals.find(filtered_dimension) == _intervals.end())
			{
				_intervals[filtered_dimension] = new std::list<right_open_interval<T> >();
			}

			_intervals[filtered_dimension]->push_back(right_open_interval<T>::make_interval(start));

			_num_intervals++;
		}

		void add_interval(const std::size_t dimension, const T start, const T finish)
		{
			std::size_t filtered_dimension = barcode_collection::shift(dimension);
			if (_intervals.find(filtered_dimension) == _intervals.end())
			{
				_intervals[filtered_dimension] = new std::list<right_open_interval<T> >();
			}

			_intervals[filtered_dimension]->push_back(right_open_interval<T>::make_interval(start, finish));

			_num_intervals++;
		}

		basic_matrix<T> get_endpoint_matrix(const T max_filtration_value) const
		{
			basic_matrix<T> endpoints(_num_intervals, 3);

			std::size_t interval_index(0);

			for (typename std::map<std::size_t, std::list<right_open_interval<T> > *>::const_iterator iter = this->_intervals.begin(); iter
					!= this->_intervals.end(); iter++)
			{
				std::list<right_open_interval<T> > * interval_set = iter->second;
				std::size_t filtered_dimension = iter->first;
				std::size_t dimension = barcode_collection::unshift(filtered_dimension);
				for (typename std::list<right_open_interval<T> >::const_iterator interval_iter = interval_set->begin(); interval_iter != interval_set->end(); interval_iter++)
				{
					endpoints.operator()(interval_index, 0) = dimension;
					endpoints.operator()(interval_index, 1) = interval_iter->start();
					if (interval_iter->is_infinite())
					{
						endpoints.operator()(interval_index, 2) = max_filtration_value;
					}
					else
					{
						endpoints.operator()(interval_index, 2) = interval_iter->finish();
					}

					interval_index++;
				}
			}

			return endpoints;
		}

		std::pair<std::vector<T>, std::vector<T> > get_startpoints(const std::size_t dimension, const bool include_infinite_intervals = false)
		{
			std::vector<T> start_points;
			std::vector<T> end_points;
			for (typename std::map<std::size_t, std::list<right_open_interval<T> > *>::const_iterator iter = this->_intervals.begin(); iter
					!= this->_intervals.end(); iter++)
			{
				std::list<right_open_interval<T> > * interval_set = iter->second;
				std::size_t filtered_dimension = iter->first;
				std::size_t iter_dimension = barcode_collection::unshift(filtered_dimension);
				if (dimension == iter_dimension)
				{
					for (typename std::list<right_open_interval<T> >::const_iterator interval_iter = interval_set->begin(); interval_iter != interval_set->end(); interval_iter++)
					{
						if (include_infinite_intervals || !interval_iter->is_infinite())
						{
							start_points.push_back(interval_iter->start());
							end_points.push_back(interval_iter->finish());
						}
					}
				}
			}

			return std::make_pair<std::vector<T>, std::vector<T> >(start_points, end_points);
		}

		friend std::ostream & operator <<(std::ostream & s, const barcode_collection<T> & collection)
		{
			for (typename std::map<std::size_t, std::list<right_open_interval<T> > *>::const_iterator iter = collection._intervals.begin(); iter
					!= collection._intervals.end(); iter++)
			{
				std::list<right_open_interval<T> > * interval_set = iter->second;
				std::size_t filtered_dimension = iter->first;
				std::size_t dimension = barcode_collection::unshift(filtered_dimension);
				for (typename std::list<right_open_interval<T> >::const_iterator interval_iter = interval_set->begin(); interval_iter != interval_set->end(); interval_iter++)
				{
					s << dimension << ": " << (*interval_iter) << std::endl;
				}
			}

			return s;
		}

	private:
		static std::size_t shift(const std::size_t dimension)
		{
			if (dimension == 0)
			{
				return barcode_collection::ZERO_REPLACEMENT;
			}
			return dimension;
		}

		static std::size_t unshift(const std::size_t dimension)
		{
			if (dimension == barcode_collection::ZERO_REPLACEMENT)
			{
				return 0;
			}
			return dimension;
		}
	};

}

#endif /* BARCODE_COLLECTION_HPP_ */
