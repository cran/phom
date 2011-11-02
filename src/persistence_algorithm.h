//============================================================================
// Name        : cph
// Author      : Andrew Tausz <atausz@stanford.edu>
// Version     : 1.0
// Copyright   : Copyright © 2011 Andrew Tausz
// Description : A basic package for persistent homology in C++
//============================================================================

#ifndef PERSISTENCE_ALGORITHM_H_
#define PERSISTENCE_ALGORITHM_H_

#include "simplex.h"
#include "simplex_stream.h"
#include "barcode_collection.h"

#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <list>
#include <functional>

namespace cph
{

	template<class B, class T>
	class persistence_algorithm
	{
	private:
		typename std::size_t _max_dimension;

		std::set<B> * _marked_simplices;
		std::map<B, std::set<B> > * _T;

	public:
		persistence_algorithm(const typename std::size_t max_dimension = 2)
			: _max_dimension(max_dimension)
		{
		}
		virtual ~persistence_algorithm()
		{
		}

		barcode_collection<T> compute_intervals(const simplex_stream<B, T> & stream)
		{
			const typename std::vector<B>::const_iterator & begin = stream.begin();
			const typename std::vector<B>::const_iterator & end = stream.end();

			this->_T = new std::map<B, std::set<B> >();
			this->_marked_simplices = new std::set<B>();

			barcode_collection<T> intervals;

			for (typename std::vector<B>::const_iterator iter = begin; iter != end; iter++)
			{
   			    B simplex = (*iter);

				if (simplex.dimension() > this->_max_dimension + 1)
				{
					continue;
				}

				std::set<B> d = this->remove_pivot_rows(simplex, stream);

				if (d.empty())
				{
					this->_marked_simplices->insert(simplex);
				}
				else
				{
					B & sigma_j = simplex;
					const B sigma_i = this->get_maximum_object(d, stream);
					typename std::size_t k = sigma_i.dimension();
					this->_T->operator[](sigma_i) = d;

					T t_i = stream.get_filtration_value(sigma_i);
					T t_j = stream.get_filtration_value(sigma_j);

					if ((t_j - t_i > 0) && (k <= this->_max_dimension))
					{
						intervals.add_interval(k, t_i, t_j);
					}
				}
			}

			for (typename std::set<B>::const_iterator iter = this->_marked_simplices->begin(); iter != this->_marked_simplices->end(); iter++)
			{
				if (this->_T->find(*iter) == this->_T->end() || this->_T->at(*iter).empty())
				{
					typename std::size_t k = (*iter).dimension();
					if (k <= this->_max_dimension)
					{
						T t = stream.get_filtration_value(*iter);
						intervals.add_interval(k, t);
					}
				}
			}

			delete (this->_T);
			delete (this->_marked_simplices);

			return intervals;
		}

	private:
		std::set<B> remove_pivot_rows(const B & simplex, const simplex_stream<B, T> & stream)
		{
			typedef std::set<B> type;
			std::set<B> d = simplex.boundary_set();

			// removed non-marked simplices from d
			std::set<B> to_be_removed;
			for (typename std::set<B>::iterator iter = d.begin(); iter != d.end(); iter++)
			{
				if (this->_marked_simplices->find(*iter) == this->_marked_simplices->end())
				{
					to_be_removed.insert(*iter);
				}
			}

			for (typename std::set<B>::iterator iter = to_be_removed.begin(); iter != to_be_removed.end(); iter++)
			{
				d.erase(*iter);
			}

			while (!d.empty())
			{
				const B & sigma_i = this->get_maximum_object(d, stream);

				if (this->_T->find(sigma_i) == this->_T->end())
				{
					break;
				}

				if (this->_T->find(sigma_i)->second.empty())
				{
					break;
				}

				if (this->get_coefficient(this->_T->at(sigma_i), sigma_i) == false)
				{
					break;
				}

				this->accumulate(d, this->_T->at(sigma_i));
			}

			return d;
		}

		inline const B & get_maximum_object(const std::set<B> & chain, const simplex_stream<B, T> & stream) const
		{
			return *std::max_element(chain.begin(), chain.end(), stream.get_filtered_comparator());
		}

		inline bool get_coefficient(const std::set<B> & chain, const B & object) const
		{
			return (chain.find(object) != chain.end());
		}

		inline void accumulate(std::set<B> & a, const B & object) const
		{
			if (a.find(object) == a.end())
			{
				a.insert(object);
			}
			else
			{
				a.erase(object);
			}
		}

		void accumulate(std::set<B> & a, const std::set<B> & b) const
		{
			for (typename std::set<B>::const_iterator iter = b.begin(); iter != b.end(); iter++)
			{
				this->accumulate(a, *iter);
			}
		}

		std::set<B> sum(const std::set<B> & a, const std::set<B> & b) const
		{
			std::set<B> result;

			std::set_symmetric_difference(a.begin(), a.end(), b.begin(), b.end(), result.begin());

			return result;
		}
	};

}

#endif /* PERSISTENCE_ALGORITHM_H_ */
