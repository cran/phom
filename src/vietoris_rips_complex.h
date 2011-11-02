//============================================================================
// Name        : cph
// Author      : Andrew Tausz <atausz@stanford.edu>
// Version     : 1.0
// Copyright   : Copyright © 2011 Andrew Tausz
// Description : A basic package for persistent homology in C++
//============================================================================

#ifndef VIETORIS_RIPS_COMPLEX_H_
#define VIETORIS_RIPS_COMPLEX_H_

#include "flag_complex.h"
#include "basic_graph.h"

namespace cph
{

	template<class T>
	class vietoris_rips_complex: public flag_complex<T>
	{

	private:

	public:
		vietoris_rips_complex(const finite_metric_space<T> & metric_space, const T & max_filtration_value, const int max_dimension) :
			flag_complex<T> (metric_space, max_filtration_value, max_dimension, metric_space.size())
		{
		}

		virtual ~vietoris_rips_complex()
		{
		}

		virtual basic_graph<T> * create_1_skeleton()
		{
			std::size_t n = this->_metric_space.size();

			basic_graph<T> * graph(new basic_graph<T> ());

			T distance;

			for (std::size_t i = 0; i < n; i++)
			{
				for (std::size_t j = i + 1; j < n; j++)
				{
					distance = this->_metric_space.distance(i, j);
					if (distance <= this->_max_filtration_value)
					{
						graph->add_edge(i, j, distance);
					}
				}
			}

			return graph;
		}

	};

}

#endif /* VIETORIS_RIPS_COMPLEX_H_ */
