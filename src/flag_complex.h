//============================================================================
// Name        : cph
// Author      : Andrew Tausz <atausz@stanford.edu>
// Version     : 1.0
// Copyright   : Copyright © 2011 Andrew Tausz
// Description : A basic package for persistent homology in C++
//============================================================================

#ifndef FLAG_COMPLEX_H_
#define FLAG_COMPLEX_H_

#include "simplex.h"
#include "simplex_stream.h"
#include "basic_graph.h"
#include "finite_metric_space.h"

#include <algorithm>

namespace cph
{

	template<class T>
	class flag_complex: public simplex_stream<simplex<std::size_t> , T>
	{

	protected:
		const finite_metric_space<T> & _metric_space;
		const T _max_filtration_value;
		const std::size_t _max_dimension;
		const std::size_t _vertex_set_size;

	public:
		flag_complex(const finite_metric_space<T> & metric_space, const T & max_filtration_value, const std::size_t max_dimension,
				const std::size_t vertex_set_size) :
			_metric_space(metric_space), _max_filtration_value(max_filtration_value), _max_dimension(max_dimension), _vertex_set_size(vertex_set_size)
		{
		}

		virtual ~flag_complex()
		{
		}

		void construct()
		{
			basic_graph<T> * graph = this->create_1_skeleton();
			this->incremental_expansion(graph, this->_max_dimension);
			this->ensure_sorted();
			delete (graph);
		}

	protected:
		virtual basic_graph<T> * create_1_skeleton() = 0;

		void incremental_expansion(basic_graph<T> * graph, const std::size_t k)
		{
			for (std::size_t u = 0; u < this->_vertex_set_size; u++)
			{
				std::set<std::size_t> * lower_neighbors = graph->lower_neighbors(u);
				this->add_cofaces(graph, k, simplex<std::size_t>::make_simplex(u), lower_neighbors, 0);
				delete (lower_neighbors);
			}
		}

		void add_cofaces(basic_graph<T> * graph, const std::size_t k, const simplex<std::size_t> & tau, const std::set<std::size_t> * N,
				const T filtration_value)
		{

			this->add_simplex(tau, filtration_value);

			if (tau.dimension() >= k)
			{
				return;
			}

			T weight(0);

			for (typename std::set<std::size_t>::iterator iterator = N->begin(); iterator != N->end(); ++iterator)
			{
				std::size_t v = (*iterator);

				simplex<std::size_t> sigma = tau.append_to(v);

				std::set<std::size_t> * M = graph->intersect_with_lower_neighbors(*N, v);

				if (sigma.dimension() == 1)
				{
					const std::size_t i = sigma[0];
					const std::size_t j = sigma[1];
					weight = graph->get_weight(i, j);
				}
				else if (sigma.dimension() > 1)
				{
					weight = filtration_value;
					for (std::size_t i = 0; i < tau.dimension() + 1; i++)
					{
						weight = std::max(weight, graph->get_weight(tau[i], v));
					}
				}

				this->add_cofaces(graph, k, sigma, M, weight);

				delete (M);
			}

		}

	};
}

#endif /* FLAG_COMPLEX_H_ */
