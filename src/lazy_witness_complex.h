//============================================================================
// Name        : cph
// Author      : Andrew Tausz <atausz@stanford.edu>
// Version     : 1.0
// Copyright   : Copyright © 2011 Andrew Tausz
// Description : A basic package for persistent homology in C++
//============================================================================

#ifndef LAZY_WITNESS_COMPLEX_HPP_
#define LAZY_WITNESS_COMPLEX_HPP_

#include <vector>
#include <algorithm>

#include "flag_complex.h"
#include "basic_graph.h"
#include "basic_matrix.h"

namespace cph
{

template<class T>
class lazy_witness_complex: public flag_complex<T>
{
protected:
	const std::vector<std::size_t> & _landmark_selection;
	const std::size_t _nu;
public:
	lazy_witness_complex(const finite_metric_space<T> & metric_space,
			const std::vector<std::size_t> & landmark_selection,
			const T & max_filtration_value,
			const int max_dimension) :
			flag_complex<T>(metric_space, max_filtration_value, max_dimension, landmark_selection.size()),
			_landmark_selection(landmark_selection), _nu(2)
	{
	}

	virtual ~lazy_witness_complex()
	{
	}

	virtual basic_graph<T> * create_1_skeleton()
	{
		std::size_t N = this->_metric_space.size();
		std::size_t L = this->_landmark_selection.size();

		basic_graph<T> * graph(new basic_graph<T>());

		/*
		 * Let N be the number of points in the metric space, and n the number of
		 * landmark points. Let D be the L x N matrix of distances between the set
		 * of landmark points, and the set of all points in the metric space.
		 *
		 * The definition of the 1-skeleton of the lazy witness complex is as follows:
		 *
		 * - If nu = 0, then define m_i = 0, otherwise define m_i to be the nu-th smallest entry
		 * in the i-th column of D.
		 * - The edge [ab] belongs to W(D, R, nu) iff there exists as witness i in {1, ..., N} such
		 * that max(D(a, i), D(b, i)) <= R + m_i
		 *
		 */

		std::vector<double> m(N);
		std::vector<std::size_t> m_temp;

		basic_matrix<double> D(L, N);

		for (std::size_t l = 0; l < L; l++)
		{
			std::size_t landmark_point = this->_landmark_selection[l];
			for (std::size_t n = 0; n < N; n++)
			{
				D(l, n) = this->_metric_space.distance(landmark_point, n);
			}
		}

		if (this->_nu > 0)
		{
			//T m_temp[L + 1];

			for (std::size_t n = 0; n < N; n++)
			{
				m_temp.clear();
				m_temp.push_back(0);
				for (std::size_t l = 0; l < L; l++)
				{
					m_temp.push_back(D(l, n));
				}
				std::sort(m_temp.begin(), m_temp.end());
				m[n] = m_temp[this->_nu];
			}

		}

		T e_ij(0), d_ij(0);
		T d_in(0), d_jn(0);

		for (std::size_t i = 0; i < L; i++)
		{
			for (std::size_t j = i + 1; j < L; j++)
			{
				e_ij = 0;
				for (std::size_t n = 0; n < N; n++)
				{
					d_in = D(i, n);
					d_jn = D(j, n);
					d_ij = (d_in > d_jn ? d_in: d_jn);
					if (d_ij < m[n])
					{
						d_ij = 0;
					}
					else
					{
						d_ij -= m[n];
					}

					if (n == 0 || (d_ij < e_ij))
					{
						e_ij = d_ij;
					}
				}

				if (e_ij < this->_max_filtration_value)
				{
					graph->add_edge(i, j, e_ij);
				}
			}
		}


		return graph;
	}
};

} /* namespace cph */
#endif /* LAZY_WITNESS_COMPLEX_HPP_ */
