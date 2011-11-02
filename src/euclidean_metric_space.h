//============================================================================
// Name        : cph
// Author      : Andrew Tausz <atausz@stanford.edu>
// Version     : 1.0
// Copyright   : Copyright © 2011 Andrew Tausz
// Description : A basic package for persistent homology in C++
//============================================================================

#ifndef EUCLIDEAN_METRIC_SPACE_H_
#define EUCLIDEAN_METRIC_SPACE_H_

#include "basic_matrix.h"
#include "utility.h"
#include "finite_metric_space.h"
#include "metrics.h"

#include <algorithm>

namespace cph
{

	template<class T>
	class euclidean_metric_space: public finite_metric_space<T>
	{

	private:
		const basic_matrix<T> * _points;
		const metric _metric_type;
		const T _p;
	public:
		euclidean_metric_space(const basic_matrix<T> * points, const metric metric_type = euclidean, const T p = 2) :
			_points(points), _metric_type(metric_type), _p(p)
		{
		}

		virtual ~euclidean_metric_space()
		{
			delete (this->_points);
		}

		const T distance(const std::size_t i, const std::size_t j) const
		{
			return this->_points->row_distance(i, j, this->_metric_type, this->_p);
		}

		const std::size_t size() const
		{
			return _points->rows();
		}

		template<class S>
		friend S & operator <<(S & s, const euclidean_metric_space<T> & value)
		{
			s << *(value._points);
			return s;
		}
	};
}

#endif /* EUCLIDEAN_METRIC_SPACE_H_ */
