/*
 * explicit_metric_space.hpp
 *
 *  Created on: Nov 7, 2011
 *      Author: atausz
 */

#ifndef EXPLICIT_METRIC_SPACE_HPP_
#define EXPLICIT_METRIC_SPACE_HPP_

#include "finite_metric_space.h"
namespace cph
{

	template<class T>
	class explicit_metric_space: public finite_metric_space<T>
	{
	public:

	private:
		const basic_matrix<T> * _distance_matrix;
	public:
		explicit_metric_space(const basic_matrix<T> * distance_matrix) :
			_distance_matrix(distance_matrix)
		{
		}

		virtual ~explicit_metric_space()
		{
			delete (this->_distance_matrix);
		}

		const T distance(const std::size_t i, const std::size_t j) const
		{
			return this->_distance_matrix->operator()(i, j);
		}

		const std::size_t size() const
		{
			return _distance_matrix->rows();
		}

		template<class S>
		friend S & operator <<(S & s, const explicit_metric_space<T> & value)
		{
			s << *(value._distance_matrix);
			return s;
		}
	};

} /* namespace cph */
#endif /* EXPLICIT_METRIC_SPACE_HPP_ */
