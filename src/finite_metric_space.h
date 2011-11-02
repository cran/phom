//============================================================================
// Name        : cph
// Author      : Andrew Tausz <atausz@stanford.edu>
// Version     : 1.0
// Copyright   : Copyright © 2011 Andrew Tausz
// Description : A basic package for persistent homology in C++
//============================================================================

#ifndef FINITE_METRIC_SPACE_H_
#define FINITE_METRIC_SPACE_H_

#include "random_utility.h"

namespace cph
{

	template<class T>
	class finite_metric_space
	{
	public:
		finite_metric_space()
		{
		}
		virtual ~finite_metric_space()
		{
		}

		virtual const T distance(const std::size_t i, const std::size_t j) const = 0;
		virtual const std::size_t size() const = 0;

		const T estimate_diameter(const std::size_t samples = 100) const
		{
			T max(0), distance(0);

			for (std::size_t s = 0; s < samples; s++)
			{
				std::size_t i = random_utility::random_integer(0, this->size() - 1);
				std::size_t j = random_utility::random_integer(0, this->size() - 1);

				distance = this->distance(i, j);
				if (distance > max)
				{
					max = distance;
				}
			}

			return max;
		}
	};

}

#endif /* FINITE_METRIC_SPACE_H_ */
