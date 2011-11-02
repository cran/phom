//============================================================================
// Name        : cph
// Author      : Andrew Tausz <atausz@stanford.edu>
// Version     : 1.0
// Copyright   : Copyright © 2011 Andrew Tausz
// Description : A basic package for persistent homology in C++
//============================================================================

#ifndef RANDOM_UTILITY_HPP_
#define RANDOM_UTILITY_HPP_

#include <stdlib.h>
#include <math.h>

namespace cph
{
	class random_utility
	{
	private:
		random_utility()
		{

		}
	public:
		virtual ~random_utility()
		{
		}

		static std::size_t random_integer(const std::size_t lowest, const std::size_t highest)
		{
			return (rand() % (highest - lowest + 1)) + lowest;
		}

	  static double random_uniform()
	  {
		return ((double) rand()) / ((double) RAND_MAX);
	  }

	  static double random_normal()
	  {
		std::size_t n(100);
		double inv_variance(12.0);
		double mean(0.5);
		double sum(0);

		for (std::size_t i(0); i < n; i++)
		  {
			sum += random_uniform();
		  }

		return sqrt(inv_variance * n) * ((sum / n) - mean);
		
	  }
	};

}

#endif /* RANDOM_UTILITY_HPP_ */
