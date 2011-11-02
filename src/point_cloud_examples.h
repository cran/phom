//============================================================================
// Name        : cph
// Author      : Andrew Tausz <atausz@stanford.edu>
// Version     : 1.0
// Copyright   : Copyright © 2011 Andrew Tausz
// Description : A basic package for persistent homology in C++
//============================================================================

#ifndef POINT_CLOUD_EXAMPLES_H_
#define POINT_CLOUD_EXAMPLES_H_

#include "basic_matrix.h"
#include "random_utility.h"
#include "euclidean_metric_space.h"

namespace cph
{

	class point_cloud_examples
	{
	public:
		point_cloud_examples();
		virtual ~point_cloud_examples();

		static euclidean_metric_space<double> create_random_sphere_points(const int n, const int d)
		{
			basic_matrix<double> * m = new basic_matrix<double> (n, d + 1);

			for (int i = 0; i < n; i++)
			{

				for (int j = 0; j < d + 1; j++)
				{
					(*m)(i, j) = random_utility::random_normal();
				}

				double row_norm = m->row_norm(i);

				for (int j = 0; j < d + 1; j++)
				{
					(*m)(i, j) /= row_norm;
				}
			}

			return euclidean_metric_space<double>(m);
		}
	};
}

#endif /* POINT_CLOUD_EXAMPLES_H_ */
