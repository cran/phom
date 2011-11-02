//============================================================================
// Name        : cph
// Author      : Andrew Tausz <atausz@stanford.edu>
// Version     : 1.0
// Copyright   : Copyright © 2011 Andrew Tausz
// Description : A basic package for persistent homology in C++
//============================================================================

#ifndef METRICS_H_
#define METRICS_H_

namespace cph
{
	enum metric
	{
		euclidean = 1, maximum = 2, manhattan = 3, canberra = 4, binary = 5, minkowski = 6
	};
}

#endif /* METRICS_H_ */
