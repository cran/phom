//============================================================================
// Name        : cph
// Author      : Andrew Tausz <atausz@stanford.edu>
// Version     : 1.0
// Copyright   : Copyright © 2011 Andrew Tausz
// Description : A basic package for persistent homology in C++
//============================================================================

#ifndef UTILITY_H_
#define UTILITY_H_

#include <algorithm>
#include <vector>
#include <set>

namespace cph
{

	template<class S, class T>
	S & operator <<(S & s, const std::vector<T> & v)
	{
		s << '{';

		for (typename std::vector<T>::const_iterator iter = v.begin(); iter != v.end(); iter++)
		{
			if (iter != v.begin())
			{
				s << ',';
			}

			s << *iter;
		}

		s << '}';

		return s;
	}

	template<class S, class T>
	S & operator <<(S & s, const std::set<T> & v)
	{
		s << '{';

		for (typename std::set<T>::const_iterator iter = v.begin(); iter != v.end(); iter++)
		{
			if (iter != v.begin())
			{
				s << ',';
			}

			s << *iter;
		}

		s << '}';

		return s;
	}

}

#endif /* UTILITY_H_ */
