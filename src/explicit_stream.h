//============================================================================
// Name        : cph
// Author      : Andrew Tausz <atausz@stanford.edu>
// Version     : 1.0
// Copyright   : Copyright © 2011 Andrew Tausz
// Description : A basic package for persistent homology in C++
//============================================================================

#ifndef EXPLICIT_STREAM_H_
#define EXPLICIT_STREAM_H_

#include "simplex_stream.h"
#include "simplex.h"

namespace cph
{

template<class V, class T>
class explicit_stream : public simplex_stream<simplex<V>, T>
{
public:
	explicit_stream(){}
	virtual ~explicit_stream(){}

	void add_simplex(const V & vertex)
	{
		simplex_stream<simplex<V>, T>::add_simplex(simplex<V>::make_simplex(vertex), T(0));
	}

	void add_simplex(const V & v1, const V & v2)
	{
		this->add_simplex(simplex<V>::make_simplex(v1, v2), 0);
	}

	void add_simplex(const V & vertex, const T & t)
	{
		this->add_simplex(simplex<V>::make_simplex(vertex), t);
	}

	void add_simplex(const V & v1, const V & v2, const T & t)
	{
		simplex_stream<simplex<V>, T>::add_simplex(simplex<V>::make_simplex(v1, v2), t);
	}
};

}

#endif /* EXPLICITSTREAM_H_ */
