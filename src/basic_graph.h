//============================================================================
// Name        : cph
// Author      : Andrew Tausz <atausz@stanford.edu>
// Version     : 1.0
// Copyright   : Copyright © 2011 Andrew Tausz
// Description : A basic package for persistent homology in C++
//============================================================================

#ifndef BASIC_GRAPH_H_
#define BASIC_GRAPH_H_

#include <map>
#include <set>

namespace cph
{

	template<class W, class V = std::size_t>
	class basic_graph
	{
	private:
		std::map<V, std::map<V, W> * > _adjacency_sets;

	public:
		basic_graph()
		{
		}
		virtual ~basic_graph()
		{
			for (typename std::map<V, std::map<V, W> * >::iterator iter = this->_adjacency_sets.begin(); iter != this->_adjacency_sets.end(); iter++)
			{
				delete iter->second;
			}
		}

		void add_edge(const V & i, const V & j, const W w)
		{
			V x = (i < j ? i : j);
			V y = (i < j ? j : i);

			if (_adjacency_sets.find(y) == _adjacency_sets.end())
			{
				_adjacency_sets[y] = new std::map<V, W>();
			}
			(_adjacency_sets[y])->operator[](x) = w;
		}

		const W get_weight(const V & i, const V & j) const
		{
			V x = (i < j ? i : j);
			V y = (i < j ? j : i);

			if (_adjacency_sets.find(y) == _adjacency_sets.end())
			{
				return W(0);
			}
			return this->_adjacency_sets.at(y)->at(x);
		}

		std::set<V> * intersect_with_lower_neighbors(const std::set<V> & set, const V & y) const
		{
			std::set<V> * result(new std::set<V>());

			if (_adjacency_sets.find(y) == _adjacency_sets.end())
			{
				return result;
			}

			const std::map<V, W> * map = _adjacency_sets.at(y);

			for (typename std::map<V, W>::const_iterator iter = map->begin(); iter != map->end(); iter++)
			{
				if (set.find((*iter).first) != set.end())
				{
					result->insert((*iter).first);
				}
			}

			return result;
		}

		std::set<V> * lower_neighbors(const V & y) const
		{
			std::set<V> * result(new std::set<V>());

			if (_adjacency_sets.find(y) == _adjacency_sets.end())
			{
				return result;
			}

			const std::map<V, W> * weighted_map = _adjacency_sets.at(y);

			typename std::map<V, W>::const_iterator iter = weighted_map->begin();
			typename std::map<V, W>::const_iterator end = weighted_map->end();

			for (; iter != end; iter++)
			{
				result->insert((*iter).first);
			}

			return result;
		}
	};

}

#endif /* BASIC_GRAPH_H_ */
