//============================================================================
// Name        : cph
// Author      : Andrew Tausz <atausz@stanford.edu>
// Version     : 1.0
// Copyright   : Copyright © 2011 Andrew Tausz
// Description : A basic package for persistent homology in C++
//============================================================================

#ifndef LANDMARK_SELECTOR_HPP_
#define LANDMARK_SELECTOR_HPP_

#include "random_utility.h"

#include <vector>
#include <algorithm>

namespace cph
{

	class landmark_selector
	{
	public:
		landmark_selector()
		{
		}
		virtual ~landmark_selector()
		{
		}

		static std::vector<std::size_t> random_landmark_selection(const std::size_t range_size, const std::size_t selection_size)
		{
			std::vector<typename std::size_t> range;
			for (std::size_t i = 0; i < range_size; i++)
			{
				range.push_back(i);
			}

			std::random_shuffle(range.begin(), range.end());

			std::size_t filtered_selection_size = (selection_size > range_size ? range_size : selection_size);

			return std::vector<typename std::size_t>(range.begin(), range.begin() + filtered_selection_size);
		}

		template<class T>
		static std::vector<std::size_t> maxmin_landmark_selection(const finite_metric_space<T> & metric_space, const std::size_t selection_size)
		{
			std::size_t initial_point = random_utility::random_integer(0, metric_space.size() - 1);

			std::vector<std::size_t> indices;

			indices.push_back(initial_point);

			/*
			 * Construct the landmark set inductively. Suppose that
			 * {l_0, ..., l_{i-1}} have been chosen as landmark points.
			 * Define the function
			 * f(z) = min{d(z, l_0), ...., d(z, l_{i-1}}
			 * and define l_i to be l_i = arg max f(z)
			 *
			 */
			while (indices.size() < selection_size)
			{
				// find point that maximizes the minimum distance to the existing landmark points
				T f_value = T(0);
				T max_f_value = T(0);
				std::size_t arg_max(0);
				for (std::size_t z = 0; z < metric_space.size(); z++)
				{
					f_value = landmark_selector::compute_min_distance(metric_space, indices, z);
					if (f_value > max_f_value)
					{
						arg_max = z;
						max_f_value = f_value;
					}
				}

				indices.push_back(arg_max);
			}

			return indices;
		}

		template<class T>
		static std::vector<std::size_t> approx_maxmin_landmark_selection(const finite_metric_space<T> & metric_space, const std::size_t selection_size,
				const std::size_t sample_size = 100)
		{
			std::size_t initial_point = random_utility::random_integer(0, metric_space.size() - 1);

			std::vector<std::size_t> indices;

			indices.push_back(initial_point);

			const std::size_t filtered_sample_size(sample_size < metric_space.size() ? sample_size : metric_space.size());

			/*
			 * Construct the landmark set inductively. Suppose that
			 * {l_0, ..., l_{i-1}} have been chosen as landmark points.
			 * Define the function
			 * f(z) = min{d(z, l_0), ...., d(z, l_{i-1}}
			 * and define l_i to be l_i = arg max f(z)
			 *
			 */
			while (indices.size() < selection_size)
			{
				// find point that maximizes the minimum distance to the existing landmark points
				T f_value = T(0);
				T max_f_value = T(0);
				std::size_t arg_max(0);
				for (std::size_t z_index = 0; z_index < filtered_sample_size; z_index++)
				{
					std::size_t z = random_utility::random_integer(0, metric_space.size() - 1);
					f_value = landmark_selector::compute_min_distance(metric_space, indices, z);
					if (f_value > max_f_value)
					{
						arg_max = z;
						max_f_value = f_value;
					}
				}

				indices.push_back(arg_max);
			}

			return indices;
		}

	private:
		template<class T>
		static inline T compute_min_distance(const finite_metric_space<T> & metric_space, const std::vector<std::size_t> & landmark_set,
				const std::size_t query_point)
		{
			T min_distance = T(0);
			T distance = T(0);
			for (typename std::vector<std::size_t>::const_iterator iter = landmark_set.begin(); iter != landmark_set.end(); iter++)
			{
				if (iter == landmark_set.begin())
				{
					min_distance = metric_space.distance(*iter, query_point);
				}
				else
				{
					distance = metric_space.distance(*iter, query_point);
					if (min_distance > distance)
					{
						min_distance = distance;
					}
				}
			}

			return min_distance;
		}
	};

} /* namespace cph */
#endif /* LANDMARK_SELECTOR_HPP_ */
