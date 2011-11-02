//============================================================================
// Name        : cph
// Author      : Andrew Tausz <atausz@stanford.edu>
// Version     : 1.0
// Copyright   : Copyright © 2011 Andrew Tausz
// Description : A basic package for persistent homology in C++
//============================================================================

#ifndef CPP_INTERFACE_HPP_
#define CPP_INTERFACE_HPP_

#include "simplex.h"
#include "utility.h"
#include "finite_metric_space.h"
#include "vietoris_rips_complex.h"
#include "persistence_algorithm.h"
#include "lazy_witness_complex.h"
#include "landmark_selector.h"
#include "basic_matrix.h"
#include "barcode_collection.h"

namespace cph
{
	template<class T>
	T estimate_diameter(const finite_metric_space<T> & metric_space);
	template<class T>
	barcode_collection<T> default_persistent_homology(const finite_metric_space<T> & metric_space, const std::size_t dimension, const T max_filtration_value);
	template<class T>
	barcode_collection<T> vr_persistent_homology(const finite_metric_space<T> & metric_space, const std::size_t dimension, const T max_filtration_value);
	template<class T>
	barcode_collection<T> lw_persistent_homology(const finite_metric_space<T> & metric_space, const std::size_t dimension, const T max_filtration_value,
			const std::size_t landmark_set_size = 50, const std::size_t maxmin_samples = 100);

	template<class T>
	T estimate_diameter(const finite_metric_space<T> & metric_space)
	{
		T estimated_diameter = metric_space.estimate_diameter();
		return estimated_diameter;
	}

	template<class T>
	barcode_collection<T> default_persistent_homology(const finite_metric_space<T> & metric_space, const std::size_t dimension, const T max_filtration_value)
	{
		if (metric_space.size() <= 100)
		{
			return cph::vr_persistent_homology(metric_space, dimension, max_filtration_value);
		}

		std::size_t landmark_set_size(0);

		if (metric_space.size() > 1000)
		{
			landmark_set_size = 100;
		}
		else
		{
			landmark_set_size = (std::size_t) (2 * std::ceil(std::sqrt(T(metric_space.size()))));
		}

		std::size_t maxmin_samples(1000);

		return cph::lw_persistent_homology(metric_space, dimension, max_filtration_value, landmark_set_size, maxmin_samples);
	}

	template<class T>
	barcode_collection<T> vr_persistent_homology(const finite_metric_space<T> & metric_space, const std::size_t dimension, const T max_filtration_value)
	{
		vietoris_rips_complex<T> complex(metric_space, max_filtration_value, dimension + 1);
		complex.construct();

		persistence_algorithm<simplex<typename std::size_t> , T> persistence(dimension);
		barcode_collection<T> intervals = persistence.compute_intervals(complex);

		return intervals;
	}

	template<class T>
	barcode_collection<T> lw_persistent_homology(const finite_metric_space<T> & metric_space, const std::size_t dimension, const T max_filtration_value,
			const std::size_t landmark_set_size = 50, const std::size_t maxmin_samples = 100)
	{
		std::vector<std::size_t> landmark_selection;

		if (metric_space.size() <= maxmin_samples)
		{
			landmark_selection = landmark_selector::maxmin_landmark_selection<double>(metric_space, landmark_set_size);
		}
		else
		{
			landmark_selection = landmark_selector::approx_maxmin_landmark_selection<double>(metric_space, landmark_set_size, maxmin_samples);
		}

		lazy_witness_complex<T> complex(metric_space, landmark_selection, max_filtration_value, dimension + 1);
		complex.construct();

		persistence_algorithm<simplex<typename std::size_t> , T> persistence(dimension);
		barcode_collection<T> intervals = persistence.compute_intervals(complex);

		return intervals;
	}

} /* namespace cph */
#endif /* CPP_INTERFACE_HPP_ */
