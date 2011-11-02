//============================================================================
// Name        : cph
// Author      : Andrew Tausz <atausz@stanford.edu>
// Version     : 1.0
// Copyright   : Copyright © 2011 Andrew Tausz
// Description : A basic package for persistent homology in C++
//============================================================================

#include "cph_R_interface.h"
#include "cph_interface.h"
#include "basic_matrix.h"
#include "euclidean_metric_space.h"
#include "explicit_metric_space.h"
#include "metrics.h"

#include <Rcpp.h>
#include <vector>

SEXP default_euclidean_phom(SEXP _matrix, SEXP _dimension, SEXP _max_filtration_value, SEXP _metric_type, SEXP _power)
{
	Rcpp::NumericMatrix X_R(_matrix);

	int dimension = Rcpp::as<int>(_dimension);
	double max_filtration_value = Rcpp::as<double>(_max_filtration_value);
	cph::basic_matrix<double> * X = new cph::basic_matrix<double>(X_R.nrow(), X_R.ncol());

	cph::metric metric_type = (cph::metric) Rcpp::as<int>(_metric_type);
	double p = Rcpp::as<double>(_power);

	for (int i(0); i < X_R.nrow(); i++)
	{
		for (int j(0); j < X_R.ncol(); j++)
		{
			(*X)(i, j) = X_R(i, j);
		}
	}

	cph::euclidean_metric_space<double> metric_space(X, metric_type, p);
	cph::barcode_collection<double> intervals = cph::default_persistent_homology(metric_space, dimension, max_filtration_value);
	cph::basic_matrix<double> endpoint_matrix(intervals.get_endpoint_matrix(max_filtration_value));
	Rcpp::NumericMatrix endpoint_matrix_R(endpoint_matrix.rows(), endpoint_matrix.columns());

	for (std::size_t i(0); i < endpoint_matrix.rows(); i++)
	{
		for (std::size_t j(0); j < endpoint_matrix.columns(); j++)
		{
			endpoint_matrix_R(i, j) = endpoint_matrix.operator()(i, j);
		}
	}

	// NB: we don't have to delete X since metric_space will delete it
	// delete (X);

	return endpoint_matrix_R;
}

SEXP default_metric_phom(SEXP _matrix, SEXP _dimension, SEXP _max_filtration_value)
{
	Rcpp::NumericMatrix X_R(_matrix);

	int dimension = Rcpp::as<int>(_dimension);
	double max_filtration_value = Rcpp::as<double>(_max_filtration_value);
	cph::basic_matrix<double> * X = new cph::basic_matrix<double>(X_R.nrow(), X_R.ncol());

	for (int i(0); i < X_R.nrow(); i++)
	{
		for (int j(0); j < X_R.ncol(); j++)
		{
			(*X)(i, j) = X_R(i, j);
		}
	}

	cph::explicit_metric_space<double> metric_space(X);
	cph::barcode_collection<double> intervals = cph::default_persistent_homology(metric_space, dimension, max_filtration_value);
	cph::basic_matrix<double> endpoint_matrix(intervals.get_endpoint_matrix(max_filtration_value));
	Rcpp::NumericMatrix endpoint_matrix_R(endpoint_matrix.rows(), endpoint_matrix.columns());

	for (std::size_t i(0); i < endpoint_matrix.rows(); i++)
	{
		for (std::size_t j(0); j < endpoint_matrix.columns(); j++)
		{
			endpoint_matrix_R(i, j) = endpoint_matrix.operator()(i, j);
		}
	}

	// NB: we don't have to delete X since metric_space will delete it
	// delete (X);

	return endpoint_matrix_R;
}

SEXP vr_euclidean_phom(SEXP _matrix, SEXP _dimension, SEXP _max_filtration_value, SEXP _metric_type, SEXP _power)
{
	Rcpp::NumericMatrix X_R(_matrix);

	int dimension = Rcpp::as<int>(_dimension);
	double max_filtration_value = Rcpp::as<double>(_max_filtration_value);
	cph::basic_matrix<double> * X = new cph::basic_matrix<double>(X_R.nrow(), X_R.ncol());

	cph::metric metric_type = (cph::metric) Rcpp::as<int>(_metric_type);
	double p = Rcpp::as<double>(_power);

	for (int i(0); i < X_R.nrow(); i++)
	{
		for (int j(0); j < X_R.ncol(); j++)
		{
			(*X)(i, j) = X_R(i, j);
		}
	}

	cph::euclidean_metric_space<double> metric_space(X, metric_type, p);
	cph::barcode_collection<double> intervals = cph::vr_persistent_homology(metric_space, dimension, max_filtration_value);
	cph::basic_matrix<double> endpoint_matrix(intervals.get_endpoint_matrix(max_filtration_value));
	Rcpp::NumericMatrix endpoint_matrix_R(endpoint_matrix.rows(), endpoint_matrix.columns());

	for (std::size_t i(0); i < endpoint_matrix.rows(); i++)
	{
		for (std::size_t j(0); j < endpoint_matrix.columns(); j++)
		{
			endpoint_matrix_R(i, j) = endpoint_matrix.operator()(i, j);
		}
	}

	// NB: we don't have to delete X since metric_space will delete it
	// delete (X);

	return endpoint_matrix_R;
}

SEXP vr_metric_phom(SEXP _matrix, SEXP _dimension, SEXP _max_filtration_value)
{
	Rcpp::NumericMatrix X_R(_matrix);

	int dimension = Rcpp::as<int>(_dimension);
	double max_filtration_value = Rcpp::as<double>(_max_filtration_value);
	cph::basic_matrix<double> * X = new cph::basic_matrix<double>(X_R.nrow(), X_R.ncol());

	for (int i(0); i < X_R.nrow(); i++)
	{
		for (int j(0); j < X_R.ncol(); j++)
		{
			(*X)(i, j) = X_R(i, j);
		}
	}

	cph::explicit_metric_space<double> metric_space(X);
	cph::barcode_collection<double> intervals = cph::vr_persistent_homology(metric_space, dimension, max_filtration_value);
	cph::basic_matrix<double> endpoint_matrix(intervals.get_endpoint_matrix(max_filtration_value));
	Rcpp::NumericMatrix endpoint_matrix_R(endpoint_matrix.rows(), endpoint_matrix.columns());

	for (std::size_t i(0); i < endpoint_matrix.rows(); i++)
	{
		for (std::size_t j(0); j < endpoint_matrix.columns(); j++)
		{
			endpoint_matrix_R(i, j) = endpoint_matrix.operator()(i, j);
		}
	}

	// NB: we don't have to delete X since metric_space will delete it
	// delete (X);

	return endpoint_matrix_R;
}

SEXP lw_euclidean_phom(SEXP _matrix, SEXP _dimension, SEXP _max_filtration_value, SEXP _landmark_set_size, SEXP _maxmin_sample_size, SEXP _metric_type, SEXP _power)
{
	Rcpp::NumericMatrix X_R(_matrix);

	int dimension = Rcpp::as<int>(_dimension);
	double max_filtration_value = Rcpp::as<double>(_max_filtration_value);
	int landmark_set_size = Rcpp::as<int>(_landmark_set_size);
	int maxmin_sample_size = Rcpp::as<int>(_maxmin_sample_size);
	cph::basic_matrix<double> * X = new cph::basic_matrix<double>(X_R.nrow(), X_R.ncol());

	cph::metric metric_type = (cph::metric) Rcpp::as<int>(_metric_type);
	double p = Rcpp::as<double>(_power);

	for (int i(0); i < X_R.nrow(); i++)
	{
		for (int j(0); j < X_R.ncol(); j++)
		{
			(*X)(i, j) = X_R(i, j);
		}
	}

	cph::euclidean_metric_space<double> metric_space(X, metric_type, p);
	cph::barcode_collection<double> intervals = cph::lw_persistent_homology(metric_space, dimension, max_filtration_value, landmark_set_size,
			maxmin_sample_size);
	cph::basic_matrix<double> endpoint_matrix(intervals.get_endpoint_matrix(max_filtration_value));
	Rcpp::NumericMatrix endpoint_matrix_R(endpoint_matrix.rows(), endpoint_matrix.columns());

	for (std::size_t i(0); i < endpoint_matrix.rows(); i++)
	{
		for (std::size_t j(0); j < endpoint_matrix.columns(); j++)
		{
			endpoint_matrix_R(i, j) = endpoint_matrix.operator()(i, j);
		}
	}

	// NB: we don't have to delete X since metric_space will delete it
	// delete (X);

	return endpoint_matrix_R;
}

SEXP lw_metric_phom(SEXP _matrix, SEXP _dimension, SEXP _max_filtration_value, SEXP _landmark_set_size, SEXP _maxmin_sample_size)
{
	Rcpp::NumericMatrix X_R(_matrix);

	int dimension = Rcpp::as<int>(_dimension);
	double max_filtration_value = Rcpp::as<double>(_max_filtration_value);
	int landmark_set_size = Rcpp::as<int>(_landmark_set_size);
	int maxmin_sample_size = Rcpp::as<int>(_maxmin_sample_size);
	cph::basic_matrix<double> * X = new cph::basic_matrix<double>(X_R.nrow(), X_R.ncol());

	for (int i(0); i < X_R.nrow(); i++)
	{
		for (int j(0); j < X_R.ncol(); j++)
		{
			(*X)(i, j) = X_R(i, j);
		}
	}

	cph::explicit_metric_space<double> metric_space(X);
	cph::barcode_collection<double> intervals = cph::lw_persistent_homology(metric_space, dimension, max_filtration_value, landmark_set_size,
			maxmin_sample_size);
	cph::basic_matrix<double> endpoint_matrix(intervals.get_endpoint_matrix(max_filtration_value));
	Rcpp::NumericMatrix endpoint_matrix_R(endpoint_matrix.rows(), endpoint_matrix.columns());

	for (std::size_t i(0); i < endpoint_matrix.rows(); i++)
	{
		for (std::size_t j(0); j < endpoint_matrix.columns(); j++)
		{
			endpoint_matrix_R(i, j) = endpoint_matrix.operator()(i, j);
		}
	}

	// NB: we don't have to delete X since metric_space will delete it
	// delete (X);

	return endpoint_matrix_R;
}
