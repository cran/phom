//============================================================================
// Name        : cph
// Author      : Andrew Tausz <atausz@stanford.edu>
// Version     : 1.0
// Copyright   : Copyright © 2011 Andrew Tausz
// Description : A basic package for persistent homology in C++
//============================================================================

#ifndef CPP_R_INTERFACE_HPP_
#define CPP_R_INTERFACE_HPP_

#include <Rcpp.h>

/*
 * note : RcppExport is an alias to `extern "C"` defined by Rcpp.
 *
 * It gives C calling convention to the rcpp_hello_world function so that
 * it can be called from .Call in R. Otherwise, the C++ compiler mangles the
 * name of the function and .Call can't find it.
 *
 * It is only useful to use RcppExport when the function is intended to be called
 * by .Call. See the thread http://thread.gmane.org/gmane.comp.lang.r.rcpp/649/focus=672
 * on Rcpp-devel for a misuse of RcppExport
 */

RcppExport SEXP default_euclidean_phom(SEXP _matrix, SEXP _dimension, SEXP _max_filtration_value, SEXP _metric_type, SEXP _power);
RcppExport SEXP default_metric_phom(SEXP _matrix, SEXP _dimension, SEXP _max_filtration_value);
RcppExport SEXP vr_euclidean_phom(SEXP _matrix, SEXP _dimension, SEXP _max_filtration_value, SEXP _metric_type, SEXP _power);
RcppExport SEXP vr_metric_phom(SEXP _matrix, SEXP _dimension, SEXP _max_filtration_value);
RcppExport SEXP lw_euclidean_phom(SEXP _matrix, SEXP _dimension, SEXP _max_filtration_value, SEXP _landmark_set_size, SEXP _maxmin_sample_size, SEXP _metric_type, SEXP _power);
RcppExport SEXP lw_metric_phom(SEXP _matrix, SEXP _dimension, SEXP _max_filtration_value, SEXP _landmark_set_size, SEXP _maxmin_sample_size);



#endif /* CPP_R_INTERFACE_HPP_ */
