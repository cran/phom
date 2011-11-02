pHom <- function(X, dimension, max_filtration_value, mode="vr", metric="euclidean", p = 2, landmark_set_size = 2 * ceiling(sqrt(length(X))), maxmin_samples = min(1000, length(X))) {
	
	modes <- c("vr", "lw")
	mode_index = pmatch(mode, modes)
	if (is.na(mode_index)) {
    		stop("Invalid mode specified.")
	}
	if (mode_index == -1) {
		stop("Ambiguous mode specified.")
	}


	metrics <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "distance_matrix")
	metric_index = pmatch(metric, metrics)
	if (is.na(metric_index)) {
    		stop("Invalid metric specified.")
	}
	if (metric_index == -1) {
		stop("Ambiguous metric specified.")
	}


	if (metric_index < 7) {
	# R^n points with given metric from 1-6
		if (mode_index == 1) {
		# VR
			out <- .Call( "vr_euclidean_phom", X, dimension, max_filtration_value, metric_index, p, PACKAGE = "phom" )
			return (out)
		} else {
		# LW
			out <- .Call( "lw_euclidean_phom", X, dimension, max_filtration_value, landmark_set_size, maxmin_samples, metric_index, p, PACKAGE = "phom" )
			return (out)
		}
	}

	# explicit distance matrix - we must have metric_index == 7
	if (mode_index == 1) {
		out <- .Call( "vr_metric_phom", X, dimension, max_filtration_value, PACKAGE = "phom" )
		return (out)
	} else {
		out <- .Call( "lw_metric_phom", X, dimension, max_filtration_value, landmark_set_size, maxmin_samples, PACKAGE = "phom" )
		return (out)
	}
}

plotPersistenceDiagram <- function(intervals, max_dim, max_f, title="Persistence Diagram") {

	plot_colors <- rainbow(max_dim + 1)
	plot_symbols <- c(1:20)

	for(i in 0:max_dim) {
		  indices <- which(intervals[, 1]==i, arr.ind=T)
		  dim_i_intervals <- intervals[indices, 2:3, drop=FALSE]

		  if (length(dim_i_intervals) > 0) {
		  	 if (i==0) {
		  	 	plot(dim_i_intervals[, 1], dim_i_intervals[, 2], pch=plot_symbols[i+1], col=plot_colors[i+1], xlim=c(0, max_f), ylim=c(0, max_f), xlab="Interval Start", ylab="Interval End", main=title)
			 } else {
		  	    points(dim_i_intervals[, 1], dim_i_intervals[, 2], pch=plot_symbols[i+1], col=plot_colors[i+1], xlim=c(0, max_f), ylim=c(0, max_f))
		  	 }
		  }
	}

	abline(a=0, b=1)
	legend(max_f * 0.8, max_f * 0.2, c(0:max_dim), cex=0.8, pch=plot_symbols[1:(max_dim+1)], col=plot_colors[1:(max_dim+1)])
}

plotBarcodeDiagram <- function(intervals, dimension, max_f, title="Barcode Diagram") {

	i <- dimension
	
	indices <- which(intervals[, 1]==i, arr.ind=T)
	dim_i_intervals <- intervals[indices, 2:3, drop=FALSE]
	K <- nrow(dim_i_intervals)

	if (K > 0) {
		
		plot(c(0,max_f),c(0,K),ann=F,bty='n',type='n',xaxt='n',yaxt='n')
		for(k in 1:K) {
			segments(x0 = dim_i_intervals[k, 1], x1 = dim_i_intervals[k, 2], y0 = k, y1 = k)
		}
		axis(1, max_f * c(0:5) / 5);
		title(main = title)
	}
}

