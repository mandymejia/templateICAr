#' Summarize a \code{"templateICA.cifti"} object
#'
#' Summary method for class \code{"templateICA.cifti"}
#'
#' @param object Object of class \code{"templateICA.cifti"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method summary templateICA.cifti
summary.templateICA.cifti <- function(object, ...) {
  x <- c(
    summary(object$subjICmean),
    object$params
  )

  class(x) <- "summary.templateICA.cifti"
  return(x)
}

#' @rdname summary.templateICA.cifti
#' @export
#'
#' @param x The template from \code{estimate_template.cifti}
#' @param ... further arguments passed to or from other methods.
#' @method print summary.templateICA.cifti
print.summary.templateICA.cifti <- function(x, ...) {
  # Get DCT output.
  dct <- x$detrend_DCT
  if (!is.null(dct)) {
    dct <- as.numeric(x$detrend_DCT)
    dct <- if (dct>1) {
      paste(dct, "DCT bases")
    } else if (dct > 0) {
      paste(dct, "DCT basis")
    } else {
      "None"
    }
  }

  cat("====TEMPLATE ICA INFO================\n")
  cat("Detrending:      ", dct, "\n")
  cat("Spatial scaling: ", x$scale, "\n")
  cat("A normalization: ", x$normA, "\n")
  cat("Variance method: ", x$var_method, "\n")
  cat("Q2 and Q2_max:   ", paste0(x$Q2, ", ", x$Q2_max), "\n")
  cat("-------------------------------------\n")
  cat("Spatial model:   ", x$spatial_model, "\n")
  cat("Dims reduced:    ", x$reduce_dim, "\n")
  cat("Maximum iters:   ", x$maxiter, "\n")
  cat("Epsilon:         ", x$epsilon, "\n")
  if (as.logical(x$spatial_model)) {
    cat("Mwall removed:   ", x$rm_mwall, "\n")
    cat("Initial Kappa:   ", x$kappa_init, "\n")
  }
  cat("\n")

  class(x) <- "summary.xifti"
  print(x)

  invisible(NULL)
}

#' @rdname summary.templateICA.cifti
#' @export
#'
#' @method print templateICA.cifti
print.templateICA.cifti <- function(x, ...) {
  print.summary.templateICA.cifti(summary(x))
}

#' Summarize a \code{"tICA"} object
#'
#' Summary method for class \code{"tICA"}
#'
#' @param object Object of class \code{"tICA"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method summary tICA
summary.tICA <- function(object, ...) {
  x <- c(
    list(nV=nrow(object$subjICmean), nL=ncol(object$subjICmean)),
    object$params
  )

  class(x) <- "summary.tICA"
  return(x)
}

#' @rdname summary.tICA
#' @export
#'
#' @param x The template from \code{estimate_template.cifti}
#' @param ... further arguments passed to or from other methods.
#' @method print summary.tICA
print.summary.tICA <- function(x, ...) {
  # Get DCT output.
  dct <- x$detrend_DCT
  if (!is.null(dct)) {
    dct <- as.numeric(x$detrend_DCT)
    dct <- if (dct>1) {
      paste(dct, "DCT bases")
    } else if (dct > 0) {
      paste(dct, "DCT basis")
    } else {
      "None"
    }
  }

  cat("====TEMPLATE INFO====================\n")
  cat("Detrending:      ", dct, "\n")
  cat("Spatial scaling: ", x$scale, "\n")
  cat("A normalization: ", x$normA, "\n")
  cat("Variance method: ", x$var_method, "\n")
  cat("Q2 and Q2_max:   ", paste0(x$Q2, ", ", x$Q2_max), "\n")
  cat("-------------------------------------\n")
  cat("Spatial model:   ", x$spatial_model, "\n")
  cat("Dims reduced:    ", x$reduce_dim, "\n")
  cat("Maximum iters:   ", x$maxiter, "\n")
  cat("Epsilon:         ", x$epsilon, "\n")
  cat("-------------------------------------\n")
  cat("# Locations:     ", x$nL, "\n")
  cat("# Template ICs:  ", x$nV, "\n")
  cat("\n")

  invisible(NULL)
}

#' @rdname summary.tICA
#' @export
#'
#' @method print tICA
print.tICA <- function(x, ...) {
  print.summary.tICA(summary(x))
}
