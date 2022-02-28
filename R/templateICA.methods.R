#' Summarize a \code{"tICA.cifti"} object
#'
#' Summary method for class \code{"tICA.cifti"}
#'
#' @param object Object of class \code{"tICA.cifti"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method summary tICA.cifti
summary.tICA.cifti <- function(object, ...) {
  x <- c(
    summary(object$subjICmean),
    object$params
  )

  class(x) <- "summary.tICA.cifti"
  return(x)
}

#' Summarize a \code{"tICA.nifti"} object
#'
#' Summary method for class \code{"tICA.nifti"}
#'
#' @param object Object of class \code{"tICA.nifti"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method summary tICA.nifti
summary.tICA.nifti <- function(object, ...) {
  x <- c(
    list(
      mask_dims=dim(object$mask),
      nV=nrow(object$subjICmean), 
      nL=ncol(object$subjICmean)
    ),
    object$params
  )

  class(x) <- "summary.tICA.nifti"
  return(x)
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

#' @rdname summary.tICA.cifti
#' @export
#'
#' @param x The result of \code{templateICA} with CIFTI data
#' @param ... further arguments passed to or from other methods.
#' @method print summary.tICA.cifti
print.summary.tICA.cifti <- function(x, ...) {
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
  cat("Variance method: ", x$tvar_method, "\n")
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

#' @rdname summary.tICA.nifti
#' @export
#'
#' @param x The template from \code{estimate_template.cifti}
#' @param ... further arguments passed to or from other methods.
#' @method print summary.tICA.nifti
print.summary.tICA.nifti <- function(x, ...) {
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
  cat("Variance method: ", x$tvar_method, "\n")
  cat("Q2 and Q2_max:   ", paste0(x$Q2, ", ", x$Q2_max), "\n")
  cat("-------------------------------------\n")
  cat("Spatial model:   ", x$spatial_model, "\n")
  cat("Dims reduced:    ", x$reduce_dim, "\n")
  cat("Maximum iters:   ", x$maxiter, "\n")
  cat("Epsilon:         ", x$epsilon, "\n")
  cat("-------------------------------------\n")
  cat("Mask dims:       ", paste0(x$mask_dims, collapse=" x "), "\n")
  cat("Vectorized dims:\n")
  cat("# Locations:     ", x$nV, "\n")
  cat("# Template ICs:  ", x$nL, "\n")
  cat("\n")

  invisible(NULL)
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
  cat("Variance method: ", x$tvar_method, "\n")
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

#' @rdname summary.tICA.cifti
#' @export
#'
#' @method print tICA.cifti
print.tICA.cifti <- function(x, ...) {
  print.summary.tICA.cifti(summary(x))
}

#' @rdname summary.tICA.nifti
#' @export
#'
#' @method print tICA.nifti
print.tICA.nifti <- function(x, ...) {
  print.summary.tICA.nifti(summary(x))
}

#' @rdname summary.tICA
#' @export
#'
#' @method print tICA
print.tICA <- function(x, ...) {
  print.summary.tICA(summary(x))
}

#' Plot template
#'
#' @param x The result of \code{templateICA} with CIFTI data
#' @param stat \code{"mean"} (default), \code{"se"}, or \code{"both"}
#' @param ... Additional arguments to \code{view_xifti}
#' @return The plot
#' @export
#' @importFrom ciftiTools view_xifti
#' @method plot tICA.cifti
plot.tICA.cifti <- function(x, stat=c("mean", "se", "both"), ...) {
  stopifnot(inherits(x, "tICA.cifti"))

  # Check `...`
  args <- list(...)
  has_title <- "title" %in% names(args)
  has_idx <- "idx" %in% names(args)
  has_fname <- "fname" %in% names(args)

  # Check `stat`
  stat <- tolower(stat)
  if (has_idx && length(args$idx)>1 && !("fname" %in% names(args))) {
    if (identical(stat, c("mean", "se", "both"))) {
      stat <- "mean"
    } else {
      stat <- match.arg(stat, c("mean", "se", "both"))
    }
    if (stat == "both") {
      if (!("fname" %in% names(args))) {
        warning(
          "For multiple `idx`, use one call to plot() ",
          "for the mean template, ",
          "and a separate one for the seiance template. ",
          "Showing the mean template now."
        )
        stat <- "mean"
      }
    }
  }
  stat <- match.arg(stat, c("mean", "se", "both"))

  # Print message saying what's happening.
  msg1 <- ifelse(has_idx,
    "Plotting the",
    "Plotting the first component's"
  )
  msg2 <- switch(stat,
    both="estimate and standard error.",
    mean="estimate.",
    se="standard error."
  )
  cat(msg1, msg2, "\n")

  # Plot
  out <- list(mean=NULL, se=NULL)
  if (stat == "both") { stat <- c("mean", "se") }
  for (ss in stat) {
    args_ss <- args
    tsfx_ss <- c(mean="", se=" (se)")[ss]
    # Handle title and idx
    if (!has_title && !has_idx) {
      c1name <- if (!is.null(x$subjICmean$meta$cifti$names)) {
        x$subjICmean$meta$cifti$names[1]
      } else {
        "First component"
      }
      args_ss$title <- paste0(c1name, tsfx_ss)
    } else if (!has_idx) {
      args_ss$title <- paste0(args_ss$title, tsfx_ss)
    }
    # Handle fname
    if (has_fname) {
      fext <- if (grepl("html$", args_ss$fname[1])) {
        "html"
      } else if (grepl("pdf$", args_ss$fname[1])) {
        "pdf"
      } else {
        "png"
      }
      args_ss$fname <- gsub(paste0(".", fext), "", args_ss$fname, fixed=TRUE)
      args_ss$fname <- paste0(args_ss$fname, "_", ss, ".", fext)
    }
    out[[ss]] <- do.call(
      view_xifti, c(list(x[[paste0("subjIC", ss)]]), args_ss)
    )
  }

  invisible(out)
}

#' Plot template
#'
#' @param x The result of \code{templateICA} with NIFTI data
#' @param ... Additional arguments
#' @return The plot
#' @export
#' @method plot tICA.nifti
plot.tICA.nifti <- function(x, ...) {
  stopifnot(inherits(x, "tICA.nifti"))
  stop("Not supported yet.")
}

#' Plot template
#'
#' @param x The result of \code{templateICA} with NIFTI data
#' @param ... Additional arguments
#' @return The plot
#' @export
#' @method plot tICA
plot.tICA <- function(x, ...) {
  stopifnot(inherits(x, "tICA"))
  stop("Not supported yet.")
}