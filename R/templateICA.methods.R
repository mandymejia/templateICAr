#' Summarize a \code{"tICA.cifti"} object
#'
#' Summary method for class \code{"tICA.cifti"}
#'
#' @param object Object of class \code{"tICA.cifti"}.
#' @param ... further arguments passed to or from other methods.
#' @return A list summarizing of the results of the templateICA analysis.
#' @export
#' @method summary tICA.cifti
summary.tICA.cifti <- function(object, ...) {
  tICA_params <- lapply(
    object$params,
    function(q) {
      if (is.null(q)) { q <- "NULL"};
      paste0(as.character(q), collapse=" ")
    }
  )

  x <- c(
    summary(object$subjICmean),
    tICA_params
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
#' @return A list summarizing of the results of the templateICA analysis.
#' @export
#' @method summary tICA.nifti
summary.tICA.nifti <- function(object, ...) {
  tICA_params <- lapply(
    object$params,
    function(q) {
      if (is.null(q)) { q <- "NULL"};
      paste0(as.character(q), collapse=" ")
    }
  )

  x <- c(
    list(
      mask_dims=dim(object$mask_nii),
      nV=nrow(object$template_mean),
      nL=ncol(object$template_mean)
    ),
    tICA_params
  )

  class(x) <- "summary.tICA.nifti"
  return(x)
}

#' Summarize a \code{"tICA.matrix"} object
#'
#' Summary method for class \code{"tICA.matrix"}
#'
#' @param object Object of class \code{"tICA.matrix"}.
#' @param ... further arguments passed to or from other methods.
#' @return A list summarizing of the results of the templateICA analysis.
#' @export
#' @method summary tICA.matrix
summary.tICA.matrix <- function(object, ...) {
  tICA_params <- lapply(
    object$params,
    function(q) {
      if (is.null(q)) { q <- "NULL"};
      paste0(as.character(q), collapse=" ")
    }
  )

  x <- c(
    list(nV=nrow(object$subjICmean), nL=ncol(object$subjICmean)),
    tICA_params
  )

  class(x) <- "summary.tICA.matrix"
  return(x)
}

#' @rdname summary.tICA.cifti
#' @export
#'
#' @param x The result of \code{templateICA} with CIFTI data
#' @param ... further arguments passed to or from other methods.
#' @return Nothing, invisibly.
#' @method print summary.tICA.cifti
print.summary.tICA.cifti <- function(x, ...) {
  cat("====TEMPLATE ICA INFO================\n")
  cat("Temporal Res.:   ", x$TR, "s.\n")
  cat("Highpass filter: ", x$hpf, "Hz\n")
  cat("Spatial scaling: ", x$scale, "\n")
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
#' @return Nothing, invisibly.
#' @method print summary.tICA.nifti
print.summary.tICA.nifti <- function(x, ...) {
  cat("====TEMPLATE INFO====================\n")
  cat("Temporal Res.:   ", x$TR, "s.\n")
  cat("Highpass filter: ", x$hpf, "Hz\n")
  cat("Spatial scaling: ", x$scale, "\n")
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

#' @rdname summary.tICA.matrix
#' @export
#'
#' @param x The template from \code{estimate_template.cifti}
#' @param ... further arguments passed to or from other methods.
#' @return Nothing, invisibly.
#' @method print summary.tICA.matrix
print.summary.tICA.matrix <- function(x, ...) {
  cat("====TEMPLATE INFO====================\n")
  cat("Temporal Res.:   ", x$TR, "s.\n")
  cat("Highpass filter: ", x$hpf, "Hz\n")
  cat("Spatial scaling: ", x$scale, "\n")
  cat("Variance method: ", x$tvar_method, "\n")
  cat("Q2 and Q2_max:   ", paste0(x$Q2, ", ", x$Q2_max), "\n")
  cat("-------------------------------------\n")
  cat("Spatial model:   ", x$spatial_model, "\n")
  cat("Dims reduced:    ", x$reduce_dim, "\n")
  cat("Maximum iters:   ", x$maxiter, "\n")
  cat("Epsilon:         ", x$epsilon, "\n")
  cat("-------------------------------------\n")
  cat("# Locations:     ", x$nV, "\n")
  cat("# Template ICs:  ", x$nL, "\n")
  cat("\n")

  invisible(NULL)
}

#' @rdname summary.tICA.cifti
#' @export
#'
#' @return Nothing, invisibly.
#' @method print tICA.cifti
print.tICA.cifti <- function(x, ...) {
  print.summary.tICA.cifti(summary(x))
}

#' @rdname summary.tICA.nifti
#' @export
#'
#' @return Nothing, invisibly.
#' @method print tICA.nifti
print.tICA.nifti <- function(x, ...) {
  print.summary.tICA.nifti(summary(x))
}

#' @rdname summary.tICA.matrix
#' @export
#'
#' @return Nothing, invisibly.
#' @method print tICA.matrix
print.tICA.matrix <- function(x, ...) {
  print.summary.tICA.matrix(summary(x))
}

#' Plot template
#'
#' @param x The result of \code{templateICA} with CIFTI data
#' @param stat \code{"mean"} (default), \code{"se"}, or \code{"both"}
#' @param ... Additional arguments to \code{view_xifti}
#' @return The plot
#' @export
#' @method plot tICA.cifti
plot.tICA.cifti <- function(x, stat=c("mean", "se", "both"), ...) {
  stopifnot(inherits(x, "tICA.cifti"))

  if (!requireNamespace("ciftiTools", quietly = TRUE)) {
    stop("Package \"ciftiTools\" needed to work with CIFTI data. Please install it.", call. = FALSE)
  }

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
      ciftiTools::view_xifti, c(list(x[[paste0("subjIC", ss)]]), args_ss)
    )
  }

  invisible(out)
}

#' Plot template
#'
#' @param x The result of \code{templateICA} with NIFTI data
#' @param stat \code{"mean"} (default), \code{"se"}
#' @param plane,n_slices,slices Anatomical plane and which slice indices to show.
#'  Default: 9 axial slices.
#' @param ... Additional arguments
#' @return The plot
#' @export
#' @method plot tICA.nifti
plot.tICA.nifti <- function(x, stat=c("mean", "se"),
  plane=c("axial", "sagittal", "coronal"), n_slices=9, slices=NULL,
  ...) {
  stopifnot(inherits(x, "tICA.nifti"))

  if (!requireNamespace("oro.nifti", quietly = TRUE)) {
    stop("Package \"oro.nifti\" needed to read NIFTI data. Please install it.", call. = FALSE)
  }

  # Check `...`
  args <- list(...)
  has_title <- "title" %in% names(args)
  has_idx <- "idx" %in% names(args)
  has_fname <- "fname" %in% names(args)

  # Check `idx`
  if (has_idx) {
    stopifnot(length(args$idx)==1)
    stopifnot(is.numeric(args$idx) && args$idx==round(args$idx))
    stopifnot(args$idx %in% seq(ncol(x$subjICmean)))
  } else {
    args$idx <- 1
  }
  idx <- args$idx; args$idx <- NULL

  # Check `stat`
  stat <- tolower(stat)
  if (has_idx && length(args$idx)>1 && !("fname" %in% names(args))) {
    if (identical(stat, c("mean", "se"))) {
      stat <- "mean"
    } else {
      stat <- match.arg(stat, c("mean", "se"))
    }
  }
  stat <- match.arg(stat, c("mean", "se"))

  # Print message saying what's happening.
  msg1 <- ifelse(has_idx,
    "Plotting the",
    "Plotting the first component's"
  )
  msg2 <- switch(stat,
    mean="estimate.",
    se="standard error."
  )
  cat(msg1, msg2, "\n")

  # Plot
  out <- list(mean=NULL, se=NULL)
  ss <- stat

  plane <- match.arg(plane, c("axial", "sagittal", "coronal"))
  args$plane <- plane
  plane_dim <- switch(plane, axial=3, coronal=2, sagittal=1)
  if (is.null(slices)) {
    if (is.null(n_slices)) { warning("Using 9 slices."); n_slices <- 9 }
    n_slices <- as.numeric(n_slices)
    if (length(n_slices) > 1) { warning("Using the first entry of `slice`."); n_slices <- n_slices[1] }
    # Pick slices that are spaced out, and with many voxels.
    mask_count <- apply(x$mask_nii, plane_dim, sum)
    ns_all <- length(mask_count)
    slices <- seq(ns_all)
    # Remove slices with zero voxels.
    slices <- slices[mask_count != 0]
    mask_count <- mask_count[mask_count != 0]
    ns_all <- length(mask_count)
    if (n_slices > length(slices)) {
      warning(
        "`n_slices` is larger than the number of non-empty slices (",
        length(slices), "). Showing all non-empty slices."
      )
      n_slices <- length(slices)
    }
    # Remove slices with few voxels.
    if (n_slices < (ns_all / 2)) {
      slices <- slices[mask_count > quantile(mask_count, .33)]
    }
    slices <- slices[round(seq(1, length(slices), length.out=n_slices))]
  } else {
    slices <- as.numeric(slices)
    stopifnot(all(slices %in% seq(dim(x$mask_nii)[plane_dim])))
  }

  tss <- x[[paste0("subjIC", ss)]]
  tss <- tss[,,,idx]

  if (plane=="axial") {
    tss <- tss[,,slices,drop=FALSE]
  } else if (plane=="coronal") {
    tss <- tss[,slices,,drop=FALSE]
  } else if (plane=="sagittal") {
    tss <- tss[slices,,,drop=FALSE]
  } else { stop() }

  args_ss <- args
  args_ss$plane <- plane
  # Handle title and idx
  if (!has_title && !has_idx) {
    c1name <- "First component"
  }
  if (has_title) { stop("Not supported yet.") }
  if (has_fname) { stop("Not supported yet. Call `pdf` or `png` beforehand, and then `dev.off`.") }
  do.call(
    oro.nifti::image,
    c(list(oro.nifti::as.nifti(tss)), args_ss)
  )
}

#' Plot template
#' 
#' This feature is not supported yet.
#'
#' @param x The result of \code{templateICA} with NIFTI data
#' @param ... Additional arguments
#' @return Nothing, because an error is raised.
#' @export
#' @method plot tICA.matrix
plot.tICA.matrix <- function(x, ...) {
  stopifnot(inherits(x, "tICA.matrix"))
  stop("Not supported yet.")
}
