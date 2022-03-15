#' Summarize a \code{"template.cifti"} object
#'
#' Summary method for class \code{"template.cifti"}
#'
#' @param object Object of class \code{"template.cifti"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method summary template.cifti
summary.template.cifti <- function(object, ...) {
  tmean <- struct_template(object$template$mean, "CIFTI", object$dat_struct, object$params)
  x <- c(
    summary(tmean),
    list(has_DR="DR" %in% names(object)),
    object$params
  )

  class(x) <- "summary.template.cifti"
  return(x)
}

#' Summarize a \code{"template.nifti"} object
#'
#' Summary method for class \code{"template.nifti"}
#'
#' @param object Object of class \code{"template.nifti"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method summary template.nifti
summary.template.nifti <- function(object, ...) {
  x <- c(
    list(
      mask_dims=dim(object$dat_struct),
      nV=nrow(object$template$mean), 
      nL=ncol(object$template$mean),
      hasDR="DR" %in% names(object)
    ),
    object$params
  )

  class(x) <- "summary.template.nifti"
  return(x)
}

#' Summarize a \code{"template.data"} object
#'
#' Summary method for class \code{"template.data"}
#'
#' @param object Object of class \code{"template.data"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method summary template.data
summary.template.data <- function(object, ...) {
  x <- c(
    list(
      nV=nrow(object$template$mean), 
      nL=ncol(object$template$mean),
      hasDR="DR" %in% names(object)
    ),
    object$params
  )

  class(x) <- "summary.template.data"
  return(x)
}

#' @rdname summary.template.cifti
#' @export
#'
#' @param x The template from \code{estimate_template.cifti}
#' @param ... further arguments passed to or from other methods.
#' @method print summary.template.cifti
print.summary.template.cifti <- function(x, ...) {
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
  cat("# Subjects:      ", x$num_subjects, "\n")
  cat("Detrending:      ", dct, "\n")
  cat("Spatial scaling: ", x$scale, "\n")
  cat("A normalization: ", x$normA, "\n")
  cat("Q2 and Q2_max:   ", paste0(x$Q2, ", ", x$Q2_max), "\n")
  cat("Pseudo retest:   ", x$pseudo_retest, "\n")
  cat("\n")

  class(x) <- "summary.xifti"
  print(x)
  invisible(NULL)
}

#' @rdname summary.template.nifti
#' @export
#'
#' @param x The template from \code{estimate_template.nifti}
#' @param ... further arguments passed to or from other methods.
#' @method print summary.template.nifti
print.summary.template.nifti <- function(x, ...) {
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
  cat("# Subjects:      ", x$num_subjects, "\n")
  cat("Detrending:      ", dct, "\n")
  cat("Spatial scaling: ", x$scale, "\n")
  cat("A normalization: ", x$normA, "\n")
  cat("Q2 and Q2_max:   ", paste0(x$Q2, ", ", x$Q2_max), "\n")
  cat("Pseudo retest:   ", x$pseudo_retest, "\n")
  cat("-------------------------------------\n")
  cat("Mask dims:       ", paste0(x$mask_dims, collapse=" x "), "\n")
  cat("Vectorized dims:\n")
  cat("# Locations:     ", x$nV, "\n")
  cat("# Template ICs:  ", x$nL, "\n")
  cat("\n")

  invisible(NULL)
}

#' @rdname summary.template.data
#' @export
#'
#' @param x The template from \code{estimate_template.cifti}
#' @param ... further arguments passed to or from other methods.
#' @method print summary.template.data
print.summary.template.data <- function(x, ...) {
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
  cat("# Subjects:      ", x$num_subjects, "\n")
  cat("Detrending:      ", dct, "\n")
  cat("Spatial scaling: ", x$scale, "\n")
  cat("A normalization: ", x$normA, "\n")
  cat("Q2 and Q2_max:   ", paste0(x$Q2, ", ", x$Q2_max), "\n")
  cat("Pseudo retest:   ", x$pseudo_retest, "\n")
  cat("-------------------------------------\n")
  cat("Dimensions:      \n")
  cat("# Locations:     ", x$nV, "\n")
  cat("# Template ICs:  ", x$nL, "\n")
  cat("\n")

  invisible(NULL)
}

#' @rdname summary.template.cifti
#' @export
#'
#' @method print template.cifti
print.template.cifti <- function(x, ...) {
  print.summary.template.cifti(summary(x))
}

#' @rdname summary.template.nifti
#' @export
#'
#' @method print template.nifti
print.template.nifti <- function(x, ...) {
  print.summary.template.nifti(summary(x))
}

#' @rdname summary.template.data
#' @export
#'
#' @method print template.data
print.template.data <- function(x, ...) {
  print.summary.template.data(summary(x))
}

#' Plot template
#'
#' @param x The template from \code{estimate_template.cifti}
#' @param stat \code{"mean"}, \code{"var"}, or \code{"both"} (default)
#' @param var_method \code{"non-negative"} (default) or \code{"unbiased"}
#' @param ... Additional arguments to \code{view_xifti}
#' @return The plot
#' @export
#' @importFrom ciftiTools view_xifti
#' @method plot template.cifti
plot.template.cifti <- function(x, stat=c("both", "mean", "var"), 
  var_method=c("non-negative", "unbiased"), ...) {
  stopifnot(inherits(x, "template.cifti"))

  var_method <- match.arg(var_method, c("non-negative", "unbiased"))

  # Check `...`
  args <- list(...)
  has_title <- "title" %in% names(args)
  has_idx <- "idx" %in% names(args)
  has_fname <- "fname" %in% names(args)

  # Check `stat`
  stat <- tolower(stat)
  if (has_idx && length(args$idx)>1 && !("fname" %in% names(args))) {
    if (identical(stat, c("both", "mean", "var"))) {
      stat <- "mean"
    } else {
      stat <- match.arg(stat, c("both", "mean", "var"))
    }
    if (stat == "both") {
      if (!("fname" %in% names(args))) {
        warning(
          "For multiple `idx`, use one call to plot() ",
          "for the mean template, ",
          "and a separate one for the variance template. ",
          "Showing the mean template now."
        )
        stat <- "mean"
      }
    }
  }
  stat <- match.arg(stat, c("both", "mean", "var"))

  # Print message saying what's happening.
  msg1 <- ifelse(has_idx,
    "Plotting the",
    "Plotting the first component's"
  )
  msg2 <- switch(stat,
    both="templates.",
    mean="mean template.",
    var="variance template."
  )
  cat(msg1, msg2, "\n")

  # Plot
  out <- list(mean=NULL, var=NULL)
  if (stat == "both") { stat <- c("mean", "var") }
  for (ss in stat) {
    ssname <- if (ss == "mean") {
      ss
    } else if (var_method=="non-negative") { 
      "varNN"
    } else {
      "varUB"
    }
    if (ss=="var" && var_method=="unbiased") { x$template[[ssname]][] <- pmax(0, x$template[[ssname]]) }
    tss <- struct_template(x$template[[ssname]], "CIFTI", x$dat_struct, x$params)
    args_ss <- args
    # Handle title and idx
    if (!has_title && !has_idx) {
      c1name <- if (!is.null(tss$meta$cifti$names)) {
        tss$meta$cifti$names[1]
      } else {
        "First component"
      }
      args_ss$title <- paste0(c1name, " (", ssname, ")")
    } else if (!has_idx) {
      args_ss$title <- paste0(args_ss$title, "(", ssname, ")")
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
      view_xifti, c(list(tss), args_ss)
    )
  }

  invisible(out)
}

#' Plot template
#'
#' @param x The template from \code{estimate_template.nifti}
#' @param ... Additional arguments
#' @return The plot
#' @export
#' @method plot template.nifti
plot.template.nifti <- function(x, ...) {
  stop("Not supported yet.")
}

#' Plot template
#'
#' @param x The template from \code{estimate_template.data}
#' @param ... Additional arguments
#' @return The plot
#' @export
#' @method plot template.data
plot.template.data <- function(x, ...) {
  stop("Not supported yet.")
}