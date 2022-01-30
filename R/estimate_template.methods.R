#' Summarize a \code{"template_cifti"} object
#'
#' Summary method for class \code{"template_cifti"}
#'
#' @param object Object of class \code{"template_cifti"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method summary template_cifti
summary.template_cifti <- function(object, ...) {

  x <- c(
    summary(object$template_mean),
    list(has_DR="DR" %in% names(object)),
    object$params
  )

  class(x) <- "summary.template_cifti"
  return(x)
}

#' @rdname summary.template_cifti
#' @export
#'
#' @param x The template from \code{estimate_template.cifti}
#' @param ... further arguments passed to or from other methods.
#' @method print summary.template_cifti
print.summary.template_cifti <- function(x, ...) {
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
  cat("Pseudo retest:   ", x$pseudo_retest, "\n")
  cat("\n")

  class(x) <- "summary.xifti"
  print(x)
}

#' @rdname summary.template_cifti
#' @export
#'
#' @method print template_cifti
print.template_cifti <- function(x, ...) {
  print.summary.template_cifti(summary(x))
}

#' Plot template
#'
#' @param x The template from \code{estimate_template.cifti}
#' @param stat \code{"mean"} (default) or \code{"var"}
#' @param ... Additional arguments to \code{view_xifti}
#' @return The plot
#' @export
#' @importFrom ciftiTools view_xifti
#' @method plot template_cifti
plot.template_cifti <- function(x, stat=c("both", "mean", "var"), ...) {
  stopifnot(inherits(x, "template_cifti"))

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
    args_ss <- args
    # Handle title and idx
    if (!has_title && !has_idx) {
      c1name <- if (!is.null(x$template_mean$meta$cifti$names)) {
        x$template_mean$meta$cifti$names[1]
      } else {
        "First component"
      }
      args_ss$title <- paste0(c1name, " (", ss, ")")
    } else if (!has_idx) {
      args_ss$title <- paste0(args_ss$title, "(", ss, ")")
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
      view_xifti, c(list(x[[paste0("template_", ss)]]), args_ss)
    )
  }

  invisible(out)
}

