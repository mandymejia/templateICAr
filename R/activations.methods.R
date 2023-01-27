#' Summarize a \code{"tICA_act.cifti"} object
#'
#' Summary method for class \code{"tICA_act.cifti"}
#'
#' @param object Object of class \code{"tICA_act.cifti"}.
#' @param ... further arguments passed to or from other methods.
#' @return A list summarizing the data and results for the activations analysis.
#' @export
#' @method summary tICA_act.cifti
summary.tICA_act.cifti <- function(object, ...) {
  act_counts <- colSums(as.matrix(object$active)==2, na.rm=TRUE)
  x <- c(
    summary(object$active),
    list(act_counts=act_counts),
    object[c("u", "alpha", "type", "method_p", "deviation")]
  )

  class(x) <- "summary.tICA_act.cifti"
  return(x)
}

#' @rdname summary.tICA_act.cifti
#' @export
#'
#' @param x The activations from \code{activations.cifti}
#' @param ... further arguments passed to or from other methods.
#' @return Nothing, invisibly.
#' @method print summary.tICA_act.cifti
print.summary.tICA_act.cifti <- function(x, ...) {

  #mapct <- paste0(" (", round(mean(x$act_counts)/sum(x$verts_per_bs)*100), "% of locations)")
  apct <- round(x$act_counts/sum(x$verts_per_bs)*100)
  pm_nice <- switch(x$method_p,
    bonferroni = "Bonferroni",
    holm = "Holm",
    hochberg = "Hochberg",
    hommel = "Hommel",
    BH = "Benjamini & Hochberg (FDR)",
    BY = "Benjamini & Yekutieli",
    fdr = "Benjamini & Hochberg (FDR)",
    none = "none"
  )

  usign <- ifelse(x$u >=0, "+", "-")
  adesc <- if (x$deviation) {
    if (x$u != 0) {
      paste("x", x$type, "mu", usign, abs(x$u))
    } else {
      paste("x", x$type, "mu")
    }
  } else {
    if (x$u != 0) {
      paste("x", x$type, x$u)
    } else {
      paste("x", x$type, "0")
    }
  }

  nMeasShow <- min(5, x$measurements)
  nMeasTriC <- ifelse(nMeasShow > 5, ", ...", "")

  cat("====ACTIVATIONS STATS================\n")
  cat("alpha:           ", x$alpha, "\n")
  cat("p-val method:    ", pm_nice, "\n")
  cat("Test:            ", adesc, "\n")
  # cat("Type:            ", x$type, "\n")
  # cat("Threshold:       ", x$u, "\n")
  # cat("Deviation:       ", x$deviation, "\n")
  cat(
    "Active Loc. (%): ", 
    paste0(paste(apct[seq(nMeasShow)], collapse=", "), nMeasTriC), "\n"
  )
  cat("\n")

  class(x) <- "summary.xifti"
  print(x)
  invisible(NULL)
}

#' @rdname summary.tICA_act.cifti
#' @export
#'
#' @return Nothing, invisibly.
#' @method print tICA_act.cifti
print.tICA_act.cifti <- function(x, ...) {
  print.summary.tICA_act.cifti(summary(x))
}

#' Plot activations
#'
#' @param x The activations from \code{activations.cifti}
#' @param stat \code{"active"} (default), \code{"pvals"}, \code{"pvals_adj"},
#'  \code{"tstats"}, or \code{"vars"}.
#' @param ... Additional arguments to \code{view_xifti}
#' @return The activations plot
#' @export
#' @method plot tICA_act.cifti
plot.tICA_act.cifti <- function(x, stat=c("active", "pvals", "pvals_adj", "tstats", "se"), ...) {
  stopifnot(inherits(x, "tICA_act.cifti"))

  if (!requireNamespace("ciftiTools", quietly = TRUE)) {
    stop("Package \"ciftiTools\" needed to read NIFTI data. Please install it.", call. = FALSE)
  }

  # Check `...`
  args <- list(...)
  has_title <- "title" %in% names(args)
  has_idx <- "idx" %in% names(args)
  has_fname <- "fname" %in% names(args)

  stat <- match.arg(stat, c("active", "pvals", "pvals_adj", "tstats", "se"))

  # Print message saying what's happening.
  msg1 <- ifelse(has_idx,
    "Plotting the",
    "Plotting the first component's"
  )
  msg2 <- switch(stat,
    active="activation maps.",
    pvals="p values.",
    pvals_adj="adjusted p values.",
    tstats="t statistics.",
    se="standard errors."
  )
  cat(msg1, msg2, "\n")

  if (stat == "active") {
    x <- x$active
  } else {
    x <- ciftiTools::newdata_xifti(x$se, as.matrix(x[[stat]]))
  }

  ss <- stat # to match `plot.template.cifti`
  args_ss <- args
  # Handle title and idx
  if (!has_title && !has_idx) {
    c1name <- if (!is.null(x$meta$cifti$names)) {
      x$meta$cifti$names[1]
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
  do.call(ciftiTools::view_xifti, c(list(x), args_ss))
}

#' Summarize a \code{"tICA_act.matrix"} object
#'
#' Summary method for class \code{"tICA_act.matrix"}
#'
#' @param object Object of class \code{"tICA_act.matrix"}.
#' @param ... further arguments passed to or from other methods.
#' @return A list summarizing the data and results for the activations analysis.
#' @export
#' @method summary tICA_act.matrix
summary.tICA_act.matrix <- function(object, ...) {
  act_counts <- colSums(as.matrix(object$active), na.rm=TRUE)
  x <- c(
    list(nV=nrow(object$active), nL=ncol(object$active)),
    list(act_counts=act_counts),
    object[c("u", "alpha", "type", "method_p", "deviation")]
  )

  class(x) <- "summary.tICA_act.matrix"
  return(x)
}

#' @rdname summary.tICA_act.matrix
#' @export
#'
#' @param x The activations from \code{activations}
#' @param ... further arguments passed to or from other methods.
#' @return Nothing, invisibly.
#' @method print summary.tICA_act.matrix
print.summary.tICA_act.matrix <- function(x, ...) {

  #mapct <- paste0(" (", round(mean(x$act_counts)/x$nV*100), "% of locations)")
  apct <- round(x$act_counts/x$nV*100)
  pm_nice <- switch(x$method_p,
    bonferroni = "Bonferroni",
    holm = "Holm",
    hochberg = "Hochberg",
    hommel = "Hommel",
    BH = "Benjamini & Hochberg (FDR)",
    BY = "Benjamini & Yekutieli",
    fdr = "Benjamini & Hochberg (FDR)",
    none = "none"
  )

  usign <- ifelse(x$u >=0, "+", "-")
  adesc <- if (x$deviation) {
    if (x$u != 0) {
      paste("x", x$type, "mu", usign, abs(x$u))
    } else {
      paste("x", x$type, "mu")
    }
  } else {
    if (x$u != 0) {
      paste("x", x$type, x$u)
    } else {
      paste("x", x$type, "0")
    }
  }

  nMeasShow <- min(5, x$measurements)
  nMeasTriC <- ifelse(nMeasShow > 5, ", ...", "")

  cat("====ACTIVATIONS STATS================\n")
  cat("alpha:           ", x$alpha, "\n")
  cat("p-val method:    ", pm_nice, "\n")
  cat("Test:            ", adesc, "\n")
  # cat("Type:            ", x$type, "\n")
  # cat("Threshold:       ", x$u, "\n")
  # cat("Deviation:       ", x$deviation, "\n")
  cat(
    "Active Loc. (%): ", 
    paste0(paste(apct[seq(nMeasShow)], collapse=", "), nMeasTriC), "\n"
  )
  cat("-------------------------------------\n")
  cat("# Locations:     ", x$nL, "\n")
  cat("# Template ICs:  ", x$nV, "\n")
  cat("\n")
  invisible(NULL)
}

#' @rdname summary.tICA_act.matrix
#' @export
#'
#' @return Nothing, invisibly.
#' @method print tICA_act.matrix
print.tICA_act.matrix <- function(x, ...) {
  print.summary.tICA_act.matrix(summary(x))
}
