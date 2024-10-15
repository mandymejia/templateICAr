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

  nC <- length(object) - 2 # list items beside "active" and "params"
  act_counts_lenient <- colSums(as.matrix(object$active)>0, na.rm=TRUE)
  act_counts_strict <- colSums(as.matrix(object$active)==nC, na.rm=TRUE)
  verts_with_data_per_bs <- vapply(
    object$active$data[!vapply(object$active$data, is.null, FALSE)],
    function(q){sum(q[,1]>-1, na.rm=TRUE)},
    0
  )
  x <- c(
    summary(object$active),
    list(nC=nC),
    list(activation_name=do.call(
      format_activation_name,
      c(object$params[c("u", "z", "type", "deviation")], list(collapse=TRUE))
    )),
    list(verts_with_data_per_bs=verts_with_data_per_bs),
    list(act_counts_lenient=act_counts_lenient),
    list(act_counts_strict=act_counts_strict),
    object$params
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

  apct_lenient <- round(x$act_counts_lenient/sum(x$verts_with_data_per_bs)*100)
  apct_strict <- round(x$act_counts_strict/sum(x$verts_with_data_per_bs)*100)
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

  cname <- ifelse(!is.null(x$z), "`c`", "`u`")

  nMeasShow <- min(5, x$measurements)
  nMeasTriC <- ifelse(x$measurements > 5, ", ...", "")

  cat("====ACTIVATIONS STATS================\n")
  cat("alpha:           ", x$alpha, "\n")
  cat("p-val method:    ", pm_nice, "\n")
  cat("Test:            ", x$activation_name, "\n")
  cat(
    "Active Loc. (%): ",
    paste0(paste(apct_lenient[seq(nMeasShow)], collapse=", "), nMeasTriC),
    ifelse(x$nC > 1, paste0("(most lenient ", cname, ")"), ""),
    "\n"
  )
  if (!all(x$act_counts_lenient == x$act_counts_strict)) {
    cat(
      "Active Loc. (%): ",
      paste0(paste(apct_strict[seq(nMeasShow)], collapse=", "), nMeasTriC),
      ifelse(x$nC > 1, paste0("(most strict ", cname, ")"), ""),
      "\n"
    )
  }
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

#' Summarize a \code{"tICA_act.nifti"} object
#'
#' Summary method for class \code{"tICA_act.nifti"}
#'
#' @param object Object of class \code{"tICA_act.nifti"}.
#' @param ... further arguments passed to or from other methods.
#' @return A list summarizing the data and results for the activations analysis.
#' @export
#' @method summary tICA_act.nifti
summary.tICA_act.nifti <- function(object, ...) {
  act_counts <- colSums(object$active, na.rm=TRUE)
  x <- c(
    summary(object$active),
    list(nV=nrow(object$active), nL=ncol(object$active)),
    list(act_counts=act_counts),
    object$params
  )

  class(x) <- "summary.tICA_act.nifti"
  return(x)
}

#' @rdname summary.tICA_act.nifti
#' @export
#'
#' @param x The activations from \code{activations}
#' @param ... further arguments passed to or from other methods.
#' @return Nothing, invisibly.
#' @method print summary.tICA_act.nifti
print.summary.tICA_act.nifti <- function(x, ...) {

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

  usign <- if (all(x$u>=0)) {
    "+"
  } else if (all(x$u<=0)) {
    "-"
  } else {
    "+/-" # not the best, but this shouldn't even happen :)
  }
  ustr <- if (length(x$z)==1) {
    ifelse(x$deviation, paste0(abs(x$z), "*z"), paste0(x$z, "*z"))
  } else if (length(x$z)>1) {
    "z"
  } else if (length(x$u)==1) {
    ifelse(x$deviation, abs(x$u), x$u)
  } else {
    "u"
  }
  adesc <- if (x$deviation) {
    if (any(x$u!=0)) {
      paste("x", x$type, "mu", usign, ustr)
    } else {
      paste("x", x$type, "mu")
    }
  } else {
    if (any(x$u!=0)) {
      paste("x", x$type, ustr)
    } else {
      paste("x", x$type, "0")
    }
  }

  nMeasShow <- min(5, x$measurements)
  nMeasTriC <- ifelse(x$measurements > nMeasShow, ", ...", "")

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
  cat("# Locations:     ", x$nV, "\n")
  cat("# Template ICs:  ", x$nL, "\n")
  cat("\n")
  invisible(NULL)
}

#' @rdname summary.tICA_act.nifti
#' @export
#'
#' @return Nothing, invisibly.
#' @method print tICA_act.nifti
print.tICA_act.nifti <- function(x, ...) {
  print.summary.tICA_act.nifti(summary(x))
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
    object$params
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

  usign <- if (all(x$u>=0)) {
    "+"
  } else if (all(x$u<=0)) {
    "-"
  } else {
    "+/-" # not the best, but this shouldn't even happen :)
  }
  ustr <- if (length(x$z)==1) {
    ifelse(x$deviation, paste0(abs(x$z), "*z"), paste0(x$z, "*z"))
  } else if (length(x$z)>1) {
    "z"
  } else if (length(x$u)==1) {
    ifelse(x$deviation, abs(x$u), x$u)
  } else {
    "u"
  }
  adesc <- if (x$deviation) {
    if (any(x$u!=0)) {
      paste("x", x$type, "mu", usign, ustr)
    } else {
      paste("x", x$type, "mu")
    }
  } else {
    if (any(x$u!=0)) {
      paste("x", x$type, ustr)
    } else {
      paste("x", x$type, "0")
    }
  }

  nMeasShow <- min(5, x$measurements)
  nMeasTriC <- ifelse(x$measurements > nMeasShow, ", ...", "")

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
  cat("# Locations:     ", x$nV, "\n")
  cat("# Template ICs:  ", x$nL, "\n")
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
