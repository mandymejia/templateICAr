#' Dual Regression
#'
#' @param dat Subject-level fMRI data (\eqn{VxT})
#' @param GICA Group-level independent components (\eqn{VxQ})
#' @param scale A logical value indicating whether the fMRI timeseries should
#'  be scaled by the image standard deviation.
#'
#' @importFrom matrixStats colVars
#' 
#' @return A list containing the subject-level independent components \strong{S} (\eqn{VxQ}), 
#'  and subject-level mixing matrix \strong{A} (\eqn{TxQ}).
#' 
#' @export
#'
dual_reg <- function(dat, GICA, scale=FALSE){

  nT <- ncol(dat) #length of timeseries
  nvox <- nrow(dat) #number of data locations
  if(nT > nvox) warning('More time points than voxels. Are you sure?')
  if(nvox != nrow(GICA)) stop('The number of voxels in dat and GICA must match')

  Q <- ncol(GICA) #number of ICs
  if(Q > nvox) warning('More ICs than voxels. Are you sure?')
  if(Q > nT) warning('More ICs than time points. Are you sure?')

  # center timeseries data across space and time (and standardize scale if scale=TRUE)
  # transpose it
  dat_ctr <- t(scale_BOLD(dat, scale=scale))

  #center each group IC over voxels
  GICA - rep(colMeans(GICA), rep.int(nvox, Q))

	#estimate A (IC timeseries)
	A <- (dat_ctr %*% GICA) %*% chol2inv(chol(crossprod(GICA)))
	#estimate S (IC maps)
	S <- solve(a=crossprod(A), b=crossprod(A, dat_ctr))

	#fix scale of spatial maps (sd=1)
	#sd_S <- sqrt(rowVars(S))
	#A <- A %*% diag(sd_S)
	#S <- diag(1/sd_S) %*% S

	#return result
	list(S = S, A = A)
}

#' Dual Regression wrapper
#' 
#' Wrapper to \code{dual_reg}. Handles NIFTI or CIFTI file input, halves the 
#' data in absense of retest data, and may remove nuisance regressors.
#' 
#' @param BOLD,BOLD2 subject-level BOLD data: a CIFTI file, NIFTI file, or data matrix.
#'  Cannot be a \code{"xifti"} object, \code{"nifti"} object, etc. If \code{BOLD2} is provided 
#'  it must be in the same format as \code{BOLD}; \code{BOLD} will be the test data and
#'  \code{BOLD2} will be the retest data. If \code{BOLD2} is not provided,
#'  \code{BOLD} will be split in half; the first half will be the test data and
#'  the second half will be the retest data.
#' 
#'  NIFTI implementation is not complete yet (will return error).
#' @param scale Logical indicating whether BOLD data should be scaled by the
#'  spatial standard deviation before template estimation.
#' @param GICA Group ICA data matrix. Cannot be a \code{"xifti"} object, \code{"nifti"} object, file, etc.
#' @param format Expected format of \code{BOLD} and \code{BOLD2}. Should be one of the following:
#'  \code{"infer"} (default), a \code{"data"} matrix, \code{"CIFTI"} file, or \code{"NIFTI"} file.
#' @param dim_expect If \code{!(format == "NIFTI")}, this is the number of rows to expect in the data.
#'  Otherwise, ... 
#' @param Q2 The number of nuisance ICs to identify. If \code{NULL}, will be estimated.
#'  Only provide \code{Q2} or \code{maxQ} but not both.
#' @param maxQ Maximum number of ICs (template+nuisance) to identify 
#'  (L <= maxQ <= T). If \code{maxQ == L}, then do not remove any nuisance regressors.
#'  Only provide \code{Q2} or \code{maxQ} but not both.
#' @param brainstructures Only applies if \code{format} is \code{"CIFTI"}. 
#'  Character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
#'  cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Can also be \code{"all"} (obtain all three brain structures).
#'  Default: \code{c("left","right")} (cortical surface only).
#' @param verbose If \code{TRUE}, display progress updates
#' 
#' @keywords internal
dual_reg2 <- function(
  BOLD, BOLD2=NULL, GICA, scale=FALSE, format=c("infer", "data", "CIFTI", "NIFTI"), 
  Q2=NULL, maxQ=NULL, dim_expect=NULL, brainstructures=c("left", "right"), verbose=FALSE){

  # TO DO: add removal of nuisance ICs.
  # If !retest, must be careful about the number of ICs since nT is halved.

  # Initialize the return value
  out <- list(
    missing = NULL,
    test = NULL,
    retest = NULL
  )

  # Retest or pseudo retest?
  retest <- !is.null(BOLD2)

  # Determine the format of `BOLD` and `BOLD2`
  format <- match.arg(format, c("infer", "data", "CIFTI", "NIFTI"))
  if (format == "infer") {
    if (is.character(BOLD)) {
      stopifnot(length(BOLD) == 1)
      format <- ifelse(
        endsWith(BOLD, ".dtseries.nii") | endsWith(BOLD, ".dscalar.nii"),
        "CIFTI", "NIFTI"
      )
    } else {
      format <- "data"
    }
    if (verbose) { cat("Inferred input format:", format, "\n") }
  }

  if (format == "NIFTI") { stop("Not completed yet.") }

  # If BOLD (and BOLD2) is a CIFTI or NIFTI, check that it exists.
  if (format %in% c("CIFTI", "NIFTI")) {
    if (retest) { 
      stopifnot(is.character(BOLD2))
      stopifnot(length(BOLD2) == 1)
    }
    BOLD_exists <- file.exists(c(BOLD, BOLD2))
    if (!any(BOLD_exists)) {
      out$missing <- c(BOLD, BOLD2)[BOLD_exists]
      return(out)
    }
  }

  # If BOLD (and BOLD2) is a CIFTI or NIFTI, read it in.
  #   CIFTI
  if (format == "CIFTI") {
    # [TO DO]: check for ciftiTools
    BOLD <- as.matrix(read_cifti(BOLD, brainstructures=brainstructures))
    if (retest) {
      BOLD2 <- as.matrix(read_cifti(BOLD2, brainstructures=brainstructures))
    }
  #   NIFTI
  } else if (format == "NIFTI") {
    # [TO DO]: check for RNifti
    BOLD <- readNIfTI(BOLD, reorient=FALSE)
    if (retest) {
      BOLD2 <- readNIfTI(BOLD2, reorient=FALSE)
    }
  }

  # Check dim
  dBOLD <- dim(BOLD)
  if (!is.null(dim_expect)) {
    dim_expect <- as.numeric(dim_expect)
    if (format == "NIFTI") {
      stop()
    } else {
      if (length(dim_expect) > 1) {
        warning("`dim_expect` should only give the number of expected rows. Using the first entry only.\n")
        dim_expect <- dim_expect[1]
      }
      if (dBOLD[1] != dim_expect) {
        stop("nrow(BOLD) is ", dBOLD[1], ", but `dim_expect` is ", dim_expect, ".")
      }
    }
  }

  # If no retest data, halve the test data
  if (!retest) {
    nT <- dBOLD[length(dBOLD)]
    part1 <- seq(round(nT/2))
    part2 <- setdiff(seq(nT), part1)
    if (format=="CIFTI") {
      BOLD2 <- BOLD[, part2, drop=FALSE]
      BOLD <- BOLD[, part1, drop=FALSE]
    } else if (format=="NIFTI") {
      BOLD2 <- BOLD[,,, part2, drop=FALSE]
      BOLD <- BOLD[,,, part1, drop=FALSE]
    }
  }

  # [TO DO]: Remove NA values?

  # Perform dual regression on test and retest data
  out$test <- dual_reg(BOLD, as.matrix(GICA), scale=scale)
  out$retest <- dual_reg(BOLD2, as.matrix(GICA), scale=scale)
  out$missing <- FALSE

  nT <- min(dim(BOLD)[length(dim(BOLD))], dim(BOLD2)[length(dim(BOLD2))])
  skip_rm_nuisIC <- (!is.null(Q2) && Q2 <= 0) || (!is.null(maxQ) && maxQ <=0)
  if (!skip_rm_nuisIC) {

    L <- ncol(GICA)

    # Check `maxQ`
    if (!is.null(maxQ)) { 
      if(round(maxQ) != maxQ || maxQ <= 0) stop('maxQ must be NULL or a round positive number.')
    } else {
      maxQ <- round(nT/2)
    }
    if (maxQ < L) {
      warning('maxQ must be at least L.  Setting maxQ=L.')
      maxQ <- L
    }
    # This is to avoid the area of the pesel objective function that spikes close 
    #   to rank(X), which often leads to nPC close to rank(X)
    if (maxQ > nT*0.75) {
      warning('maxQ too high, setting to 75% of T.')
      maxQ <- round(nT*0.75)
    }

    # Scale
    BOLD <- scale_BOLD(BOLD, scale=scale)
    BOLD2 <- scale_BOLD(BOLD2, scale=scale)

    # Estimate and deal with nuisance ICs
    if (maxQ > L) {
      BOLD <- rm_nuisIC(BOLD, DR=out$test, Q2=Q2, Q2_max=maxQ-L, verbose=verbose)
      BOLD2 <- rm_nuisIC(BOLD2, DR=out$retest, Q2=Q2, Q2_max=maxQ-L, verbose=verbose)
    }

    # Do DR again.
    out$test_preclean <- out$test$S # [TO DO]: Rename?
    out$test <- dual_reg(BOLD, as.matrix(GICA), scale=scale)
    out$retest_preclean <- out$retest$S
    out$retest <- dual_reg(BOLD2, as.matrix(GICA), scale=scale)
  }

  out$test <- out$test$S
  out$retest <- out$retest$S
  out
}