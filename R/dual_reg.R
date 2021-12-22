#' Dual Regression
#'
#' @param BOLD Subject-level fMRI data matrix (\eqn{V \times T})
#' @param GICA Group-level independent components (\eqn{V \times Q})
#' @param center_rows,center_cols Center BOLD data across rows (each data location's time series) or columns (each time point's image)? Default: \code{TRUE} for both.
#' @param center_Gcols Center GICA across columns (each ICA)? Default: \code{TRUE}.
#' @param scale A logical value indicating whether the fMRI timeseries should be scaled by the image standard deviation.
#' @param detrend_DCT Detrend the data? This is the number of DCT bases to use for detrending. If \code{0} (default), do not detrend.
#' @param normA Scale each IC timeseries (column of \eqn{A}) in the dual regression 
#'  estimates? Default: \code{FALSE}. (The opposite scaling will be applied to \eqn{S}
#'  such that the product \eqn{A \times S} remains the same).
#'
#' @importFrom matrixStats colVars
#' 
#' @return A list containing the subject-level independent components \strong{S} (\eqn{V \times Q}), 
#'  and subject-level mixing matrix \strong{A} (\eqn{TxQ}).
#' 
#' @export
#'
dual_reg <- function(
  BOLD, GICA, 
  center_rows=TRUE, center_cols=TRUE, scale=FALSE, detrend_DCT=0, 
  center_Gcols=TRUE, normA=FALSE){

  stopifnot(is.matrix(BOLD))
  stopifnot(is.matrix(GICA))
  stopifnot(is.logical(center_rows) && length(center_rows)==1)
  stopifnot(is.logical(center_cols) && length(center_cols)==1)
  stopifnot(is.logical(scale) && length(scale)==1)
  stopifnot(is.logical(normA) && length(normA)==1)

  nT <- ncol(BOLD) #length of timeseries
  nV <- nrow(BOLD) #number of data locations
  if(nT > nV) warning('More time points than voxels. Are you sure?')
  if(nV != nrow(GICA)) {
    stop('The number of voxels in dat (', nV, ') and GICA (', nrow(GICA), ') must match')
  }

  nQ <- ncol(GICA) #number of ICs
  if(nQ > nV) warning('More ICs than voxels. Are you sure?')
  if(nQ > nT) warning('More ICs than time points. Are you sure?')

  # Center timeseries data across space and time if `center`
  # Standardize data scale if `scale`
  # Transpose it
  BOLD <- t(norm_BOLD(
    BOLD, center_rows=center_rows, center_cols=center_cols, 
    scale=scale, detrend_DCT=detrend_DCT
  ))

  # Center each group IC over voxels
  if (center_Gcols) { GICA - rep(colMeans(GICA), rep.int(nV, nQ)) }

	# Estimate A (IC timeseries)
	A <- (BOLD %*% GICA) %*% chol2inv(chol(crossprod(GICA)))
	if (normA) { A <- scale(A) }

	# Estimate S (IC maps)
	S <- solve(a=crossprod(A), b=crossprod(A, BOLD))

	#return result
	list(S = S, A = A)
}

#' Dual Regression wrapper
#' 
#' Wrapper to \code{dual_reg}. Can handle NIFTI or CIFTI file input, halving the 
#' data in absense of retest data, and removing nuisance regressors.
#' 
#' Used to compute dual regression for the purpose of template estimation, and 
#'  to obtain initial subject-level IC estimates for template ICA. In the latter
#'  case, \code{GICA} should be the template mean. 
#' 
#' @param BOLD,BOLD2 Subject-level fMRI data in one of the following formats: 
#'  a CIFTI file path, a \code{"xifti"}, a NIFTI file path, a \code{"nifti"}, or
#'  a \eqn{V \times T} numeric matrix.
#' 
#'  If \code{BOLD2} is provided it must be in the same format as \code{BOLD}; 
#'  \code{BOLD} will be the test data and \code{BOLD2} will be the retest data. 
#'  If \code{BOLD2} is not provided, \code{BOLD} will be split in half; 
#'  the first half will be the test data and the second half will be the retest data.
#' @param GICA Group ICA maps in a format compatible with \code{BOLD}. Can also be a
#'  (vectorized) numeric matrix (\eqn{V \times Q}) no matter the format of \code{BOLD}.
#' @param center_rows,center_cols Center BOLD data across rows (each data location's time series) or columns (each time point's image)? Default: \code{TRUE} for both.
#' @param center_Gcols Center GICA across columns (each ICA)? Default: \code{TRUE}.
#' @param scale A logical value indicating whether the fMRI timeseries should be scaled by the image standard deviation.
#' @param detrend_DCT Detrend the data? This is the number of DCT bases to use for detrending. If \code{0} (default), do not detrend.
#' @param normA Scale each IC timeseries (column of \eqn{A}) in the dual regression 
#'  estimates? Default: \code{FALSE}. (The opposite scaling will be applied to \eqn{S}
#'  such that the product \eqn{A \times S} remains the same).
#' @param format Expected format of \code{BOLD} and \code{BOLD2}. Should be one of the following:
#'  \code{"infer"} (default), a \code{"CIFTI"} file path, a \code{"xifti"} object, a
#'  \code{"NIFTI"} file path, a \code{"nifti"} object, or a \code{"data"} matrix.
#' @param brainstructures Only applies if \code{format} is \code{"CIFTI"}. 
#'  Character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
#'  cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Can also be \code{"all"} (obtain all three brain structures).
#'  Default: \code{c("left","right")} (cortical surface only).
#' @param mask Required if and only if \code{BOLD} is a NIFTI file path or
#'  \code{"nifti"} object. This is a brain map formatted as a binary array of the same 
#'  size as the fMRI data, with \code{TRUE} corresponding to in-mask voxels.
#' @param Q2,maxQ Obtain dual regression estimates after denoising? Denoising is based on modeling and
#'  removing nuisance ICs. It may result in a cleaner estimate for smaller datasets, but it may be unnecessary (and time-consuming) for larger datasets.
#'  If both arguments are \code{NULL}, denoising will be performed, with the number of nuisance 
#'  ICs estimated for \code{BOLD} and \code{BOLD2} separately. Otherwise, specify one or the other:
#'  use \code{Q2} to specify the number of nuisance ICs, or \code{maxQ} to specify the number of
#'  total ICs (group + nuisance, or \eqn{Q + Q2}). Set either to zero to skip denoising.
#'  Default: \code{Q2==0} (do not denoise).
#'  
#'  The valid inputs are \eqn{Q <= (Q+Q2) = maxQ <= T}, where \eqn{Q} is the number
#'  of group ICs and \eqn{T} is the number of timepoints in each fMRI scan. 
#' @param verbose Display progress updates? Default: \code{TRUE}.
#' 
#' @keywords internal
dual_reg2 <- function(
  BOLD, BOLD2=NULL, 
  GICA, 
  center_rows=TRUE, center_cols=TRUE, scale=TRUE, detrend_DCT=0, 
  center_Gcols=TRUE, normA=FALSE,
  format=c("infer", "CIFTI", "xifti", "NIFTI", "nifti", "data"),
  brainstructures=c("left", "right"), mask=NULL, 
  Q2=0, maxQ=NULL, 
  verbose=TRUE){

  # Check arguments ------------------------------------------
  stopifnot(is.logical(center_rows) && length(center_rows)==1)
  stopifnot(is.logical(center_cols) && length(center_cols)==1)
  stopifnot(is.logical(scale) && length(scale)==1)
  stopifnot(is.numeric(detrend_DCT) && length(detrend_DCT)==1)
  stopifnot(detrend_DCT >=0 && detrend_DCT==round(detrend_DCT))
  stopifnot(is.logical(normA) && length(normA)==1)
  if (!is.null(Q2) && !is.null(maxQ)) { stop("Specify one of `Q2` or `maxQ`.") }

  retest <- !is.null(BOLD2) # Retest or pseudo retest?

  # Initialize the return value.
  out <- list(test = NULL, retest = NULL)

  # Determine the format of `BOLD` and `BOLD2`.
  format <- match.arg(format, c("infer", "CIFTI", "xifti", "NIFTI", "nifti", "data"))
  if (format == "infer") { format <- infer_BOLD_format(BOLD) }
  if (format %in% c("CIFTI", "NIFTI") && length(BOLD) > 1) { stop("`BOLD` should be length one.") }
  if (retest) {
    format2 <- infer_BOLD_format(BOLD2)
    if (format2 != format) {
      stop("`BOLD` format is ", format, ", but `BOLD2` format is ", format2, ".")
    }
    if (format %in% c("CIFTI", "NIFTI") && length(BOLD2) > 1) { stop("`BOLD2` should be length one.") }
  }
  FORMAT <- switch(format,
    CIFTI = "CIFTI",
    xifti = "CIFTI",
    NIFTI = "NIFTI",
    nifti = "NIFTI",
    data = "DATA"
  )

  # If BOLD (and BOLD2) is a CIFTI or NIFTI file, check that the file paths exist.
  stopifnot(all(file.exists(c(BOLD, BOLD2))))

  # Get `GICA` as a numeric data matrix or array.
  if (FORMAT == "CIFTI") {
    if (is.character(GICA)) { GICA <- ciftiTools::read_xifti(GICA, brainstructures=brainstructures) }
    if (is.xifti(GICA)) { 
      xii1 <- select_xifti(GICA, 1) # for formatting output
      GICA <- as.matrix(GICA)
    }
    stopifnot(is.matrix(GICA))
  } else if (FORMAT == "NIFTI") {
    if (is.character(GICA)) { GICA <- oro.nifti::readNIfTI(GICA, reorient=FALSE) }
    stopifnot(length(dim(GICA)) > 1)
  } else {
    stopifnot(is.matrix(GICA))
  }
  nQ <- dim(GICA)[length(dim(GICA))]

  # Get `mask` as a logical array.
  #   Check `GICA` and `mask` dimensions match.
  #   Vectorize `GICA`.
  if (FORMAT == "NIFTI") {
    if (is.null(mask)) { stop("`mask` is required.") }
    if (is.character(mask)) { mask <- oro.nifti::readNIfTI(mask, reorient=FALSE) }
    if (dim(mask)[length(dim(mask))] == 1) { mask <- array(mask, dim=dim(mask)[length(dim(mask))-1]) }
    if (is.numeric(mask)) { 
      cat("Coercing `mask` to a logical array with `as.logical`.\n")
      mask[] <- as.logical(mask)
    }
    nI <- dim(mask); nV <- sum(mask)
    stopifnot(length(dim(GICA)) %in% c(2, length(nI)+1))
    if (length(dim(GICA)) == length(nI)+1) {
      if (length(dim(GICA)) != 2) {
        stopifnot(all(dim(GICA)[length(dim(GICA))-1] == nI))
      }
      if (all(dim(GICA)[length(dim(GICA))-1] == nI)) {
        GICA <- matrix(GICA[rep(mask, nQ)], ncol=nQ)
        stopifnot(nrow(GICA) == nV)  
      }
    }
  } else {
    nI <- nV <- nrow(GICA)
  }

  # Center `GICA` columns.
  if (center_Gcols) { GICA - rep(colMeans(GICA), rep.int(nV, nQ)) }

  # Get BOLD (and BOLD2) as a data matrix or array. 
  if (format == "CIFTI") {
    BOLD <- ciftiTools::read_cifti(BOLD, brainstructures=brainstructures)
    if (is.xifti(BOLD)) { BOLD <- as.matrix(BOLD) }
    stopifnot(is.matrix(BOLD))
    if (retest) {
      BOLD2 <- ciftiTools::read_cifti(BOLD2, brainstructures=brainstructures)
      if (is.xifti(BOLD2)) { BOLD <- as.matrix(BOLD2) }
      stopifnot(is.matrix(BOLD2))
    }
  } else if (format == "NIFTI") {
    BOLD <- oro.nifti::readNIfTI(BOLD, reorient=FALSE)
    stopifnot(length(dim(BOLD)) > 1)
    if (retest) {
      BOLD2 <- oro.nifti::readNIfTI(BOLD2, reorient=FALSE)
      stopifnot(length(dim(BOLD)) > 1)
    }
  }
  dBOLD <- dim(BOLD)
  ldB <- length(dim(BOLD))
  nT <- dim(BOLD)[ldB]

  # Check BOLD (and BOLD2) dimensions correspond with `GICA` and `mask`.
  stopifnot(ldB-1 == length(nI))
  stopifnot(all(dBOLD[seq(ldB-1)] == nI))
  if (retest) { stopifnot(all(dim(BOLD2)[seq(ldB-1)] == nI)) }

  # Vectorize BOLD (and BOLD2).
  if (FORMAT=="NIFTI") {
    BOLD <- matrix(BOLD[rep(mask, dBOLD[ldB])], ncol=nT)
    stopifnot(nrow(BOLD) == nV)
    if (retest) { 
      BOLD2 <- matrix(BOLD2[rep(mask, dBOLD[ldB])], ncol=nT)
      stopifnot(nrow(BOLD2) == nV)
    }
  }

  # If no retest data, halve the test data.
  if (!retest) {
    part1 <- seq(round(nT/2))
    part2 <- setdiff(seq(nT), part1)
    BOLD2 <- BOLD[, part2, drop=FALSE]
    BOLD <- BOLD[, part1, drop=FALSE]
  }

  # [TO DO]: Remove NA values?

  # Normalize BOLD (and BOLD2) again
  BOLD <- norm_BOLD(
    BOLD, center_rows=center_rows, center_cols=center_cols, 
    scale=scale, detrend_DCT=detrend_DCT
  )
  BOLD2 <- norm_BOLD(
    BOLD2, center_rows=center_rows, center_cols=center_cols, 
    scale=scale, detrend_DCT=detrend_DCT
  )

  # Perform dual regression on test and retest data
  out$test <- dual_reg(
    BOLD, GICA, center_rows=FALSE, center_cols=FALSE, 
    scale=FALSE, center_Gcols=FALSE, detrend_DCT=0
  )$S
  out$retest <- dual_reg(
    BOLD2, GICA, center_rows=FALSE, center_cols=FALSE, 
    scale=FALSE, center_Gcols=FALSE, detrend_DCT=0
  )$S

  # If !retest, must be careful about the number of ICs since nT is halved.
  nT2 <- min(ncol(BOLD), ncol(BOLD2))
  denoise <- !((!is.null(Q2) && Q2 <= 0) || (!is.null(maxQ) && maxQ <=0))

  if (!denoise) { return(out) }

  # Check `maxQ`: nQ <= maxQ <= nT2*.75
  maxQ <- maxQ_check(maxQ, L=nQ, T=nT2)

  # Estimate and deal with nuisance ICs
  if (maxQ > nQ) {
    BOLD <- rm_nuisIC(BOLD, DR=out$test, Q2=Q2, Q2_max=maxQ-nQ, verbose=verbose)
    BOLD2 <- rm_nuisIC(BOLD2, DR=out$retest, Q2=Q2, Q2_max=maxQ-nQ, verbose=verbose)
  }

  # Center and scale BOLD (and BOLD2) again, but do not detrend again.
  BOLD <- norm_BOLD(
    BOLD, center_rows=center_rows, center_cols=center_cols, 
    scale=scale, detrend_DCT=FALSE
  )
  BOLD2 <- norm_BOLD(
    BOLD2, center_rows=center_rows, center_cols=center_cols, 
    scale=scale, detrend_DCT=FALSE
  )

  # Do DR again.
  out$test_preclean <- out$test
  out$test <- dual_reg(
    BOLD, GICA, center_rows=FALSE, center_cols=FALSE, 
    scale=FALSE, detrend_DCT=0
  )$S
  out$retest_preclean <- out$retest
  out$retest <- dual_reg(
    BOLD2, GICA, center_rows=FALSE, center_cols=FALSE, 
    scale=FALSE, detrend_DCT=0
  )$S

  out
}