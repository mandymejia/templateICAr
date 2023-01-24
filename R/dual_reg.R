#' Dual Regression
#'
#' @param BOLD Subject-level fMRI data matrix (\eqn{V \times T}). Rows will be
#'  centered.
#' @param GICA Group-level independent components (\eqn{V \times Q})
#' @inheritParams center_Bcols_Param
#' @inheritParams scale_Param
#' @param scale_sm_xifti,scale_sm_FWHM Only applies if \code{scale=="local"} and
#'  \code{BOLD} represents CIFTI-format data. To smooth the standard deviation
#'  estimates used for local scaling, provide a \code{"xifti"} object with data
#'  locations in alignment with \code{BOLD}, as well as the smoothing FWHM
#'  (default: \code{2}). If no \code{"xifti"} object is provided (default), do
#'  not smooth.
#' @inheritParams detrend_DCT_Param
#' @inheritParams normA_Param
#'
#' @importFrom fMRItools colCenter
#' @importFrom matrixStats colVars
#'
#' @return A list containing
#'  the subject-level independent components \strong{S} (\eqn{V \times Q}),
#'  and subject-level mixing matrix \strong{A} (\eqn{TxQ}).
#'
#' @export
#' @examples
#' nT <- 30
#' nV <- 400
#' nQ <- 7
#' mU <- matrix(rnorm(nV*nQ), nrow=nV)
#' mS <- mU %*% diag(seq(nQ, 1)) %*% matrix(rnorm(nQ*nT), nrow=nQ)
#' BOLD <- mS + rnorm(nV*nT, sd=.05)
#' GICA <- mU
#' dual_reg(BOLD=BOLD, GICA=mU, scale="local")
#'
dual_reg <- function(
  BOLD, GICA,
  scale=c("global", "local", "none"), scale_sm_xifti=NULL, scale_sm_FWHM=2,
  detrend_DCT=0, center_Bcols=FALSE, normA=FALSE){

  stopifnot(is.matrix(BOLD))
  stopifnot(is.matrix(GICA))
  if (is.null(scale) || isFALSE(scale)) { scale <- "none" }
  if (isTRUE(scale)) {
    warning(
      "Setting `scale='global'`. Use `'global'` or `'local'` ",
      "instead of `TRUE`, which has been deprecated."
    )
    scale <- "global"
  }
  scale <- match.arg(scale, c("global", "local", "none"))
  if (!is.null(scale_sm_xifti)) { stopifnot(ciftiTools::is.xifti(scale_sm_xifti)) }
  stopifnot(is.numeric(scale_sm_FWHM) && length(scale_sm_FWHM)==1)
  stopifnot(is.logical(normA) && length(normA)==1)

  if (any(is.na(BOLD))) { stop("`NA` values in `BOLD` not supported with DR.") }
  if (any(is.na(GICA))) { stop("`NA` values in `GICA` not supported with DR.") }

  nV <- nrow(BOLD) #number of data locations
  nT <- ncol(BOLD) #length of timeseries
  if(nV < nT) warning('More time points than voxels. Are you sure?')
  if(nV != nrow(GICA)) {
    stop('The number of voxels in dat (', nV, ') and GICA (', nrow(GICA), ') must match')
  }

  nQ <- ncol(GICA) #number of ICs
  if(nQ > nV) warning('More ICs than voxels. Are you sure?')
  if(nQ > nT) warning('More ICs than time points. Are you sure?')

  # Center each voxel timecourse. Do not center the image at each timepoint.
  # Standardize scale if `scale`, and detrend if `detrend_DCT`.
  # Transpose it: now `BOLD` is TxV.
  BOLD <- t(norm_BOLD(
    BOLD, center_rows=TRUE, center_cols=center_Bcols,
    scale=scale, scale_sm_xifti=scale_sm_xifti, scale_sm_FWHM=scale_sm_FWHM,
    detrend_DCT=detrend_DCT
  ))

  # Center each group IC across space. (Used to be a function argument.)
  center_Gcols <- TRUE
  if (center_Gcols) { GICA <- fMRItools::colCenter(GICA) }

  # Estimate A (IC timeseries).
  # We need to center `BOLD` across space because the linear model has no intercept.
  A <- ((BOLD - rowMeans(BOLD, na.rm=TRUE)) %*% GICA) %*% chol2inv(chol(crossprod(GICA)))

  # Center each subject IC timecourse across time.
  # (Redundant. Since BOLD is column-centered, A is already column-centered.)
  # A <- fMRItools::colCenter(A)

  # Normalize each subject IC timecourse if `normA`.
  if (normA) { A <- scale(A) }

  # Estimate S (IC maps).
  # Don't worry about the intercept: `BOLD` and `A` are centered across time.
  S <- solve(a=crossprod(A), b=crossprod(A, BOLD))

  #return result
  list(S = S, A = A)
}

#' Dual Regression wrapper
#'
#' Wrapper to \code{dual_reg} used by `estimate_template`. The format of `BOLD`
#'  (and `BOLD2`) must be provided, and `GICA` must be vectorized if applicable.
#'
#' @param BOLD,BOLD2 Subject-level fMRI data in one of the following formats:
#'  a CIFTI file path, a \code{"xifti"} object, a NIFTI file path, a \code{"nifti"} object, or
#'  \eqn{V \times T} numeric matrices, where \eqn{V} is the number of data locations and
#'  \eqn{T} is the number of timepoints.
#'
#'  If \code{BOLD2} is provided it must be in the same format as \code{BOLD};
#'  \code{BOLD} will be the test data and \code{BOLD2} will be the retest data.
#'
#'  If \code{BOLD2} is not provided, \code{BOLD} will be split in half;
#'  the first half will be the test data and the second half will be the retest data.
#' @param GICA Group ICA maps as a (vectorized) numeric matrix
#'  (\eqn{V \times Q}). Its columns will be centered.
#' @param keepA Keep the resulting \strong{A} matrices, or only return the \strong{S} matrices
#'  (default)?
#' @inheritParams scale_Param
#' @param scale_sm_surfL,scale_sm_surfR,scale_sm_FWHM Only applies if
#'  \code{scale=="local"} and \code{BOLD} represents CIFTI-format data. To
#'  smooth the standard deviation estimates used for local scaling, provide the
#'  surface geometries along which to smooth as GIFTI geometry files or
#'  \code{"surf"} objects, as well as the smoothing FWHM (default: \code{2}).
#'
#'  If \code{scale_sm_FWHM==0}, no smoothing of the local standard deviation
#'  estimates will be performed.
#'
#'  If \code{scale_sm_FWHM>0} but \code{scale_sm_surfL} and
#'  \code{scale_sm_surfR} are not provided, the default inflated surfaces from
#'  the HCP will be used.
#'
#'  To create a \code{"surf"} object from data, see
#'  \code{\link[ciftiTools]{make_surf}}. The surfaces must be in the same
#'  resolution as the \code{BOLD} data.
#' @inheritParams detrend_DCT_Param
#' @inheritParams center_Bcols_Param
#' @inheritParams normA_Param
#' @param brainstructures Only applies if the entries of \code{BOLD} are CIFTI file paths.
#'  Character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
#'  cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Can also be \code{"all"} (obtain all three brain structures).
#'  Default: \code{c("left","right")} (cortical surface only).
#' @param mask Required if and only if the entries of \code{BOLD} are NIFTI file paths or
#'  \code{"nifti"} objects. This is a brain map formatted as a binary array of the same
#'  size as the fMRI data, with \code{TRUE} corresponding to in-mask voxels.
#' @param format Expected format of \code{BOLD} and \code{BOLD2}. Should be one
#'  of the following: a \code{"CIFTI"} file path, a \code{"xifti"} object, a
#'  \code{"NIFTI"} file path, a \code{"nifti"} object, or a \code{"data"} matrix.
#' @param Q2,Q2_max Obtain dual regression estimates after denoising? Denoising is
#'  based on modeling and removing nuisance ICs. It may result in a cleaner
#'  estimate for smaller datasets, but it may be unnecessary (and time-consuming)
#'  for larger datasets.
#'
#'  Set \code{Q2} to control denoising: use a positive integer to specify the
#'  number of nuisance ICs, \code{NULL} to have the number of nuisance ICs
#'  estimated by PESEL, or zero (default) to skip denoising.
#'
#'  If \code{is.null(Q2)}, use \code{Q2_max} to specify the maximum number of
#'  nuisance ICs that should be estimated by PESEL. \code{Q2_max} must be less
#'  than \eqn{T * .75 - Q} where \eqn{T} is the minimum number of timepoints in
#'  each fMRI scan and \eqn{Q} is the number of group ICs. If \code{NULL}
#'  (default), \code{Q2_max} will be set to \eqn{T * .50 - Q}, rounded.
#' @inheritParams varTol_Param
#' @param maskTol Tolerance for number of locations masked out due to low
#'  variance or missing values. If more than this many locations are masked out,
#'  this subject is skipped without calculating dual regression. \code{maskTol}
#'  can be specified either as a proportion of the number of locations (between
#'  zero and one), or as a number of locations (integers greater than one).
#'  Default: \code{.1}, i.e. up to 10\% of locations can be masked out.
#'
#'  If \code{BOLD2} is provided, masks are calculated for each scan and then
#'  the intersection of the masks is used.
#' @param verbose Display progress updates? Default: \code{TRUE}.
#'
#' @return The dual regression \strong{S} matrices, or both the \strong{S}
#'  and \strong{A} matrices if \code{keepA}, or \code{NULL} if dual
#'  regression was skipped due to too many masked data locations.
#'
#' @keywords internal
dual_reg2 <- function(
  BOLD, BOLD2=NULL, 
  format=c("CIFTI", "xifti", "GIFTI", "gifti", "NIFTI", "nifti", "RDS", "data"),
  GICA,
  keepA=FALSE,
  scale=c("global", "local", "none"),
  scale_sm_surfL=NULL, scale_sm_surfR=NULL, scale_sm_FWHM=2,
  detrend_DCT=0,
  center_Bcols=FALSE, normA=FALSE,
  Q2=0, Q2_max=NULL, NA_limit=.1,
  brainstructures=c("left", "right"), mask=NULL,
  varTol=1e-6, maskTol=.1,
  verbose=TRUE){

  if (verbose) { extime <- Sys.time() }

  keepA <- as.logical(keepA); stopifnot(length(keepA)==1)
  # No other arg checks: check them before calling this function.

  # For `"xifti"` data for handling the medial wall and smoothing.
  xii1 <- NULL

  # Prepare output.
  out <- list(test = NULL, retest = NULL)

  # Load helper variables.
  retest <- !is.null(BOLD2)
  format <- match.arg(format, c("CIFTI", "xifti", "GIFTI", "gifti", "NIFTI", "nifti", "RDS", "data"))
  FORMAT <- get_FORMAT(format)
  nQ <- ncol(GICA)
  nI <- nV <- nrow(GICA)

  check_req_ifti_pkg(FORMAT)

  # Get `BOLD` (and `BOLD2`) as a data matrix or array.  -----------------------
  if (verbose) { cat("\tReading and formatting data...") }
  if (FORMAT == "CIFTI") {
    if (is.character(BOLD)) { BOLD <- ciftiTools::read_cifti(BOLD, brainstructures=brainstructures) }
    if (ciftiTools::is.xifti(BOLD)) {
      if (scale == "local") {
        xii1 <- ciftiTools::convert_xifti(ciftiTools::select_xifti(BOLD, 1), "dscalar") * 0
      }
      BOLD <- as.matrix(BOLD)
    }
    stopifnot(is.matrix(BOLD))
    if (retest) {
      if (is.character(BOLD2)) { BOLD2 <- ciftiTools::read_cifti(BOLD2, brainstructures=brainstructures) }
      if (ciftiTools::is.xifti(BOLD2)) { BOLD2 <- as.matrix(BOLD2) }
      stopifnot(is.matrix(BOLD2))
    }
    nI <- nV <- nrow(GICA)
  } else if (FORMAT == "GIFTI") {
    if (is.character(BOLD)) { BOLD <- gifti::readgii(BOLD) }
    stopifnot(gifti::is.gifti(BOLD))
    ghemi <- BOLD$file_meta["AnatomicalStructurePrimary"]
    if (!(ghemi %in% c("CortexLeft", "CortexRight"))) {
      stop("AnatomicalStructurePrimary metadata missing or invalid for GICA.")
    }
    ghemi <- switch(ghemi, CortexLeft="left", CortexRight="right")
    if (scale == "local") {
      if (ghemi == "left") {
        xii1 <- ciftiTools::select_xifti(ciftiTools::as.xifti(cortexL=do.call(cbind, BOLD$data)), 1) * 0
      } else if (ghemi == "right") {
        xii1 <- ciftiTools::select_xifti(ciftiTools::as.xifti(cortexR=do.call(cbind, BOLD$data)), 1) * 0
      } else { stop() }
      xii1$meta$cifti$intent <- 3006
    }
    BOLD <- do.call(cbind, BOLD$data)

    stopifnot(is.matrix(BOLD))
    if (retest) {
      if (is.character(BOLD2)) { BOLD2 <- gifti::readgii(BOLD2) }
      if (inherits(BOLD2, "gifti")) { BOLD2 <- do.call(cbind, BOLD2$data) }
      stopifnot(is.matrix(BOLD2))
    }
    nI <- nV <- nrow(GICA)
  } else if (FORMAT == "NIFTI") {
    if (is.character(BOLD)) { BOLD <- RNifti::readNifti(BOLD) }
    stopifnot(length(dim(BOLD)) > 1)
    if (retest) {
      if (is.character(BOLD2)) { BOLD2 <- RNifti::readNifti(BOLD2) }
      stopifnot(length(dim(BOLD2)) > 1)
    }
    stopifnot(!is.null(mask))
    nI <- dim(drop(mask)); nV <- sum(mask)
  } else if (FORMAT == "MATRIX") {
    if (is.character(BOLD)) { BOLD <- readRDS(BOLD) }
    stopifnot(is.matrix(BOLD))
    if (retest) {
      if (is.character(BOLD2)) { BOLD2 <- readRDS(BOLD2) }
      stopifnot(is.matrix(BOLD2))
    }
    nI <- nV <- nrow(GICA)
  } else { stop() }
  dBOLD <- dim(BOLD)
  ldB <- length(dim(BOLD))
  nT <- dim(BOLD)[ldB]

  # If `retest`, ensure that spatial dimensions of `BOLD2` match with `BOLD`.
  if (retest) {
    stopifnot(length(dim(BOLD)) == length(dim(BOLD2)))
    stopifnot(all(dBOLD[seq(ldB-1)] == dim(BOLD2)[seq(ldB-1)]))
  }

  # Check BOLD (and BOLD2) dimensions correspond with `GICA` and `mask`.
  if(!(ldB-1 == length(nI))) { stop("`GICA` and BOLD spatial dimensions do not match.") }
  if(!all(dBOLD[seq(ldB-1)] == nI)) { stop("`GICA` and BOLD spatial dimensions do not match.") }

  # Vectorize `BOLD` (and `BOLD2`). --------------------------------------------
  if (FORMAT=="NIFTI") {
    BOLD <- matrix(BOLD[rep(as.logical(mask), dBOLD[ldB])], ncol=nT)
    stopifnot(nrow(BOLD) == nV)
    if (retest) {
      BOLD2 <- matrix(BOLD2[rep(as.logical(mask), dBOLD[ldB])], ncol=nT)
      stopifnot(nrow(BOLD2) == nV)
    }
  }

  # Check for missing values. --------------------------------------------------
  nV0 <- nV # not used
  mask <- make_mask(BOLD, varTol=varTol)
  if (retest) { mask <- mask & make_mask(BOLD2, varTol=varTol) }
  use_mask <- !all(mask)
  if (use_mask) {
    # Coerce `maskTol` to number of locations.
    stopifnot(is.numeric(maskTol) && length(maskTol)==1 && maskTol >= 0)
    if (maskTol < 1) { maskTol <- maskTol * nV }
    # Skip this scan if `maskTol` is surpassed.
    if (sum(!mask) > maskTol) { return(NULL) }
    # Mask out the locations.
    BOLD <- BOLD[mask,,drop=FALSE]
    GICA <- GICA[mask,,drop=FALSE]
    if (retest) { BOLD2 <- BOLD2[mask,,drop=FALSE] }
    if (!is.null(xii1)) {
      xiitmp <- as.matrix(xii1)
      xiitmp[!mask,] <- NA
      xii1 <- ciftiTools::move_to_mwall(ciftiTools::newdata_xifti(xii1, xiitmp))
    }
    nV <- nrow(BOLD)

    # [TO DO]: replace with fMRIscrub::unmask_mat(..., mask_dim=2)
    # For later
    unmask <- function(S, mask) {
      S2 <- matrix(NA, nrow=nrow(S), ncol=length(mask))
      S2[,mask] <- S
      S2
    }
  }

  if (!is.null(xii1) && scale=="local" && scale_sm_FWHM > 0) {
    xii1 <- ciftiTools::add_surf(xii1, surfL=scale_sm_surfL, surfR=scale_sm_surfR)
  }

  # Helper functions
  this_norm_BOLD <- function(B){ norm_BOLD(
    B, center_rows=TRUE, center_cols=center_Bcols,
    scale=scale, scale_sm_xifti=xii1, scale_sm_FWHM=scale_sm_FWHM,
    detrend_DCT=detrend_DCT
  ) }

  dual_reg_yesNorm <- function(B){ dual_reg(
    B, GICA,
    scale=scale, scale_sm_xifti=xii1, scale_sm_FWHM=scale_sm_FWHM,
    detrend_DCT=detrend_DCT, center_Bcols=center_Bcols, normA=normA
  ) }

  dual_reg_noNorm <- function(B){ dual_reg(
    B, GICA,
    scale="none", detrend_DCT=0, center_Bcols=FALSE, normA=normA
  ) }

  # Get the first dual regression results. -------------------------------------
  if (verbose) { cat(" Computing DR...") }
  # If using pseudo-retest data, compute DR on the halves of `BOLD`.
  # Do this before normalizating `BOLD` so to avoid normalizing twice.
  if (!retest) {
    part1 <- seq(round(nT/2))
    part2 <- setdiff(seq(nT), part1)
    out$test <- dual_reg_yesNorm(BOLD[, part1, drop=FALSE])
    out$retest <- dual_reg_yesNorm(BOLD[, part2, drop=FALSE])
  } else {
    # If retest, normalize `BOLD` and `BOLD2`, and then compute DR.
    BOLD <- this_norm_BOLD(BOLD)
    BOLD2 <- this_norm_BOLD(BOLD2)
    # (No need to normalize again.)
    out$test <- dual_reg_noNorm(BOLD)
    out$retest <- dual_reg_noNorm(BOLD2)
  }

  # Return these DR results if denoising is not needed. ------------------------
  if ((!is.null(Q2) && Q2==0) || (!is.null(Q2_max) && Q2_max==0)) {
    if (!keepA) {
      out$test$A <- NULL
      out$retest$A <- NULL
    }
    if (use_mask) {
      out$test$S <- unmask(out$test$S, mask)
      out$retest$S <- unmask(out$retest$S, mask)
    }
    if (verbose) { cat(" Done!\n") }
    if (verbose) { print(Sys.time() - extime) }
    return(out)
  }

  # Estimate and deal with nuisance ICs. ---------------------------------------
  if (verbose) { cat(" Denoising...") }
  # If !retest, we prefer to estimate nuisance ICs across the full scan
  # and then halve it after.
  if (!retest) {
    BOLD <- this_norm_BOLD(BOLD)
    BOLD_DR <- dual_reg_noNorm(BOLD)
    BOLD <- rm_nuisIC(BOLD, DR=BOLD_DR, Q2=Q2, Q2_max=Q2_max, verbose=verbose)
    rm(BOLD_DR)
    BOLD2 <- BOLD[, part2, drop=FALSE]
    BOLD <- BOLD[, part1, drop=FALSE]
  } else {
    BOLD <- rm_nuisIC(BOLD, DR=out$test, Q2=Q2, Q2_max=Q2_max, verbose=verbose)
    BOLD2 <- rm_nuisIC(BOLD2, DR=out$retest, Q2=Q2, Q2_max=Q2_max, verbose=verbose)
  }

  # Center and scale `BOLD` and `BOLD2` (again), but do not detrend again. -----
  BOLD <- norm_BOLD(
    BOLD, center_rows=TRUE, center_cols=center_Bcols,
    scale=scale, scale_sm_xifti=xii1, scale_sm_FWHM=scale_sm_FWHM,
    detrend_DCT=0
  )
  BOLD2 <- norm_BOLD(
    BOLD2, center_rows=TRUE, center_cols=center_Bcols,
    scale=scale, scale_sm_xifti=xii1, scale_sm_FWHM=scale_sm_FWHM,
    detrend_DCT=0
  )

  # Do DR again. ---------------------------------------------------------------
  if (verbose) { cat(" Computing DR again...") }

  out$test_preclean <- out$test
  out$test <- dual_reg_noNorm(BOLD)
  out$retest_preclean <- out$retest
  out$retest <- dual_reg_noNorm(BOLD2)

  if (!keepA) {
    out$test_preclean$A <- NULL
    out$test$A <- NULL
    out$retest_preclean$A <- NULL
    out$retest$A <- NULL
  }

  if (use_mask) {
    out$test_preclean$S <- unmask(out$test_preclean$S, mask)
    out$retest_preclean$S <- unmask(out$retest_preclean$S, mask)
    out$test$S <- unmask(out$test$S, mask)
    out$retest$S <- unmask(out$retest$S, mask)
  }

  if (verbose) { cat(" Done!\n") }
  if (verbose) { print(Sys.time() - extime) }
  out
}
