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
#' @param GICA_parc_table Is the GICA actually a parcellation? If so, provide
#'  the parcellation table here. Default: \code{NULL}.
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
#' @inheritParams TR_param
#' @inheritParams hpf_param
#' @inheritParams GSR_Param
#' @param brainstructures Only applies if the entries of \code{BOLD} are CIFTI file paths.
#'  Character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
#'  cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Can also be \code{"all"} (obtain all three brain structures).
#'  Default: \code{c("all")}.
#' @param resamp_res Only applies if the entries of \code{BOLD} are CIFTI file paths.
#'  Resample the data upon reading it in? Default: \code{NULL} (no resampling).
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
#' @importFrom fMRItools dual_reg
#'
#' @keywords internal
dual_reg2 <- function(
  BOLD, BOLD2=NULL,
  format=c("CIFTI", "xifti", "GIFTI", "gifti", "NIFTI", "nifti", "RDS", "data"),
  GICA, GICA_parc_table=NULL,
  mask=NULL,
  keepA=FALSE,
  scale=c("local", "global", "none"),
  scale_sm_surfL=NULL, scale_sm_surfR=NULL, scale_sm_FWHM=2,
  TR=NULL, hpf=.01,
  GSR=FALSE,
  Q2=0, Q2_max=NULL, 
  NA_limit=.1,
  brainstructures="all", resamp_res=NULL,
  varTol=1e-6, maskTol=.1,
  verbose=TRUE){

  if (verbose) { extime <- Sys.time() }

  keepA <- as.logical(keepA); stopifnot(length(keepA)==1)
  scale <- match.arg(scale, c("local", "global", "none"))
  # No other arg checks: check them before calling this function.

  # For `"xifti"` data for handling the medial wall and smoothing.
  xii1 <- NULL

  # Prepare output.
  out <- list(test = NULL, retest = NULL)

  # Load helper variables.
  retest <- !is.null(BOLD2)
  format <- match.arg(format, c("CIFTI", "xifti", "GIFTI", "gifti", "NIFTI", "nifti", "RDS", "data"))
  FORMAT <- get_FORMAT(format)
  check_req_ifti_pkg(FORMAT)

  GICA_parc <- !is.null(GICA_parc_table)
  nQ <- if (GICA_parc) { nrow(GICA_parc_table) } else { ncol(GICA) }

  if (is.null(mask)) {
    nI <- nV <- nrow(GICA)
  } else if (FORMAT=="NIFTI") {
    nI <- dim(drop(mask))
    nV <- sum(mask)
  } else {
    nI <- length(mask); nV <- sum(mask)
  }

  # Get `BOLD` (and `BOLD2`) as a data matrix or array.  -----------------------
  if (verbose) { cat("\tReading in data... ") }
  if (FORMAT == "CIFTI") {
    if (is.character(BOLD)) { BOLD <- ciftiTools::read_cifti(BOLD, brainstructures=brainstructures, resamp_res=resamp_res) }
    if (ciftiTools::is.xifti(BOLD)) {
      if (scale == "local") {
        xii1 <- ciftiTools::convert_xifti(ciftiTools::select_xifti(BOLD, 1), "dscalar") * 0
      }
      BOLD <- as.matrix(BOLD)
    }
    stopifnot(is.matrix(BOLD))
    if (retest) {
      if (is.character(BOLD2)) { BOLD2 <- ciftiTools::read_cifti(BOLD2, brainstructures=brainstructures, resamp_res=resamp_res) }
      if (ciftiTools::is.xifti(BOLD2)) { BOLD2 <- as.matrix(BOLD2) }
      stopifnot(is.matrix(BOLD2))
    }
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
  } else if (!is.null(mask)) {
    # Mask out the locations.
    BOLD <- BOLD[mask,,drop=FALSE]
    if (!is.null(xii1)) {
      xiitmp <- as.matrix(xii1)
      xiitmp[!mask,] <- NA
      xii1 <- ciftiTools::move_to_mwall(ciftiTools::newdata_xifti(xii1, xiitmp))
    }
    nV <- nrow(BOLD)
    if (retest) {
      BOLD2 <- BOLD2[mask,,drop=FALSE]
      stopifnot(nrow(BOLD2)==nV)
    }
  }

  # Check for missing values. --------------------------------------------------
  nV0 <- nV # not used
  mask2 <- make_mask(BOLD, varTol=varTol)
  if (retest) { mask2 <- mask2 & make_mask(BOLD2, varTol=varTol) }
  use_mask2 <- !all(mask2)
  if (use_mask2) {
    # Coerce `maskTol` to number of locations.
    stopifnot(is.numeric(maskTol) && length(maskTol)==1 && maskTol >= 0)
    if (maskTol < 1) { maskTol <- maskTol * nV }
    # Skip this scan if `maskTol` is surpassed.
    if (sum(!mask2) > maskTol) { return(NULL) }
    # Mask out the locations.
    BOLD <- BOLD[mask2,,drop=FALSE]
    GICA <- GICA[mask2,,drop=FALSE]
    if (retest) { BOLD2 <- BOLD2[mask2,,drop=FALSE] }
    if (!is.null(xii1)) {
      xiitmp <- as.matrix(xii1)
      xiitmp[!mask2,] <- NA
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
    unmask_vec <- function(vec, mask) {
      vec2 <- rep(NA, length(mask))
      vec2[mask] <- vec
      vec2
    }
  }

  if (!is.null(xii1) && scale=="local" && scale_sm_FWHM > 0) {
    xii1 <- ciftiTools::add_surf(xii1, surfL=scale_sm_surfL, surfR=scale_sm_surfR)
  }

  # Helper functions
  this_norm_BOLD <- function(B){ norm_BOLD(
    B, center_rows=TRUE, center_cols=GSR,
    scale=scale, scale_sm_xifti=xii1, scale_sm_FWHM=scale_sm_FWHM,
    TR=TR, hpf=hpf
  ) }

  DR_FUN <- if (GICA_parc) {
    function(GICA, ...) { fMRItools::dual_reg_parc(parc=GICA, ...) }
  } else {
    # Do twice to get timecourse estimate w/ subject maps, rather than w/ GICA (`A2`)
    function(GICA, parc_vals, ...) {
      out <- fMRItools::dual_reg(GICA=GICA, ...)
      GICA <- t(out$S)
      out$A2 <- fMRItools::dual_reg(GICA=GICA, ...)$A
      out
    }
  }

  dual_reg_yesNorm <- function(B){ DR_FUN(
    B, GICA=GICA, parc_vals=GICA_parc_table$Key,
    scale=scale, scale_sm_xifti=xii1, scale_sm_FWHM=scale_sm_FWHM,
    TR=TR, hpf=hpf, GSR=GSR
  ) }

  dual_reg_noNorm <- function(B){ DR_FUN(
    B, GICA=GICA, parc_vals=GICA_parc_table$Key,
    scale="none", hpf=0, GSR=FALSE
  ) }

  # Get the first dual regression results. -------------------------------------
  if (verbose) { cat("\n\tDual regression... ") }
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

  BOLDss <- list(
    test = if (!retest) { BOLD[, part1, drop=FALSE] } else { BOLD },
    retest = if (!retest) { BOLD[, part2, drop=FALSE] } else { BOLD2 }
  )
  BOLDss$test_preclean <- BOLDss$test
  BOLDss$retest_preclean <- BOLDss$retest

  # Return these DR results if denoising is not needed. ------------------------
  if ((!is.null(Q2) && Q2==0) || (!is.null(Q2_max) && Q2_max==0)) {
    for (sess in c("test", "retest")) {
      out[[sess]]$sigma_sq <- colSums((out[[sess]]$A %*% out[[sess]]$S - t(BOLDss[[sess]]))^2)/nT # part inside colSums() is TxV
      if (use_mask2) { out[[sess]]$sigma_sq <- unmask_vec(out[[sess]]$sigma_sq, mask2) }
      if (!keepA) { out[[sess]]$A <- NULL; out[[sess]]$A2 <- NULL }
      if (use_mask2) { out[[sess]]$S <- unmask(out[[sess]]$S, mask2) }
    }

    if (verbose) { cat(" Done!\n") }
    if (verbose) { print(Sys.time() - extime) }
    return(out)
  }

  # Estimate and deal with nuisance ICs. ---------------------------------------
  if (verbose) { cat(" Denoising... ") }
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
    BOLD, center_rows=TRUE, center_cols=GSR,
    scale=scale, scale_sm_xifti=xii1, scale_sm_FWHM=scale_sm_FWHM,
    TR=TR, hpf=0
  )
  BOLD2 <- norm_BOLD(
    BOLD2, center_rows=TRUE, center_cols=GSR,
    scale=scale, scale_sm_xifti=xii1, scale_sm_FWHM=scale_sm_FWHM,
    TR=TR, hpf=0
  )

  # Do DR again. ---------------------------------------------------------------
  if (verbose) { cat("\n\tDual regression again... ") }

  BOLDss <- list(test = BOLD, retest = BOLD2)
  BOLDss$test_preclean <- BOLDss$test
  BOLDss$retest_preclean <- BOLDss$retest

  out$test_preclean <- out$test
  out$test <- dual_reg_noNorm(BOLD)
  out$retest_preclean <- out$retest
  out$retest <- dual_reg_noNorm(BOLD2)

  for (sess in c("test", "retest", "test_preclean", "retest_preclean")) {
    out[[sess]]$sigma_sq <- colSums((out[[sess]]$A %*% out[[sess]]$S - t(BOLDss[[sess]]))^2)/nT # part inside colSums() is TxV
    if (use_mask2) { out[[sess]]$sigma_sq <- unmask_vec(out[[sess]]$sigma_sq, mask2) }
    if (!keepA) { out[[sess]]$A <- NULL; out[[sess]]$A2 <- NULL }
    if (use_mask2) { out[[sess]]$S <- unmask(out[[sess]]$S, mask2) }
  }

  if (verbose) { cat(" Done!\n") }
  if (verbose) { print(Sys.time() - extime) }
  out
}
