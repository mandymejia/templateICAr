#' MIGP
#'
#' Melodic's Incremental Group PCA
#'
#' https://doi.org/10.1016/j.neuroimage.2014.07.051
#'
#' Appendix A
#'
#' @param dat The files or data to compute PCA for, in list format, with each
#'  subject as a separate list entry.
#' @param datProcFUN Processes each entry of `dat` into a \eqn{T \times V}
#'  matrix (timepoints by locations) which has been column-centered (each
#'  location's timecourse is mean zero) and variance normalized.
#' @param checkColCentered Check that each subject's data is column centered
#'  after processing it with \code{datProcFUN}? Default: \code{TRUE}.
#' @param nM The size of the "internal" PCA space. Must be larger than the
#'  typical number of timepoints in each individual dataset, \eqn{T}. If
#'  \code{NULL} (default), set to \code{T*2}.
#' @param nP The number of final group PCA components to obtain. If
#'  \code{NULL} (default), will be set to \code{nM}.
#' @param initW The initial \eqn{D \times V^t} matrix, for example from a previous
#'  call to MIGP. This argument could be used to iteratively build an MIGP result
#'  across many fMRI scans. Default: \code{NULL} (no prior data).
#' @param verbose Occasional updates? Default: \code{TRUE}.
#' @param ... Additional arguments to \code{datProcFUN}
#' 
#' @return The dimension-reduced data
#' 
#' @export
MIGP <- function(dat, datProcFUN, checkColCentered=TRUE, nM=NULL, nP=NULL, initW=NULL, verbose=TRUE, ...){
  # Arg checks -----------------------------------------------------------------
  if (is.character(dat)) { message("Treating each file as a separate subject."); dat <- as.list(dat) }
  stopifnot(is.list(dat))
  stopifnot(is.function(datProcFUN))

  nN <- length(dat)
  cat('Number of subjects:            ', nN, "\n")

  # First subject --------------------------------------------------------------
  # Read in and process.
  dn <- dat[[1]]
  dn <- datProcFUN(dn, ...)
  nT <- nrow(dn)
  nV <- ncol(dn)
  cat('Number of data locations:      ', nV, "\n")
  if (verbose) { cat("Subject 1: ") }
  if (verbose) { cat(nT, " timepoints.\n") }

  # Checks
  stopifnot(all(!is.na(dn)))
  if (nT > nV) { warning(
    "Data should be TxV after processing, ",
    "but for the first scan there are more rows than columns."
  )}
  if (checkColCentered) {
    if (max(colMeans(dn)) > 1e-8) {
      stop("First subject's data columns are not demeaned.")
    }
  }

  # Set nM and nP
  if (is.null(nM)) {
    nM <- nT * 2
  } else {
    stopifnot(is.numeric(nM))
    stopifnot(length(nM)==1)
    if (nM <= nT) {warning(
      "`nM` should be larger than the ",
      "typical number of timepoints per subject."
    ) }
  }
  if (is.null(nP)) {
    nP <- nM
  } else {
    stopifnot(is.numeric(nP))
    stopifnot(length(nP)==1)
    stopifnot(nP > 1)
  }

  # Initialize W, the running estimate of the final group-average eigenvectors
  W <- if (is.null(initW)) { dn } else { rbind(initW, dn) }
  
  # All other subjects ---------------------------------------------------------
  for (nn in seq(2, nN)) {
    if (verbose) { cat(paste0("Subject ", nn, ": ")) }
    # Read in and process.
    dn <- dat[[nn]]
    dn <- datProcFUN(dn, ...)
    if (verbose) { cat(nrow(dn), " timepoints.") }

    # Checks
    stopifnot(all(!is.na(dn)))
    if (ncol(dn) != nV) {
      stop("Subject ", nn, " has ", ncol(dn), " locations (", nV, " expected).")
    }
    if (checkColCentered) {
      if (max(colMeans(dn)) > 1e-8) {
        stop("Data columns for subject ", nn, " are not demeaned.")
      }
    }

    # Concatenate
    W <- rbind(W, dn)

    # PCA
    if (nrow(W) > nM) {
      if (verbose) { cat(paste0(" Doing PCA.\n")) }
      W <- W - rowMeans(W)
      Wu <- svd(tcrossprod(W), nu=nM, nv=0)$u
      # U'W = U'UDV' = DV'
      W <- crossprod(Wu, W)[seq(nM),,drop=FALSE]
    } else {
      if (verbose) { cat("\n") }
    }
    cat("\t", nn, "~", dim(W)[1], ",", dim(W)[2], "\n")
  }

  # Move this into Group ICA code
  # Possible [TO DO]:
  # https://www.humanconnectome.org/study/hcp-young-adult/document/mound-and-moat-effect
  # Fit a Wishart ('pure noise') eigenspectrum to the tail of the estimated 
  # eigenspectrum and adjust the estimated eigenspectrum by subtracting the 
  # noise part. This process is referred to as ROW ('roll-off Wishart').

  W[seq(min(nP, nrow(W))),,drop=FALSE]
}

#' CIFTI data processing function for MIGP
#'
#' @param dat The CIFTI data 
#' @param brainstructures Character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
#'  cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Can also be \code{"all"} (obtain all three brain structures).
#'  Default: \code{c("all")}.
#' @param resamp_res The target resolution for resampling (number of cortical
#'  surface vertices per hemisphere).
#' @param GSR,scale,scale_sm_FWHM Center BOLD columns, scale by the
#'  standard deviation, and detrend voxel timecourses? See 
#'  \code{\link{norm_BOLD}}. Normalization is applied separately to each scan.
#'  Defaults: Center BOLD columns, scale by the local standard deviation, but
#'  do not detrend.
#' @inheritParams TR_param
#' @inheritParams hpf_param
#' 
#'  Note that elsewhere in \code{templateICAr} global scaling is used, but 
#'  to match the MELODIC/MIGP default local scaling is used here.
#' 
#' @keywords internal
datProcFUN.cifti <- function(
  dat, brainstructures="all", resamp_res=NULL,
  GSR=FALSE, 
  scale=c("local", "global", "none"), scale_sm_FWHM=2,
  TR=NULL, hpf=.01){

  # Simple argument checks.
  if (is.null(scale) || isFALSE(scale)) { scale <- "none" }
  if (isTRUE(scale)) { 
    warning(
      "Setting `scale='local'`. Use `'local'` or `'global'` ",
      "instead of `TRUE`, which has been deprecated."
    )
    scale <- "local"
  }
  scale <- match.arg(scale, c("local", "global", "none"))
  stopifnot(is.numeric(scale_sm_FWHM) && length(scale_sm_FWHM)==1)
  stopifnot(is.logical(GSR) && length(GSR)==1)

  brainstructures <- match.arg(
    brainstructures,
    c("left", "right", "subcortical", "all"),
    several.ok=TRUE
  )

  # Read in data.
  if (is.character(dat)) {
    dat <- lapply(
      dat, read_cifti,
      brainstructures=brainstructures, resamp_res=resamp_res
    )
  }
  stopifnot(all(vapply(dat, ciftiTools::is.xifti, messages=FALSE, FALSE)))

  # Normalize each scan (keep in `"xifti"` format for `merge_xifti` next).
  dat <- lapply(dat, function(x){
    ciftiTools::newdata_xifti(x, norm_BOLD(
      as.matrix(x), center_cols=GSR, 
      scale=scale, 
      scale_sm_xifti=ciftiTools::select_xifti(x, 1), scale_sm_FWHM=scale_sm_FWHM,
      TR=TR, hpf=hpf
    ))
  })

  # Concatenate and convert to matrix.
  # `merge_xifti` will check that voxels and vertices align;
  #   i.e. that the resolutions are the same.
  dat <- as.matrix(merge_xifti(xifti_list=dat))

  t(dat) # Return TxV matrix
}
