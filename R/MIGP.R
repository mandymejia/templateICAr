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
#' @param verbose Occasional updates? Default: \code{TRUE}.
#' @param ... Additional arguments to \code{datProcFUN}
MIGP <- function(dat, datProcFUN, checkColCentered=TRUE, nM=NULL, nP=NULL, verbose=TRUE, ...){
  # Arg checks -----------------------------------------------------------------
  stopifnot(is.list(dat))
  stopifnot(is.function(datProcFUN))

  nN <- length(dat)

  # First subject --------------------------------------------------------------
  # Read in and process.
  dn <- dat[[1]]
  dn <- datProcFUN(dn, ...)
  nT <- nrow(dn)
  nV <- ncol(dn)
  cat('Number of subjects:            ', nN, "\n")
  cat('Number of data locations:      ', nV, "\n")
  if (verbose) { cat("Subject 1: ") }
  if (verbose) { cat(nT, " timepoints.\n") }

  # Checks
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
  W <- dn # TxV

  # All other subjects ---------------------------------------------------------
  for (nn in seq(2, nN)) {
    if (verbose) { cat(paste0("Subject ", nn, ": ")) }
    # Read in and process.
    dn <- dat[[nn]]
    dn <- datProcFUN(dn, ...)
    if (verbose) { cat(nrow(dn), " timepoints.") }

    # Checks
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
      W <- crossprod(Wu, W)[seq(nM),,drop=FALSE]
    } else {
      if (verbose) { cat("\n") }
    }
    cat("\t", nn, "~", dim(W)[1], ",", dim(W)[2], "\n")
  }

  W[seq(min(nP, nrow(W))),,drop=FALSE]
}

#' CIFTI data processing function for MIGP
#'
datProcFUN.cifti <- function(
  dat, brainstructures=c("left", "right"), resamp_res=NULL,
  center_Bcols=FALSE, scale=TRUE, detrend_DCT=0){

  # Simple argument checks.
  stopifnot(is.logical(scale) && length(scale)==1)
  if (isFALSE(detrend_DCT)) { detrend_DCT <- 0 }
  stopifnot(is.numeric(detrend_DCT) && length(detrend_DCT)==1)
  stopifnot(detrend_DCT >=0 && detrend_DCT==round(detrend_DCT))
  stopifnot(is.logical(center_Bcols) && length(center_Bcols)==1)

  brainstructures <- match.arg(
    brainstructures,
    c("left", "right", "subcortical", "all"),
    several.ok=TRUE
  )

  # Read in data.
  if (is.character(dat)) {
    dat <- lapply(
      dat, read_xifti,
      brainstructures=brainstructures, resamp_res=resamp_res
    )
  }
  stopifnot(all(vapply(dat, is.xifti, messages=FALSE, FALSE)))

  # Normalize each scan (keep in `"xifti"` format for `merge_xifti` next).
  dat <- lapply(dat, function(x){
    newdata_xifti(x, norm_BOLD(
      as.matrix(x), center_cols=center_Bcols, scale=scale, detrend_DCT=detrend_DCT
    ))
  })

  # Concatenate and convert to matrix.
  # `merge_xifti` will check that voxels and vertices align;
  #   i.e. that the resolutions are the same.
  dat <- as.matrix(merge_xifti(xifti_list=dat))

  t(dat) # Return TxV matrix
}
