#' Estimate template from DR
#'
#' Estimate variance decomposition and templates from DR estimates.
#'
#' @param DR the test/retest dual regression estimates, as an array with
#'  dimensions \eqn{M \times N \times (L \times V)}, where \eqn{M} is the number
#'  of visits (2), \eqn{N} is the number of subjects, \eqn{L} is the number of
#'  IC networks, and \eqn{V} is the number of data locations.
#'
#'  (\eqn{L} and \eqn{V} are collapsed because they are treated equivalently
#'  in the context of calculating the variance decomposition and templates).
#' @param LV A length-two integer vector giving the dimensions \eqn{L} and
#'  \eqn{V} to reshape the result. Default: \code{NULL} (do not reshape the
#'  result).
#'
#' @return List of two elements: the templates and the variance decomposition.
#'
#'  There are two version of the variance template: \code{varUB} gives the
#'  unbiased variance estimate, and \code{varNN} gives the upwardly-biased
#'  non-negative variance estimate. Values in \code{varUB} will need to be
#'  clamped above zero before using in \code{templateICA}.
#'
#' @importFrom fMRItools var_decomp
#' @export
estimate_template_from_DR <- function(
  DR, LV=NULL){

  # Check arguments.
  stopifnot(length(dim(DR)) == 3)
  nM <- dim(DR)[1]  # visits
  nN <- dim(DR)[2]  # subjects
  nLV <- dim(DR)[3] # locations & networks
  if (!is.null(LV)) {
    stopifnot(is.numeric(nLV) && all(nLV > 0) && all(nLV == round(nLV)))
    stopifnot(prod(LV) == nLV)
  }

  # Variance decomposition
  vd <- var_decomp(DR)

  # Template calculation
  # Below true for M==2. Double check correct for M > 3? (Not used currently.)
  MSB_divM <- (vd$SSB / (nN-1)) / nM
  MSE_divM <- (vd$SSR / ((nM-1)*(nN-1))) / nM
  template <- list(
    mean = vd$grand_mean,
    varUB = MSB_divM - MSE_divM,
    varNN = MSB_divM
  )

  # Format `vd`
  vd$nM <- vd$nS <- vd$grand_mean <- NULL # Get rid of redundant entries

  ## Format `template`: clamp var est above zero.
  # template$varUB <- pmax(0, template$varUB)

  # Format as matrix.
  if (!is.null(LV)) {
    template <- lapply(template, function(x){ matrix(x, nrow=LV[1], ncol=LV[2]) })
    vd <- lapply(vd, function(x){ matrix(x, nrow=LV[1], ncol=LV[2]) })
  }

  # Return
  list(template=template, var_decomp=vd)
}

#' Estimate template from DR estimates (when there are two measurements)
#'
#' Legacy version of \code{\link{estimate_template_from_DR}}
#'
#' @param DR1,DR2 the test and retest dual regression estimates (\eqn{N \times L \times V})
#'
#' @return List of two elements: the mean and variance templates
#' @keywords internal
estimate_template_from_DR_two <- function(DR1, DR2){

  # Check arguments.
  stopifnot(length(dim(DR1)) == length(dim(DR2)))
  stopifnot(all(dim(DR1) == dim(DR2)))
  N <- dim(DR1)[1]

  template <- list(mean=NULL, var=NULL)

  # Mean.
  template$mean <- t(colMeans(DR1 + DR2, na.rm=TRUE) / 2)

  # Variance.
  SSB <- 2 * colSums(((DR1 + DR2)/2 - rep(t(template$mean), each=N))^2, na.rm=TRUE)
  template$var_nn <- t(SSB / (N-1)) / 2 # MSB/2
  # Unbiased.
  # 1. Fastest method.
  var_noise <- t( (1/2) * apply(DR1 - DR2, c(2,3), var, na.rm=TRUE) )
  template$var_ub <- template$var_nn - var_noise/2

  # # 2. Previous, equivalent calculation.
  # var_tot1 <- apply(DR1, c(2,3), var, na.rm=TRUE)
  # var_tot2 <- apply(DR2, c(2,3), var, na.rm=TRUE)
  # var_tot <- t((var_tot1 + var_tot2)/2)
  # # noise (within-subject) variance
  # DR_diff <- DR1 - DR2;
  # var_noise <- t((1/2)*apply(DR_diff, c(2,3), var, na.rm=TRUE))
  # # signal (between-subject) variance
  # template$var <- var_tot - var_noise
  #
  # # 3. Another equivalent calculation.
  # template$var <- t(apply(
  #   abind::abind(DR1, DR2, along=1),
  #   seq(2, 3),
  #   function(q){ cov(q[seq(N)], q[seq(N+1, 2*N)], use="complete.obs") }
  # ))

  # Make negative estimates equal to zero.
  template$var_ub[template$var_ub < 0] <- 0

  template
}

#' Estimate FC template
#'
#' @param FC0 The FC estimates from \code{\link{estimate_template}}.
#' @importFrom matrixStats colVars
#' @keywords internal
estimate_template_FC <- function(FC0){

  nL <- dim(FC0)[3]
  stopifnot(nL == dim(FC0)[4])

  FC1 <- FC0[1,,,]; FC2 <- FC0[2,,,]
  mean_FC <- (colMeans(FC1, na.rm=TRUE) + colMeans(FC2, na.rm=TRUE)) / 2
  var_FC_tot  <- (apply(FC1, c(2, 3), var, na.rm=TRUE) + apply(FC2, c(2, 3), var, na.rm=TRUE))/2
  var_FC_within  <- 1/2*(apply(FC1-FC2, c(2, 3), var, na.rm=TRUE))
  var_FC_between <- var_FC_tot - var_FC_within
  var_FC_between[var_FC_between < 0] <- NA

  nu_est <- estimate_nu_matrix(var_FC_between, mean_FC)
  nu_est1 <- quantile(nu_est[upper.tri(nu_est, diag=TRUE)], 0.01, na.rm = TRUE)

  psi <- mean_FC*(nu_est1 - nL - 1)

  list(nu = nu_est1, psi = psi)
}

#' Estimate template
#'
#' Estimate template for Template ICA based on fMRI data
#'
#' All fMRI data (entries in \code{BOLD} and \code{BOLD2}, and \code{GICA}) must
#'  be in the same spatial resolution.
#'
#' @param BOLD,BOLD2 Vector of subject-level fMRI data in one of the following
#'  formats: CIFTI file paths, \code{"xifti"} objects, GIFTI file paths,
#'  \code{"gifti"} objects, NIFTI file paths, \code{"nifti"} objects,
#'  or \eqn{V \times T} numeric matrices, where \eqn{V} is the number of data
#'  locations and \eqn{T} is the number of timepoints.
#
#  If GIFTI or \code{"gifti"}, the input can also be a length two list,
#  where the first list element is a length \eqn{N} vector for the left
#  hemisphere and the second list element is a length \eqn{N} vector for the
#  right hemisphere.
#'
#'  If \code{BOLD2} is provided it must be in the same format as \code{BOLD};
#'  \code{BOLD} will be the test data and \code{BOLD2} will be the retest data.
#'  \code{BOLD2} should be the same length as \code{BOLD} and have the same
#'  subjects in the same order. If \code{BOLD2} is not provided, \code{BOLD}
#'  will be split in half; the first half will be the test data and the second
#'  half will be the retest data.
#' @param GICA Group ICA maps in a format compatible with \code{BOLD}. Can also
#'  be a (vectorized) numeric matrix (\eqn{V \times Q}) no matter the format of
#'  \code{BOLD}. Its columns will be centered.
#' @param inds Numeric indices of the group ICs to include in the template. If
#'  \code{NULL}, use all group ICs (default).
#'
#'  If \code{inds} is provided, the ICs not included will be removed after
#'  calculating dual regression, not before. This is because removing the ICs
#'  prior to dual regression would leave unmodelled signals in the data, which
#'  could bias the templates.
#' @inheritParams scale_Param
#' @param scale_sm_surfL,scale_sm_surfR,scale_sm_FWHM Only applies if
#'  \code{scale=="local"} and \code{BOLD} represents surface data (CIFTI or
#'  GIFTI). To smooth the standard deviation estimates used for local scaling,
#'  provide the surface geometries along which to smooth as GIFTI geometry files
#'  or \code{"surf"} objects, as well as the smoothing FWHM (default: \code{2}).
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
#' @param brainstructures Only applies if the entries of \code{BOLD} are CIFTI
#'  file paths. This is a character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
#'  cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Can also be \code{"all"} (obtain all three brain structures).
#'  Default: \code{c("left","right")} (cortical surface only).
#' @param mask Required if and only if the entries of \code{BOLD} are NIFTI
#'  file paths or \code{"nifti"} objects. This is a brain map formatted as a
#'  binary array of the same spatial dimensions as the fMRI data, with
#'  \code{TRUE} corresponding to in-mask voxels.
#' @param keep_DR Keep the DR estimates? If \code{FALSE} (default), do not save
#'  the DR estimates and only return the templates. If \code{TRUE}, the DR
#'  estimates are returned too. If a single file path, save the DR estimates as
#'  an RDS file at that location rather than returning them.
#   [TO DO] If a list of two vectors of file paths with the same lengths as
#   \code{BOLD}, save the DR estimates as individual files at these locations in
#   the appropriate format (CIFTI, NIFTI, or RDS files, depending on \code{BOLD}).
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
#' @param FC Include the functional connectivity template? Default: \code{FALSE}
#'  (not fully supported yet.)
#' @param varTol Tolerance for variance of each data location. For each scan,
#'  locations which do not meet this threshold are masked out of the analysis.
#'  Default: \code{1e-6}. Variance is calculated on the original data, before
#'  any normalization.
#' @param maskTol For computing the dual regression results for each subject:
#'  tolerance for number of locations masked out due to low
#'  variance or missing values. If more than this many locations are masked out,
#'  a subject is skipped without calculating dual regression. \code{maskTol}
#'  can be specified either as a proportion of the number of locations (between
#'  zero and one), or as a number of locations (integers greater than one).
#'  Default: \code{.1}, i.e. up to 10 percent of locations can be masked out.
#'
#'  If \code{BOLD2} is provided, masks are calculated for both scans and then
#'  the intersection of the masks is used, for each subject.
#' @param missingTol For computing the variance decomposition across all subjects:
#'  tolerance for number of subjects masked out due to low variance or missing
#'  values at a given location. If more than this many subjects are masked out,
#'  the location's value will be \code{NA} in the templates. \code{missingTol}
#'  can be specified either as a proportion of the number of locations (between
#'  zero and one), or as a number of locations (integers greater than one).
#'  Default: \code{.1}, i.e. up to 10 percent of subjects can be masked out
#'  at a given location.
#' @param usePar,wb_path Parallelize the DR computations over subjects? Default:
#'  \code{FALSE}. Can be the number of cores to use or \code{TRUE}, which will
#'  use the number on the PC minus two. If the input data is in CIFTI format, the
#'  \code{wb_path} must also be provided.
#' @param verbose Display progress updates? Default: \code{TRUE}.
#'
#' @importFrom stats cov quantile
#' @importFrom fMRItools is_1 is_integer is_posNum colCenter unmask_mat infer_format_ifti_vec
#' @importFrom abind abind
#'
#' @return A list: the \code{template} and \code{var_decomp} with entries in
#'  matrix format; the \code{mask} of locations without template values due to
#'  too many low variance or missing values; the function \code{params} such as
#'  the type of scaling and detrending performed; the {dat_struct} which can be
#'  used to convert \code{template} and \code{var_decomp} to \code{"xifti"} or
#'  \code{"nifti"} objects if the \code{BOLD} format was CIFTI or NIFTI data;
#'  and \code{DR} if \code{isTRUE(keep_DR)}.
#'
#'  Use \code{summary} to print a description of the template results, and
#'  for CIFTI-format data use \code{plot} to plot the template mean and variance
#'  estimates. Use \code{\link{export_template}} to save the templates to
#'  individual RDS, CIFTI, or NIFTI files (depending on the \code{BOLD} format).
#' @export
#'
#' @examples
#' nT <- 30
#' nV <- 400
#' nQ <- 7
#' mU <- matrix(rnorm(nV*nQ), nrow=nV)
#' mS <- mU %*% diag(seq(nQ, 1)) %*% matrix(rnorm(nQ*nT), nrow=nQ)
#' BOLD <- list(B1=mS, B2=mS, B3=mS)
#' BOLD <- lapply(BOLD, function(x){x + rnorm(nV*nT, sd=.05)})
#' GICA <- mU
#' estimate_template(BOLD=BOLD, GICA=mU)
#'
#' \dontrun{
#'  estimate_template(
#'    run1_cifti_fnames, run2_cifti_fnames,
#'    gICA_cifti_fname, brainstructures="all",
#'    scale="local", detrend_DCT=7, Q2=NULL, varTol=10
#'  )
#' }
estimate_template <- function(
  BOLD, BOLD2=NULL,
  GICA, inds=NULL,
  scale=c("global", "local", "none"),
  scale_sm_surfL=NULL, scale_sm_surfR=NULL, scale_sm_FWHM=2,
  detrend_DCT=0,
  center_Bcols=FALSE, normA=FALSE,
  Q2=0, Q2_max=NULL,
  brainstructures=c("left","right"), mask=NULL,
  keep_DR=FALSE,
  FC=FALSE,
  varTol=1e-6, maskTol=.1, missingTol=.1,
  usePar=FALSE, wb_path=NULL,
  verbose=TRUE) {

  # Check arguments ------------------------------------------------------------

  # Simple argument checks.
  if (is.null(scale) || isFALSE(scale)) { scale <- "none" }
  if (isTRUE(scale)) {
    warning(
      "Setting `scale='global'`. Use `'global'` or `'local'` ",
      "instead of `TRUE`, which has been deprecated."
    )
    scale <- "global"
  }
  scale <- match.arg(scale, c("global", "local", "none"))
  stopifnot(fMRItools::is_1(scale_sm_FWHM, "numeric"))
  if (isFALSE(detrend_DCT)) { detrend_DCT <- 0 }
  stopifnot(fMRItools::is_integer(detrend_DCT, nneg=TRUE))
  stopifnot(fMRItools::is_1(center_Bcols, "logical"))
  stopifnot(fMRItools::is_1(normA, "logical"))
  if (!is.null(Q2)) { # Q2_max checked later.
    stopifnot(fMRItools::is_integer(Q2) && (Q2 >= 0))
  }
  stopifnot(fMRItools::is_1(FC, "logical"))
  if (isTRUE(FC)) { warning("FC template is still under development.") }
  stopifnot(fMRItools::is_1(varTol, "numeric"))
  if (varTol < 0) { cat("Setting `varTol=0`."); varTol <- 0 }
  stopifnot(fMRItools::is_posNum(maskTol))
  stopifnot(fMRItools::is_posNum(missingTol))
  stopifnot(fMRItools::is_1(verbose, "logical"))
  real_retest <- !is.null(BOLD2)

  # `keep_DR`
  if (is.logical(keep_DR)) {
    stopifnot(length(keep_DR)==1)
  } else {
    if (is.character(keep_DR)) {
      stopifnot(length(keep_DR)==1)
      if (!dir.exists(dirname(keep_DR))) { stop('Directory part of `keep_DR` does not exist.') }
      if (!endsWith(keep_DR, ".rds")) { keep_DR <- paste0(keep_DR, ".rds") }
    } else if (is.list(keep_DR)) {
      stop("Not supported: `keep_DR` must be `TRUE`, `FALSE`, or a single file path.")
      # [TO DO]
      # if (length(keep_DR) != 2) {
      #   stop("If `keep_DR` is a list it must have two entries, each being a vector of file paths the same length as `BOLD`.")
      # }
      # if (length(keep_DR[[1]]) != nN || length(keep_DR[[2]]) != nN) {
      #   stop("If `keep_DR` is a list it must have two entries, each being a vector of file paths the same length as `BOLD`.")
      # }
      # if (!all(dir.exists(dirname(do.call(c, keep_DR))))) { stop('At least one directory part of `keep_DR` does not exist.') }
    }
  }

  # `usePar`
  if (!isFALSE(usePar)) {
    parPkgs <- c("parallel", "doParallel")
    parPkgs_missing <- !vapply(parPkgs, function(x){requireNamespace(x, quietly=TRUE)}, FALSE)
    if (any(parPkgs_missing)) {
      if (all(parPkgs_missing)) {
        stop(paste0(
          "Packages `parallel` and `doParallel` needed ",
          "for `usePar` to loop over subjects. Please install them."), call.=FALSE
        )
      } else {
        stop(paste0(
          "Package `", parPkgs[parPkgs_missing], "` needed ",
          "for `usePar` to loop over subjects. Please install it."), call.=FALSE
        )
      }
    }

    cores <- parallel::detectCores()
    if (isTRUE(usePar)) { nCores <- cores[1] - 2 } else { nCores <- usePar; usePar <- TRUE }
    if (nCores < 2) {
      usePar <- FALSE
    } else {
      cluster <- parallel::makeCluster(nCores, outfile="")
      doParallel::registerDoParallel(cluster)
    }
  }

  # `BOLD` and `BOLD2` ---------------------------------------------------------
  # Determine the format of `BOLD` and `BOLD2`.
  format <- infer_format_ifti_vec(BOLD)[1]
  if (real_retest) {
    format2 <- infer_format_ifti_vec(BOLD2)[1]
    if (format2 != format) {
      stop("`BOLD` format is ", format, ", but `BOLD2` format is ", format2, ".")
    }
  }
  FORMAT <- get_FORMAT(format)
  FORMAT_extn <- switch(FORMAT,
    CIFTI=".dscalar.nii",
    GIFTI=".func.gii",
    NIFTI=".nii",
    MATRIX=".rds"
  )
  nN <- length(BOLD)

  check_req_ifti_pkg(FORMAT)

  # Ensure `BOLD2` is the same length.
  if (real_retest) {
    if (length(BOLD) != length(BOLD2)) {
      stop("`BOLD2` represents corresponding retest data for `BOLD`, so it must have the same length as `BOLD`.")
    }
  }

  # If BOLD (and BOLD2) is a CIFTI, GIFTI, NIFTI, or RDS file, check that the file paths exist.
  if (format %in% c("CIFTI", "GIFTI", "NIFTI", "RDS")) {
    missing_BOLD <- !file.exists(BOLD)
    if (all(missing_BOLD)) stop('The files in `BOLD` do not exist.')
    if (real_retest) {
      missing_BOLD2 <- !file.exists(BOLD2)
      if (all(missing_BOLD2)) stop('The files in `BOLD2` do not exist.')
      # determine pairs with at least one missing scan
      missing_BOLD <- missing_BOLD | missing_BOLD2
      rm(missing_BOLD2)
      if (all(missing_BOLD)) stop('Files in `BOLD` and/or `BOLD2` are missing such that no complete pair of data exists.')
    }
    if (any(missing_BOLD)) {
      if (real_retest) {
        warning('There are ', missing_BOLD, ' pairs of `BOLD` and `BOLD2` with at least one non-existent scan. These pairs will be excluded from template estimation.')
        BOLD <- BOLD[!missing_BOLD]
        BOLD2 <- BOLD[!missing_BOLD]
      } else {
        warning('There are ', missing_BOLD, ' scans in `BOLD` that do not exist. These scans will be excluded from template estimation.')
        BOLD <- BOLD[!missing_BOLD]
      }
    }
  }

  # Check `scale_sm_FWHM`
  if (scale_sm_FWHM !=0 && FORMAT %in% c("NIFTI", "MATRIX")) {
    if (scale_sm_FWHM==2) {
      cat("Setting `scale_sm_FWHM == 0`.\n")
    } else {
      if (FORMAT == "NIFTI") {
        warning( "Setting `scale_sm_FWHM == 0` (Scale smoothing not available for volumetric data.).\n")
      } else {
        warning( "Setting `scale_sm_FWHM == 0` (Scale smoothing not available for data matrices: use CIFTI/GIFTI files.).\n")
      }
    }
    scale_sm_FWHM <- 0
  }

  # [TO DO]: Mysteriously, MATRIX FORMAT is not working with parallel.
  if (usePar && format=="RDS") { stop("Parallel computation not working with RDS file input. Working on this!") }

  # `GICA` ---------------------------------------------------------------------
  # Conver `GICA` to a numeric data matrix or array.
  if (FORMAT == "CIFTI") {
    if (is.character(GICA)) { GICA <- ciftiTools::read_cifti(GICA, brainstructures=brainstructures) }
    if (ciftiTools::is.xifti(GICA, messages=FALSE)) {
      xii1 <- ciftiTools::select_xifti(GICA, 1) # for formatting output
      GICA <- as.matrix(GICA)
    } else {
      # Get `xii1` from first data entry.
      xii1 <- BOLD[[1]]
      if (is.character(xii1)) {
        xii1 <- ciftiTools::read_cifti(xii1, brainstructures=brainstructures, idx=1)
      }
      xii1 <- ciftiTools::convert_xifti(ciftiTools::select_xifti(xii1, 1), "dscalar")
    }
    stopifnot(is.matrix(GICA))
  } else if (FORMAT == "GIFTI") {
    if (is.character(GICA)) { GICA <- gifti::readgii(GICA) }
    ghemi <- GICA$file_meta["AnatomicalStructurePrimary"]
    if (!(ghemi %in% c("CortexLeft", "CortexRight"))) {
      stop("AnatomicalStructurePrimary metadata missing or invalid for GICA.")
    }
    ghemi <- switch(ghemi, CortexLeft="left", CortexRight="right")
    GICA <- do.call(cbind, GICA$data)
  } else if (FORMAT == "NIFTI") {
    if (is.character(GICA)) { GICA <- RNifti::readNifti(GICA) }
    stopifnot(length(dim(GICA)) > 1)
  } else if (FORMAT == "MATRIX") {
    if (is.character(GICA)) { GICA <- readRDS(GICA) }
    stopifnot(is.matrix(GICA))
  }
  nQ <- dim(GICA)[length(dim(GICA))]

  # `inds`.
  if (!is.null(inds)) {
    if (!all(inds %in% seq(nQ))) stop('Invalid entries in inds argument.')
    nL <- length(inds)
  } else {
    inds <- seq(nQ)
    nL <- nQ
  }

  # [TO DO]: NA in GICA?

  # `mask` ---------------------------------------------------------------------
  # Get `mask` as a logical array.
  # Check `GICA` and `mask` dimensions match.
  # Append NIFTI header from GICA to `mask`.
  # Vectorize `GICA`.
  if (FORMAT == "NIFTI") {
    if (is.null(mask)) { stop("`mask` is required.") }
    if (is.character(mask)) { mask <- RNifti::readNifti(mask); mask <- array(as.logical(mask), dim=dim(mask)) }
    if (dim(mask)[length(dim(mask))] == 1) { mask <- array(mask, dim=dim(mask)[length(dim(mask))-1]) }
    if (is.numeric(mask)) {
      cat("Coercing `mask` to a logical array.\n")
      mask <- array(as.logical(mask), dim=dim(mask))
    }
    nI <- dim(mask); nV <- sum(mask)
    stopifnot(length(dim(GICA)) %in% c(2, length(nI)+1))
    if (length(dim(GICA)) == length(nI)+1) {
      if (length(dim(GICA)) != 2) {
        stopifnot(all(dim(GICA)[seq(length(dim(GICA))-1)] == nI))
      }
      # Append NIFTI header.
      mask <- RNifti::asNifti(array(mask, dim=c(dim(mask), 1)), reference=GICA)
      # Vectorize `GICA`.
      if (all(dim(GICA)[seq(length(dim(GICA))-1)] == nI)) {
        GICA <- matrix(GICA[rep(as.logical(mask), nQ)], ncol=nQ)
        stopifnot(nrow(GICA) == nV)
      }
    }
  } else {
    if (!is.null(mask)) {
      warning("Ignoring `mask`, which is only applicable to NIFTI data.")
      mask <- NULL
    }
    nI <- nV <- nrow(GICA)
  }

  # Center `GICA` columns.
  center_Gcols <- TRUE
  if (center_Gcols) { GICA <- fMRItools::colCenter(GICA) }

  # Print summary of data ------------------------------------------------------
  format2 <- if (format == "data") { "numeric matrix" } else { format }
  if (verbose) {
    cat("Data input format:             ", format2, "\n")
    cat('Number of data locations:      ', nV, "\n")
    if (FORMAT == "NIFTI") {
      cat("Unmasked dimensions:           ", paste(nI, collapse=" x "), "\n")
    }
    if (FORMAT == "GIFTI") {
      cat("Cortex Hemisphere:             ", ghemi, "\n")
    }
    cat('Number of original group ICs:  ', nQ, "\n")
    cat('Number of template ICs:        ', nL, "\n")
    cat('Number of training subjects:   ', nN, "\n")
  }

  # Process each scan ----------------------------------------------------------
  nM <- 2

  if (usePar) {
    # Check `foreach`.
    if (!requireNamespace("foreach", quietly = TRUE)) {
      stop(
        "Package \"foreach\" needed to parallel loop over scans. Please install it.",
        call. = FALSE
      )
    }

    # Loop over subjects.
    `%dopar%` <- foreach::`%dopar%`
    q <- foreach::foreach(ii = seq(nN)) %dopar% {
      if (FORMAT=="CIFTI" || FORMAT=="GIFTI") {
        # Load the workbench.
        if (is.null(wb_path)) {
          stop("`wb_path` is required for parallel computation.")
        }
        requireNamespace("ciftiTools")
        ciftiTools::ciftiTools.setOption("wb_path", wb_path)
      }

      # Initialize output.
      out <- list(DR=array(NA, dim=c(nM, 1, nL, nV)))
      if (FC) { out$FC <- array(NA, dim=c(nM, 1, nL, nL)) }

      # Dual regression.
      if(verbose) { cat(paste0(
        '\nReading and analyzing data for subject ', ii,' of ', nN, ".\n"
      )) }
      if (real_retest) { B2 <- BOLD2[[ii]] } else { B2 <- NULL }
      DR_ii <- dual_reg2(
        BOLD[[ii]], BOLD2=B2,
        format=format,
        GICA=GICA,
        keepA=FC,
        center_Bcols=center_Bcols,
        scale=scale,
        scale_sm_surfL=scale_sm_surfL, scale_sm_surfR=scale_sm_surfR,
        scale_sm_FWHM=scale_sm_FWHM,
        detrend_DCT=detrend_DCT,
        normA=normA,
        Q2=Q2, Q2_max=Q2_max,
        brainstructures=brainstructures, mask=mask,
        varTol=varTol, maskTol=maskTol,
        verbose=verbose
      )

      # Add results if this subject was not skipped.
      # (Subjects are skipped if too many locations are masked out.)
      if (!is.null(DR_ii)) {
        out$DR[1,,,] <- DR_ii$test$S[inds,]
        out$DR[2,,,] <- DR_ii$retest$S[inds,]
        if(FC) {
          out$FC[1,,,] <- cov(DR_ii$test$A[,inds])
          out$FC[2,,,] <- cov(DR_ii$retest$A[,inds])
        }
      } else {
        if(verbose) { cat(paste0(
          '\nSubject ', ii,' of ', nN, " was skipped (too many masked locations).\n"
        )) }
      }
      out
    }

    # Aggregate.
    DR0 <- abind::abind(lapply(q, `[[`, "DR"), along=2)
    if (FC) { FC0 <- abind::abind(lapply(q, `[[`, "FC"), along=2) }

    doParallel::stopImplicitCluster()

  } else {
    # Initialize output.
    DR0 <- array(NA, dim=c(nM, nN, nL, nV)) # measurements by subjects by components by locations
    if(FC) FC0 <- array(NA, dim=c(nM, nN, nL, nL)) # for functional connectivity template

    for (ii in seq(nN)) {
      if(verbose) { cat(paste0(
        '\nReading and analyzing data for subject ', ii,' of ', nN, '.\n'
      )) }
      if (real_retest) { B2 <- BOLD2[[ii]] } else { B2 <- NULL }

      DR_ii <- dual_reg2(
        BOLD[[ii]], BOLD2=B2,
        format=format,
        GICA=GICA,
        keepA=FC,
        center_Bcols=center_Bcols,
        scale=scale,
        scale_sm_surfL=scale_sm_surfL, scale_sm_surfR=scale_sm_surfR,
        scale_sm_FWHM=scale_sm_FWHM,
        detrend_DCT=detrend_DCT,
        normA=normA,
        Q2=Q2, Q2_max=Q2_max,
        brainstructures=brainstructures, mask=mask,
        varTol=varTol, maskTol=maskTol,
        verbose=verbose
      )

      # Add results if this subject was not skipped.
      # (Subjects are skipped if too many locations are masked out.)
      if (!is.null(DR_ii)) {
        DR0[1,ii,,] <- DR_ii$test$S[inds,]
        DR0[2,ii,,] <- DR_ii$retest$S[inds,]
        if(FC) {
          FC0[1,ii,,] <- cov(DR_ii$test$A[,inds])
          FC0[2,ii,,] <- cov(DR_ii$retest$A[,inds])
        }
      } else {
        if(verbose) { cat(paste0(
          '\tSubject ', ii,' was skipped (too many masked locations).\n'
        )) }
      }
    }
    rm(DR_ii)
  }

  # Aggregate results, and compute templates. ----------------------------------

  # Mask out locations for which too many subjects do not have data
  nVm <- nV # Number of locations after data masking.
  if (missingTol < 1) { missingTol <- missingTol * nN }
  # Assume is.na(DR0[mm,,,]) is the same for all mm
  # Assume is.na(DR0[,,cc,]) is the same for all cc
  NA_counts <- colSums(is.na(DR0[1,,1,]))
  maskAll <- NA_counts < missingTol
  use_mask <- !all(maskAll)
  if (use_mask) {
    if (all(!maskAll)) { stop(
      "No locations meet the minimum number of subjects with data (",
      round(missingTol, digits=3), "). Check `maskTol` and `missingTol`."
    ) }
    DR0 <- DR0[,,,maskAll,drop=FALSE]
    nVm <- sum(maskAll)
  }
  # Note that `NA` values may still exist in `DR0`.

  # Vectorize components and locations
  DR0 <- array(DR0, dim=c(nM, nN, nL*nVm))

  if (verbose) { cat("\nCalculating template.\n") }
  # Estimate the mean and variance templates.
  # Also obtain the variance decomposition.
  x <- estimate_template_from_DR(DR0, c(nL, nVm))
  template <- lapply(x$template, t)
  var_decomp <- lapply(x$var_decomp, t)
  rm(x)

  # Unmask the data matrices.
  if (use_mask) {
    template <- lapply(template, fMRItools::unmask_mat, mask=maskAll)
    var_decomp <- lapply(var_decomp, fMRItools::unmask_mat, mask=maskAll)
  }

  # Estimate FC template
  if(FC){ template$FC <- estimate_template_FC(FC0) }

  # Format result ---------------------------------------------------
  # Keep DR
  if (!isFALSE(keep_DR)) {
    DR0 <- array(DR0, dim=c(nM, nN, nL, nVm)) # Undo vectorize
    if (use_mask) {
      DR0temp <- array(NA, dim=c(dim(DR0)[seq(3)], length(maskAll)))
      DR0temp[,,,maskAll] <- DR0
      DR0 <- DR0temp
    }
    if (is.character(keep_DR)) {
      saveRDS(DR0, keep_DR)
      keep_DR <- FALSE # no longer need it.
    } else if (!isTRUE(keep_DR)) {
      warning("`keep_DR` should be `TRUE`, `FALSE`, or a file path. Using `FALSE`.")
      keep_DR <- FALSE
    }
  }
  if (!keep_DR) { rm(DR0) }

  # Params, formatted as length-one character vectors to be compatible with
  #   putting in "xifti" metadata
  indsp <- if (all(seq(nQ) %in% inds)) { paste("all", nQ) } else { inds }
  tparams <- list(
    FC=FC,
    num_subjects=nN, num_visits=nM,
    inds=indsp, center_Bcols=center_Bcols,
    scale=scale, scale_sm_FWHM=scale_sm_FWHM,
    detrend_DCT=detrend_DCT, normA=normA,
    Q2=Q2, Q2_max=Q2_max,
    brainstructures=brainstructures,
    varTol=varTol, maskTol=maskTol, missingTol=missingTol,
    pseudo_retest=!real_retest
  )
  tparams <- lapply(
    tparams,
    function(x) {
      if (is.null(x)) { x <- "NULL"};
      paste0(as.character(x), collapse=" ")
    }
  )

  # If masking was performed, and if verbose, report the NA count in templates.
  # If no masking was performed, remove unnecessary metadata.
  if (use_mask) {
    if (verbose) {
      cat("Number of template locations with missing values:", sum(!maskAll), "\n")
    }
  } else {
    maskAll <- NULL
  }

  dat_struct <- switch(FORMAT,
    CIFTI = ciftiTools::newdata_xifti(xii1, 0),
    GIFTI = list(hemisphere=ghemi),
    NIFTI = mask,
    MATRIX = NULL
  )

  # Results list.
  result <- list(
    template=template,
    var_decomp=var_decomp,
    mask=maskAll,
    params=tparams,
    dat_struct=dat_struct
  )

  # Add DR if applicable.
  if (keep_DR) { result$DR <- DR0 }

  # Return results.
  class(result) <- paste0("template.", tolower(FORMAT))
  result
}

#' @rdname estimate_template
#'
estimate_template.cifti <- function(
  BOLD, BOLD2=NULL,
  GICA, inds=NULL,
  scale=c("global", "local", "none"),
  scale_sm_surfL=NULL, scale_sm_surfR=NULL, scale_sm_FWHM=2,
  detrend_DCT=0,
  center_Bcols=FALSE, normA=FALSE,
  Q2=0, Q2_max=NULL,
  brainstructures=c("left","right"),
  keep_DR=FALSE,
  FC=FALSE,
  varTol=1e-6, maskTol=.1, missingTol=.1,
  usePar=FALSE, wb_path=NULL,
  verbose=TRUE) {

  estimate_template(
    BOLD=BOLD, BOLD2=BOLD2,
    GICA=GICA, inds=inds,
    scale=scale, scale_sm_surfL=scale_sm_surfL, scale_sm_surfR=scale_sm_surfR,
    scale_sm_FWHM=scale_sm_FWHM,
    detrend_DCT=detrend_DCT,
    center_Bcols=center_Bcols, normA=normA,
    Q2=Q2, Q2_max=Q2_max,
    brainstructures=brainstructures,
    keep_DR=keep_DR,
    FC=FC,
    varTol=varTol, maskTol=maskTol, missingTol=missingTol,
    usePar=usePar, wb_path=wb_path,
    verbose=verbose
  )
}

#' @rdname estimate_template
#'
estimate_template.gifti <- function(
  BOLD, BOLD2=NULL,
  GICA, inds=NULL,
  scale=c("global", "local", "none"),
  scale_sm_surfL=NULL, scale_sm_surfR=NULL, scale_sm_FWHM=2,
  detrend_DCT=0,
  center_Bcols=FALSE, normA=FALSE,
  Q2=0, Q2_max=NULL,
  brainstructures=c("left","right"),
  keep_DR=FALSE,
  FC=FALSE,
  varTol=1e-6, maskTol=.1, missingTol=.1,
  usePar=FALSE, wb_path=NULL,
  verbose=TRUE) {

  estimate_template(
    BOLD=BOLD, BOLD2=BOLD2,
    GICA=GICA, inds=inds,
    scale=scale, scale_sm_surfL=scale_sm_surfL, scale_sm_surfR=scale_sm_surfR,
    scale_sm_FWHM=scale_sm_FWHM,
    detrend_DCT=detrend_DCT,
    center_Bcols=center_Bcols, normA=normA,
    Q2=Q2, Q2_max=Q2_max,
    brainstructures=brainstructures,
    keep_DR=keep_DR,
    FC=FC,
    varTol=varTol, maskTol=maskTol, missingTol=missingTol,
    usePar=usePar, wb_path=wb_path,
    verbose=verbose
  )
}

#' @rdname estimate_template
#'
estimate_template.nifti <- function(
  BOLD, BOLD2=NULL,
  GICA, inds=NULL,
  scale=c("global", "local", "none"),
  detrend_DCT=0,
  center_Bcols=FALSE, normA=FALSE,
  Q2=0, Q2_max=NULL,
  mask=NULL,
  keep_DR=FALSE,
  FC=FALSE,
  varTol=1e-6, maskTol=.1, missingTol=.1,
  usePar=FALSE, wb_path=NULL,
  verbose=TRUE) {

  estimate_template(
    BOLD=BOLD, BOLD2=BOLD2,
    GICA=GICA, inds=inds,
    scale=scale,
    detrend_DCT=detrend_DCT,
    center_Bcols=center_Bcols, normA=normA,
    Q2=Q2, Q2_max=Q2_max,
    mask=mask,
    keep_DR=keep_DR,
    FC=FC,
    varTol=varTol, maskTol=maskTol, missingTol=missingTol,
    usePar=usePar, wb_path=wb_path,
    verbose=verbose
  )
}

