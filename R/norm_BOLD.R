#' Normalize BOLD data
#'
#' Center the data across space and/or time, detrend, and scale, in that order.
#'  For dual regression, row centering is required and column centering is not
#'  recommended. Scaling and detrending depend on the user preference.
#'
#' @template superseded-by-BayesBrainMap
#' @param BOLD fMRI numeric data matrix (\eqn{V \times T})
#' @param center_rows,center_cols Center BOLD data across rows (each data
#'  location's time series) or columns (each time point's image)? Default:
#'  \code{TRUE} for row centering, and \code{FALSE} for column centering.
#' @param scale \code{"global"} (default), \code{"local"}, or \code{"none"}.
#'  Global scaling will divide the entire data matrix by the mean image standard
#'  deviation (\code{mean(sqrt(rowVars(BOLD)))}). Local scaling will divide each
#'  data location's time series by its estimated standard deviation.
#' @param scale_sm_xifti,scale_sm_FWHM Only applies if \code{scale=="local"} and
#'  \code{BOLD} represents CIFTI-format data. To smooth the standard deviation
#'  estimates used for local scaling, provide a \code{"xifti"} object with data
#'  locations in alignment with \code{BOLD}, as well as the smoothing FWHM
#'  (default: \code{2}). If no \code{"xifti"} object is provided (default), do
#'  not smooth.
#' @param scale_sm_xifti_mask For local scaling with smoothing, the data must
#'  be unmasked to be mapped back to the surface. So if the data are masked,
#'  provide the mask here.
#' @param TR The temporal resolution of the data, i.e. the time between volumes,
#'  in seconds. \code{TR} is required for detrending with \code{hpf}.
#' @param hpf,lpf The frequencies at which to apply temporal filtering to the
#'  data during pre-processing, in Hertz. Set either to \code{NULL} to disable.
#'  Default: \code{0.01} Hz highpass filter, and \code{NULL} for the lowpass 
#'  filter (disabled). Filtering is accomplished by nuisance regression of
#'  discrete cosine transform (DCT) bases.
#' 
#'  The highpass filter serves to detrend the data, since low-frequency 
#'  variance is associated with noise. The lowpass filter removes high-frequency
#'  variance, which is also thought to be from non-neuronal noise.
#' 
#'  Note the \code{TR} argument is required for temporal filtering. If
#'  \code{TR} is not provided, \code{hpf} and \code{lpf} will be ignored.
#'
#' @return Normalized BOLD data matrix (\eqn{V \times T})
#'
#' @export
#'
#' @importFrom fMRItools nuisance_regression temporal_filter
#'
norm_BOLD <- function(
  BOLD, center_rows=TRUE, center_cols=FALSE,
  scale=c("local", "global", "none"), scale_sm_xifti=NULL, scale_sm_FWHM=2,
  scale_sm_xifti_mask=NULL,
  TR=NULL, hpf=.01, lpf=NULL){

  nT <- ncol(BOLD)
  nV <- nrow(BOLD)
  if (nT > nV) { warning('More time points than voxels. Are you sure?') }

  stopifnot(is.logical(center_rows) && length(center_rows)==1)
  stopifnot(is.logical(center_cols) && length(center_cols)==1)
  if (is.null(scale) || isFALSE(scale)) { scale <- "none" }
  if (isTRUE(scale)) {
    warning(
      "Setting `scale='global'`. Use `'global'` or `'local'` ",
      "instead of `TRUE`, which has been deprecated."
    )
    scale <- "global"
  }
  scale <- match.arg(scale, c("local", "global", "none"))
  if (!is.null(scale_sm_xifti)) {
    if (!requireNamespace("ciftiTools", quietly = TRUE)) {
      stop("Package \"ciftiTools\" needed to work with CIFTI data. Please install it.", call. = FALSE)
    }
    stopifnot(ciftiTools::is.xifti(scale_sm_xifti))
    if (!is.null(scale_sm_xifti_mask)) {
      stopifnot(is.vector(scale_sm_xifti_mask) && is.logical(scale_sm_xifti_mask))
      stopifnot(sum(scale_sm_xifti_mask) == nV)
    }
  }
  stopifnot(is.numeric(scale_sm_FWHM) && length(scale_sm_FWHM)==1)
  if (length(hpf)==1 && hpf==0) { hpf <- NULL }
  if (length(lpf)==1 && lpf==Inf) { lpf <- NULL }
  if (is.null(TR)) {
    if (!is.null(hpf)) {
      if (hpf==.01) {
        message("Setting `hpf=NULL` because `TR` was not provided. Either provide `TR` or set `hpf=NULL` to disable this message.")
        hpf <- NULL
      } else {
        stop("Cannot apply `hpf` because `TR` was not provided. Either provide `TR` or set `hpf=NULL`.")
      }
    }
    if (!is.null(lpf)) {
      stop("Cannot apply `lpf` because `TR` was not provided. Either provide `TR` or set `lpf=NULL`.")
    }
  } else {
    stopifnot(is_posNum(TR))
    stopifnot(is.null(hpf) || is_posNum(hpf))
    stopifnot(is.null(lpf) || is_posNum(lpf))
  }

  # Center.
  if (center_rows || center_cols) {
    # `BOLD` is transposed twice.
    # Center each voxel time series (across time).
    if (center_rows) {
      BOLD <- t(BOLD - rowMeans(BOLD, na.rm=TRUE))
    } else {
      BOLD <- t(BOLD)
    }
    # Center each image (across space).
    if (center_cols) {
      BOLD <- t(BOLD - rowMeans(BOLD, na.rm=TRUE))
    } else {
      BOLD <- t(BOLD)
    }
  }

  # Apply the temporal filter.
  # [NOTE]: If `center_cols`, columns won't be exactly centered anymore after the filter.
  if (!is.null(hpf) || !is.null(lpf)) {
    if (!center_rows) { voxMeans <- rowMeans(BOLD, na.rm=TRUE) }
    dct <- fMRItools::temporal_filter(
      X=ncol(BOLD), TR=TR, hpf=hpf, lpf=lpf, method="DCT", verbose=FALSE
    ) # [TO DO] carry over verbose arg?
    BOLD <- nuisance_regression(BOLD, cbind(1, dct))
    if (!center_rows) { BOLD <- BOLD + voxMeans }
  }

  # Scale.
  # Get scale at each location.
  if (scale != "none") { sig <- sqrt(rowVars(BOLD, na.rm=TRUE)) }
  # Global scaling: take mean scale across all locations, and use that.
  if (scale == "global") {
    sig <- mean(sig, na.rm=TRUE)
    if (sig < 1e-8) {
      warning("Estimated scale is near zero. Skipping scaling.")
    } else {
      # Apply global scaling.
      BOLD <- BOLD / sig
    }
  # Local scaling: use estimate of scale at each location.
  } else if (scale == "local") {
    # Smooth estimates, if applicable.
    if (!is.null(scale_sm_xifti) && (scale_sm_FWHM != 0)) {
      # Check `scale_sm_xifti` is valid.
      is_masked <- !is.null(scale_sm_xifti_mask)
      # Un-mask, if applicable.
      if (is_masked) {
        sig <- c(unmask_mat(as.matrix(sig), scale_sm_xifti_mask))
        nV <- length(sig)
      }
      if (nV != nrow(scale_sm_xifti)) {
        stop("`scale_sm_xifti` not compatible with `BOLD`: different spatial dimensions.")
      }
      if (!is.null(scale_sm_xifti$meta$cifti$intent) && scale_sm_xifti$meta$cifti$intent == 3007) {
        scale_sm_xifti <- ciftiTools::convert_xifti(scale_sm_xifti, "dscalar")
      }
      # Compute and smooth the SD.
      sig <- ciftiTools::newdata_xifti(ciftiTools::select_xifti(scale_sm_xifti, 1), sig)
      sig <- ciftiTools::move_to_mwall(sig, NA)
      if (!is.null(sig$data$subcort)) {
        sub_mask <- !is.na(sig$data$subcort[,1])
        sig$data$subcort <- sig$data$subcort[sub_mask,,drop=FALSE]
        sig$meta$subcort$labels <- sig$meta$subcort$labels[sub_mask]
        sig$meta$subcort$mask[sig$meta$subcort$mask][!sub_mask] <- FALSE
      }
      sig <- ciftiTools::smooth_xifti(sig, surf_FWHM=scale_sm_FWHM, vol_FWHM=scale_sm_FWHM)
      sig <- c(as.matrix(sig))
    }
    # Apply local scaling.
    BOLD <- BOLD / sig
  }

  BOLD
}