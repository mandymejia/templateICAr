#' Multiple regression for parcel data
#'
#' @template superseded-by-BayesBrainMap
#' @param BOLD Subject-level fMRI data matrix (\eqn{V \times T}). Rows will be
#'  centered.
#' @param parc The parcellation as an integer vector.
#' @param parc_vals The parcel values (keys) in desired order, e.g.
#'  \code{sort(unique(parc))}.
#' @param GSR Center BOLD across columns (each image)? This
#'  is equivalent to performing global signal regression. Default:
#'  \code{FALSE}.
#' @param scale \code{"local"} (default), \code{"global"}, or \code{"none"}.
#'  Local scaling will divide each data location's time series by its estimated
#'  standard deviation. Global scaling will divide the entire data matrix by the
#'  mean image standard deviation (\code{mean(sqrt(rowVars(BOLD)))}).
#' @param scale_sm_xifti,scale_sm_FWHM Only applies if \code{scale=="local"} and
#'  \code{BOLD} represents CIFTI-format data. To smooth the standard deviation
#'  estimates used for local scaling, provide a \code{"xifti"} object with data
#'  locations in alignment with \code{BOLD}, as well as the smoothing FWHM
#'  (default: \code{2}). If no \code{"xifti"} object is provided (default), do
#'  not smooth.
#' @param TR The temporal resolution of the data, i.e. the time between volumes,
#'  in seconds. \code{TR} is required for detrending with \code{hpf}.
#' @param hpf The frequency at which to apply a highpass filter to the data
#'  during pre-processing, in Hertz. Default: \code{0.01} Hertz. Set to \code{0}
#'  to disable the highpass filter.
#'
#'
#'  The highpass filter serves to detrend the data, since low-frequency
#'  variance is associated with noise. Highpass filtering is accomplished by
#'  nuisance regression of discrete cosine transform (DCT) bases.
#'
#'  Note the \code{TR} argument is required for highpass filtering. If
#'  \code{TR} is not provided, \code{hpf} will be ignored.
#'
#' @return A list containing
#'  the subject-level independent components \strong{S} (\eqn{Q \times V}),
#'  and subject-level mixing matrix \strong{A} (\eqn{TxQ}).
#'
#' @importFrom matrixStats rowMedians
#' @export
#'
dual_reg_parc <- function(
  BOLD, parc, parc_vals,
  scale=c("local", "global", "none"), scale_sm_xifti=NULL, scale_sm_FWHM=2,
  TR=NULL, hpf=.01,
  GSR=FALSE){

  stopifnot(is.matrix(BOLD))
  stopifnot(is.numeric(parc))
  parc <- as.matrix(parc)
  if (is.null(scale) || isFALSE(scale)) { scale <- "none" }
  if (isTRUE(scale)) {
    warning(
      "Setting `scale='global'`. Use `'global'` or `'local'` ",
      "instead of `TRUE`, which has been deprecated."
    )
    scale <- "global"
  }
  scale <- match.arg(scale, c("local", "global", "none"))
  if (!is.null(scale_sm_xifti)) { stopifnot(ciftiTools::is.xifti(scale_sm_xifti)) }
  stopifnot(is.numeric(scale_sm_FWHM) && length(scale_sm_FWHM)==1)

  if (any(is.na(BOLD))) { stop("`NA` values in `BOLD` not supported with DR.") }
  if (any(is.na(parc))) { stop("`NA` values in `parc` not supported with DR.") }

  nV <- nrow(BOLD) #number of data locations
  nT <- ncol(BOLD) #length of timeseries
  if(nV < nT) warning('More time points than voxels. Are you sure?')
  if(nV != nrow(parc)) {
    stop('The number of voxels in dat (', nV, ') and parc (', nrow(parc), ') must match')
  }

  stopifnot(all(unique(c(parc)) %in% parc_vals))
  stopifnot(all(parc_vals %in% parc))
  nQ <- length(parc_vals) #number of parcels
  if(nQ > nV) warning('More parcels than voxels. Are you sure?')
  if(nQ > nT) warning('More parcels than time points. Are you sure?')

  # Center each voxel timecourse. Do not center the image at each timepoint.
  # Standardize scale if `scale`, and detrend if `hpf>0`.
  # Transpose it: now `BOLD` is TxV.
  BOLD <- t(norm_BOLD(
    BOLD, center_rows=TRUE, center_cols=GSR,
    scale=scale, scale_sm_xifti=scale_sm_xifti, scale_sm_FWHM=scale_sm_FWHM,
    TR=TR, hpf=hpf
  ))

  # Estimate A (parcel timeseries).
  # Neccesary to temporarily center BOLD (like in standard dual regression since there is no intecept)
  BOLD_rowMeans <- rowMeans(BOLD, na.rm=TRUE)
  A <- matrix(NA, nrow=nT, ncol=nQ)
  for (qq in seq(nQ)) {
    BOLD_qq <- BOLD[,parc==parc_vals[qq],drop=FALSE]
    A[,qq] <- matrixStats::rowMedians(BOLD_qq - BOLD_rowMeans)
  }
  rm(BOLD_rowMeans)

  # Normalize each subject parcel timecourse. (Used to be a function argument.)
  normA <- TRUE
  if (normA) { A <- scale(A) }

  # Check rank of `A`.
  A_rank <- qr(A)$rank
  if (A_rank < ncol(A)) {
    warning(
      "DR has estimated an `A` matrix that has ", ncol(A), " columns, but its rank is ", A_rank, ". ",
      "An `A` matrix that is not full rank can occur when the number of group ICs approaches the number of volumes in the subject data. ",
      "This problem can be avoided by using a group ICA with fewer components, ",
      "or by providing more volumes of data. ",
      "Continuing, but an error may occur in further calculations."
    )
  }

  # Estimate S (IC maps). VxQ
  S <- solve(a=crossprod(A), b=crossprod(A, BOLD))

  #return result
  list(S = S, A = A)
}