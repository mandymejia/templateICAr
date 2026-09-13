#' Dual Regression
#'
#' @template superseded-by-BayesBrainMap
#' @param BOLD Subject-level fMRI data matrix (\eqn{V \times T}). Rows will be
#'  centered.
#' @param GICA Group-level independent components (\eqn{V \times Q})
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
#' @param hpf,lpf The frequencies at which to apply a highpass filter or lowpass
#'  filter to the data during pre-processing, in Hertz. Set either to 
#'  \code{NULL} to disable filtering. Default: \code{0.01} Hertz for the 
#'  highpass filter, and \code{NULL} for the lowpass filter.
#'
#'  The highpass filter serves to detrend the data, since low-frequency
#'  variance is associated with noise. Highpass filtering is accomplished by
#'  nuisance regression of discrete cosine transform (DCT) bases.
#' 
#'  The lowpass filter removes high-frequency variance also thought to be
#'  associated with non-neuronal noise. 
#'
#'  Note the \code{TR} argument is required for temporal filtering. If
#'  \code{TR} is not provided, \code{hpf} and \code{lpf} will be ignored.
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
  scale=c("local", "global", "none"), scale_sm_xifti=NULL, scale_sm_FWHM=2,
  TR=NULL, hpf=.01, lpf=NULL,
  GSR=FALSE){

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
  scale <- match.arg(scale, c("local", "global", "none"))
  if (!is.null(scale_sm_xifti)) { stopifnot(ciftiTools::is.xifti(scale_sm_xifti)) }
  stopifnot(is.numeric(scale_sm_FWHM) && length(scale_sm_FWHM)==1)

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

  # Center each voxel timecourse. Do not center the image at each timepoint unless GSR = TRUE.
  # Standardize scale if `scale`, and detrend if `hpf>0`.
  # Transpose it: now `BOLD` is TxV.
  BOLD <- t(norm_BOLD(
    BOLD, center_rows=TRUE, center_cols=GSR,
    scale=scale, scale_sm_xifti=scale_sm_xifti, scale_sm_FWHM=scale_sm_FWHM,
    TR=TR, hpf=hpf, lpf=lpf
  ))

  # Center each group IC across space. (Used to be a function argument.)
  GICA <- colCenter(GICA)

  # Estimate A (IC timeseries).
  # We need to center `BOLD` across space because the linear model has no intercept.
  A <- ((BOLD - rowMeans(BOLD, na.rm=TRUE)) %*% GICA) %*% chol2inv(chol(crossprod(GICA)))

  # Center each subject IC timecourse across time.
  # (Redundant. Since BOLD is column-centered, A is already column-centered.)
  # A <- colCenter(A)

  # Normalize each subject IC timecourse to constrain the ICA. (Used to be a function argument.)
  A <- scale(A)

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

  # Estimate S (IC maps).
  # Don't worry about the intercept: `BOLD` and `A` are centered across time.
  S <- solve(a=crossprod(A), b=crossprod(A, BOLD))

  # Re-estimate A (IC timeseries) based on the subject-level IC maps
  # We need to center `BOLD` across space because the linear model has no intercept.
  S_ctr <- colCenter(t(S))
  A2 <- ((BOLD - rowMeans(BOLD, na.rm=TRUE)) %*% S_ctr) %*% chol2inv(chol(crossprod(S_ctr)))

  #return result
  list(S = S, A = A, A2 = A2)
}