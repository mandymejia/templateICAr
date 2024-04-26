#' TR
#' 
#' @param TR The temporal resolution of the data, i.e. the time between volumes,
#'  in seconds. \code{TR} is required for detrending with \code{hpf}.
#' 
#' @name TR_param
#' @keywords internal
NULL

#' hpf
#' 
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
#' @name hpf_param
#' @keywords internal
NULL

#' varTol
#' 
#' @param varTol Tolerance for variance of each data location. For each scan,
#'  locations which do not meet this threshold are masked out of the analysis.
#'  Default: \code{1e-6}. Variance is calculated on the original data, before
#'  any normalization. Set to \code{0} to avoid removing locations due to 
#'  low variance.
#' 
#' @name varTol_Param
#' @keywords internal
NULL

#' scale
#' 
#' @param scale \code{"local"} (default), \code{"global"}, or \code{"none"}.
#'  Local scaling will divide each data location's time series by its estimated 
#'  standard deviation. Global scaling will divide the entire data matrix by the 
#'  mean image standard deviation (\code{mean(sqrt(rowVars(BOLD)))}). 
#' 
#' @name scale_Param
#' @keywords internal
NULL

#' GSR
#' 
#' @param GSR Center BOLD across columns (each image)? This
#'  is equivalent to performing global signal regression. Default: 
#'  \code{FALSE}. 
#' 
#' @name GSR_Param
#' @keywords internal
NULL