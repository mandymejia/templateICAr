#' detrend_DCT
#' 
#' @param detrend_DCT Detrend the data? This is an integer number of DCT bases
#'  to use for detrending. If \code{0} (default), do not detrend.
#' 
#' @name detrend_DCT_Param
#' @keywords internal
NULL

#' normA
#' 
#' @param normA Scale each IC timeseries (column of \eqn{A}) in the dual
#'  regression estimates? Default: \code{FALSE} (not recommended). Note that the
#'  product \eqn{A \times S} remains the same with either option.
#' 
#' @name normA_Param
#' @keywords internal
NULL

#' varTol
#' 
#' @param varTol Tolerance for variance of each data location. For each scan,
#'  locations which do not meet this threshold are masked out of the analysis.
#'  Default: \code{1e-6}. Variance is calculated on the original data, before
#'  any normalization.
#' 
#' @name varTol_Param
#' @keywords internal
NULL

#' scale_sm_FWHM
#' 
#' @param scale_sm_FWHM Only applies if \code{scale=="local"} and if \code{BOLD}
#'  represents surface data. To smooth the standard deviation estimates used for
#'  local scaling, provide the smoothing FWHM (default: \code{2}). As the FWHM
#'  becomes large, local scaling approaches global scaling. 
#' 
#' @name scale_sm_FWHM_Param
#' @keywords internal
NULL

#' scale
#' 
#' @param scale \code{"global"} (default), \code{"local"}, or \code{"none"}.
#'  Global scaling will divide the entire data matrix by the mean image standard 
#'  deviation (\code{mean(sqrt(rowVars(BOLD)))}). Local scaling will divide each
#'  data location's time series by its estimated standard deviation. 
#' 
#' @name scale_Param
#' @keywords internal
NULL

#' center_Bcols
#' 
#' @param center_Bcols Center BOLD across columns (each image)? Default: 
#'  \code{FALSE} (not recommended).
#' 
#' @name center_Bcols_Param
#' @keywords internal
NULL