#' Infer fMRI data format for a vector of inputs
#' 
#' Vectorized version of \code{infer_format_ifti}
#' 
#' Raises an error if the elements of \code{BOLD} do not share the same format.
#'
#' @param BOLD The vector of fMRI data, expected to be of one format
#' @param verbose Print the format? Default: \code{FALSE}.
#' @return A length-two vector. The first element indicates the format:
#'  \code{"CIFTI"} file path, \code{"xifti"} object,
#'  \code{"GIFTI"} file path, \code{"gifti"} object,  
#'  \code{"NIFTI"} file path, \code{"nifti"} object, 
#'  \code{"RDS"} file path, or \code{"data"}. The second element indicates 
#'  the sub-format if relevant; i.e. the type of CIFTI or GIFTI file/object.
#' 
#' @keywords internal
infer_format_ifti_vec <- function(BOLD, verbose=FALSE){
  BOLD <- as.list(BOLD)
  Bformat <- lapply(BOLD, infer_format_ifti, verbose=verbose)
  Bformat <- unique(Bformat)
  if (length(Bformat)>1) {
    stop(paste(
      "The formats are not identical: ", 
      paste(lapply(Bformat, paste, collapse=", "), collapse="; ")
    ))
  }
  Bformat[[1]]
}