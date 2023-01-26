#' Resample CIFTI template
#' 
#' Resample a CIFTI template to a new spatial resolution.
#' 
#' @param x The \code{"template.cifti"} object.
#' @param resamp_res The new resampling resolution.
#' @param verbose Give occasional updates? Default: \code{FALSE}.
#' 
#' @return The resampled \code{"template.cifti"} object.
#' @export
#' 
resample_template <- function(x, resamp_res, verbose=FALSE){
  stopifnot(inherits(x, "template.cifti"))

  # Resample the data.
  if (verbose) { cat("Resampling templates.\n") }
  x$template <- lapply(x$template, function(y){
    as.matrix(ciftiTools::resample_xifti(ciftiTools::newdata_xifti(x$dat_struct, y), resamp_res=resamp_res, verbose=FALSE))
  })
  if (verbose) { cat("Resampling variance decomposition.\n") }
  x$var_decomp <- lapply(x$var_decomp, function(y){
    as.matrix(ciftiTools::resample_xifti(ciftiTools::newdata_xifti(x$dat_struct, y), resamp_res=resamp_res, verbose=FALSE))
  })

  if (verbose) { cat("Formatting new template object.\n") }
  # Replace `NaN` values with NA values.
  x$template <- lapply(x$template, function(y){y[] <- ifelse(is.nan(y), NA, y)})
  x$var_decomp <- lapply(x$var_decomp, function(y){y[] <- ifelse(is.nan(y), NA, y)})

  # Get new `dat_struct` and mask.
  x$dat_struct <- ciftiTools::resample_xifti(x$dat_struct, resamp_res=resamp_res)
  x$mask <- !is.na(x$template$mean[,1])

  x
}