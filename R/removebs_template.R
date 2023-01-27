#' Remove brain structure from CIFTI template
#' 
#' Remove a brain structure from a CIFTI template
#' 
#' @param x The \code{"template.cifti"} object.
#' @param remove \code{"cortex_left"}, \code{"cortex_right"}, and/or \code{"subcortical"}.
#' 
#' @keywords internal
removebs_template <- function(x, remove=NULL){
  stopifnot(inherits(x, "template.cifti"))
  remove <- match.arg(remove, c("cortex_left", "cortex_right", "subcortical"), several.ok=TRUE)

  # Remove brain structure(s) from data.
  x$template <- lapply(x$template, function(y){
    as.matrix(ciftiTools::remove_xifti(ciftiTools::newdata_xifti(x$dat_struct, y), remove=remove))
  })
  x$var_decomp <- lapply(x$var_decomp, function(y){
    as.matrix(ciftiTools::remove_xifti(ciftiTools::newdata_xifti(x$dat_struct, y), remove=remove))
  })

  # Get new `dat_struct` and mask.
  x$dat_struct <- ciftiTools::remove_xifti(x$dat_struct, remove=remove)
  x$mask <- !is.na(x$template$mean[,1])

  x
}