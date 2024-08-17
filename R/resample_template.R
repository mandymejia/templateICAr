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

  if (!requireNamespace("ciftiTools", quietly = TRUE)) {
    stop("Package \"ciftiTools\" needed to work with CIFTI data. Please install it.", call. = FALSE)
  }

  tm_struct_mask <- !(names(x$template) %in% c("FC", "FC_Chol")) # mean, varUB, varNN

  # Resample the data.
  if (verbose) { cat("Resampling templates.\n") }
  x$template[tm_struct_mask] <- lapply(x$template[tm_struct_mask], function(y){
    as.matrix(ciftiTools::resample_xifti(ciftiTools::newdata_xifti(x$dat_struct, y), resamp_res=resamp_res, verbose=FALSE))
  })

  if (verbose) { cat("Resampling variance decomposition.\n") }
  x$var_decomp <- lapply(x$var_decomp, function(y){
    as.matrix(ciftiTools::resample_xifti(ciftiTools::newdata_xifti(x$dat_struct, y), resamp_res=resamp_res, verbose=FALSE))
  })
  if (!is.null(x$sigma_sq0)) {
    x$sigma_sq0 <- c(as.matrix(ciftiTools::resample_xifti(ciftiTools::newdata_xifti(x$dat_struct, x$sigma_sq), resamp_res=resamp_res, verbose=FALSE)))
  }

  if (verbose) { cat("Formatting new template object.\n") }
  # Replace `NaN` values with NA values.
  x$template[tm_struct_mask] <- lapply(x$template[tm_struct_mask], function(y){y[] <- ifelse(is.nan(y), NA, y)})
  x$var_decomp <- lapply(x$var_decomp, function(y){y[] <- ifelse(is.nan(y), NA, y)})
  if (!is.null(x$sigma_sq0)) {
    x$sigma_sq0 <- ifelse(is.nan(x$sigma_sq0), NA, x$sigma_sq0)
  }

  # Get new `dat_struct` and mask.
  x$dat_struct <- ciftiTools::resample_xifti(x$dat_struct, resamp_res=resamp_res)
  x$mask <- !is.na(x$template$mean[,1])

  x
}
