#' Resample CIFTI template
#'
#' Resample a CIFTI template to a new spatial resolution.
#'
#' @param x The \code{"template.cifti"} object.
#' @param resamp_res The new resampling resolution. If NULL, do not perform resampling.
#' @param brainstructures The brainstructures to retain.
#' @param verbose Give occasional updates? Default: \code{FALSE}.
#'
#' @return The resampled \code{"template.cifti"} object.
#' @export
#'
resample_template <- function(x, resamp_res, brainstructures, verbose=FALSE){
  stopifnot(inherits(x, "template.cifti"))

  if (!requireNamespace("ciftiTools", quietly = TRUE)) {
    stop("Package \"ciftiTools\" needed to work with CIFTI data. Please install it.", call. = FALSE)
  }

  tm_struct_mask <- !(names(x$template) %in% c("FC", "FC_Chol"))
  xii <- x$dat_struct #example data structure (i.e. xifti object)

  ### Identify brainstructures to remove from template.

  bs_template <- x$params$brainstructures
  if(bs_template == 'all') bs_template <- c('left','right','subcortical')
  if(!all(brainstructures %in% bs_template)) stop('At least one brainstructure does not exist in the template. Adjust brainstructures argument to be a subset of template brainstructures.')
  bs_template_rm <- bs_template[!(bs_template %in% brainstructures)]
  if(length(bs_template_rm)==0) bs_template_rm <- NULL
  if (verbose & !is.null(bs_template_rm)) cat("\t Removing extra brainstructures from template: ", paste(bs_template_rm, collapse=', '), "\n")
  if(!is.null(bs_template_rm)){
    bs_template_rm[bs_template_rm=='left'] <- 'cortex_left'
    bs_template_rm[bs_template_rm=='right'] <- 'cortex_right'
  }

  ### Resample the data & remove extraneous brainstructures.

  #function to resample the matrix-format data and remove extraneous brainstructures
  fn_resamp_remove <- function(y){
    y_xii <- ciftiTools::newdata_xifti(xii, y) #temporarily convert from matrix to xifti
    if(!is.null(resamp_res)) y_xii <- ciftiTools::resample_xifti(y_xii, resamp_res=resamp_res, verbose=FALSE) # resample
    if(!is.null(bs_template_rm)) y_xii <- ciftiTools::remove_xifti(y_xii, bs_template_rm) # remove extra brainstructures
    as.matrix(y_xii)
  }

  if (verbose) { cat("\t Resampling template mean and variance ... ") }
  x$template[tm_struct_mask] <- lapply(x$template[tm_struct_mask], fn_resamp_remove)
  if (verbose) { cat("variance decomposition ... ") }
  x$var_decomp <- lapply(x$var_decomp, fn_resamp_remove)
  if (verbose) { cat("... and sigma_sq0.\n") }
  x$sigma_sq0 <- fn_resamp_remove(x$sigma_sq0)

  # Replace `NaN` values with NA values.
  x$template[tm_struct_mask] <- lapply(x$template[tm_struct_mask], function(y){y[] <- ifelse(is.nan(y), NA, y)})
  x$var_decomp <- lapply(x$var_decomp, function(y){y[] <- ifelse(is.nan(y), NA, y)})
  x$sigma_sq0 <- ifelse(is.nan(x$sigma_sq0), NA, x$sigma_sq0)

  ### Update `dat_struct` and mask.

  if (verbose) { cat("\t Formatting new template object.\n") }

  if(!is.null(resamp_res))  xii <- ciftiTools::resample_xifti(xii, resamp_res=resamp_res) #optionally resample
  if(!is.null(bs_template_rm)) xii <- ciftiTools::remove_xifti(xii, bs_template_rm) #remove extra brainstructures
  x$dat_struct <- xii
  x$mask <- !is.na(x$template$mean[,1])

  x
}
