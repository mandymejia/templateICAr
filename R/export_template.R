#' Apply data structure to templates
#'
#' @param template The template
#' @param FORMAT "CIFTI", "GIFTI", "NIFTI", or "DATA"
#' @param dat_struct The data structure
#' @param mask_input The input mask
#' @param params The params
#'
#' @keywords internal
struct_template <- function(template, FORMAT, mask_input, params, dat_struct, GICA_parc_table){
  # Un-apply the input mask.
  if (!is.null(mask_input)) {
    if (FORMAT=="NIFTI") {
      template <- fMRItools::unvec_vol(template, drop(mask_input))
    } else {
      template <- fMRItools::unmask_mat(template, mask_input)
    }
  }

  # Add metadata.
  if (FORMAT == "CIFTI") {
    if (!requireNamespace("ciftiTools", quietly = TRUE)) {
      stop("Package \"ciftiTools\" needed to work with CIFTI data. Please install it.", call. = FALSE)
    }
    template <- ciftiTools::newdata_xifti(dat_struct, template)
    template$meta$cifti$names <- if (!is.null(GICA_parc_table)) {
      rownames(GICA_parc_table)
    } else {
      paste("IC", params$inds)
    }

  } else if (FORMAT == "GIFTI") {
    if (!requireNamespace("ciftiTools", quietly = TRUE)) {
      stop("Package \"ciftiTools\" needed to work with GIFTI data. Please install it.", call. = FALSE)
    }

    template <- ciftiTools:::as.metric_gifti(
      template, hemisphere=dat_struct$hemisphere
    )

  } else if (FORMAT == "NIFTI") {
    template <- RNifti::asNifti(template, reference=mask_input)
  }

  template
}

#' Export template
#'
#' Export the templates (mean and variance) as separate files for
#'  visualization or processing outside of \code{templateICAr}.
#'
#' @param x The result of \code{estimate_template}
#' @param out_fname Use \code{NULL} (default) to just return the template
#'  objects directly. Otherwise, use a character vector of length 3 or 4 of file
#'  path(s) to save the output to:
#'  the mean template, the variance template, the variance decomposition, and
#'  the FC template if present, in that order. If one file name is provided,
#'  it will be appended with
#'  \code{"_mean.[file_ext]"} for the template mean map,
#'  \code{"_var.[file_ext]"} for the template variance map,
#'  \code{"_varDecomp.rds"} for the variance decomposition, and
#'  \code{"_FC.rds"} where \code{[file_ext]}
#'  will be \code{"dscalar.nii"} for CIFTI input, \code{"nii"} for NIFTI input,
#'  and \code{"rds"} for data input.
#' @param var_method \code{"non-negative"} (default) or \code{"unbiased"}
#'
#' @return If \code{is.null(out_fname)}, the templates in data matrix,
#'  \code{"xifti"}, or \code{"nifti"} format, to match the format of the
#'  original BOLD data. Otherwise, the paths to the new files specified by
#'  \code{out_fname}. If template includes functional connectivity components,
#'  the FC template and its mean and variance will be included.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'  tm <- estimate_template(cii1_fnames, cii2_fnames, gICA_fname)
#'  export_template(tm, out_fname="my_template", var_method="unbiased")
#' }
#'
export_template <- function(x, out_fname=NULL, var_method=c("non-negative", "unbiased")){

  # Check template format.
  FORMAT <- class(x)[grepl("template", class(x))]
  if (length(FORMAT) != 1) { stop("Not a template.") }
  FORMAT <- switch(FORMAT,
    template.cifti = "CIFTI",
    template.gifti = "GIFTI",
    template.nifti = "NIFTI",
    template.matrix = "MATRIX"
  )
  FORMAT_extn <- switch(FORMAT, CIFTI=".dscalar.nii", GIFTI=".func.gii", NIFTI=".nii", MATRIX=".rds")

  var_method <- match.arg(var_method, c("non-negative", "unbiased"))
  var_name <- switch(var_method, `non-negative`="varNN", unbiased="varUB")

  x$template$varUB[] <- pmax(0, x$template$varUB)

  FC <- "FC" %in% names(x$template)

  # Fix for old version
  if (FORMAT == "CIFTI" && !is.null(x$dat_struct$meta$subcort$labels)) {
    if (!requireNamespace("ciftiTools", quietly = TRUE)) {
      stop("Package \"ciftiTools\" needed to export template. Please install it.", call. = FALSE)
    }
    sub_levs <- levels(x$dat_struct$meta$subcort$labels)
    if (length(sub_levs) != length(ciftiTools::substructure_table()$ciftiTools_Name)) {
      x$dat_struct$meta$subcort$labels <- factor(
        x$dat_struct$meta$subcort$labels,
        levels = ciftiTools::substructure_table()$ciftiTools_Name
      )
      stopifnot(ciftiTools:::is.subcort_labs(x$dat_struct$meta$subcort$labels))
    }
  }

  # `out_fname` ----------------------------------------------------------------
  if (!is.null(out_fname)) {
    out_fname <- as.character(out_fname)
    if (!all(dir.exists(dirname(out_fname)))) { stop('Directory part of `out_fname` does not exist.') }
    if (length(out_fname) == 1) {
      if (!endsWith(out_fname, FORMAT_extn)) { out_fname <- paste0(out_fname, FORMAT_extn) }
      out_fname <- c(
        gsub(paste0(FORMAT_extn, "$"), paste0("_mean", FORMAT_extn), out_fname),
        gsub(paste0(FORMAT_extn, "$"), paste0("_var", FORMAT_extn), out_fname),
        gsub(paste0(FORMAT_extn, "$"), paste0("_varDecomp.rds"), out_fname),
        gsub(paste0(FORMAT_extn, "$"), paste0("_FC.rds"), out_fname)
      )
      if (!FC) { out_fname <- out_fname[seq(3)] }
    } else if (length(out_fname) == 3 + as.numeric(FC)) {
      if (!all(endsWith(out_fname[seq(2)], FORMAT_extn))) {
        out_fname[seq(2)] <- paste0(out_fname[seq(2)], FORMAT_extn)
      }
      if (!endsWith(out_fname[3], ".rds")) {
        out_fname[3] <- paste0(out_fname[3], ".rds")
      }
      if (FC && !endsWith(out_fname[4], ".rds")) {
        out_fname[4] <- paste0(out_fname[4], ".rds")
      }
    } else {
      stop(
        "`out_fname` should be a length 1 or 3/4 character vector giving the ",
        "names for:\n\tThe mean template,\n\tthe variance template,",
        "\n\tthe variance decomposition, and \n\tthe FC template.\n"
      )
    }
  }

  tm_struct_mask <- !(names(x$template) %in% c("FC", "FC_Chol"))
  x$template[tm_struct_mask] <- lapply(
    x$template[tm_struct_mask], struct_template,
    FORMAT, x$mask_input, x$params, x$dat_struct, x$GICA_parc_table
  )

  # Select the chosen variance decomposition.
  x$template <- list(
    mean = x$template$mean,
    var = x$template[[var_name]],
    FC = x$template$FC
    #, FC_Chol=FC_Chol # [TO DO]: want to export anything in `FC_Chol`?
  )

  #compute mean and variance of FC
  if(FC){
    Q <- nrow(x$template$FC$psi)
    FC_mean <- x$template$FC$psi/(x$template$FC$nu - Q - 1)
    FC_var <- FC_mean*0
    for(q1 in 1:Q){
      for(q2 in 1:Q){
        FC_var[q1,q2] <- IW_var(x$template$FC$nu, Q, FC_mean[q1,q2], FC_mean[q1,q1], FC_mean[q2,q2])
      }
    }
    x$template$FC$mean <- FC_mean
    x$template$FC$var <- FC_var
  }

  # Add params to `"xifti"` metadata; resample it.
  if (FORMAT == "CIFTI") {
    x$params <- lapply(
      x$params,
      function(q) {
        if (is.null(q)) { q <- "NULL"};
        paste0(as.character(q), collapse=" ")
      }
    )
    x$template$mean$meta$cifti$misc <- c(list(template="mean"), x$params)
    x$template$var$meta$cifti$misc <- c(list(template="var"), x$params)
  }

  # Save
  if (!is.null(out_fname)) {
    if (FORMAT == "CIFTI") {
      ciftiTools::write_cifti(x$template$mean, out_fname[1])
      ciftiTools::write_cifti(x$template$var, out_fname[2])
    } else if (FORMAT == "GIFTI") {
      if (!requireNamespace("gifti", quietly = TRUE)) {
        stop("Package \"gifti\" needed to write NIFTI data. Please install it.", call. = FALSE)
      }
      gifti::writegii(x$template$mean, out_fname[1])
      gifti::writegii(x$template$var, out_fname[2])
    } else if (FORMAT == "NIFTI") {
      if (!requireNamespace("RNifti", quietly = TRUE)) {
        stop("Package \"RNifti\" needed to write NIFTI data. Please install it.", call. = FALSE)
      }
      RNifti::writeNifti(x$template$mean, out_fname[1])
      RNifti::writeNifti(x$template$var, out_fname[2])
    } else {
      saveRDS(x$template$mean, out_fname[1])
      saveRDS(x$template$var, out_fname[2])
    }
    saveRDS(x$var_decomp, out_fname[3])
    if (FC) { saveRDS(x$template$FC, out_fname[4]) }
  }

  if (is.null(out_fname)) {
    return(x$template)
  } else {
    return(invisible(out_fname))
  }
}



