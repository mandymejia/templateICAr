#' Apply data structure to templates
#' 
#' @param template The template
#' @param FORMAT "CIFTI", "NIFTI", or "DATA"
#' @param dat_struct The data structure
#' @param params The params
#' 
#' @keywords internal
struct_template <- function(template, FORMAT, dat_struct, params){
  if (FORMAT == "CIFTI") {
    template <- newdata_xifti(dat_struct, template)
    if (params$inds == paste("all", ncol(template))) {
      template$meta$cifti$names <- paste(
        "IC", seq(ncol(template))
      )
    } else {
      template$meta$cifti$names <- paste(
        "IC", strsplit(params$inds, " ")[[1]]
      )
    }
  } else if (FORMAT == "NIFTI") {
    RNifti::asNifti(
      unmask_subcortex(template, dat_struct, fill=NA)
    )
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
#'  objects directly. Otherwise, use a length-3 character vector of file 
#'  path(s) to save the output to:
#'  the mean template, the variance template, and the variance decomposition,
#'  in that order. If one file name is provided, it will be appended with
#'  \code{"_mean.[file_ext]"} for the template mean map,
#'  \code{"_var.[file_ext]"} for the template variance map, and
#'  \code{"_varDecomp.rds"} for the variance decomposition, where \code{[file_ext]}
#'  will be \code{"dscalar.nii"} for CIFTI input, \code{"nii"} for NIFTI input,
#'  and \code{"rds"} for data input.
#' @param var_method \code{"non-negative"} (default) or \code{"unbiased"}
#' 
#' @return If \code{is.null(out_fname)}, the templates in data matrix, 
#'  \code{"xifti"}, or \code{"nifti"} format, to match the format of the
#'  original BOLD data. Otherwise, the paths to the new files specified by
#'  \code{out_fname}.
#'
#' @export
export_template <- function(x, out_fname=NULL, var_method=c("non-negative", "unbiased")){
  
  # Check template format.
  FORMAT <- class(x)[grepl("template", class(x))]
  if (length(FORMAT) != 1) { stop("Not a template.") }
  FORMAT <- switch(FORMAT,
    template.cifti = "CIFTI",
    template.nifti = "NIFTI",
    template.data = "DATA"
  )
  FORMAT_extn <- switch(FORMAT, CIFTI=".dscalar.nii", NIFTI=".nii", DATA=".rds")

  var_method <- match.arg(var_method, c("non-negative", "unbiased"))
  var_name <- switch(var_method, `non-negative`="varNN", unbiased="varUB")

  x$template$varUB[] <- pmax(0, x$template$varUB)

  # `out_fname` ----------------------------------------------------------------
  if (!is.null(out_fname)) {
    out_fname <- as.character(out_fname)
    if (!all(dir.exists(dirname(out_fname)))) { stop('Directory part of `out_fname` does not exist.') }
    if (length(out_fname) == 1) {
      if (!endsWith(out_fname, FORMAT_extn)) { out_fname <- paste0(out_fname, FORMAT_extn) }
      out_fname <- c(
        gsub(FORMAT_extn, paste0("_mean", FORMAT_extn), out_fname),
        gsub(FORMAT_extn, paste0("_var", FORMAT_extn), out_fname),
        gsub(FORMAT_extn, paste0("_varDecomp.rds"), out_fname)
      )
    } else if (length(out_fname) == 3) {
      if (!all(endsWith(out_fname[seq(2)], FORMAT_extn))) {
        out_fname[seq(2)] <- paste0(out_fname[seq(2)], FORMAT_extn)
      }
      if (!endsWith(out_fname[3], ".rds")) {
        out_fname[3] <- paste0(out_fname[3], ".rds")
      }
    } else {
      stop(
        "`out_fname` should be a length 1 or 3 character vector giving the ",
        "names for:\n\tThe mean template,\n\tThe variance template,",
        "\n\tand the variance decomposition.\n"
      )
    }
  }

  x$template <- lapply(x$template, struct_template,
    FORMAT, x$dat_struct, x$params
  )

  # Select the chosen variance decomposition.
  x$template <- list(
    mean = x$template$mean,
    var = x$template[[var_name]]
  )

  # Add params to `"xifti"` metadata; resample it.
  if (FORMAT == "CIFTI") {
    x$template$mean$meta$cifti$misc <- c(list(template="mean"), x$params)
    x$template$var$meta$cifti$misc <- c(list(template="var"), x$params)
  }

  # Save
  if (!is.null(out_fname)) {
    if (FORMAT == "CIFTI") {
      write_cifti(x$template$mean, out_fname[1])
      write_cifti(x$template$var, out_fname[2])
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
  }

  if (is.null(out_fname)) {
    return(x$template)
  } else {
    return(invisible(out_fname))
  }
}