#' Infer fMRI data format
#'
#' @param BOLD The fMRI data
#' @param verbose Print the format? Default: \code{FALSE}.
#' @return The format: \code{"CIFTI"} file path, \code{"xifti"} object,
#'  \code{"GIFTI"} file path, \code{"gifti"} object,  
#'  \code{"NIFTI"} file path, \code{"nifti"} object, 
#'  or \code{"data"}.
#' @keywords internal
infer_BOLD_format <- function(BOLD, verbose=FALSE){

  Bformat <- NULL

  # Character vector: CIFTI, GIFTI, or NIFTI
  if (is.character(BOLD)) {
    Bformat <- ifelse(
      endsWith(BOLD, ".dtseries.nii") | endsWith(BOLD, ".dscalar.nii"),
      "CIFTI",
      ifelse(endsWith(BOLD, ".gii"), "GIFTI", "NIFTI")
    )
    if (length(unique(Bformat)) == 1) {
      Bformat <- Bformat[1]
    } else {
      if (all(endsWith(BOLD, "nii") | endsWith(BOLD, "gii"))) {
        stop("BOLD format seems to be a mix. Use the same format or correct the file names.")
      } else {
        stop("BOLD list seems to include unexpected file types. Use a correct format or fix the names.")
      }
    }

  # The xifti, gifti, or niftiImage data object
  } else if (inherits(BOLD, "xifti")) {
    Bformat <- "xifti"
  } else if (inherits(BOLD, "gifti")) {
    Bformat <- "gifti"
  } else if (inherits(BOLD, "niftiImage")) {
    Bformat <- "nifti"

  # Non-character vector: xifti, gifti, or niftiImage data object
  } else if (inherits(BOLD[[1]], "xifti")) {
    if (all(vapply(BOLD, inherits, what="xifti", FALSE))) {
      Bformat <- "xifti"
    } else {
      stop("BOLD format seems to be a mix of `xifti` files and something else. Use the same format for all.")
    }
  } else if (inherits(BOLD[[1]], "gifti")) {
    if (all(vapply(BOLD, inherits, what="gifti", FALSE))) {
      Bformat <- "gifti"
    } else {
      stop("BOLD format seems to be a mix of `gifti` files and something else. Use the same format for all.")
    }
  } else if (inherits(BOLD[[1]], "niftiImage")) {
    if (all(vapply(BOLD, inherits, what="niftiImage", FALSE))) {
      Bformat <- "nifti"
    } else {
      stop("BOLD format seems to be a mix of `nifti` files and something else. Use the same format for all.")
    }

  # List: GIFTI right and left
  } else if (is.list(BOLD)) {
    if ((length(BOLD) == 2) && length(BOLD[[1]]) == length(BOLD[[2]])) {
      if (is.character(BOLD[[1]]) && is.character(BOLD[[2]])) {
        if (all(endsWith(BOLD[[1]], "gii")) && all(endsWith(BOLD[[2]], "gii"))) {
          Bformat <- "GIFTI2"
        }
      }
    } else if (all(vapply(do.call(c, BOLD), inherits, what="gifti", FALSE))) {
      Bformat <- "gifti2"
    }
  }

  # Data matrix
  if (is.null(Bformat)) {
    if (is.numeric(BOLD) || (is.list(BOLD) && is.numeric(BOLD[[1]]))) {
      if (!is.list(BOLD)) { BOLD <- list(BOLD) }
      BOLD_dims <- lapply(BOLD, dim)
      BOLD_dims_lens <- sort(unique(vapply(BOLD_dims, length, 0)))
      if (length(BOLD_dims_lens) > 1) {
        stop("BOLD data have inconsistent dimension lengths. fMRI data should be provided as matrices, not vectors or arrays.")
      } else if (BOLD_dims_lens==4) {
        Bformat <- "nifti" # 4D array: treat as a "nifti"
      } else if (BOLD_dims_lens!=2) {
        stop("BOLD data should be provided as matrices, not vectors or arrays.")
      } else {
        Bformat <- "data"
      }
    }
  }

  if (is.null(Bformat)) {
    stop("Could not infer BOLD format.")
  }

  if (verbose) { cat("Inferred input format:", Bformat, "\n") }
  Bformat
}