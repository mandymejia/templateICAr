#' Check for INLA
#'
#' @return \code{NULL}, invisibly
#'
#' @keywords internal
INLA_check <- function(){
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop(
      "Package \"INLA\" needed to for spatial modeling. ",
      "Please install it at https://www.r-inla.org/download-install.",
      call. = FALSE
    )
  }
}

