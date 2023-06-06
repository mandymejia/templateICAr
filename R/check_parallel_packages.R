#' Check for required parallel packages
#' 
#' Check that required packages for parallel computation are available. If not,
#'  stop execution with a helpful error message.
#' 
#' @return \code{NULL}, invisibly
#' @keywords internal
check_parallel_packages <- function() {
  need_par_pkg <- FALSE
  if (!requireNamespace("foreach", quietly = TRUE)) {
    need_par_pkg <- TRUE
    message("Package \"foreach\" needed for parallel computation. Please install it.")
  }
  if (!requireNamespace("parallel", quietly = TRUE)) {
    need_par_pkg <- TRUE
    message("Package \"parallel\" needed for parallel computation. Please install it.")
  }
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    need_par_pkg <- TRUE
    message("Package \"doParallel\" needed for parallel computation. Please install it.")
  }
  if (need_par_pkg) { stop("Install required packages (see above), or disable parallel computation.") }
  invisible(NULL)
}