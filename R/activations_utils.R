#' Combine additive terms in string
#'
#' Combine two terms with "+" or "-" in a string
#'
#' @param a,b The two terms. Each is a length-one character vector.
#'
#' @return The result as a length-one character vector.
#' @keywords internal
add_str <- function(a,b){

  stopifnot(is_1(a, "character"))
  stopifnot(is.character(b))

  if (a=="0") {
    ifelse(b=="0", "0", as.character(b))
  } else {
    ifelse(b=="0", as.character(a),
       ifelse(
         grepl("^-", b),
         paste(a, "-", gsub("^-", "", b)),
         paste(a, "+", b)
       )
    )
  }
}

#' Format activation name
#'
#' Format the name of an activation, given parameters of the statistical test.
#'
#' @param u,z,type,deviation See \code{\link{activations}}.
#' @param collapse If multiple \code{u} or \code{z} were provided, should just
#'  one name be returned? Default: \code{FALSE}
#'
#' @keywords internal
#'
format_activation_name <- function(u, z, type, deviation, collapse=FALSE){

  stopifnot(is_1(collapse, "logical"))

  # `z` will have been converted to `u`, so `u` will always not be NULL.
  use_z <- !is.null(z)

  use_abs <- type == "abs >"

  # left-hand side.
  if (use_abs) {
    LHS <- "|x|"
    type <- ">"
  } else {
    LHS <- "x"
  }

  # right-hand side, build from right to left.
  RHS <- if (length(z)==1 || (length(z)>1 && !collapse)) {
    ifelse(z==0, "0", paste0(z, "*z"))
  } else if (length(z)>1) {
    "c*z (multiple `c` evaluated)"
  } else if (length(u)==1 || (length(u)>1 && !collapse)) {
    as.character(u)
  } else if (length(u) > 1) {
    "u (multiple `u` evaluated)"
  } else {
    "0"
  }

  RHS <- if (deviation) { add_str("mu", RHS) } else { add_str("0", RHS) }

  # all together now.
  # [Note]: could optionally gsub("1*", "", result, fixed=TRUE)
  paste(LHS, type, RHS)
}
