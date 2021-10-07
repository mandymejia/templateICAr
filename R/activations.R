#' Activations of (s)tICA
#'
#' Identify areas of activation in each independent component map
#'
#' @param result Fitted stICA or tICA model object (of class stICA or tICA)
#' @param u Activation threshold, default = 0
#' @param alpha Significance level for joint PPM, default = 0.1
#' @param type Type of region.  Default is '>' (positive excursion region).
#' @param method_p If result is type tICA, the type of multiple comparisons correction to use for p-values, or NULL for no correction.  See \code{help(p.adjust)}.
#' @param verbose If \code{TRUE}, display progress of algorithm. Default: \code{FALSE}.
#' @param which.ICs Indices of ICs for which to identify activations.  If NULL, use all ICs.
#' @param deviation If \code{TRUE}. identify significant deviations from the template mean, rather than significant areas of engagement
#'
#' @return A list containing activation maps for each IC and the joint and marginal PPMs for each IC.
#'
#' @export
#'
#' @importFrom excursions excursions
#' @importFrom stats pnorm p.adjust
#'
activations <- function(result, u=0, alpha=0.01, type=">", method_p='BH', verbose=FALSE, which.ICs=NULL, deviation=FALSE){

  if(!(type %in% c('>','<','!='))) stop("type must be one of: '>', '<', '!='")
  if(alpha <= 0 | alpha >= 1) stop('alpha must be between 0 and 1')
  if(!(class(result) %in% c('stICA','tICA'))) stop("result must be of class stICA or tICA")

  L <- ncol(result$A)
  if(is.null(which.ICs)) which.ICs <- 1:L
  if(min((which.ICs) %in% (1:L))==0) stop('Invalid entries in which.ICs')

  template_mean <- result$template_mean
  template_var <- result$template_var

  if (inherits(result, "stICA")) {

    if(verbose) cat('Determining areas of activations based on joint posterior distribution of latent fields\n')

    nvox <- nrow(result$subjICmean)
    L <- ncol(result$subjICmean)

    #identify areas of activation in each IC
    active <- jointPPM <- marginalPPM <- vars <- matrix(NA, nrow=nvox, ncol=L)

     for(q in which.ICs){
      if(verbose) cat(paste0('.. IC ',q,' (',which(which.ICs==q),' of ',length(which.ICs),') \n'))
      inds_q <- (1:nvox) + (q-1)*nvox
      if(deviation){
        Dinv_mu_s <-  (as.vector(result$subjICmean) - as.vector(template_mean) - u)/as.vector(sqrt(template_var))
      } else {
        Dinv_mu_s <- (as.vector(result$subjICmean) - u)/as.vector(sqrt(template_var))
      }

      if(q==which.ICs[1]) {
        #we scale mu by D^(-1) to use Omega for precision (actual precision of s|y is D^(-1) * Omega * D^(-1) )
        #we subtract u first since rescaling by D^(-1) would affect u too
        #save rho from first time running excursions, pass into excursions for other ICs
        tmp <- system.time(res_q <- excursions(alpha = alpha, mu = Dinv_mu_s, Q = result$Omega, type = type, u = 0, ind = inds_q))
        if(verbose) print(tmp)
        rho <- res_q$rho
      } else {
        tmp <- system.time(res_q <- excursions(alpha = alpha, mu = Dinv_mu_s, Q = result$Omega, type = type, u = 0, ind = inds_q, rho=rho))
        if(verbose) print(tmp)
      }
      active[,q] <- res_q$E[inds_q]
      jointPPM[,q] <- res_q$F[inds_q]
      marginalPPM[,q] <- res_q$rho[inds_q]
      vars[,q] <- res_q$vars[inds_q]
    }

    result <- list(active=active, jointPPM=jointPPM, marginalPPM=marginalPPM, vars=vars, u = u, alpha = alpha, type = type, deviation=deviation)
  }

  if (inherits(result, "tICA")) {

    if(verbose) cat('Determining areas of activations based on hypothesis testing at each location\n')

    nvox <- nrow(result$subjICmean)
    L <- ncol(result$subjICmean)

    if(deviation){
      t_stat <- as.matrix((result$subjICmean - template_mean - u)/sqrt(result$subjICvar))
    } else {
      t_stat <- as.matrix((result$subjICmean - u)/sqrt(result$subjICvar))
    }

    if(type=='>') pvals <- 1-pnorm(t_stat)
    if(type=='<') pvals <- pnorm(t_stat)
    if(type=='!=') pvals <- 2*(1-pnorm(abs(t_stat)))

    if(verbose) cat(paste0('Correcting for multiple comparisons with method ', method_p))

    pvals_adj <- active <- matrix(NA, nrow=nvox, ncol=L)
    if(is.null(method_p)) method_p <- 'none'
    for(q in which.ICs){
      pvals_adj[,q] <- p.adjust(pvals[,q], method=method_p)
      active[,q] <- (pvals_adj[,q] < alpha)
    }

    result <- list(
      active = active, pvals = pvals, pvals_adj = pvals_adj, tstats = t_stat,
      vars = result$subjICvar, u = u, alpha = alpha, method_p =
      method_p, deviation=deviation
    )
  }

  return(result)
}

#' Activations of (s)tICA
#'
#' Identify areas of activation in each independent component map
#'
#' @param result Result of templateICA.cifti model call
#' @param spatial_model Should spatial model result be used, if available?  If FALSE, will use standard template ICA result. If NULL, use spatial model result if available.
#' @param u Activation threshold, default = 0
#' @param alpha Significance level for joint PPM, default = 0.1
#' @param type Type of region.  Default is '>' (positive excursion region).
#' @param method_p Type of multiple comparisons correction to use for p-values for standard template ICA, or NULL for no correction.  See \code{help(p.adjust)}.
#' @param verbose If \code{TRUE}, display progress of algorithm. Default: \code{FALSE}.
#' @param which.ICs Indices of ICs for which to identify activations.  If NULL, use all ICs.
#' @param deviation If \code{TRUE}. identify significant deviations from the template mean, rather than significant areas of engagement
#'
#' @return A list containing activation maps for each IC and the joint and marginal PPMs for each IC.
#'
#' @importFrom ciftiTools newdata_xifti transform_xifti
#' @export
#'
#'
activations.cifti <- function(result, spatial_model=NULL, u=0, alpha=0.01, type=">", method_p='BH', verbose=FALSE, which.ICs=NULL, deviation=FALSE){

  if( (!inherits(result, "tICA")) && (!inherits(result, "stICA")) ) stop("result must be of class stICA or tICA")

  # Select stICA or tICA
  if (is.null(spatial_model)) { spatial_model <- inherits(result, "stICA") }
  if (isTRUE(spatial_model) && inherits(result, "tICA")) {
    warning(
      'spatial_model set to TRUE but class of model result is tICA. ', 
      'Setting spatial_model = FALSE, performing inference using standard ', 
      'template ICA.'
    )
    spatial_model <- FALSE
  }
  if (isFALSE(spatial_model) && inherits(result, "stICA")) {
    result <- result$result_tICA
  } 

  #run activations function
  activations_result <- activations(
    result, u=u, alpha=alpha, type=type, method_p=method_p, 
    verbose=verbose, which.ICs=which.ICs, deviation=deviation
  )

  activations_result$active <- newdata_xifti(
    result$subjICmean, as.numeric(activations_result$active)
  )
  activations_result$active <- convert_xifti(activations_result$active, "dlabel", colors="red")
  for (ii in seq(ncol(activations_result$active))) {
    rownames(activations_result$active$meta$cifti$labels[[ii]]) <- c("Inactive", "Active")
  }

  class(activations_result) <- "tICA_activations"

  activations_result
}

#' Summarize a \code{"tICA_activations"} object
#'
#' Summary method for class \code{"tICA_activations"}
#'
#' @param object Object of class \code{"tICA_activations"}. 
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method summary tICA_activations
summary.tICA_activations <- function(object, ...) {

  act_counts <- colSums(as.matrix(object$active), na.rm=TRUE)

  x <- c(
    summary(object$vars),
    list(act_counts=act_counts),
    object[c("u", "alpha", "method_p", "deviation")]
  )

  class(x) <- "summary.tICA_activations"
  return(x)
}

#' @rdname summary.tICA_activations
#' @export
#' 
#' @param x The activations from \code{activations.cifti}
#' @param ... further arguments passed to or from other methods.
#' @method print summary.tICA_activations
print.summary.tICA_activations <- function(x, ...) {
  
  cat("====ACTIVATIONS STATS================\n")
  cat("Threshold:       ", x$u, "\n")
  cat("alpha:           ", x$alpha, "\n")
  pm_nice <- switch(x$method_p,
    bonferroni = "Bonferroni",
    holm = "Holm",
    hochberg = "Hochberg",
    hommel = "Hommel",
    BH = "Benjamini & Hochberg (FDR)",
    fdr = "Benjamini & Hochberg (FDR)",
    by = "Benjamini & Yekutieli",
    none = "none"
  )
  cat("p-val method:    ", pm_nice, "\n")
  cat("Deviation:       ", x$deviation, "\n")
  # [TO DO]: add activation counts
  cat("\n")

  class(x) <- "summary.xifti"
  print(x) 
}

#' @rdname summary.tICA_activations
#' @export
#' 
#' @method print tICA_activations
print.tICA_activations <- function(x, ...) {
  print.summary.tICA_activations(summary(x))
}

#' Plot activations
#' 
#' @param x The activations from \code{activations.cifti}
#' @param stat \code{"active"} (default), \code{"pvals"}, \code{"pvals_adj"},
#'  \code{"tstats"}, or \code{"vars"}.
#' @param ... Additional arguments to \code{view_xifti}
#' @return The plot
#' @export
#' @importFrom ciftiTools view_xifti
#' @method plot tICA_activations
plot.tICA_activations <- function(x, stat=c("active", "pvals", "pvals_adj", "tstats", "vars"), ...) {
  stopifnot(inherits(x, "tICA_activations"))
  stat <- match.arg(stat, c("active", "pvals", "pvals_adj", "tstats", "vars"))
  if (stat == "active") {
    x <- x$active
  } else {
    x <- newdata_xifti(x$vars, as.matrix(x[[stat]]))
  }
  view_xifti(x, ...)
}