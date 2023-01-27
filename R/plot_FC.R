#' Plot FC template
#' 
#' @param mat The FC template matrix from \code{estimate_template.cifti}
#' @param colFUN A \code{scale_fill} function. If \code{NULL}, use this default:
#'  \code{function(limits=c(0,1), ...){ viridis::scale_fill_viridis(option="inferno", limits=limits, ...) }}
#' @param title The plot title. Default: \code{"FC template"}.
#' @param legTitle The legend title. Default: \code{"FC"}.
#' @param group_divs,group_cols Split the FC matrix into groups of contiguous
#'  rows/columns? Use \code{group_divs} to indicate the index of the first 
#'  element in each group. For example, if the groups are 1-8, 9-15, and 16 to
#'  25, set \code{group_divs=c(1, 9, 16). Groups will be indicated by a color
#'  bar along the left and top sides of the plot. Use \code{group_cols} to 
#'  change the colors. By default, \code{group_divs} is \code{NULL} (no groups).
#' @param uppertri_means Use the upper triangle of the FC matrix to show the
#'  mean of each group indicated by \code{group_divs}? Default: \code{TRUE}.
#' @param divColor Color of group dividers and outline of the FC matrix. 
#'  Default: \code{"black"}.
#' @param lim Length-two numeric vector indicating the limits for the color bar.
#'  Values beyond \code{lim} will be clipped. If \code{NULL} (default), the
#'  limits will be set to the largest magnitude value in \code{mat} (no
#'  clipping).
#' @param diagVal On-diagonal values will be set to this value. 
#'  (\code{uppertri_means} are calculated before \code{diagVal} is used.)
#' @return The plot
#' @export
#' @importFrom RColorBrewer brewer.pal
#' @importFrom viridis scale_fill_viridis
#' @importFrom ggcorrplot ggcorrplot
#' @importFrom ggplot2 coord_equal guide_colorbar labs theme element_blank margin annotate geom_rect aes
#' @method plot template.cifti
plot.FC_template <- function(
  mat,
  colFUN=NULL, 
  title="FC template", legTitle="FC", 
  group_divs=NULL, group_cols=RColorBrewer::brewer.pal(8, "Set2"), 
  uppertri_means=TRUE, divColor="black", lim=NULL, diagVal=1
  ){

  nD <- nrow(mat)
  stopifnot(nD==ncol(mat))

  dw <- .004 * nD

  use_groups <- !is.null(group_divs)
  if (use_groups) {
    stopifnot(all(diff(group_divs) > 0))
    stopifnot(group_divs[1] > 0)
    if (group_divs[1] != 1) { group_divs <- c(1, group_divs) }
    stopifnot(group_divs[length(group_divs)] < nD+2)
    if (group_divs[length(group_divs)] != nD+1) { group_divs <- c(group_divs, nD+1) }
    group_divs2 <- group_divs[seq(2, length(group_divs)-1)]
    gg_pdv <- data.frame(xmin=group_divs2-.5-(dw/2), xmax=group_divs2-.5+(dw/2))
    gg_pdv$Var1 <- Var2 <- 0

    the_cols <- as.numeric(cut(seq(nD)+1, c(-Inf, group_divs2, Inf)))
  }

  if (is.null(colFUN)) {
    colFUN <- function(limits=c(0,1), ...){ 
      viridis::scale_fill_viridis(option="inferno", limits=limits, ...)
    }
  }

  # Take network-network means in upper tri
  if (uppertri_means && use_groups) {
    mat2 <- mat
    for (jj in seq(length(group_divs)-1)) {
      for (kk in seq(length(group_divs)-1)) {
        mat2[seq(group_divs[jj], group_divs[jj+1]-1), nD+1 - seq(group_divs[kk], group_divs[kk+1]-1)] <- mean(
          mat2[seq(group_divs[jj], group_divs[jj+1]-1), nD+1 - seq(group_divs[kk], group_divs[kk+1]-1)], na.rm=TRUE
        )
      }
    }
  }
  if (!is.null(diagVal)) {
    mat[is.na(mat)] <- diagVal
  }
  if (uppertri_means && use_groups) {
    mat[seq(nD), rev(seq(nD))][lower.tri(mat)] <- mat2[seq(nD), rev(seq(nD))][lower.tri(mat)]
  }

  # Limits
  if (!is.null(lim)) { 
    if (length(lim)==1) { lim <- c(-lim, lim) }
    stopifnot(length(lim)==2)
    mat[] <- pmax(lim[1], pmin(lim[2], mat))
  } else {
    lim <- max(abs(mat[]))
    lim <- c(-lim, lim)
  }

  lim_expand <- ifelse(use_groups, nD*.1, 0)

  plt <- ggcorrplot(mat, outline.color = "#00000000", title=title, digits=12) + 
    #scale_y_discrete(labels=cor_mat_ylabs) +
    coord_equal(xlim=.5+c(-1-lim_expand-dw, nD+dw), ylim=.5+c(-dw, nD+1+lim_expand+dw), expand=FALSE) +
    colFUN(
      c(lim[1], lim[2]), guide=guide_colorbar(ticks.colour = divColor, ticks.linewidth = 1), 
      labels = function(x){gsub("0.", ".", x, fixed=TRUE)}, na.value="red"
    ) +
    labs(fill=legTitle) +
    theme(
      panel.grid.major = element_blank(),
      axis.text.y = element_text(margin=margin(r=10)),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "bottom"#, legend.text.align = 1
    ) +
    annotate(
      "raster", x=seq(nD), y=rev(seq(nD)), fill=divColor
    ) +
    # brain region label color bar
    annotate(
      "rect", xmin=.5-lim_expand-dw, xmax=.5, ymin=.5+nD-seq(nD), ymax=.5+nD-seq(nD)+1,
      fill=group_cols[the_cols]
    ) +
    annotate(
      "rect", ymin=.5+nD, ymax=.5+nD+lim_expand+dw, xmin=seq(nD)-.5, xmax=seq(nD)+.5,
      fill=group_cols[the_cols]
    )

  if (use_groups) {
    plt <- plt +
      # dividers
      geom_rect(aes(xmin=xmin, xmax=xmax, ymin=0, ymax=.5+nD+lim_expand+dw), data=gg_pdv, fill=divColor) +
      geom_rect(aes(ymin=nD+1-xmin, ymax=nD+1-xmax, xmin=.5-lim_expand-dw, xmax=.5+nD), data=gg_pdv, fill=divColor)
  }

  plt <- plt + 
    # four edges
    annotate("rect", xmin=.5-dw, xmax=.5,       ymin=.5-dw,            ymax=.5+nD+lim_expand+dw, fill=divColor) +
    annotate("rect", xmin=.5+nD, xmax=.5+nD+dw, ymin=.5-dw,            ymax=.5+nD+lim_expand+dw, fill=divColor) + 
    annotate("rect", ymin=.5-dw, ymax=.5,       xmin=.5-lim_expand-dw, xmax=.5+nD+dw,            fill=divColor) +
    annotate("rect", ymin=.5+nD, ymax=.5+nD+dw, xmin=.5-lim_expand-dw, xmax=.5+nD+dw,            fill=divColor)

  if (use_groups) {
    plt <- plt + 
      # separate brain region label color bar
      annotate("rect", ymin=-dw, ymax=nD+dw+lim_expand, xmin=.5-dw-3, xmax=.5-dw, fill="white") +
      annotate("rect", ymin=nD+.5+dw, ymax=nD+.5+dw+3, xmin=-.5-dw, xmax=nD+.5+dw, fill="white")
  }

  plt
}