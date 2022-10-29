#' Plot FC template
#'
#' @param x The FC template
#' @param zlim Color limits
#' @param title Plot title
#' @param cols Colors
#' @param break_by Color breaks
#' @param cors \code{TRUE}
#' @export
plot_FC <- function(x, zlim=c(-1,1), title=NULL, cols=c('darkblue','turquoise','white','pink','red'), break_by=0.5, cor=TRUE){

  require(grDevices)

  #set color scale and breaks
  breaks <- seq(zlim[1], zlim[2], length.out=100)
  levs <- seq(zlim[1], zlim[2], break_by)
  palfun <- grDevices::colorRampPalette(cols)
  pal <- palfun(100-1)

  #make plot
  layout(matrix(c(1,2,0,3), nrow=2, ncol=2), widths=c(5,1), heights=c(1.2,5))
  #title
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.25, title, cex = 2, col = "black")
  # matrix
  if(cor) diag(x) <- NA
  size <- ncol(x)
  x[x >= zlim[2]] <- (zlim[2]-0.0001) #truncate above
  x[x <= zlim[1]] <- (zlim[1]+0.0001) #trunate below
  par(mar=c(1,2,0,1))
  image(seq(size), seq(size), t(x[size:1,]), col=pal, breaks=breaks-1e-8, xaxt="n", yaxt="n", ylab="", xlab="")
  abline(h=seq(size)-0.5, v=seq(size)-0.5)
  # color scale
  par(mar=c(1,1,0,3))
  image.scale(mat, col=pal, breaks=breaks-1e-8, axis.pos=4, add.axis=FALSE)
  axis(4,at=levs, las=2)
  abline(h=levs)
}

image.scale <- function(z, zlim, col = heat.colors(12),
                        breaks, axis.pos=1, add.axis=TRUE, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  if(axis.pos %in% c(1,3)){ylim<-c(0,1); xlim<-range(breaks)}
  if(axis.pos %in% c(2,4)){ylim<-range(breaks); xlim<-c(0,1)}
  plot(1,1,t="n",ylim=ylim, xlim=xlim, axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i", ...)
  for(i in seq(poly)){
    if(axis.pos %in% c(1,3)){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(axis.pos %in% c(2,4)){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
  box()
  if(add.axis) {axis(axis.pos)}
}
