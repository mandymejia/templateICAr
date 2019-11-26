#' Create INLA mesh and observation weight matrix based on a binary brain mask
#'
#' @param mask Brain mask (matrix of 0/1 or TRUE/FALSE)
#'
#' @return List containing INLA mesh, observation weight matrix (A) for translating between mesh locations and original data locations, the brain mask used to create the mesh, and the number of original and mesh data locations
#' @export
#' @importFrom INLA inla.nonconvex.hull inla.mesh.2d inla.spde.make.A inla.spde2.matern
#'
make_mesh <- function(mask){

  # Check only 0s and 1s
  values <- sort(unique(as.numeric(mask)))
  if(min(values %in% 0:1) == FALSE) stop("Mask should be composed of only 0s and 1s")

  xy.in <- which(mask==1, arr.ind=TRUE)[,2:1]
  boundary <- inla.nonconvex.hull(xy.in, resolution = 100)
  mesh <- inla.mesh.2d(loc = xy.in, boundary = boundary, max.edge = c(2, 4))
  spde <- inla.spde2.matern(mesh, alpha=2)
  Amat <- inla.spde.make.A(mesh, loc=xy.in)

  result <- list(mesh=mesh, A=Amat, spde=spde, mask=mask, n.mask = sum(mask), n.mesh = mesh$n)
  class(result) <- 'templateICA_mesh'
  return(result)

}
