#' Create INLA mesh and observation weight matrix based on a binary brain mask
#'
#' @param mask Brain mask (matrix of 0/1 or TRUE/FALSE) (only supply mask OR vertices and faces)
#' @param vertices Matrix of vertices (only supply vertices and faces OR mask)
#' @param faces Matrix of faces (only supply vertices and faces OR mask)
#' @param keep Logical or 0/1 vector indicating vertices to keep, if using a submesh
#'
#' @return List containing INLA mesh, observation weight matrix (A) for translating between mesh locations and original data locations, the brain mask used to create the mesh, and the number of original and mesh data locations
#' @export
#' @importFrom INLA inla.nonconvex.hull inla.mesh.2d inla.spde.make.A inla.spde2.matern inla.mesh.create
#' @importFrom excursions submesh.mesh
#'
#' @note This function requires the \code{INLA} package, which is not a CRAN package. See \url{http://www.r-inla.org/download} for easy installation instructions.
make_mesh <- function(mask=NULL, vertices=NULL, faces=NULL, keep=NULL){

  hasmask <- !is.null(mask)
  hasverts <- !is.null(vertices)
  hasfaces <- !is.null(faces)

  if(hasmask & (hasverts | hasfaces)) stop('Must supply EITHER mask OR vertices and faces')

  if(hasmask){

    # Check only 0s and 1s
    values <- sort(unique(as.numeric(mask)))
    if(min(values %in% 0:1) == FALSE) stop("Mask should be composed of only 0s and 1s")

    xy.in <- which(mask==1, arr.ind=TRUE)[,2:1]
    boundary <- inla.nonconvex.hull(xy.in, resolution = 100)
    mesh <- inla.mesh.2d(loc = xy.in, boundary = boundary, max.edge = c(2, 4))
    Amat <- inla.spde.make.A(mesh, loc=xy.in)
    n.mask = sum(mask)

  } else {

    if(!(hasverts & hasfaces)) stop('Must supply vertices AND faces')

    # Check index of faces
    if(min(faces) == 0){ faces <- faces + 1 }

    # Construct INLA mesh
    mesh <- inla.mesh.create(loc = as.matrix(vertices), tv = as.matrix(faces))
    if(!is.null(keep)) mesh <- submesh.mesh(keep, mesh)

    mask <- NULL
    n.mask <- NULL
    Amat <- diag(rep(1, mesh$n))

  }

  spde <- inla.spde2.matern(mesh)

  result <- list(mesh=mesh, A=Amat, spde=spde, mask=mask, n.mask = n.mask, n.mesh = mesh$n)
  class(result) <- 'templateICA_mesh'
  return(result)

}
