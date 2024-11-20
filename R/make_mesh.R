#' Make INLA mesh from \code{"surf"} object
#' 
#' Create INLA mesh and observation weight matrix based on a \code{"surf"} object
#'
#' @param surf Object of class \code{"surf"}. See
#'  \code{\link[ciftiTools]{make_surf}} and \code{\link[ciftiTools]{is.surf}}.
#' @param inds_data Subset of vertices to include in analysis, e.g. non-medial 
#'  wall locations.
#' @param inds_mesh Subset of vertices to retain in mesh, e.g. non-medial wall 
#'  locations. Must be a superset of \code{inds_data}.
#'
#' @return List containing INLA mesh, observation weight matrix \strong{A} for 
#'  translating between mesh locations and original data locations, the brain 
#'  mask used to create the mesh, and the number of original and mesh data 
#'  locations.
#' 
#' @export
#' 
# @importFrom INLA inla.spde2.matern inla.mesh.create
#' @importFrom Matrix Diagonal
#'
make_mesh <- function(surf=NULL, inds_data=NULL, inds_mesh=NULL){

  INLA_check()

  if (!requireNamespace("ciftiTools", quietly = TRUE)) {
    stop("Package \"ciftiTools\" needed to work with CIFTI data. Please install it.", call. = FALSE)
  }

  if (!requireNamespace("excursions", quietly = TRUE)) {
    stop("Package \"excursions\" needed for `make_mesh`. Please install.", call. = FALSE)
  }

  stopifnot(ciftiTools::is.surf(surf))

  #if inds_mesh is NULL, keep all current vertices in the mesh
  nmesh_orig <- nrow(surf$vertices)
  if (is.null(inds_mesh)) { inds_mesh <- seq(nmesh_orig) }

  #check that inds_data is a subset of inds_mesh
  if (!is.null(inds_data)) {
    if(any(!(inds_data %in% seq(nmesh_orig)))) stop(paste0('`inds_data` should contain only indices from 1 to ', nmesh_orig))
    if(any(!(inds_data %in% inds_mesh))) stop('`inds_data` must be a subset of `inds_mesh`')
  }
  if (any(!(inds_mesh %in% seq(nmesh_orig)))) stop(paste0('`inds_mesh` should contain only indices from 1 to ', nmesh_orig))

  # 1. Construct INLA mesh
  # check locs, should be seq(nmesh_orig)
  mesh <- INLA::inla.mesh.create(
    loc = as.matrix(surf$vertices), 
    tv = as.matrix(surf$faces)
  )

  # 2. Use submesh.mesh to exclude vertices not in inds_mesh
  nmesh_new <- length(inds_mesh)
  keep <- (seq(nmesh_orig)) %in% inds_mesh
  mesh <- excursions::submesh.mesh(keep, mesh) #check locs, should be 1:nmesh_new
  mesh$idx$loc <- mesh$idx$loc[!is.na(mesh$idx$loc)]

  # 3. Record which mesh locations are data locations & adjust Amat
  Amat <- Diagonal(nmesh_new, x=1)
  if (!is.null(inds_data)) {
    inds_data_mesh <- which(inds_mesh %in% inds_data)
    mesh$idx$loc <- mesh$idx$loc[inds_data_mesh] #remove masked-out vertices from vector of data locations in mesh$idx$loc
    Amat <- Amat[inds_data_mesh,] #data projection matrix, project to only vertices in inds_data
  }

  spde <- INLA::inla.spde2.matern(mesh)

  result <- list(
    mesh=mesh, A=Amat, spde=spde, n.mesh = mesh$n, 
    inds_data = inds_data, inds_mesh = inds_mesh
  )
  class(result) <- 'templateICA_mesh'
  result
}

#' Make 2D INLA mesh
#' 
#' Create INLA mesh and observation weight matrix based on a binary brain mask
#'
#' This function requires the \code{INLA} package, which is not a CRAN 
#'  package. See \url{https://www.r-inla.org/download-install} for easy 
#'  installation instructions.
#' 
#' @param mask Brain mask (matrix of 0 and 1 or \code{TRUE} and \code{FALSE}). 
#'
# @importFrom INLA inla.nonconvex.hull inla.mesh.2d inla.spde.make.A inla.spde2.matern
#' 
#' @return List containing INLA mesh, observation weight matrix \strong{A} for 
#'  translating between mesh locations and original data locations, the brain 
#'  mask used to create the mesh, and the number of original and mesh data 
#'  locations.
#' 
#' @export
make_mesh_2D <- function(mask){

  INLA_check()

  mask[] <- as.logical(mask)

  # Check only 0s and 1s
  if (!all(mask %in% c(TRUE, FALSE))) { 
    stop("Mask should be composed of only `FALSE` and `TRUE`")
  }

  xy.in <- which(mask, arr.ind=TRUE)[,2:1]
  boundary <- INLA::inla.nonconvex.hull(xy.in, resolution = 100)
  mesh <- INLA::inla.mesh.2d(loc = xy.in, boundary = boundary, max.edge = c(2, 4))
  Amat <- INLA::inla.spde.make.A(mesh, loc=xy.in)
  n.mask <- sum(mask)

  spde <- INLA::inla.spde2.matern(mesh)

  result <- list(
    mesh=mesh, A=Amat, spde=spde, mask=mask, n.mask = n.mask, n.mesh = mesh$n
  )
  class(result) <- 'templateICA_mesh_2D'
  result
}
