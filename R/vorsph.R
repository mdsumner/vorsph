

#' Longitude latitude (and z) to Cartesian XYZ
#'
#' @param llh matrix of longitude, latitude, and optionally z (height above surface)
#' @param rad radius of sphere, assumed Earthish
#' @param exag exaggeration factor for height values
#'
#' @return matrix of XYZ coordinates (3 columns)
#' @export
#'
#' @examples
#' ll2cart(cbind(c(-180, 0, 90, 180), c(-90, -42, 0, 45, 90)))
ll2cart <- function (llh, rad = 6378137, exag = 1)
{
  if (ncol(llh) == 2L) llh <- cbind(llh, 0)
  cosLat = cos(llh[, 2] * pi/180)
  sinLat = sin(llh[, 2] * pi/180)
  cosLon = cos(llh[, 1] * pi/180)
  sinLon = sin(llh[, 1] * pi/180)
  rad <- (exag * llh[, 3] + rad)
  x = rad * cosLat * cosLon
  y = rad * cosLat * sinLon
  z = rad * sinLat
  cbind(x, y, z)
}

#' Spherical Voronoi
#'
#' Generate a Voronoi mesh on spherical coordinates
#'
#' An input matrix of XYZ coordinates is triangulated using convex hull (Delaunay on the sphere),
#' coordinates are assumed to lie on the surface of a sphere.
#'
#' A list is returned with `$edges` a matrix of index pairs into `$centroids`, a matrix
#' of XYZ coordinates that are the centres of the triangles of the convex hull of `xyz`.
#' @param xyz matrix of columns longitude, latitude, and optionally z
#'
#' @return list with matrix elements `edges` and `centroids`.
#' @export
#' @importFrom geometry convhulln
#' @examples
#' ## does the hard work (convex hull on xyz points is Delaunay triangulation)
#' library(geometry)
#' ## random longlat coordinates on the (geo)sphere
#' pts <- geosphere::randomCoordinates(8000)
#' library(rgl)
#' vs <- voronoi_spherical(ll2cart(pts))
#' clear3d()
#' segments3d(vs$centroids[t(vs$edges), ])
#'
#' #pts <- do.call(cbind, maps::map(plot = FALSE)[c("x", "y")])
#' #pts <- pts[!is.na(pts[,1]), ]
#' #pts <- pts[sort(sample(1:nrow(pts), 8000)), ]
voronoi_spherical <- function(xyz) {
  ## the index of triangles
  idx <- geometry::convhulln(xyz)
  edge <- triangles_to_edge(idx)
  ## each = 3 is the triangles each edge belongs to, we then extract the triangle pairs each belongs to by its edge_id (using pmin/pmax sort to identify them)
  trs <- matrix(rep(1:nrow(idx), each = 3)[order(edge_ids(edge))],  ncol = 2L, byrow = TRUE)
  list(edges = trs, centroids =  centres(xyz[t(idx), ]))
}
