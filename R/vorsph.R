

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
#' A list is returned with `triangles` the index of triangles created from `xyz`,
#' `$edges` a matrix of index pairs into `$centroids`, a matrix
#' of XYZ coordinates that are the centres of the triangles of the convex hull of `xyz`.
#' @param xyz matrix of columns longitude, latitude, and optionally z
#'
#' @return list with matrix elements `triangles`, `edges` and `centroids`.
#' @export
#' @importFrom geometry convhulln
#' @examples
#' ## does the hard work (convex hull on xyz points is Delaunay triangulation)
#' library(geometry)
#' ## random longlat coordinates on the (geo)sphere
#' pts <- geosphere::randomCoordinates(5000)
#' library(rgl)
#' vs <- voronoi_spherical(ll2cart(pts))
#' clear3d()
#' segments3d(vs$centroids[t(vs$edges), ], lit = FALSE, lwd = 2)
#' spheres3d(0, 0, 0, 6310000, color = "lightgrey")
#' ## add the triangles
#' wire3d(tmesh3d(t(ll2cart(pts, rad = 6350000)), t(vs$triangles),
#'   homogeneous = FALSE), color = grey(0.3), alpha = .8)
#'
#'
#' ## regular coordinates
#'
#' ptsr <- regular_coords(row = 12)
#' vs <- voronoi_spherical(ll2cart(ptsr))
#' clear3d()
#' segments3d(vs$centroids[t(vs$edges), ], lit = FALSE, lwd = 2)
#'
#' ptss <- regular_coords(row = 1e5, prob = 0.000002)
#' vs <- voronoi_spherical(ll2cart(ptss))
#' clear3d()
#' segments3d(vs$centroids[t(vs$edges), ], lit = FALSE)
#'
#'
#' #pts <- do.call(cbind, maps::map(plot = FALSE)[c("x", "y")])
#' #pts <- pts[!is.na(pts[,1]), ]
#' #pts <- pts[sort(sample(1:nrow(pts), 8000)), ]
voronoi_spherical <- function(xyz) {
  if (ncol(xyz) < 3) {
    stop("input 'xyz` must be a matrix of spherical coordinates, \n perhaps use 'll2cart()' first to create from longitude,latitude values")
  }
  ## the index of triangles
  idx <- geometry::convhulln(xyz)
  edge <- triangles_to_edge(idx)
  ## each = 3 is the triangles each edge belongs to, we then extract the triangle pairs each belongs to by its edge_id (using pmin/pmax sort to identify them)
  trs <- matrix(rep(1:nrow(idx), each = 3)[order(edge_ids(edge))],  ncol = 2L, byrow = TRUE)
  list(triangles = idx, edges = trs, centroids =  centres(xyz[t(idx), ]))
}

