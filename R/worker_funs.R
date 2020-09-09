## triangle centres in xyz (these are the end points of each Voronoi cell edge)
centres <- function(x) cbind(colMeans(matrix(x[,1], 3)), colMeans(matrix(x[,2], 3)), colMeans(matrix(x[,3], 3)))


## convert to edges
triangles_to_edge <- function(x) {
  out <- matrix(t(cbind(cbind(x[,1], x[,2]),
                        cbind(x[,2], x[,3]),
                        cbind(x[,3], x[,1]))),
                ncol = 2L, byrow = TRUE)
  colnames(out) <- c(".vx1", ".vx2")
  out
}
edge_ids <- function(x) {
  as.integer(factor(paste(pmin(x[,1L], x[,2L]), pmax(x[,1L], x[,2L]), sep = "-")))
}
