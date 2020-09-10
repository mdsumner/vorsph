initbin <- function(NUMROWS = 2160) {
  ## TODO options for lon-lat sub-sets
  latbin <- (((seq(NUMROWS) - 1) + 0.5) * 180 / NUMROWS ) - 90
  ## this will overflow at 2^31-1
  #numbin <- as.integer(2 * NUMROWS * cos(latbin * pi/180) + 0.5)
  numbin <- trunc(2 * NUMROWS * cos(latbin * pi/180) + 0.5)
  basebin <- cumsum(c(1L, numbin[-length(numbin)]))
  totbins = basebin[NUMROWS] + numbin[NUMROWS] - 1
  list(latbin = latbin, numbin = numbin, basebin = basebin, totbins = totbins)
}

bin2lonlat <- function(bins, nrows) {
  row <- seq_len(nrows) - 1
  latbin = ((row + 0.5)*180.0/nrows) - 90.0;
  numbin <- trunc((2*nrows*cos(latbin*pi/180.0) + 0.5))
  basebin <- c(1L, utils::head(cumsum(numbin) + 1L, -1L))
  totbins <- utils::tail(basebin, 1) + utils::tail(numbin, 1) - 1
  index <- findInterval(bins, basebin)
  lat <- latbin[index]
  lon <- 360.0*(bins - basebin[index] + 0.5)/numbin[index] - 180.0;
  cbind(lon, lat)
}


#' Generate regular longitude latitude coordinates
#'
#' Input is a number of `rows`, default is quite low.
#'
#' `rows` argument is how many distinct latitude bands.
#'
#' Code is taken from the croc package, which bundles logic from the Ocean Color
#' NASA website. https://github.com/sosoc/croc
#'
#' @param rows Number of rows of coordinates see Details
#' @param force set to `TRUE` to ignore restraint and create a lot of points
#' @param prob defaults to 1, otherwise set to a fraction to subsample the returned points
#'
#' @return matrix of longitude,latitude
#' @export
#'
#' @examples
#' plot(regular_coords(10))
#'
#'
regular_coords <- function(rows = 24, force = FALSE, prob = 1) {
  binit <- initbin(rows)
  if (prob <= 0) stop("'prob' must be between 0, and 1")
  if ((binit$totbins * prob) > 1e6 && !force) {
    message(sprintf("'rows = %i' is A LOT of points (%0.f), are you sure?\ Use 'prob' to get a smaller fraction, or 'force = TRUE' to ignore restraint.",
                 as.integer(rows),
                 round(binit$totbins * prob)))

    stop()
  }
  if (prob < 1 && prob > 0) {
    bins <- sort(sample(binit$totbins, max(c(1, binit$totbins * prob))))
  } else {
    bins <- seq_len(binit$totbins)
  }
  bin2lonlat(bins, rows)
}
