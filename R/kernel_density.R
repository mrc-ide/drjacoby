#------------------------------------------------
# drop points into bins in 2D
#' @noRd
bin_2d <- function(x, y, x_breaks, y_breaks) {
  
  # check inputs
  assert_numeric(x)
  assert_vector(x)
  assert_numeric(y)
  assert_vector(y)
  assert_same_length(x, y)
  assert_numeric(x_breaks)
  assert_vector(x_breaks)
  assert_numeric(y_breaks)
  assert_vector(y_breaks)
  
  # get number of breaks in each dimension
  nx <- length(x_breaks)
  ny <- length(y_breaks)
  
  # create table of binned values
  tb <- table(findInterval(x, x_breaks), findInterval(y, y_breaks))
  
  # convert to dataframe and force numeric
  df <- as.data.frame(tb, stringsAsFactors = FALSE)
  names(df) <- c("x", "y", "count")
  df$x <- as.numeric(df$x)
  df$y <- as.numeric(df$y)
  
  # subset to within breaks range
  df <- subset(df, x > 0 & x < nx & y > 0 & y < ny)
  
  # fill in matrix
  z <- matrix(0, ny-1, nx-1)
  z[cbind(df$y, df$x)] <- df$count
  
  # return output as list
  ret <- list(x_mids = midpoints(x_breaks),
              y_mids = midpoints(y_breaks),
              z = z)
  return(ret)
}

#------------------------------------------------
# produce a smooth surface using 2D kernel density smoothing
#' @importFrom stats dist dnorm median
#' @noRd
kernel_smooth <- function(x, y, breaks_x, breaks_y, lambda = NULL) {
  
  # check inputs
  assert_numeric(x)
  assert_numeric(y)
  assert_same_length(x, y)
  assert_numeric(breaks_x)
  assert_numeric(breaks_y)
  if (!is.null(lambda)) {
    assert_single_pos(lambda, zero_allowed = FALSE)
  }
  
  # set default value of lambda
  if (is.null(lambda)) {
    lambda <- 0.2*median(dist(cbind(x, y)))
  }
  
  # get properties of cells in each dimension
  cells_x <- length(breaks_x) - 1
  cells_y <- length(breaks_y) - 1
  centre_x <- mean(breaks_x)
  centre_y <- mean(breaks_y)
  cell_size_x <- diff(breaks_x[1:2])
  cell_size_y <- diff(breaks_y[1:2])
  
  # bin x/y values in two dimensions and check that at least one value in
  # chosen region
  surface_raw <- bin_2d(x, y, breaks_x, breaks_y)$z
  if (all(surface_raw == 0)) {
    stop('chosen window contains no posterior draws')
  }
  
  # temporarily add guard rail to surface to avoid Fourier series bleeding round
  # edges
  rail_size_x <- cells_x
  rail_size_y <- cells_y
  rail_mat_x <- matrix(0, cells_x, rail_size_y)
  rail_mat_y <- matrix(0, rail_size_y, cells_x + 2*rail_size_x)
  
  surface_normalised <- surface_raw/sum(surface_raw)
  surface_normalised <- cbind(rail_mat_x, surface_normalised, rail_mat_x)
  surface_normalised <- rbind(rail_mat_y, surface_normalised, rail_mat_y)
  
  # calculate Fourier transform of posterior surface
  f1 = fftwtools::fftw2d(surface_normalised)
  
  # produce surface over which kernel will be calculated. This surface wraps
  # around in both x and y (i.e. the kernel is actually defined over a torus)
  kernel_x <- cell_size_x * c(0:floor(ncol(surface_normalised)/2), floor((ncol(surface_normalised) - 1)/2):1)
  kernel_y <- cell_size_y * c(0:floor(nrow(surface_normalised)/2), floor((nrow(surface_normalised) - 1)/2):1)
  kernel_x_mat <- outer(rep(1,length(kernel_y)), kernel_x)
  kernel_y_mat <- outer(kernel_y, rep(1,length(kernel_x)))
  kernel_s_mat <- sqrt(kernel_x_mat^2 + kernel_y_mat^2)
  
  # define kernel
  kernel <- dnorm(kernel_s_mat, sd = lambda)
  f2 = fftwtools::fftw2d(kernel)
  
  # combine Fourier transformed surfaces and take inverse. f4 is the surface of
  # interest
  f3 = f1*f2
  f4 = Re(fftwtools::fftw2d(f3, inverse = T))/length(surface_normalised)
  
  # remove guard rail
  f4 <- f4[,(rail_size_x+1):(ncol(f4)-rail_size_x)]
  f4 <- f4[(rail_size_y+1):(nrow(f4)-rail_size_y),]
  
  # return surface
  return(f4)
}
