#' Random generator for the Bernoulli distribution.
#'
#' @param n     number of samples to generate
#' @param prob  event probability
#' @return an \code{integer} vector of 0 (non-event) and 1 (event)
rbern <- function(n, prob) {
  rbinom(n, size=1, prob)
}

#' Random generator for the Dirichlet distribution.
#'
#' @param n      number of vectors to generate
#' @param alpha  vector of parameters of the Dirichlet distribution
#' @return a \code{matrix} in which each row vector is Dirichlet distributed
rdirichlet <- function(n, alpha) {
  l <- length(alpha);
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE);
  sm <- x %*% rep(1, l);
  
  x / as.vector(sm)
}

#' Density function for Gamma-Poisson distribution.
#'
#' Data follow the Poisson distribution parameterized by a mean parameter
#' that follows a gamma distribution.
#' 
#' @param x    vector of x values
#' @param a    shape parameter for gamma distribution on mean parameter
#' @param b    rate parameter for gamma distribution on mean parameter
#' @param log  whether to return the density in log scale
#' @return density values
dgammapois <- function(x, a, b=1, log=FALSE) {
  f <- a*log(b) + lgamma(a + x) - lgamma(a) - (a + x)*log(b + 1) - lgamma(x + 1);
  if (log) {
    f
  } else {
    exp(f)
  }
}

#' Generate a random binary vector.
#'
#' @param n      size of binary vector
#' @param prob   event probability
#'               (not accounting for minimum event constraint)
#' @param e.min  minimum number of events
#' @return an \code{integer} vector of 0 and 1
rbvec <- function(n, prob, e.min=0) {
  e <- min(rpois(1, n * prob) + e.min, n);
  idx <- sample(1:n, e);
  x <- integer(n);
  x[idx] <- 1;
  x
}

