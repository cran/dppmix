#' @export
estimate.dppmix_mcmc <- function(object, pars, ...) {
  i <- which_min_error(object);

  if (missing(pars)) {
    pars <- c("z", "K");
  }

  object[[i]][pars]
}

# Find mcmc sample that minimizes error
which_min_error <- function(mcmc) {
  m <- length(mcmc);
  N <- length(mcmc[[1]]$z);

  Hmatrix <- matrix(0, N, N);
  tmp <- array(0, c(N, N, m));
  for (i in 1:m) {
    tmp1 <- mcmc[[i]]$z %*% t(mcmc[[i]]$z);
    # look for squares: 1*1, 2*2, 3*3, ...
    id <- sqrt(tmp1) %in% 1:N;
    tmp1[id] <- 1;
    tmp1[-id] <- 0;
    tmp[,,i] <- tmp1;
    Hmatrix <- Hmatrix + tmp[,,i];
  }
  Hmatrix <- Hmatrix / m;

  d <- rep(0, m);
  for (i in 1:m) {
    d[i] <- sum((Hmatrix - tmp[,,i])^2);
  }

  which.min(d)
}

