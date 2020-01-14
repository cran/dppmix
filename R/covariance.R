# Apply covariance function on all pairs of rows of X
# C[x, x] = 1  since x - x = 0
cov_matrix_rows <- function(X, hparams) {
  K <- nrow(X);
  C <- matrix(1.0, K, K);

  # TODO Optimize in Rcpp
  if (K > 1) {
    # matrix is symmetric
    for (i in 1:(K-1)) {
      for (j in (i+1):K) {
        cij <- exp(-sum( similarity_kernel(X[i,], X[j,], hparams) ));
        C[i,j] <- C[j,i] <- cij;
      }
    }
  }

  C
}

# Covariance function for a row vector against each row of a matrix.
cov_f_rvec_mat <- function(rvec, mat, hparams) {
  # replicate row vector by row-wise to match dimension of mat
  rmat <- matrix(rvec, nrow(mat), ncol(mat), byrow=TRUE);
  matrix(exp(-rowSums( similarity_kernel(rmat, mat, hparams) )), nrow=1)
}

# Apply covariance function on all pairs of columns of X
# C[x, x] = 1  since x - x = 0
cov_matrix_cols <- function(X, hparams) {
  K <- ncol(X);
  C <- matrix(1.0, K, K);

  # TODO Optimize in Rcpp
  if (K > 1) {
    # matrix is symmetric
    for (i in 1:(K-1)) {
      for (j in (i+1):K) {
        cij <- exp(-sum( similarity_kernel(X[,i], X[,j], hparams) ));
        C[i,j] <- C[j,i] <- cij;
      }
    }
  }

  C
}

# Covariance function for a column vector against each column of a matrix.
cov_f_cvec_mat <- function(cvec, mat, hparams) {
  # replicate column vector will automaticlaly be replicated column-wise to
  # match dimension of the matrix
  matrix(exp(-colSums( similarity_kernel(cvec, mat, hparams) )), nrow=1)
}

# Similiarty kernel with the covariance function
similarity_kernel <- function(x, y, hparams) {
  ((x - y) / hparams$theta)^2
}
