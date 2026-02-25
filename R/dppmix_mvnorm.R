#' Fit a determinantal point process multivariate normal mixture model.
#' 
#' Discover clusters in multidimensional data using a multivariate normal 
#' mixture model with a determinantal point process prior.
#'
#' A determinantal point process (DPP) prior is a repulsive prior.
#' Compare to mixture models using independent priors, a DPP mixutre model
#' will often discover a parsimonious set of mixture components (clusters).
#'
#' Model fitting is done by sampling parameters from the posterior
#' distribution using a reversible jump Markov chain Monte Carlo sampling
#' approach.
#' 
#' Given \eqn{X = [x_i]}, where each \eqn{x_i} is a D-dimensional real vector,
#' we seek the posterior distribution the latent variable \eqn{z = [z_i]}, where
#' each \eqn{z_i} is an integer representing cluster membership.
#'
#' \deqn{ x_i \mid z_i  \sim  Normal(\mu_k, \Sigma_k) }
#' \deqn{ z_i           \sim  Categorical(w) }
#' \deqn{ w             \sim  Dirichlet([\delta ... \delta]) }
#' \deqn{ \mu_k         \sim  DPP(C) }
#'
#' where \eqn{C} is the covariance function that evaluates the distances among the
#' data points:
#'
#' \deqn{ C(x_1, x_2) = exp( - \sum_d \frac{ (x_1 - x_2)^2 }{ \theta^2 } ) }
#'
#' We also define \eqn{\Sigma_k = E_k \Lambda_k E_k^\top}, where \eqn{E_k} is an
#' orthonormal matrix whose column represents eigenvectors.
#' We further assume that \eqn{E_k = E} is fixed across all cluster components 
#' so that \eqn{E} can be estimated as the eigenvectors of the covariance matrix of
#' the data matrix \eqn{X}. Finally, we put a prior on the entries of the
#' \eqn{\Lambda_k} diagonal matrix:
#'
#' \deqn{ \lambda_{kd}^{-1}  \sim  Gamma( a_0, b_0 ) }
#'
#' Hence, the hyperameters of the model include:
#' \code{delta, a0, b0, theta}, as well as sampling hyperparameter
#' \code{sigma_pro_mu}, which controls the spread of the Gaussian
#' proposal distribution for the random-walk Metropolis-Hastings update of
#' the \eqn{\mu} parameter.
#'
#' The parameters (and their dimensions) in the model include:
#' \code{K}, \code{z (N x 1)}, \code{w (K x 1)}, \code{lambda (K x J)},
#' \code{mu (K x J)}, \code{Sigma (J x J x K)}.
#' If any parameter is fixed, then \code{K} must be fixed as well.
#'
#' @param X        \code{N x J} data matrix of \code{N} observations and
#'                 \code{J} features
#' @param hparams  a list of hyperparameter values:
#'                 \code{delta, a0, b0, theta, sigma_prop_mu}
#' @param store    a vector of character strings specifying additional vars of
#'                 interest; a value of \code{NA} indicates that
#'                 samples of all parameters in the model will be stored
#' @param control  a list of control parameters:
#'                 \code{niter, burnin, thin}
#' @param fixed    a list of fixed parameter values
#' @param verbose  whether to emit verbose message
#'
#' @return a \code{dppmix_mcmc} object containing posterior samples of
#'         the parameters
#'
#' @references Yanxun Xu, Peter Mueller, Donatello Telesca.
#'             Bayesian Inference for Latent Biologic Structure with
#'             Determinantal Point Processes.
#'             Biometrics. 2016;72(3):955-64.
#'
#' @examples
#' set.seed(1)
#' ns <- c(3, 3)
#' means <- list(c(-6, -3), c(0, 4))
#' d <- rmvnorm_clusters(ns, means)
#'
#' mcmc <- dppmix_mvnorm(d$X, verbose=FALSE)
#' res <- estimate(mcmc)
#' table(d$cl, res$z)
#'
#' @import stats mvtnorm
#' @export
dppmix_mvnorm <- function(X, hparams=NULL, store=NULL, control=NULL, fixed=NULL, verbose=TRUE) {
  stime <- proc.time();

  hparams.default <- list(
    a0 = 2,
    b0 = 2,
    delta = 1,
    sigma_pro_mu = 0.4,
    theta = 5
  );

  hparams <- c(hparams, hparams.default);

  control.default <- list(
    niter = 100,
    burnin = 50,
    thin = 2
  );

  control <- c(control, control.default);

  E <- dppmix_mvnorm_compute_E(X);

  smpl <- dppmix_mvnorm_initialize(X, E, fixed);

  if (is.null(store)) {
    store <- c("z", "K");
  } else if (length(store) == 1 && is.na(store)) {
    store <- c(names(smpl), "Sigma");
  }

  if (control$burnin <= 0) {
    mcmc <- list(smpl);
  } else {
    mcmc <- list();
  }

  Sigma <- dppmix_mvnorm_update_Sigma(smpl$lambda, smpl$K, E);

  i <- 1; 
  for (iter in 2:control$niter) {
    if (is.null(fixed$z)) {
      smpl <- dppmix_mvnorm_update_z(smpl$K, smpl$mu, smpl$w, smpl$lambda, Sigma, X, E);
      # mu may also be updated due to empty components
    }

    if (is.null(fixed$w)) {
      smpl$w <- dppmix_mvnorm_update_w(smpl$z, hparams);
    }

    if (is.null(fixed$lambda)) {
      smpl$lambda <- dppmix_mvnorm_update_lambda(smpl$K, smpl$z, smpl$mu, X, E, hparams);
    }

    if (is.null(fixed$Sigma)) {
      Sigma <- dppmix_mvnorm_update_Sigma(smpl$lambda, smpl$K, E);
    }

    if (is.null(fixed$mu)) {
      smpl$mu <- dppmix_mvnorm_update_mu(smpl$K, smpl$z, smpl$mu, smpl$lambda, Sigma, X, hparams);
    }

    if (is.null(fixed$K)) {
      smpl <- dppmix_mvnorm_update_K(smpl$K, smpl$z, smpl$mu, smpl$w, smpl$lambda, Sigma, X, E, hparams)
      # all other parameters will change if K changes
    }

    Sigma <- smpl$Sigma;

    if (iter > control$burnin && iter %% control$thin == 0) {
      mcmc[[i]] <- smpl[store];
      i <- i + 1;
    }

    if (verbose) {
      message("iter ", iter, "  K = ", smpl$K)
    }
  }

  if (verbose) {
    message("Elapsed time: ", proc.time()[3] - stime[3])
  }

  structure(mcmc, class="dppmix_mcmc")
}

# Initialize the MCMC chain with an estimate.
dppmix_mvnorm_initialize <- function(X, E, fixed) {
  N <- nrow(X);
  J <- ncol(X);

  # z is N x 1
  # mu is K x J
  # w is K x 1
  # lambda is K x J

  # infer fixed values of K from fixed values of parameters
  if (is.null(fixed$K)) {
    if (!is.null(fixed$mu)) {
      fixed$K <- nrow(fixed$mu);
    } else if (!is.null(fixed$w)) {
      fixed$K <- length(fixed$w);
    } else if (!is.null(fixed$lambda)) {
      fixed$K <- nrow(fixed$lambda);
    } else if (!is.null(fixed$z)) {
      fixed$K <- length(unique(fixed$z));
    }
  }

  # K is variable but K <= N
  if (is.null(fixed$K)) {
    K <- ceiling(log(nrow(X)));
  } else {
    K <- fixed$K;
  }

  if (! all(c("z", "mu", "w") %in% names(fixed)) ) {
    # initialize clusters with K-means
    cl <- kmeans(X, K);
  }

  if (is.null(fixed$z)) {
    z <- cl$cluster;
  } else {
    z <- fixed$z;
  }

  if (is.null(fixed$mu)) {
    mu <- cl$centers;
  } else {
    mu <- fixed$mu;
  }

  if (is.null(fixed$w)) {
    w <- cl$size / sum(cl$size);
  } else {
    w <- fixed$w;
  }

  if (is.null(fixed$lambda)) {
    lambda <- matrix(1.0, K, J);
  } else {
    lambda <- fixed$lambda;
  }

  list(
    K = K,
    z = z,
    mu = mu,
    w = w,
    lambda = lambda
  )
}

# Compute eigenvectors of the data covariance matrix.
# These eigenvectors are used to re-construct Sigma.
dppmix_mvnorm_compute_E <- function(X) {
  X.c <- scale(X, center=TRUE, scale=FALSE);
  svd(X.c)$v
}

# Update z.
#
# Re-sample z from complete conditional posterior:
# p(z_i = k \mid ...) = \frac{1}{Z} w_k Normal(y_i; \mu_k, \Sigma_k)
# where numeric calculation of the normalization constant Z is trivial.
#
# Update can cause K to change; in turn, mu can also change.
# lambda and w can also change in size, but they will be updated imminently.
# Sigma can also change, but we will update it after updating lambda.
dppmix_mvnorm_update_z <- function(K, mu, w, lambda, Sigma, X, E) {
  N <- nrow(X);
  ks <- 1:K;

  # compute unnormlized density
  # ds is N x K
  # bottleneck
  ds <- matrix(
    unlist(lapply(ks,
      function(k) {
        w[k] * likelihood(X, mu[k, ], Sigma[,  , k])
      }
    )),
    nrow = nrow(X)
  );

  # re-sample the assignment of each data point to component
  # bottleneck
  z <- apply(ds, 1,
    function(d) {
      # `sample` will normalize d to sum to 1
      sample(ks, 1, prob=d)
    }
  );

  # due to sampling, some components may now be empty:
  # no data points assigned
  components <- sort(unique(z));

  K.new <- length(components);

  # remove empty components
  mu <- mu[ks %in% components, , drop=FALSE];

  # re-code the index s.t. empty components are skipped
  for (k in 1:K.new) {
    z[z == components[k]] <- k;
  }

  list(K=K.new, z=z, mu=mu)
}

# Update Sigma.
#
# Re-construct Sigma from eigenvalues (lambda) and eigenvectors (fixed E),
# for each possible K.
# Sigma = E diag(\lambda) E^\top
dppmix_mvnorm_update_Sigma <- function(lambda, K, E) {
  J <- nrow(E);
  Sigma <- array(0, c(J, J, K));
  for (k in 1:K) {
    Sigma[, , k] <- E %*% diag(lambda[k, ], ncol=J, nrow=J) %*% t(E);
  }

  Sigma
}

# Update w.
#
# Re-sample from the complete conditional posterior.
# w \sim Dirichlet(\delta + m)
# where m is the vector of component sizes:
#       m_k = \sum_{i=1}^N I(z_i = k)
dppmix_mvnorm_update_w <- function(z, hparams) {
  # z must be have components 1 ... K
  # s.t. table(z) has size K
  m <- as.numeric(table(z));
  rdirichlet(1, hparams$delta + m)
}

# Compute Sk matrix (J x J) by definition
dppmix_mvnorm_compute_Sk_v1 <- function(X, mu, idx, k) {
  J <- ncol(X);
  Sk <- matrix(0, J, J);
  for (i in 1:length(idx)) {
    d <- X[idx[i],] - mu[k,];
    Sk <- Sk + (d %*% t(d));
  }

  Sk
}

# Compute Sk matrix (J x J)
dppmix_mvnorm_compute_Sk_v2 <- function(X, mu, idx, k) {
  # vector mu is automatically expanded columnwise
  # ensure that x stays as a matrix
  d <- t(X[idx, , drop=FALSE]) - mu[k, ];
  d %*% t(d)
}

dppmix_mvnorm_compute_Sk <- dppmix_mvnorm_compute_Sk_v2;

# Update lambda (eigenvalues of Sigma).
#
# Prior: 
# \lambda_{kd}^{-1} \sim Gamma( \frac{a}{2}, \frac{b}{2} )
#
# Complete conditional posterior:
# \lambda_{kd}^{-1} \sim 
#   Gamma( \frac{a + m_k}{2}, \frac{b + e_d^\top S_k e_d}{2} )
# where e_d is the d-th eigenvector (column of E)
#       S_k is a J x J distance matrix:
#       S_k = \sum_{i: z_i = k} (y_i - \mu_k) (y_i - \mu_k)^\top
dppmix_mvnorm_update_lambda <- function(K, z, mu, X, E, hparams) {
  J <- ncol(X);

  lambda <- matrix(0, K, J);
  for (k in 1:K) {
    idx <- which(z == k);

    # Consider: handle empty components
    # if (length(idx) == 0) next;
    
    Sk <- dppmix_mvnorm_compute_Sk(X, mu, idx, k);

    aa <- hparams$a0 + 1 + length(idx);
    bb <- hparams$b0 + diag(t(E) %*% Sk %*% E);
    lambda[k, ] <- 1 / rgamma(J, aa, bb);
  }

  lambda
}

# Update mu.
# 
# p(\mu_k \mid ...) \propto
#   det(C_{\mu_1 ... \mu_K})
#   \prod_{s_i = k} Normal( x_i ; \mu_k, \Sigma_k )
#
# After applying the Schur determinant identity, we get
#
# p(\mu_k \mid ...) \propto
#   ( C(\mu_k, \mu_k} - b C_{\mu_{-k}}^{-1} b^\top )
#   \prod_{s_i = k} Normal( x_i ; \mu_k, \Sigma_k )
# where
#   b = C(\mu_k, \mu_{-k})
#   \mu_{-k} = (\mu_j}_{j /neq k}
# 
# Since the normalization constant is not computable, we use
# a Metropolis-Hastings (random-walk) step using a multivariate
# normal proposal. Since this proposal distribution is symmetric,
# the forward and reverse proposals cancel out in the MH ratio.
dppmix_mvnorm_update_mu <- function(K, z, mu, lambda, Sigma, X, hparams) {
  J <- ncol(X);

  # covariance matrix for multivariate normal proposal distribution
  S.pro <- hparams$sigma_pro_mu^2 * diag(J);

  # covariance matrix for DPP
  # NB If a pair of mu elements are too close,
  #    then C would be computationally singular and non-invertible!
  C <- cov_matrix_rows(mu, hparams);

	# repeat until we get a valid mu
  repeat {

  if (K == 1) {

    mk <- nrow(X);

    # K == 1 => C_{\mu} = [C_{\mu_1, \mu_1}] = [1]
    # thus, det(C) = det([1]) = 1
    # and p(\mu_1 \mid ...) = \prod_i Normal( x_i ; \mu_1, \Sigma_1 )
    mu.k.old <- mu[1, , drop=FALSE];
    mu.k.new <- mvtnorm::rmvnorm(1, mu.k.old, S.pro);

    S.inv <- solve(Sigma[,,1]);
    D.old <- X - matrix(mu.k.old, mk, J, byrow=TRUE);
    D.new <- X - matrix(mu.k.new, mk, J, byrow=TRUE);
    # terms not in the exponent cancel out in the ratio
    ll.old <- -0.5 * sum(apply(D.old, 1, function(d) t(d) %*% S.inv %*% d));
    ll.new <- -0.5 * sum(apply(D.new, 1, function(d) t(d) %*% S.inv %*% d));
    ratio <- exp(ll.new - ll.old);

    if (runif(1) < ratio) {
      mu[1,] <- mu.k.new;
    }

  } else {

    for (k in 1:K) {
      # If K == 1, C is 1 x 1 and C[-k,-k] is 0 x 0, causing crash
      C.sinv <- solve(C[-k,-k]);

      mu.k.old <- mu[k, , drop=FALSE];

      # sample new mu_k
      # bottleneck
      mu.k.new <- mvtnorm::rmvnorm(1, mu.k.old, S.pro);

      idx <- which(z == k);
      mk <- length(idx);

      # calculate ratio of the likelihood of the new mu_k to the old mu_k
      S.inv <- solve(Sigma[,,k]);
      D.old <- X[idx,] - matrix(mu.k.old, mk, J, byrow=TRUE);
      D.new <- X[idx,] - matrix(mu.k.new, mk, J, byrow=TRUE);
      # terms not in the exponent cancel out in the ratio
      ll.old <- -0.5 * sum(apply(D.old, 1, function(d) t(d) %*% S.inv %*% d));
      ll.new <- -0.5 * sum(apply(D.new, 1, function(d) t(d) %*% S.inv %*% d));
      ratio.like <- exp(ll.new - ll.old);

      # b = C(\mu_k, \mu_{-k})
      b.old <- C[k, -k, drop=FALSE];
      # apply covariance function on new mu_k and mu for other components
      b.new <- cov_f_rvec_mat(mu.k.new, mu[-k, , drop=FALSE], hparams);

      # Since covariance function C is defined with only the similarity kernel,
      # C(\mu_k, \mu_k) = 1
      ratio <- (1 - b.new %*% C.sinv %*% t(b.new)) / (1 - b.old %*% C.sinv %*% t(b.old)) * ratio.like;

      # probabilistic step
      if (runif(1) < ratio) {
        mu[k,] <- mu.k.new;
      }
    }

  }  # if (K == 1)

	if (mu_is_valid(mu, hparams)) {
		break;
  }

  } # repeat

  mu
}

# Log likelihood function by definition.
# 
# p(x_i \mid z_i = k) = Normal( x_i; \mu_k, \Sigma_k)
# p(X \mid z) = \prod_i Normal( x_i; \mu_{z_i}, \Sigma_{z_i})
dppmix_mvnorm_llikelihood_v1 <- function(X, z, mu, Sigma) {
  ll <- 0;
  for (i in 1:length(z)) {
    zi <- z[i];
    ll <- ll + likelihood(X[i,], mu[zi,], Sigma[,,zi], log=TRUE);
  }

  sum(ll)
}

# Log likelihood function.
dppmix_mvnorm_llikelihood_v2 <- function(X, z, mu, Sigma) {
  K <- nrow(mu);
  sum(unlist(lapply(1:K,
    function(k) {
      likelihood(X[z == k,], mu[k,], Sigma[,,k], log=TRUE)
    }
  )))
}

# bottleneck
dppmix_mvnorm_llikelihood <- dppmix_mvnorm_llikelihood_v2;

dppmix_mvnorm_split_move <- function(K, z, mu, w, lambda, Sigma, X, E, hparams) {
  N <- nrow(X);
  J <- ncol(X);

  C.old <- cov_matrix_rows(mu, hparams);
  ll.old <- dppmix_mvnorm_llikelihood(X, z, mu, Sigma);

  K.new <- K + 1;
  k1 <- K.new - 1;
  k2 <- K.new;

	# repeat until we get a valid mu
	repeat {

  # randomly determine the terms of the divorce settlement
  alpha <- rbeta(1, 1, 1);  # E[alpha] = 0.5
  beta <- rbeta(J, 1, 1);   # E[beta_d] = 0.5

  # choose a component to split with uniform probability
  k <- sample(1:K, 1);
  k.idx <- which(z == k);
  mk <- length(k.idx);

  # NB the random sign of r is not described in Xu, Mueller, Telesca 2015
  r <- rbeta(J, 2, 2) * sample(c(-1, 1), J, replace=TRUE);

  # \mu_1^{new} = \mu - sqrt(w_2^{new} / w_1^{new}) \sum_d k_{1d}^{1/2} r_d e_d
  # \mu_2^{new} = \mu + sqrt(w_1^{new} / w_2^{new}) \sum_d k_{1d}^{1/2} r_d e_d
  # ensure that lambda[k,] remains a 1 x J row vector
  # dd will be a 1 x J row vector
  dd <- (sqrt(lambda[k, , drop=FALSE]) * r) %*% E;
  mu1 <- mu[k,] - sqrt( (1-alpha)/alpha ) * dd;
  mu2 <- mu[k,] + sqrt( alpha/(1-alpha) ) * dd;

  mu.new <- rbind(mu[-k, , drop=FALSE], mu1, mu2);

	if (mu_is_valid(mu.new, hparams)) {
		break;
	}

	} # repeat

	# update other variables

  # w_1^{new} = w_1 \alpha
  # w_2^{new} = w_1 (1 - \alpha)
  w1 <- w[k] * alpha;
  w2 <- w[k] * (1 - alpha);

  # \lambda_{1d}^{new} = \beta_d * (1 - r_d^2) * w_1 / w_1^{new} k_{1d}
  # \lambda_{2d}^{new} = (1 - \beta_d) * (1 - r_d^2) * w_1 / w_1^{new} k_{2d}
  lambda1 <- beta * (1 - r^2) / alpha * lambda[k,];
  lambda2 <- (1 - beta) * (1 - r^2) / (1 - alpha) * lambda[k,];

  # delete split component and append new components
  w.new <- c(w[-k], w1, w2);
  lambda.new = rbind(lambda[-k, , drop=FALSE], lambda1, lambda2);
  Sigma.new <- dppmix_mvnorm_update_Sigma(lambda.new, K.new, E);
  C.new <- cov_matrix_rows(mu.new, hparams);

  # allocate data points in split component into new components
  # based on complete conditional posterior
  p1 <- w1 * likelihood(X[k.idx,], mu1, Sigma.new[,,k1]);
  p2 <- w2 * likelihood(X[k.idx,], mu2, Sigma.new[,,k2]);
  # probability of assigning data point to first new component
  prob <- ifelse(p1 + p2 > 0, p1 / (p1 + p2), 0);

  k1.idx <- rbern(mk, prob);
  k1.idx <- k1.idx == 1;
  k2.idx <- !k1.idx;

  # update index
  z.new <- z;
  # delete selected component: shift components with greater index toward 1
  after.k <- z > k;
  z.new[after.k] <- z[after.k] - 1;
  # update index of new components
  z.new[k.idx[k1.idx]] <- k1;
  z.new[k.idx[k2.idx]] <- k2;

  # sizes of new components
  n1 <- sum(k1.idx);
  n2 <- mk - n1;

  ll.new <- dppmix_mvnorm_llikelihood(X, z.new, mu.new, Sigma.new);
  lratio_like <- ll.new - ll.old;

  # \lambda_{kd^{-1} \sim Gamma( a_0, b_0 )
  a <- hparams$a0;
  b <- hparams$b0;
  delta <- hparams$delta;

  lratio_prior <-
    # log p(\tilde{w})
    (delta - 1 + n1) * log(w1) + (delta - 1 + n2) * log(w2) - 
    # log p(w)
    ( (delta - 1 + mk) * log(w[k]) + lbeta(delta, K * delta) ) +
    # log p(\tilde{\mu})
    ldet(C.new) - 
    # log p(\mu)
    ldet(C.old) + 
    # log p(\tilde{\lambda}) - log p(\lambda)
    -J * lgamma(a) +
    sum(
      log(lambda1)*(1 - a) + log(lambda2)*(1 - a) + log(b)*(a) - 
        log(lambda[k,])*(1 - a)
    ) -
    b * sum(1/lambda1 + 1/lambda2 - 1/lambda[k,]);

  # proposal distributions for auxilary split parameters
  lratio_propose <- - dbeta(alpha, 1, 1, log=TRUE) - sum(dbeta(beta, 1, 1, log=TRUE)) -
    sum(dbeta(abs(r), 2, 2, log=TRUE));

  # ratio of proposing to combine vs. split
  # q_{K+1,d} / q_{K,u} = 1
  # since q_{K+1,d} = q_{K,u} = 0.5
  
  # ratio of prob. of combining (j_1, j_2) vs. prob. of splitting j
  # q_{K+1,c}(j_1, j_2) / q_{K,s}(j) = 1 / (K + 1)
  # where q_{K+1,c}(j_1, j_2) = 1 / (K*(K + 1))
  #       because pairs are sampled uniformly for combination,
  #       and q_{K,s}(j) = 1 / K
  #       because a component is sampled uniformly for splitting.
  lratio_move <- - log(K + 1);

  # proposal distribution to data point allocation to components
  lratio_alloc <- - sum(log(c(prob[k1.idx], 1 - prob[k2.idx])));

  ldet_jacobian <- (3*J + 1)*log(w[k]) - 1.5*J*log(w1*w2) + 
    sum(1.5*log(lambda[k,]) + log(1 - r^2));

  # extra factor (K + 1) in the denominator is related to the use of densities with
  # respect to a unit rate Poisson process
  # (see equation 3.8 in Xu, Mueller, Telesca 2015)
  ratio <- exp(lratio_like + lratio_prior - log(K + 1) + lratio_move + 
               lratio_alloc + lratio_propose + ldet_jacobian);

  list(ratio=ratio, K=K.new, z=z.new, mu=mu.new, w=w.new, lambda=lambda.new, Sigma=Sigma.new)
}

dppmix_mvnorm_combine_move <- function(K, z, mu, w, lambda, Sigma, X, E, hparams) {
  N <- nrow(X);
  J <- ncol(X);

  C.old <- cov_matrix_rows(mu, hparams);
  ll.old <- dppmix_mvnorm_llikelihood(X, z, mu, Sigma);

  K.new <- K - 1;

  # uniformly sample 2 components to combine
  kk <- sort(sample(1:K, 2));
  k1 <- kk[1];
  k2 <- kk[2];
  idx1 <- z == k1;
  idx2 <- z == k2;

  # sizes of components being combined
  n1 <- sum(idx1);
  n2 <- sum(idx2);

  # parameters of combined cluster
  wc <- w[k1] + w[k2];
  muc <- (w[k1]*mu[k1,] + w[k2]*mu[k2,]) / wc;
  alpha <- w[k1] / wc;
  lambdac <- alpha*lambda[k1,] + (1-alpha)*lambda[k2,] + 
    alpha*(1-alpha)*(mu[k1,]-mu[k2,])^2;

  # These probabilities were not mentioend in Xu, Mueller, Telesca 2015
  p11 <- w[k1] * likelihood(X[idx1,], mu[k1,], Sigma[,,k1]);
  p12 <- w[k2] * likelihood(X[idx1,], mu[k2,], Sigma[,,k2]);
  p21 <- w[k1] * likelihood(X[idx2,], mu[k1,], Sigma[,,k1]);
  p22 <- w[k2] * likelihood(X[idx2,], mu[k2,], Sigma[,,k2]);
  prob1 <- ifelse(p11 + p12 > 0, p11 / (p11 + p12), 0);
  prob2 <- ifelse(p21 + p22 > 0, p22 / (p21 + p22), 0);

  # NB Parameters u and beta were not mentioned for the combined move 
  #    in the Xu, Muller, Telesca 2015
  u <- abs(solve( t(sqrt(lambdac)*t(E)), (muc-mu[k1,])/sqrt(w[k2]/w[k1]) ));
  # NB Why does u not appear in the definition of beta?
  tt <- lambda[k1,] / lambda[k2,] * alpha/(1-alpha);
  beta <- tt / (tt + 1);

  # delete selected components and append combined component at the end
  w.new <- c(w[-kk], wc);
  mu.new <- rbind(mu[-kk, , drop=FALSE], muc);
  lambda.new <- rbind(lambda[-kk, , drop=FALSE], lambdac);
  Sigma.new <- dppmix_mvnorm_update_Sigma(lambda.new, K.new, E);
  C.new <- cov_matrix_rows(mu.new, hparams);

  # update index
  z.new <- z;
  after.k1 <- z > k1 & z <= k2;
  z.new[after.k1] <- z[after.k1] - 1;
  after.k2 <- z > k2;
  z.new[after.k2] <- z[after.k2] - 2;
  z.new[z == k1 | z == k2] <- K.new;
  
  ll.new <- dppmix_mvnorm_llikelihood(X, z.new, mu.new, Sigma.new);
  lratio_like <- ll.new - ll.old;

  lratio_alloc <- sum(c(log(prob1), log(prob2)));

  # ratio of prob. of splitting j vs. prob. of combining (j_1, j_2)
  # q_{K,s}(j) / q_{K+1,c}(j_1, j_2) = K + 1
  # where q_{K+1,c}(j_1, j_2) = 1 / (K*(K + 1))
  #       because pairs are sampled uniformly for combination,
  #       and q_{K,s}(j) = 1 / K
  #       because a component is sampled uniformly for splitting.
  lratio_move <- log(K + 1);

  # \lambda_{kd^{-1} \sim Gamma( a_0, b_0 )
  a <- hparams$a0;
  b <- hparams$b0;
  delta <- hparams$delta;

  # order: + new - old
  lratio_prior <-
    # log p(\tilde{w})
    (delta - 1 + n1 + n2)*log(wc) + lbeta(delta, K.new*delta) -
    # log p(w)
    ( (delta - 1 + n1)*log(w[k1]) + (delta - 1 + n2)*log(w[k2]) ) +
    # log p(\tilde{\mu})
    ldet(C.new) - 
    # log p(\mu)
    ldet(C.old) +
    # log p(\tilde{\lambda}) - log p(\lambda)
    -J * lgamma(a) +
    sum(
      log(lambdac)*(1 - a) - log(lambda[k1,])*(1 - a) -
      log(lambda[k2,])*(1 - a) - log(b)*(a)
    ) -
    b * sum(1/lambdac - 1/lambda[k1,] - 1/lambda[k2,]);
  
  lratio_propose <- dbeta(alpha, 1, 1, log=TRUE) + 
    sum(dbeta(beta, 1, 1, log=TRUE)) + sum(dbeta(abs(u), 2, 2, log=TRUE));
  
  ldet_jacobian <- -(3*J + 1)*log(wc) + 1.5*J*log(w[k1]*w[k2]) - 
    sum(log(abs(lambda[kk,]^1.5 * (1 - u^2))));

  ratio <- exp(lratio_like + lratio_prior + log(K) + lratio_move + 
               lratio_alloc + lratio_propose + ldet_jacobian);

  list(ratio=ratio, K=K.new, z=z.new, mu=mu.new, w=w.new, lambda=lambda.new, Sigma=Sigma.new)
}

# Update K using reversible jump.
# All other parameters can change too.
dppmix_mvnorm_update_K <- function(K, z, mu, w, lambda, Sigma, X, E, hparams) {
  K.min <- 1;
  m <- NULL;

  # split and combine moves have equal probability, subject to constraint on K
  if (runif(1) < 0.5) {
    # consider splitting a component into two
    m <- dppmix_mvnorm_split_move(K, z, mu, w, lambda, Sigma, X, E, hparams);

    # move with probability min(1, ratio)
    if (runif(1) >= m$ratio) {
      m <- NULL;
    }

  } else if (K > K.min) {
    # consider combining two components
    m <- dppmix_mvnorm_combine_move(K, z, mu, w, lambda, Sigma, X, E, hparams);

    # move with probability min(1, ratio)
    if (runif(1) >= m$ratio) {
      m <- NULL;
    }
  }

  if (!is.null(m)) {
    K <- m$K;
    z <- m$z;
    mu <- m$mu;
    w <- m$w;
    lambda <- m$lambda;
    Sigma <- m$Sigma;
  }
  
  list(K=K, z=z, mu=mu, w=w, lambda=lambda, Sigma=Sigma)
}

likelihood <- function(X, mu, Sigma, log=FALSE) {
  if (length(mu) > 1) {
    mvtnorm::dmvnorm(X, mu, Sigma, log=log)
  } else {
    dnorm(X, mu, Sigma, log=log)
  }
}

ldet <- function(x) {
  as.numeric(determinant(x)$modulus)
}

# check that no pair of values are mu are too close, causing C to be singular
mu_is_valid <- function(mu, hparams) {
  # check if covariance matrix is invertible
  C <- cov_matrix_rows(mu, hparams);
  if (inherits(try(solve(C)), "try-error")) {
		message("mu = ")
		print(as.numeric(mu))
		FALSE
	} else {
		TRUE
	}
}
