#' Generate random multivarate clusters
#'
#' @param ns  number of data points in each cluster
#' @param means  centers of each cluster
#' @return list containing matrix \code{X} and labels \code{cl}
#' @examples
#' ns <- c(5, 8, 7)
#' means <- list(c(-6, 1), c(-1, -1), c(0, 4))
#' d <- rmvnorm_clusters(ns, means)
#' @export
rmvnorm_clusters <- function(ns, means) {
  cl <- unlist(
    mapply(
      function(id, cl) rep(id, cl),
      1:length(ns),
      ns
    )
  );

  X <- do.call(rbind,
    mapply(
      function(n, mu) mvtnorm::rmvnorm(n, mu),
      ns, means, SIMPLIFY=FALSE
    )
  );

  idx <- sample(1:length(cl));
  X <- X[idx, ];
  cl <- cl[idx];

  list(X=X, cl=cl)
}
