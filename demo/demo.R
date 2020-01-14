library(devtools)

load_all()


# Simulate data ##############################################################

set.seed(1)

ns <- c(15, 8, 7);
means <- list(c(-6, 1), c(-1, -1), c(0, 4));
#means <- list(c(0, 0), c(0, 0), c(0, 0));

d <- rmvnorm_clusters(ns, means);
X <- d$X;
cl <- d$cl;



# Set parameters #############################################################

# hyperparameters
hparams <- list(
  a0 = 2,
  b0 = 2,
  delta = 1,
  sigma_pro_mu = 0.2,
  theta = 5
);

# MCMC control parameters
control <- list(
  niter = 100,
  burnin = 50,
  thin = 2
);


# Fit model ##################################################################

mcmc <- dppmix_mvnorm(X, hparams, control=control);
res <- estimate(mcmc);
print(res)

# Plot results ###############################################################

plot(X[,1], X[,2], pch=cl)
plot(X[,1], X[,2], pch=res$z)

# confusion matrix
# label swap is expected
table(cl, res$z)

Ks <- unlist(lapply(mcmc, function(x) x$K));
hist(Ks, xlim=c(0, nrow(X)))
print(table(Ks))


# Write output ###############################################################

dir.create("out", showWarnings=FALSE)
saveRDS(mcmc, file="out/mcmc_dpp-mvnorm.rds");

cluster <- data.frame(index = 1:nrow(X), truth = cl, cluster = res$z);
write.table(cluster, file = "out/clusters.tsv", sep="\t", row.names=FALSE, quote=FALSE)

