

# generate data from a sparse matrix
# first compute covariance matrix
S = matrix(0.7, nrow = 5, ncol = 5)
for (i in 1:5){
  for (j in 1:5){
    S[i, j] = S[i, j]^abs(i - j)
  }
}

# generate 100 x 5 matrix with rows drawn from iid N_p(0, S)
Z = matrix(rnorm(100*5), nrow = 100, ncol = 5)
out = eigen(S, symmetric = TRUE)
S.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)
X = Z %*% S.sqrt

# calculate sample covariance
(nrow(X) - 1)/nrow(X)*cov(X)

# elastic-net type penalty (use CV for optimal lambda and alpha)
expect_error(CVglasso(X), NA)
expect_warning(CVglasso(X), NA)

# lasso penalty (lam = 0.1)
expect_error(CVglasso(X, lam = 0.1), NA)
expect_warning(CVglasso(X, lam = 0.1), NA)

expect_error(CVglasso(S = S, lam = 0.1), NA)
expect_warning(CVglasso(S = S, lam = 0.1), NA)

# parallel CV
expect_error(CVglasso(X, cores = 2), NA)
expect_warning(CVglasso(X, cores = 2), NA)

# adjmaxit
expect_error(CVglasso(X, adjmaxit = 2), NA)
expect_warning(CVglasso(X, adjmaxit = 2), NA)

# parallel adjmaxit
expect_error(CVglasso(X, adjmaxit = 2, cores = 2), NA)
expect_warning(CVglasso(X, adjmaxit = 2, cores = 2), NA)

# path
expect_error(CVglasso(X, path = TRUE), NA)
expect_warning(CVglasso(X, path = TRUE), NA)

expect_error(CVglasso(S = S, path = TRUE), NA)
expect_warning(CVglasso(S = S, path = TRUE), NA)
