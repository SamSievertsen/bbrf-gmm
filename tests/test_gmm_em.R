# =============================================================================
# test_gmm_em.R
# Simulation-based verification of the from-scratch GMM/EM implementation
#
# Logic: if the EM is correct, then on data simulated from KNOWN Gaussian
# components it must (a) recover the true component means and weights to within
# Monte-Carlo tolerance, and (b) have BIC prefer the true number of components.
# We test both the multivariate and univariate kernels. This is the unit test
# that licenses trusting the model on the real BBRF data
# =============================================================================

source(file.path("src", "gmm_em.R"))

set.seed(123)

## 1. Simulate from 3 known, moderately separated 3-D Gaussians ##

#1.1 Ground-truth parameters
true_K <- 3L
true_weights <- c(0.5, 0.3, 0.2)
true_mu <- rbind(
  c(0, 0, 0),
  c(4, 4, 0),
  c(-3, 2, 5))

true_Sigma <- list(
  diag(c(1.0, 1.5, 1.0)),
  matrix(c(1.0, 0.5, 0.0,
           0.5, 1.0, 0.0,
           0.0, 0.0, 1.0), 3, 3),
  diag(c(1.5, 1.0, 2.0)))

#1.2 Draw n observations according to the mixing weights
n <- 1500L
comp <- sample.int(true_K, n, replace = TRUE, prob = true_weights)

#1.3 Cholesky-based sampler (from scratch: standard normals -> correlated normals)
rmvnorm_chol <- function(n, mu, Sigma) {
  d <- length(mu)
  U <- chol(Sigma)
  Z <- matrix(rnorm(n * d), nrow = n, ncol = d)
  sweep(Z %*% U, 2L, mu, "+")
}

X <- matrix(NA_real_, nrow = n, ncol = 3L)

for (k in seq_len(true_K)) {
  idx <- which(comp == k)
  X[idx, ] <- rmvnorm_chol(length(idx), true_mu[k, ], true_Sigma[[k]])
}


## 2. Fit and check parameter recovery at the true K ##

fit <- gmm_em_fit(X, K = true_K, n_restarts = 30L, seed = 1)

#2.1 Greedily match estimated components to true components by nearest mean
match_components <- function(est_mu, true_mu) {
  K <- nrow(true_mu)
  used <- integer(0)
  mapping <- integer(K)
  
  for (j in seq_len(K)) {
    dists <- apply(est_mu, 1L, function(m) sqrt(sum((m - true_mu[j, ])^2)))
    dists[used] <- Inf
    pick <- which.min(dists)
    mapping[j] <- pick
    used <- c(used, pick)
  }
  mapping
}

map <- match_components(fit$mu, true_mu)
mu_err <- max(abs(fit$mu[map, ] - true_mu))
w_err <- max(abs(fit$weights[map] - true_weights))

cat("\n=== Parameter recovery (K = 3) ===\n")
cat(sprintf("Converged: %s in %d iterations\n", fit$converged, fit$iterations))
cat(sprintf("Max abs error in component MEANS:   %.4f\n", mu_err))
cat(sprintf("Max abs error in mixing WEIGHTS:    %.4f\n", w_err))


## 3. Check BIC selects the true number of components ##

sel <- gmm_select_k(X, k_range = 1:6, n_restarts = 20L, seed = 7)
cat("\n=== BIC model selection ===\n")
print(sel$bic_table, row.names = FALSE)
cat(sprintf("BIC-selected K: %d (true K: %d)\n", sel$best_k, true_K))


## 4. Univariate kernel: recover a known 2-component 1-D mixture ##

set.seed(42)
x1 <- c(rnorm(800, mean = -2, sd = 1), rnorm(400, mean = 3, sd = 0.7))
fit1d <- gmm_em_fit_1d(x1, K = 2L, n_restarts = 30L, seed = 5)
mu1d_sorted <- sort(fit1d$mu)
cat("\n=== Univariate recovery (true means -2 and 3) ===\n")
cat(sprintf("Recovered means: %.3f, %.3f\n", mu1d_sorted[1], mu1d_sorted[2]))


## 5. Overall pass/fail gate ##

pass_mu <- mu_err < 0.30
pass_w  <- w_err < 0.06
pass_k  <- sel$best_k == true_K
pass_1d <- abs(mu1d_sorted[1] - (-2)) < 0.25 && abs(mu1d_sorted[2] - 3) < 0.25

cat("\n=== TEST SUMMARY ===\n")
cat(sprintf(" [%s] multivariate means recovered\n", ifelse(pass_mu, "PASS", "FAIL")))
cat(sprintf(" [%s] mixing weights recovered\n", ifelse(pass_w, "PASS", "FAIL")))
cat(sprintf(" [%s] BIC selected true K\n", ifelse(pass_k, "PASS", "FAIL")))
cat(sprintf(" [%s] univariate means recovered\n", ifelse(pass_1d, "PASS", "FAIL")))

if (all(c(pass_mu, pass_w, pass_k, pass_1d))) {
  cat("\nALL TESTS PASSED\n")
} else {
  stop("One or more EM verification tests FAILED.")
}
