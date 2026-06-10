# =============================================================================
# gmm_em.R
# From-scratch Gaussian Mixture Model (GMM) via Expectation-Maximization (EM)
#
# Author: Sam A. Sievertsen, Oregon Health & Science University
#
# Purpose:
#   Implement GMM/EM entirely from first principles (no mclust / mixtools /
#   flexmix etc.), as required by the ML course final project. We create the
#   Gaussian log-density, the E-step responsibilities, the M-step parameter
#   updates, the observed-data log-likelihood, model selection via BIC, and
#   multi-restart initialization. Only base R + the `stats` linear-algebra
#   primitives (chol, forwardsolve) are used; these are standard numerical
#   building tools.
#
# Design notes:
#   - We work on the log scale and use the log-sum-exp trick throughout so that
#     responsibilities and the log-likelihood are numerically stable even when
#     component densities aren't
#   - Each M-step adds a small ridge (lambda * I) to every component covariance.
#     With ~800 person-days nested in ~32 people and 6 correlated features, a
#     component can otherwise collapse onto a near-singular covariance; the ridge
#     guarantees positive-definiteness and well-posed inverses
#   - Features are assumed standardized (z-scored) BEFORE fitting, since our raw
#     inputs are on very different scales (0-100 sliders, proportions, ratios).
#     Standardization is the caller's job (see utils.R), not the models
# =============================================================================


## 0. Small numerical helpers ##

#0.1 Log-sum-exp over a numeric vector (stable normalization constant)
log_sum_exp <- function(v) {
  m <- max(v)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(v - m)))
}

#0.2 Multivariate Gaussian LOG-density for every row of X, via Cholesky
#    Returns a length-n vector of log N(x_i | mu, Sigma). Implemented by hand:
#    log N = -0.5 * (d*log(2*pi) + log|Sigma| + mahalanobis(x, mu, Sigma)).
#    We solve through the Cholesky factor rather than inverting Sigma directly,
#    which is both faster and more numerically stable
dmvnorm_log <- function(X, mu, Sigma) {
  d <- length(mu)

  #0.2.1 Cholesky: R's chol() returns upper-triangular U with t(U) %*% U = Sigma
  U <- tryCatch(chol(Sigma), error = function(e) NULL)
  if (is.null(U)) {
    
    #0.2.1.1 Defensive fallback to nudge toward positive-definite, then retry once
    Sigma <- Sigma + diag(1e-6, d)
    U <- chol(Sigma)
  }

  #0.2.2 log|Sigma| = 2 * sum(log(diag(U)))
  log_det <- 2 * sum(log(diag(U)))

  #0.2.3 Center the data, then solve t(U) %*% z = (x - mu)' for each row.
  #      The squared column norms of z are exactly the Mahalanobis distances
  Xc <- sweep(X, 2L, mu, "-")
  z <- forwardsolve(t(U), t(Xc))          # d x n (lower-triangular solve)
  maha <- colSums(z^2)                    # length n

  -0.5 * (d * log(2 * pi) + log_det + maha)
}

#0.3 Univariate Gaussian LOG-density (vectorized); used by the 1-D EM
dnorm_log <- function(x, mu, sigma2) {
  -0.5 * (log(2 * pi) + log(sigma2) + (x - mu)^2 / sigma2)
}


## 1. Initialization (k-means++ style seeding) ##

#1.1 Choose K initial centers via k-means++ (spreads seeds out -> better optima)
#    Implemented from scratch: first center uniform-random, each subsequent center
#    drawn with probability proportional to squared distance from the nearest
#    chosen center. Returns a K x d matrix of seed means
kmeanspp_init <- function(X, K) {
  n <- nrow(X)
  centers <- matrix(NA_real_, nrow = K, ncol = ncol(X))

  #1.1.1 First center: a uniformly random observation
  idx <- sample.int(n, 1L)
  centers[1L, ] <- X[idx, ]

  #1.1.2 Remaining centers: distance-weighted sampling
  if (K >= 2L) {
    d2 <- colSums((t(X) - centers[1L, ])^2)
    for (k in 2:K) {
      probs <- d2 / sum(d2)
      idx <- sample.int(n, 1L, prob = probs)
      centers[k, ] <- X[idx, ]
      
      #1.1.2.1 Update nearest-center squared distance for every point
      d2_new <- colSums((t(X) - centers[k, ])^2)
      d2 <- pmin(d2, d2_new)
    }
  }
  centers
}


## 2. Multivariate GMM via EM (single fit) ##

#2.1 Fit one multivariate GMM with a single random start. Returns parameters,
#    responsibilities, the converged log-likelihood, BIC, and convergence info.
#    Arguments:
#      X : n x d numeric matrix (standardized features, complete cases)
#      K : number of mixture components
#      max_iter : EM iteration cap
#      tol : convergence tolerance on the change in log-likelihood
#      ridge : covariance regularization (lambda) added to each Sigma_k
gmm_em_fit_once <- function(X, K, max_iter = 500L, tol = 1e-6, ridge = 1e-4) {
  n <- nrow(X); d <- ncol(X)

  #2.1.1 Initialize means (k-means++), covariances (global cov), weights (uniform)
  mu <- kmeanspp_init(X, K)
  Sigma_global <- stats::cov(X) + diag(ridge, d)
  Sigma <- replicate(K, Sigma_global, simplify = FALSE)
  weights <- rep(1 / K, K)

  log_lik <- -Inf
  converged <- FALSE

  for (iter in seq_len(max_iter)) {

    #2.1.2 E-step: build the n x K matrix of log joint densities, then normalize
    #      each row with log-sum-exp to obtain responsibilities and the row's
    #      contribution to the observed-data log-likelihood.
    log_dens <- matrix(NA_real_, nrow = n, ncol = K)
    for (k in seq_len(K)) {
      log_dens[, k] <- log(weights[k]) + dmvnorm_log(X, mu[k, ], Sigma[[k]])
    }
    
    row_norm <- apply(log_dens, 1L, log_sum_exp)
    resp <- exp(log_dens - row_norm)          # n x K responsibilities
    new_log_lik <- sum(row_norm)

    #2.1.3 Convergence check on the log-likelihood
    if (is.finite(new_log_lik) && abs(new_log_lik - log_lik) < tol) {
      log_lik <- new_log_lik
      converged <- TRUE
      break
    }
    
    log_lik <- new_log_lik

    #2.1.4 M-step: effective counts, mixing weights, means, regularized covariances
    Nk <- colSums(resp)
    weights <- Nk / n
    for (k in seq_len(K)) {
      
      #2.1.4.1 Weighted mean
      mu[k, ] <- colSums(resp[, k] * X) / Nk[k]
      
      #2.1.4.2 Weighted covariance + ridge (guarantees positive-definiteness)
      Xc <- sweep(X, 2L, mu[k, ], "-")
      Sigma_k <- crossprod(Xc * resp[, k], Xc) / Nk[k]
      Sigma[[k]] <- Sigma_k + diag(ridge, d)
    }
  }

  #2.1.5 Free parameters for BIC: (K-1) weights + K*d means + K * d(d+1)/2 covs
  n_params <- (K - 1L) + K * d + K * (d * (d + 1L) / 2L)
  bic <- -2 * log_lik + n_params * log(n)
  aic <- -2 * log_lik + 2 * n_params

  list(
    K = K, weights = weights, mu = mu, Sigma = Sigma,
    resp = resp, log_lik = log_lik, bic = bic, aic = aic,
    n_params = n_params, iterations = iter, converged = converged)
}

#2.2 Fit a multivariate GMM with multiple random restarts; keep the best
#    (highest log-likelihood) solution. EM is only guaranteed to find a local
#    optimum, so restarts are essential for a trustworthy fit.
gmm_em_fit <- function(X, K, n_restarts = 20L, max_iter = 500L,
                       tol = 1e-6, ridge = 1e-4, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- as.matrix(X)

  best <- NULL
  
  for (r in seq_len(n_restarts)) {
    fit <- tryCatch(
      gmm_em_fit_once(X, K, max_iter = max_iter, tol = tol, ridge = ridge),
      error = function(e) NULL)
    
    if (is.null(fit)) next
    
    if (is.null(best) || (is.finite(fit$log_lik) && fit$log_lik > best$log_lik)) {
      best <- fit
    }
  }
  
  if (is.null(best)) stop("gmm_em_fit(): all restarts failed for K = ", K, ".")
  best$n_restarts <- n_restarts
  best
}


## 3. Univariate GMM via EM with 1-D stepping stone + per-feature EDA ##

#3.1 Fit a one-dimensional GMM with multiple restarts. This is both the
#    unit-testable kernel that the multivariate code generalizes, and a useful
#    EDA tool for asking whether a single feature is itself multimodal across
#    person-days (e.g., is commission-error rate bimodal?)
gmm_em_fit_1d <- function(x, K, n_restarts = 20L, max_iter = 500L,
                          tol = 1e-6, ridge = 1e-4, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  x <- as.numeric(x)
  n <- length(x)

  fit_one <- function() {
    
    #3.1.1 Init: random data points as means, global variance, uniform weights
    mu <- sample(x, K)
    sigma2 <- rep(stats::var(x) + ridge, K)
    weights <- rep(1 / K, K)
    log_lik <- -Inf; converged <- FALSE

    for (iter in seq_len(max_iter)) {
      
      #3.1.2 E-step
      log_dens <- sapply(seq_len(K), function(k) {
        log(weights[k]) + dnorm_log(x, mu[k], sigma2[k])
      })
      
      row_norm <- apply(log_dens, 1L, log_sum_exp)
      resp <- exp(log_dens - row_norm)
      new_log_lik <- sum(row_norm)

      if (is.finite(new_log_lik) && abs(new_log_lik - log_lik) < tol) {
        log_lik <- new_log_lik; converged <- TRUE; break
      }
      
      log_lik <- new_log_lik

      #3.1.3 M-step (variance floored by ridge to avoid degenerate spikes)
      Nk <- colSums(resp)
      weights <- Nk / n
      mu <- colSums(resp * x) / Nk
      sigma2 <- sapply(seq_len(K), function(k) {
        sum(resp[, k] * (x - mu[k])^2) / Nk[k] + ridge
      })
    }

    #3.1.4 Weights + means + variances
    n_params <- (K - 1L) + K + K            
    bic <- -2 * log_lik + n_params * log(n)
    list(K = K, weights = weights, mu = mu, sigma2 = sigma2, resp = resp,
         log_lik = log_lik, bic = bic, n_params = n_params,
         iterations = iter, converged = converged)
  }

  best <- NULL
  
  for (r in seq_len(n_restarts)) {
    fit <- tryCatch(fit_one(), error = function(e) NULL)
    if (is.null(fit)) next
    if (is.null(best) || fit$log_lik > best$log_lik) best <- fit
  }
  
  if (is.null(best)) stop("gmm_em_fit_1d(): all restarts failed for K = ", K, ".")
  best
}


## 4. Model selection over K (BIC) ##

#4.1 Fit the multivariate GMM across a range of K and return a tidy BIC/AIC
#    table plus the full fit objects. Selection is by minimum BIC, but the
#    caller should always inspect the BIC curve and the stability/interpretability
#    of the chosen solution rather than blindly trusting the minimum.
gmm_select_k <- function(X, k_range = 1:6, n_restarts = 20L, max_iter = 500L,
                         tol = 1e-6, ridge = 1e-4, seed = NULL) {
  X <- as.matrix(X)
  fits <- vector("list", length(k_range))
  names(fits) <- as.character(k_range)

  rows <- vector("list", length(k_range))
  for (i in seq_along(k_range)) {
    K <- k_range[i]
    
    #4.1.1 Stagger the seed per K so restarts are reproducible yet not identical
    seed_k <- if (is.null(seed)) NULL else seed + K
    fit <- gmm_em_fit(X, K, n_restarts = n_restarts, max_iter = max_iter,
                      tol = tol, ridge = ridge, seed = seed_k)
    fits[[i]] <- fit
    rows[[i]] <- data.frame(
      K = K, log_lik = fit$log_lik, n_params = fit$n_params,
      bic = fit$bic, aic = fit$aic, converged = fit$converged,
      iterations = fit$iterations)
  }

  bic_table <- do.call(rbind, rows)
  best_k <- bic_table$K[which.min(bic_table$bic)]

  list(bic_table = bic_table, fits = fits, best_k = best_k,
       best_fit = fits[[as.character(best_k)]])
}


## 5. Post-processing helpers (what we DO with the fitted states) ##

#5.1 Hard (MAP) state assignment for each observation from a fitted GMM
gmm_map_assign <- function(fit) {
  apply(fit$resp, 1L, which.max)
}

#5.2 Maximum responsibility per observation (a soft-assignment confidence)
gmm_max_resp <- function(fit) {
  apply(fit$resp, 1L, max)
}

#5.3 Person x state cross-tabulation (the critical artifact check):
#    if any state is dominated by a single person's days, that "state" is really
#    a person-identity artifact rather than a shared psychophysiological state.
#    Returns counts and the within-state share of the most-represented person.
gmm_person_state_check <- function(fit, person_id) {
  state <- gmm_map_assign(fit)
  tab <- table(person_id, state)
  state_totals <- colSums(tab)
  max_share <- apply(tab, 2L, function(col) max(col) / sum(col))
  list(
    table = tab,
    state_totals = state_totals,
    
    #5.3.1 share per state; high => artifact
    max_person_share = max_share)
}
