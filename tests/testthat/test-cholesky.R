# ============================================================================
# test-cholesky.R — testthat tests for cholesky-transform.R
#
# Sections:
#   A. Algebraic correctness    (A Sigma A' = V_sand, L factors)
#   B. Mean preservation        (exact to machine precision)
#   C. Variance recovery        (finite-sample tolerance)
#   D. Wald CI math             (width, symmetry, significance, labels)
#   E. Shrinkage direction      (PI > 10 => A_dd < 1, PI ~ 1 => A_dd ~ 1)
#   F. Input validation         (classes, levels, dimensions)
#   G. PD handling              (ridge, nearPD, singularity)
#   H. Edge cases               (P=1, identity transform, multiple levels)
#   I. Print method             (output, sections, corrupted, digits)
#   J. Internal helpers         (.make_stan_param_names, .verify, .ensure_pd)
#   K. Stress tests             (large P=10, small M=10, degenerate draws)
#   L. Comparison table / DER   (DER, PI, width ratio diagnostics)
# ============================================================================


# ---- Mock object constructors -----------------------------------------------

#' Create a mock hbb_fit object for testing
#'
#' Generates draws from MVN(theta_true, Sigma_true) and wraps in
#' a CmdStanR-like fit interface.
#'
#' @param P Number of covariates per margin.
#' @param M Number of MCMC draws.
#' @param seed Random seed.
#' @param theta_true Optional true mean vector of length 2P+1.
#' @param Sigma_true Optional true covariance matrix (2P+1 x 2P+1).
#' @return An S3 object of class "hbb_fit".
create_mock_fit <- function(P = 3, M = 2000, seed = 42,
                            theta_true = NULL, Sigma_true = NULL) {
  D <- 2L * P + 1L
  param_names <- c(
    paste0("alpha[", seq_len(P), "]"),
    paste0("beta[", seq_len(P), "]"),
    "log_kappa"
  )

  set.seed(seed)
  if (is.null(theta_true)) theta_true <- rnorm(D, 0, 0.3)
  if (is.null(Sigma_true)) Sigma_true <- diag(D) * 0.5 + 0.1

  draws <- MASS::mvrnorm(M, theta_true, Sigma_true)
  colnames(draws) <- param_names

  mock_cmdstan <- list(
    draws = function(variables = NULL, format = "matrix") {
      if (!is.null(variables)) draws[, variables, drop = FALSE]
      else draws
    }
  )

  structure(
    list(fit = mock_cmdstan, model_type = "weighted", hbb_data = list(P = P)),
    class = "hbb_fit"
  )
}

#' Create a mock hbb_sandwich object for testing
#'
#' @param P Number of covariates per margin.
#' @param Sigma_MCMC MCMC covariance matrix. If NULL, uses default.
#' @param V_sand Sandwich variance. If NULL, uses default (much smaller).
#' @param H_obs_inv Inverse Hessian. If NULL, uses default.
#' @return An S3 object of class "hbb_sandwich".
create_mock_sandwich <- function(P = 3, Sigma_MCMC = NULL, V_sand = NULL,
                                  H_obs_inv = NULL) {
  D <- 2L * P + 1L
  if (is.null(Sigma_MCMC)) Sigma_MCMC <- diag(D) * 0.5 + 0.1
  if (is.null(V_sand))     V_sand     <- diag(D) * 0.01 + 0.001
  if (is.null(H_obs_inv))  H_obs_inv  <- diag(D) * 0.005

  plabs <- c(
    paste0("alpha_", c("intercept", paste0("x", seq_len(P - 1)))),
    paste0("beta_",  c("intercept", paste0("x", seq_len(P - 1)))),
    "log_kappa"
  )

  structure(
    list(
      V_sand       = V_sand,
      H_obs        = solve(H_obs_inv),
      H_obs_inv    = H_obs_inv,
      J_cluster    = V_sand,
      Sigma_MCMC   = Sigma_MCMC,
      DER          = diag(V_sand) / diag(H_obs_inv),
      param_labels = plabs,
      D            = D,
      P            = P,
      N            = 1000L,
      model_type   = "weighted"
    ),
    class = "hbb_sandwich"
  )
}

#' Suppress cli messages during testing
quietly <- function(expr) {
  suppressMessages(expr)
}


# ==========================================================================
# Section A: Algebraic correctness — A Sigma_MCMC A' = V_sand
# ==========================================================================

test_that("A. covariance identity holds exactly (diagonal case)", {
  fit <- create_mock_fit(P = 3, M = 2000)
  sw  <- create_mock_sandwich(P = 3)
  res <- quietly(cholesky_correct(fit, sw))

  # Exact matrix identity (not statistical): A %*% Sigma_MCMC %*% t(A) = V_sand
  A_Sigma_At <- res$A %*% sw$Sigma_MCMC %*% t(res$A)
  expect_equal(A_Sigma_At, sw$V_sand, tolerance = 1e-10)
  expect_true(res$verification$A_pass)
})

test_that("A. covariance identity with dense off-diagonal structure", {
  P <- 3
  D <- 2L * P + 1L
  set.seed(99)
  R1 <- matrix(rnorm(D * D), D, D)
  Sigma_MCMC <- crossprod(R1) + diag(D)
  R2 <- matrix(rnorm(D * D), D, D)
  V_sand <- crossprod(R2) * 0.01 + diag(D) * 0.01

  fit <- create_mock_fit(P = P, M = 3000, seed = 77)
  sw  <- create_mock_sandwich(P = P, Sigma_MCMC = Sigma_MCMC, V_sand = V_sand)
  res <- quietly(cholesky_correct(fit, sw))

  A_Sigma_At <- res$A %*% Sigma_MCMC %*% t(res$A)
  expect_equal(A_Sigma_At, V_sand, tolerance = 1e-10)
})

test_that("A. Cholesky factors are lower triangular", {
  fit <- create_mock_fit(P = 3, M = 500)
  sw  <- create_mock_sandwich(P = 3)
  res <- quietly(cholesky_correct(fit, sw))

  # Upper triangle (excluding diagonal) must be zero
  expect_true(all(res$L_MCMC[upper.tri(res$L_MCMC)] == 0))
  expect_true(all(res$L_sand[upper.tri(res$L_sand)] == 0))

  # Cholesky reconstruction: L L' = original matrix
  expect_equal(res$L_MCMC %*% t(res$L_MCMC), sw$Sigma_MCMC, tolerance = 1e-12)
  expect_equal(res$L_sand %*% t(res$L_sand), sw$V_sand, tolerance = 1e-12)
})

test_that("A. return class and dimensions", {
  P <- 3
  D <- 2 * P + 1
  M <- 1000

  fit <- create_mock_fit(P = P, M = M)
  sw  <- create_mock_sandwich(P = P)
  res <- quietly(cholesky_correct(fit, sw))

  expect_s3_class(res, "hbb_cholesky")
  expect_equal(nrow(res$theta_corrected), M)
  expect_equal(ncol(res$theta_corrected), D)
  expect_equal(length(res$theta_hat), D)
  expect_equal(nrow(res$A), D)
  expect_equal(ncol(res$A), D)
  expect_equal(dim(res$L_MCMC), c(D, D))
  expect_equal(dim(res$L_sand), c(D, D))
  expect_equal(res$D, D)
  expect_equal(res$M, M)
  expect_equal(res$P, P)
})


# ==========================================================================
# Section B: Mean preservation — exact to machine precision
# ==========================================================================

test_that("B. mean preservation is exact", {
  fit <- create_mock_fit(P = 3, M = 4000)
  sw  <- create_mock_sandwich(P = 3)
  res <- quietly(cholesky_correct(fit, sw))

  mean_corrected <- colMeans(res$theta_corrected)
  expect_equal(unname(mean_corrected), unname(res$theta_hat), tolerance = 1e-12)
  expect_true(res$verification$mean_pass)
  expect_true(res$verification$mean_preservation_error < 1e-10)
})

test_that("B. mean preservation with large P and high PI", {
  P <- 5
  D <- 2L * P + 1L
  M <- 3000

  # High prior domination: Sigma_MCMC >> V_sand
  Sigma_MCMC <- diag(D) * 2.0 + 0.3
  V_sand <- diag(D) * 0.005 + 0.0005

  fit <- create_mock_fit(P = P, M = M, Sigma_true = Sigma_MCMC, seed = 7)
  sw  <- create_mock_sandwich(P = P, Sigma_MCMC = Sigma_MCMC, V_sand = V_sand)
  res <- quietly(cholesky_correct(fit, sw))

  expect_true(res$verification$mean_pass)
  mean_corrected <- colMeans(res$theta_corrected)
  expect_true(max(abs(unname(mean_corrected) - unname(res$theta_hat))) < 1e-10)
})


# ==========================================================================
# Section C: Variance recovery — finite-sample tolerance
# ==========================================================================

test_that("C. corrected sample variance approximates V_sand diagonals", {
  fit <- create_mock_fit(P = 3, M = 10000, seed = 123)
  sw  <- create_mock_sandwich(P = 3)
  res <- quietly(cholesky_correct(fit, sw))

  var_corrected <- apply(res$theta_corrected, 2, var)
  var_target    <- diag(sw$V_sand)
  rel_diff <- abs(var_corrected - var_target) / var_target
  expect_true(all(rel_diff < 0.05))
})

test_that("C. corrected full covariance within absolute tolerance", {
  fit <- create_mock_fit(P = 3, M = 10000, seed = 123)
  sw  <- create_mock_sandwich(P = 3)
  res <- quietly(cholesky_correct(fit, sw))

  Sigma_corr <- cov(res$theta_corrected)
  scale <- max(diag(sw$V_sand))
  max_abs_err <- max(abs(Sigma_corr - sw$V_sand))
  expect_true(max_abs_err / scale < 0.10,
              label = paste("scaled abs err:", signif(max_abs_err / scale, 3)))
})

test_that("C. off-diagonal covariance recovery for non-diagonal case", {
  set.seed(22)
  P <- 3
  D <- 2 * P + 1
  M <- 5000

  # Generate full PD matrices with off-diagonal structure
  L1 <- matrix(rnorm(D * D), D, D) * 0.5
  Sigma_MCMC <- L1 %*% t(L1) + diag(D) * 0.01
  L2 <- matrix(rnorm(D * D), D, D) * 0.05
  V_sand <- L2 %*% t(L2) + diag(D) * 0.001

  fit <- create_mock_fit(P = P, M = M, Sigma_true = Sigma_MCMC, seed = 22)
  sw  <- create_mock_sandwich(P = P, Sigma_MCMC = Sigma_MCMC, V_sand = V_sand)
  res <- quietly(cholesky_correct(fit, sw))

  expect_true(res$verification$A_pass)
  Cov_corrected <- cov(res$theta_corrected)
  off_diag_error <- max(abs(Cov_corrected[upper.tri(Cov_corrected)] -
                             sw$V_sand[upper.tri(sw$V_sand)])) /
    max(abs(sw$V_sand[upper.tri(sw$V_sand)]))
  expect_true(off_diag_error < 0.15,
    info = paste("Off-diagonal relative error:", round(off_diag_error, 4)))
})


# ==========================================================================
# Section D: Wald CI math
# ==========================================================================

test_that("D. Wald CI width equals 2 * z * se", {
  theta_hat <- c(alpha_1 = -0.324, beta_1 = 0.090, log_kappa = 1.92)
  V_diag <- c(0.0027, 0.00035, 0.0041)
  V_sand <- diag(V_diag)

  for (level in c(0.90, 0.95, 0.99)) {
    wald <- compute_wald_ci(theta_hat, V_sand, level = level)
    z_crit <- qnorm(1 - (1 - level) / 2)
    expected_width <- 2 * z_crit * sqrt(V_diag)
    expect_equal(wald$ci_width, expected_width, tolerance = 1e-12,
      info = paste("Level:", level))
  }
})

test_that("D. Wald CI is symmetric around theta_hat", {
  theta_hat <- c(a = 1.5, b = -0.5, c = 0.0)
  V_sand <- diag(c(0.01, 0.04, 0.09))

  wald <- compute_wald_ci(theta_hat, V_sand, level = 0.95)
  midpoints <- (wald$ci_lo + wald$ci_hi) / 2
  expect_equal(midpoints, as.numeric(theta_hat), tolerance = 1e-14)
})

test_that("D. Wald CI significance consistent with p-value", {
  set.seed(20)
  D <- 10
  theta_hat <- c(rep(0.001, 5), rep(5, 5))
  V_sand <- diag(D) * 0.1

  wald <- compute_wald_ci(theta_hat, V_sand, level = 0.95)

  for (i in seq_len(D)) {
    if (wald$significant[i]) {
      expect_true(wald$p_value[i] < 0.05)
    } else {
      expect_true(wald$p_value[i] >= 0.05)
    }
  }
})

test_that("D. Wald CI uses param_labels when provided", {
  theta_hat <- c(1.0, 2.0)
  V_sand <- diag(c(0.1, 0.2))
  labels <- c("poverty_ext", "poverty_int")

  wald <- compute_wald_ci(theta_hat, V_sand, param_labels = labels)
  expect_equal(wald$parameter, labels)
})

test_that("D. Wald CI uses names from theta_hat when param_labels NULL", {
  theta_hat <- c(alpha = 1.0, beta = 2.0)
  V_sand <- diag(c(0.1, 0.2))
  wald <- compute_wald_ci(theta_hat, V_sand)
  expect_equal(wald$parameter, c("alpha", "beta"))
})

test_that("D. Wald CI generates default labels when none available", {
  theta_hat <- c(1.0, 2.0, 3.0)
  names(theta_hat) <- NULL
  V_sand <- diag(c(0.1, 0.2, 0.3))
  wald <- compute_wald_ci(theta_hat, V_sand)
  expect_equal(wald$parameter, c("param_1", "param_2", "param_3"))
})

test_that("D. Wald CI nested intervals across levels", {
  theta_hat <- c(1.0, -0.5)
  V_sand <- diag(c(0.04, 0.09))

  w90 <- compute_wald_ci(theta_hat, V_sand, level = 0.90)
  w95 <- compute_wald_ci(theta_hat, V_sand, level = 0.95)
  w99 <- compute_wald_ci(theta_hat, V_sand, level = 0.99)

  for (i in seq_along(theta_hat)) {
    expect_true(w90$ci_lo[i] >= w95$ci_lo[i])
    expect_true(w90$ci_hi[i] <= w95$ci_hi[i])
    expect_true(w95$ci_lo[i] >= w99$ci_lo[i])
    expect_true(w95$ci_hi[i] <= w99$ci_hi[i])
  }
})

test_that("D. Wald CI has correct output columns", {
  set.seed(10)
  D <- 5
  theta_hat <- rnorm(D)
  V_sand <- diag(runif(D, 0.01, 0.5))

  wald <- compute_wald_ci(theta_hat, V_sand, level = 0.95)

  expect_equal(nrow(wald), D)
  expected_cols <- c("parameter", "post_mean", "se", "z_stat", "p_value",
                     "ci_lo", "ci_hi", "ci_width", "significant")
  expect_true(all(expected_cols %in% names(wald)))

  # z_stat = theta_hat / se
  se <- sqrt(diag(V_sand))
  expect_equal(wald$z_stat, theta_hat / se, tolerance = 1e-14)
  expect_equal(wald$se, se, tolerance = 1e-14)
})

test_that("D. Wald CI handles non-positive diagonal gracefully", {
  D <- 5
  set.seed(55)
  theta_hat <- rnorm(D)
  V_sand <- diag(D)
  V_sand[3, 3] <- 0  # zero variance for param 3

  wald <- suppressMessages(compute_wald_ci(theta_hat, V_sand))

  # SE clamped to 0, so CI width = 0
  expect_equal(wald$se[3], 0)
  expect_equal(wald$ci_width[3], 0)
})


# ==========================================================================
# Section E: Shrinkage direction
# ==========================================================================

test_that("E. high PI yields A diagonal < 1 (shrinkage)", {
  P <- 3
  D <- 2 * P + 1

  Sigma_MCMC <- diag(D) * 5.0       # Large MCMC variance
  V_sand     <- diag(D) * 0.01      # Small sandwich variance
  H_obs_inv  <- diag(D) * 0.05      # PI = 5.0/0.05 = 100

  fit <- create_mock_fit(P = P, M = 2000, Sigma_true = Sigma_MCMC, seed = 55)
  sw  <- create_mock_sandwich(
    P = P, Sigma_MCMC = Sigma_MCMC, V_sand = V_sand, H_obs_inv = H_obs_inv
  )
  res <- quietly(cholesky_correct(fit, sw))

  a_diag <- diag(res$A)
  expect_true(all(a_diag < 1),
    info = paste("A diag:", paste(round(a_diag, 4), collapse = ", ")))

  PI <- res$comparison_table$prior_inflation
  expect_true(all(PI > 10))

  expect_true(all(res$comparison_table$width_ratio < 1))
})

test_that("E. low PI (V_sand > Sigma_MCMC) yields A diagonal > 1 (expansion)", {
  P <- 3
  D <- 2L * P + 1L

  Sigma_MCMC <- diag(D) * 0.001
  V_sand     <- diag(D) * 0.05

  fit <- create_mock_fit(P = P, M = 2000, seed = 404)
  sw  <- create_mock_sandwich(P = P, Sigma_MCMC = Sigma_MCMC, V_sand = V_sand)
  res <- quietly(cholesky_correct(fit, sw))

  a_diag <- diag(res$A)
  expect_true(all(a_diag > 1),
    info = paste("A diag:", paste(signif(a_diag, 3), collapse = ", ")))
})

test_that("E. PI ~ 1 produces A diagonal near 1", {
  P <- 2
  D <- 2 * P + 1

  Sigma_MCMC <- diag(D) * 0.01
  V_sand     <- diag(D) * 0.01
  H_obs_inv  <- diag(D) * 0.01

  fit <- create_mock_fit(P = P, M = 2000, Sigma_true = Sigma_MCMC, seed = 77)
  sw  <- create_mock_sandwich(
    P = P, Sigma_MCMC = Sigma_MCMC, V_sand = V_sand, H_obs_inv = H_obs_inv
  )
  res <- quietly(cholesky_correct(fit, sw))

  a_diag <- diag(res$A)
  expect_true(all(abs(a_diag - 1) < 0.01),
    info = paste("A diag:", paste(round(a_diag, 4), collapse = ", ")))
})


# ==========================================================================
# Section F: Input validation
# ==========================================================================

test_that("F. rejects non-hbb_fit object", {
  sw <- create_mock_sandwich(P = 3)
  expect_error(cholesky_correct("not_a_fit", sw), "hbb_fit")
  expect_error(cholesky_correct(list(a = 1), sw), "hbb_fit")
})

test_that("F. rejects non-hbb_sandwich object", {
  fit <- create_mock_fit(P = 3)
  expect_error(cholesky_correct(fit, list(a = 1)), "hbb_sandwich")
})

test_that("F. rejects sandwich with missing fields", {
  fit <- create_mock_fit(P = 3)
  bad_sw <- structure(list(V_sand = diag(7), D = 7L), class = "hbb_sandwich")
  expect_error(cholesky_correct(fit, bad_sw), "missing")
})

test_that("F. rejects mismatched P between fit and sandwich", {
  fit <- create_mock_fit(P = 3)
  sw  <- create_mock_sandwich(P = 4)
  expect_error(cholesky_correct(fit, sw), "mismatch")
})

test_that("F. rejects invalid level values", {
  fit <- create_mock_fit(P = 3)
  sw  <- create_mock_sandwich(P = 3)
  expect_error(cholesky_correct(fit, sw, level = 0), "\\(0, 1\\)")
  expect_error(cholesky_correct(fit, sw, level = 1), "\\(0, 1\\)")
  expect_error(cholesky_correct(fit, sw, level = -0.5), "\\(0, 1\\)")
  expect_error(cholesky_correct(fit, sw, level = "a"), "single numeric")
  expect_error(cholesky_correct(fit, sw, level = c(0.9, 0.95)), "single numeric")
  expect_error(cholesky_correct(fit, sw, level = NA), "single numeric")
})

test_that("F. rejects D vs 2P+1 mismatch in sandwich", {
  sw <- create_mock_sandwich(P = 2)
  sw$D <- 99  # inconsistent: P=2 => D should be 5
  fit <- create_mock_fit(P = 2)
  expect_error(cholesky_correct(fit, sw), "mismatch")
})

test_that("F. detects NULL fit$fit", {
  fit <- create_mock_fit()
  fit$fit <- NULL
  sw <- create_mock_sandwich()
  expect_error(quietly(cholesky_correct(fit, sw)), "NULL")
})

test_that("F. compute_wald_ci rejects invalid inputs", {
  expect_error(compute_wald_ci("not numeric", diag(2)), "numeric")
  expect_error(compute_wald_ci(c(1, 2), "not matrix"), "matrix")
  expect_error(compute_wald_ci(c(1, 2), diag(3)), "2 x 2")
  expect_error(compute_wald_ci(c(1, 2), diag(2), level = 2), "\\(0, 1\\)")
  expect_error(
    compute_wald_ci(c(1, 2), diag(2), param_labels = c("a", "b", "c")),
    "length 2"
  )
})


# ==========================================================================
# Section G: PD handling
# ==========================================================================

test_that("G. .ensure_cholesky_pd passes through already-PD matrix", {
  mat <- diag(5) * 2 + 0.1
  result <- .ensure_cholesky_pd(mat, "test")
  expect_equal(result$mat, mat)
  expect_false(result$corrected)
  expect_true(result$min_eig > 0)
})

test_that("G. .ensure_cholesky_pd ridge method fixes non-PD matrix", {
  D <- 3
  mat <- diag(D)
  mat[1, 1] <- -0.01  # Make non-PD

  result <- suppressMessages(.ensure_cholesky_pd(mat, "test_mat", method = "ridge"))

  expect_true(result$corrected)
  expect_true(!is.null(result$details))
  expect_true(grepl("ridge", tolower(result$details)))
  # Result must be PD (Cholesky should succeed)
  expect_no_error(chol(result$mat))
  eig <- eigen(result$mat, symmetric = TRUE)$values
  expect_true(all(eig > 0))
})

test_that("G. .ensure_cholesky_pd nearpd method fixes non-PD matrix", {
  D <- 3
  mat <- diag(c(1, 0.5, -0.05))

  result <- suppressMessages(.ensure_cholesky_pd(mat, "test_mat", method = "nearpd"))

  expect_true(result$corrected)
  expect_true(!is.null(result$details))
  expect_true(grepl("nearPD", result$details))
  eig <- eigen(result$mat, symmetric = TRUE)$values
  expect_true(all(eig > 0))
})

test_that("G. .ensure_cholesky_pd enforces symmetry", {
  D <- 3
  mat <- diag(D)
  mat[1, 2] <- 0.5
  mat[2, 1] <- 0.50001  # slightly asymmetric
  result <- .ensure_cholesky_pd(mat, "asym")
  expect_equal(result$mat, t(result$mat))  # exact symmetry
})

test_that("G. cholesky_correct handles non-PD Sigma_MCMC gracefully", {
  P <- 2
  D <- 2 * P + 1

  # Create Sigma_MCMC with a negative eigenvalue
  set.seed(88)
  U <- qr.Q(qr(matrix(rnorm(D^2), D, D)))
  evals <- c(rep(1, D - 1), -0.01)  # last eigenvalue is negative
  Sigma_nonPD <- U %*% diag(evals) %*% t(U)

  fit <- create_mock_fit(P = P, M = 1000, Sigma_true = diag(D) * 0.01, seed = 88)
  sw  <- create_mock_sandwich(P = P, Sigma_MCMC = Sigma_nonPD)

  # Should apply PD correction and succeed
  res <- suppressMessages(cholesky_correct(fit, sw))
  expect_s3_class(res, "hbb_cholesky")
  expect_true(res$pd_corrections$Sigma_MCMC$corrected)
})

test_that("G. cholesky_correct handles non-PD V_sand gracefully", {
  P <- 3
  D <- 2L * P + 1L

  # Create V_sand with a clear negative eigenvalue
  set.seed(888)
  U <- qr.Q(qr(matrix(rnorm(D^2), D, D)))
  evals <- c(rep(0.01, D - 1), -0.001)
  V_nonPD <- U %*% diag(evals) %*% t(U)

  fit <- create_mock_fit(P = P, M = 500, seed = 888)
  sw  <- create_mock_sandwich(P = P, V_sand = V_nonPD)

  res <- suppressMessages(cholesky_correct(fit, sw))
  expect_s3_class(res, "hbb_cholesky")
  expect_true(res$pd_corrections$V_sand$corrected)
})


# ==========================================================================
# Section H: Edge cases
# ==========================================================================

test_that("H. minimal P=1 (D=3) works correctly", {
  P <- 1
  D <- 3
  M <- 2000

  Sigma_MCMC <- matrix(c(0.5, 0.1, 0.05,
                          0.1, 0.3, 0.02,
                          0.05, 0.02, 0.2), 3, 3)
  V_sand <- matrix(c(0.01, 0.002, 0.001,
                      0.002, 0.008, 0.0005,
                      0.001, 0.0005, 0.015), 3, 3)

  fit <- create_mock_fit(P = P, M = M, Sigma_true = Sigma_MCMC, seed = 11)
  sw  <- create_mock_sandwich(P = P, Sigma_MCMC = Sigma_MCMC, V_sand = V_sand)
  res <- quietly(cholesky_correct(fit, sw))

  expect_s3_class(res, "hbb_cholesky")
  expect_equal(res$D, 3)
  expect_equal(res$P, 1)
  expect_true(res$verification$mean_pass)
  expect_true(res$verification$A_pass)
})

test_that("H. identity transform when V_sand = Sigma_MCMC", {
  P <- 3
  D <- 2L * P + 1L
  Sigma <- diag(D) * 0.5 + 0.1

  fit <- create_mock_fit(P = P, M = 2000, seed = 700)
  sw  <- create_mock_sandwich(P = P, Sigma_MCMC = Sigma, V_sand = Sigma)
  res <- quietly(cholesky_correct(fit, sw))

  # A should be identity
  expect_equal(res$A, diag(D), tolerance = 1e-10)

  # Corrected draws should equal original draws
  param_names <- c(paste0("alpha[", 1:P, "]"), paste0("beta[", 1:P, "]"), "log_kappa")
  theta_original <- fit$fit$draws(variables = param_names, format = "matrix")
  expect_equal(unname(res$theta_corrected), unname(theta_original), tolerance = 1e-10)
})

test_that("H. level argument is stored and affects intervals", {
  fit <- create_mock_fit(P = 3, M = 1000)
  sw  <- create_mock_sandwich(P = 3)

  r90 <- quietly(cholesky_correct(fit, sw, level = 0.90))
  r95 <- quietly(cholesky_correct(fit, sw, level = 0.95))
  r99 <- quietly(cholesky_correct(fit, sw, level = 0.99))

  expect_equal(r90$level, 0.90)
  expect_equal(r95$level, 0.95)
  expect_equal(r99$level, 0.99)

  # Wald widths should increase with level
  expect_true(all(r90$comparison_table$wald_width < r95$comparison_table$wald_width))
  expect_true(all(r95$comparison_table$wald_width < r99$comparison_table$wald_width))

  # Corrected draws should be identical across levels (only intervals differ)
  expect_equal(r90$theta_corrected, r95$theta_corrected)
  expect_equal(r95$theta_corrected, r99$theta_corrected)
})

test_that("H. corrected CIs converge to Wald CIs as M increases", {
  set.seed(444)
  P <- 2
  D <- 2 * P + 1
  M <- 10000

  Sigma_MCMC <- diag(D) * 0.5
  V_sand <- diag(D) * 0.01

  fit <- create_mock_fit(P = P, M = M, Sigma_true = Sigma_MCMC, seed = 444)
  sw  <- create_mock_sandwich(P = P, Sigma_MCMC = Sigma_MCMC, V_sand = V_sand)
  res <- quietly(cholesky_correct(fit, sw, level = 0.95))

  tbl <- res$comparison_table
  width_agreement <- abs(tbl$corrected_width - tbl$wald_width) / tbl$wald_width
  expect_true(max(width_agreement) < 0.10,
    info = paste("Max width disagreement:", max(width_agreement)))
})

test_that("H. deterministic given identical inputs", {
  fit <- create_mock_fit(P = 2, M = 1000, seed = 42)
  sw  <- create_mock_sandwich(P = 2)

  r1 <- quietly(cholesky_correct(fit, sw))
  r2 <- quietly(cholesky_correct(fit, sw))

  expect_identical(r1$theta_corrected, r2$theta_corrected)
  expect_identical(r1$A, r2$A)
  expect_identical(r1$theta_hat, r2$theta_hat)
})


# ==========================================================================
# Section I: Print method
# ==========================================================================

test_that("I. print runs without error and produces output", {
  fit <- create_mock_fit(P = 2, M = 500)
  sw  <- create_mock_sandwich(P = 2)
  res <- quietly(cholesky_correct(fit, sw))

  out <- capture.output(print(res))
  expect_true(length(out) > 0)
})

test_that("I. print returns object invisibly", {
  fit <- create_mock_fit(P = 2, M = 500)
  sw  <- create_mock_sandwich(P = 2)
  res <- quietly(cholesky_correct(fit, sw))

  out <- capture.output(ret <- print(res))
  expect_identical(ret, res)
})

test_that("I. print output contains key section headers", {
  fit <- create_mock_fit(P = 2, M = 500)
  sw  <- create_mock_sandwich(P = 2)
  res <- quietly(cholesky_correct(fit, sw))

  out <- capture.output(print(res))
  full <- paste(out, collapse = "\n")
  expect_true(grepl("Williams-Savitsky", full, ignore.case = TRUE))
  expect_true(grepl("Transformation Matrix A", full))
  expect_true(grepl("Verification", full))
  expect_true(grepl("DER", full))
})

test_that("I. print shows SHRINK for A < 1 and inflate for A > 1", {
  P <- 2
  D <- 2 * P + 1

  # Large PI => SHRINK
  Sigma_MCMC_large <- diag(D) * 10.0
  V_sand_small <- diag(D) * 0.01

  fit <- create_mock_fit(P = P, M = 500, Sigma_true = Sigma_MCMC_large, seed = 33)
  sw  <- create_mock_sandwich(
    P = P, Sigma_MCMC = Sigma_MCMC_large, V_sand = V_sand_small
  )
  res <- quietly(cholesky_correct(fit, sw))
  output <- paste(capture.output(print(res)), collapse = "\n")
  expect_true(grepl("SHRINK", output))

  # Small PI => inflate
  Sigma_MCMC_small <- diag(D) * 0.01
  V_sand_large <- diag(D) * 0.02

  fit2 <- create_mock_fit(P = P, M = 500, Sigma_true = Sigma_MCMC_small, seed = 34)
  sw2  <- create_mock_sandwich(
    P = P, Sigma_MCMC = Sigma_MCMC_small, V_sand = V_sand_large
  )
  res2 <- quietly(cholesky_correct(fit2, sw2))
  output2 <- paste(capture.output(print(res2)), collapse = "\n")
  expect_true(grepl("inflate", output2))
})

test_that("I. print respects digits argument", {
  fit <- create_mock_fit(P = 2, M = 500)
  sw  <- create_mock_sandwich(P = 2)
  res <- quietly(cholesky_correct(fit, sw))

  out2 <- capture.output(print(res, digits = 2))
  out6 <- capture.output(print(res, digits = 6))
  expect_false(identical(out2, out6))
})

test_that("I. print handles corrupted object gracefully", {
  corrupted <- structure(
    list(
      D = 5,
      P = 2,
      M = NULL,
      level = 0.95,
      comparison_table = NULL,
      verification = NULL,
      pd_corrections = NULL
    ),
    class = "hbb_cholesky"
  )
  # Should not error; tryCatch inside print handles it
  expect_no_error(capture.output(print(corrupted)))
})

test_that("I. print shows PD corrections when present", {
  P <- 2
  D <- 2L * P + 1L

  # Create clearly non-PD V_sand
  set.seed(199)
  U <- qr.Q(qr(matrix(rnorm(D^2), D, D)))
  evals <- c(rep(0.01, D - 1), -0.001)
  V_nonPD <- U %*% diag(evals) %*% t(U)

  fit <- create_mock_fit(P = P, M = 500, seed = 199)
  sw  <- create_mock_sandwich(P = P, V_sand = V_nonPD)

  res <- suppressMessages(cholesky_correct(fit, sw))
  output <- paste(capture.output(print(res)), collapse = "\n")
  expect_true(
    grepl("PD|ridge|correction|nearPD", output, ignore.case = TRUE),
    info = "Print output should mention PD corrections"
  )
})


# ==========================================================================
# Section J: Internal helpers
# ==========================================================================

test_that("J. .make_stan_param_names generates correct names", {
  names1 <- .make_stan_param_names(1)
  expect_equal(names1, c("alpha[1]", "beta[1]", "log_kappa"))

  names3 <- .make_stan_param_names(3)
  expect_equal(names3, c("alpha[1]", "alpha[2]", "alpha[3]",
                          "beta[1]", "beta[2]", "beta[3]",
                          "log_kappa"))
  expect_equal(length(names3), 7)

  names5 <- .make_stan_param_names(5)
  expect_equal(length(names5), 11)
})

test_that("J. .verify_cholesky_transform produces correct structure", {
  D <- 3
  M <- 500
  theta_hat <- rep(0, D)
  A <- diag(D) * 0.5
  Sigma_MCMC <- diag(D)
  V_sand <- diag(D) * 0.25

  set.seed(100)
  theta_draws <- MASS::mvrnorm(M, theta_hat, Sigma_MCMC)
  centered <- sweep(theta_draws, 2, theta_hat, "-")
  theta_corrected <- sweep(centered %*% t(A), 2, theta_hat, "+")

  v <- .verify_cholesky_transform(theta_corrected, theta_hat, A, Sigma_MCMC, V_sand)

  expect_true(is.list(v))
  expect_true("mean_preservation_error" %in% names(v))
  expect_true("variance_relative_error" %in% names(v))
  expect_true("A_verification_error" %in% names(v))
  expect_true("mean_pass" %in% names(v))
  expect_true("variance_pass" %in% names(v))
  expect_true("A_pass" %in% names(v))

  # A * Sigma * A' = 0.5*I * I * 0.5*I = 0.25*I = V_sand
  expect_true(v$A_pass)
})

test_that("J. .validate_level rejects boundary and invalid values", {
  expect_error(.validate_level(0), regexp = "strictly")
  expect_error(.validate_level(1), regexp = "strictly")
  expect_error(.validate_level(-0.5), regexp = "strictly")
  expect_error(.validate_level(2), regexp = "strictly")
  expect_error(.validate_level("x"), regexp = "single numeric")
  expect_error(.validate_level(NA_real_), regexp = "single numeric")
  expect_error(.validate_level(c(0.5, 0.9)), regexp = "single numeric")

  # Valid levels should not error
  expect_no_error(.validate_level(0.50))
  expect_no_error(.validate_level(0.95))
  expect_no_error(.validate_level(0.001))
  expect_no_error(.validate_level(0.999))
})

test_that("J. .safe_chol_lower returns lower-triangular matrix", {
  mat <- diag(3) * 2 + 0.3
  L <- .safe_chol_lower(mat, "test")

  expect_true(all(L[upper.tri(L)] == 0))
  expect_equal(L %*% t(L), mat, tolerance = 1e-12)
})

test_that("J. .safe_chol_lower errors on non-PD matrix", {
  bad_mat <- diag(c(1, -1, 1))
  expect_error(.safe_chol_lower(bad_mat, "test"), "Cholesky")
})


# ==========================================================================
# Section K: Stress tests
# ==========================================================================

test_that("K. large P=10 (D=21) works correctly", {
  P <- 10
  D <- 2L * P + 1L
  M <- 3000

  set.seed(999)
  R <- matrix(rnorm(D * D), D, D)
  Sigma_MCMC <- crossprod(R) / D + diag(D) * 0.5
  R2 <- matrix(rnorm(D * D), D, D)
  V_sand <- crossprod(R2) / (D * 100) + diag(D) * 0.01
  H_obs_inv <- diag(D) * 0.005

  theta_true <- rnorm(D, 0, 0.2)
  fit <- create_mock_fit(P = P, M = M, seed = 999,
                          Sigma_true = Sigma_MCMC, theta_true = theta_true)
  sw  <- create_mock_sandwich(
    P = P, Sigma_MCMC = Sigma_MCMC, V_sand = V_sand, H_obs_inv = H_obs_inv
  )

  res <- quietly(cholesky_correct(fit, sw))

  # Algebraic identity
  A_Sigma_At <- res$A %*% Sigma_MCMC %*% t(res$A)
  expect_equal(A_Sigma_At, V_sand, tolerance = 1e-10)

  # Mean preservation
  expect_true(max(abs(colMeans(res$theta_corrected) - res$theta_hat)) < 1e-10)

  # Dimensions
  expect_equal(res$D, D)
  expect_equal(res$P, P)
  expect_equal(dim(res$theta_corrected), c(M, D))
  expect_equal(nrow(res$comparison_table), D)
})

test_that("K. very small M=10 does not crash", {
  fit <- create_mock_fit(P = 2, M = 10, seed = 300)
  sw  <- create_mock_sandwich(P = 2)

  res <- quietly(cholesky_correct(fit, sw))
  expect_s3_class(res, "hbb_cholesky")
  expect_equal(res$M, 10)
})

test_that("K. degenerate draws (near-constant) don't crash", {
  P <- 2
  D <- 2L * P + 1L
  param_names <- c(paste0("alpha[", seq_len(P), "]"),
                   paste0("beta[", seq_len(P), "]"),
                   "log_kappa")

  # Create draws where all values are nearly identical
  set.seed(501)
  M <- 50
  base <- rnorm(D, 0, 0.01)
  draws <- matrix(rep(base, M), nrow = M, ncol = D, byrow = TRUE)
  draws <- draws + matrix(rnorm(M * D, 0, 1e-10), M, D)
  colnames(draws) <- param_names

  mock_cmdstan <- list(
    draws = function(variables = NULL, format = "matrix") {
      if (!is.null(variables)) draws[, variables, drop = FALSE]
      else draws
    }
  )
  fit <- structure(
    list(fit = mock_cmdstan, model_type = "weighted", hbb_data = list(P = P)),
    class = "hbb_fit"
  )

  sw <- create_mock_sandwich(P = P)

  # Draws have near-zero variance but sandwich Sigma_MCMC is standard PD.
  # The function should complete without crashing, even though the
  # finite-sample variance check will show poor recovery.
  res <- suppressMessages(cholesky_correct(fit, sw))
  expect_s3_class(res, "hbb_cholesky")
  expect_equal(res$M, M)
  # Algebraic identity (A * Sigma * A' = V_sand) should still pass
  # because it depends only on Sigma_MCMC from the sandwich, not the draws
  expect_true(res$verification$A_pass)
})

test_that("K. extreme Wald CI levels work", {
  theta_hat <- c(0.5, -0.3)
  V_sand <- diag(c(0.01, 0.04))

  # Very narrow level
  w_narrow <- compute_wald_ci(theta_hat, V_sand, level = 0.001)
  expect_equal(nrow(w_narrow), 2L)
  expect_true(all(w_narrow$ci_width > 0))
  expect_true(all(w_narrow$ci_width < 0.01))

  # Very wide level
  w_wide <- compute_wald_ci(theta_hat, V_sand, level = 0.999)
  expect_equal(nrow(w_wide), 2L)
  expect_true(all(w_wide$ci_width > 0.5))
})


# ==========================================================================
# Section L: Comparison table and DER diagnostics
# ==========================================================================

test_that("L. comparison_table has all 16 columns", {
  fit <- create_mock_fit(P = 2, M = 500)
  sw  <- create_mock_sandwich(P = 2)
  res <- quietly(cholesky_correct(fit, sw))

  tbl <- res$comparison_table
  expect_s3_class(tbl, "data.frame")

  required_cols <- c("parameter", "post_mean", "naive_lo", "naive_hi",
                     "naive_width", "corrected_lo", "corrected_hi",
                     "corrected_width", "wald_lo", "wald_hi",
                     "wald_width", "width_ratio", "DER",
                     "DER_vs_MCMC", "prior_inflation", "sqrt_DER")
  for (col in required_cols) {
    expect_true(col %in% names(tbl), info = paste("Missing column:", col))
  }

  D <- 2 * 2 + 1
  expect_equal(nrow(tbl), D)
})

test_that("L. DER and PI computed correctly from known diagonals", {
  P <- 2
  D <- 2 * P + 1

  V_sand_diag <- c(0.01, 0.02, 0.03, 0.04, 0.05)
  Sigma_MCMC_diag <- c(0.5, 1.0, 1.5, 2.0, 2.5)
  H_obs_inv_diag <- c(0.005, 0.01, 0.015, 0.02, 0.025)

  V_sand <- diag(V_sand_diag)
  Sigma_MCMC <- diag(Sigma_MCMC_diag)
  H_obs_inv <- diag(H_obs_inv_diag)

  fit <- create_mock_fit(P = P, M = 1000, Sigma_true = Sigma_MCMC, seed = 55)
  sw  <- create_mock_sandwich(
    P = P, Sigma_MCMC = Sigma_MCMC, V_sand = V_sand, H_obs_inv = H_obs_inv
  )
  res <- quietly(cholesky_correct(fit, sw))
  tbl <- res$comparison_table

  # DER = V_sand / H_obs_inv
  expected_DER <- V_sand_diag / H_obs_inv_diag
  expect_equal(tbl$DER, expected_DER, tolerance = 1e-12)

  # PI = Sigma_MCMC / H_obs_inv
  expected_PI <- Sigma_MCMC_diag / H_obs_inv_diag
  expect_equal(tbl$prior_inflation, expected_PI, tolerance = 1e-12)

  # DER_vs_MCMC = V_sand / Sigma_MCMC
  expected_DER_MCMC <- V_sand_diag / Sigma_MCMC_diag
  expect_equal(tbl$DER_vs_MCMC, expected_DER_MCMC, tolerance = 1e-12)

  # sqrt_DER
  expect_equal(tbl$sqrt_DER, sqrt(expected_DER), tolerance = 1e-12)
})

test_that("L. width_ratio < 1 when V_sand << Sigma_MCMC", {
  P <- 3
  D <- 2L * P + 1L
  Sigma_MCMC <- diag(D) * 10
  V_sand     <- diag(D) * 0.01

  fit <- create_mock_fit(P = P, M = 4000, seed = 505)
  sw  <- create_mock_sandwich(P = P, Sigma_MCMC = Sigma_MCMC, V_sand = V_sand)
  res <- quietly(cholesky_correct(fit, sw))

  expect_true(all(res$comparison_table$width_ratio < 1))
})

test_that("L. param_labels fallback when missing from sandwich", {
  sw <- create_mock_sandwich(P = 2)
  sw$param_labels <- NULL

  fit <- create_mock_fit(P = 2)

  # param_labels is no longer a required field; the function falls back
  # to CmdStanR names gracefully
  res <- suppressMessages(cholesky_correct(fit, sw))
  expect_equal(length(res$comparison_table$parameter), res$D)
})

test_that("L. param_labels fallback when wrong length", {
  sw <- create_mock_sandwich(P = 2)
  sw$param_labels <- c("a", "b")  # wrong length (should be 5)

  fit <- create_mock_fit(P = 2)

  res <- suppressMessages(cholesky_correct(fit, sw))
  expect_equal(length(res$comparison_table$parameter), res$D)
})
