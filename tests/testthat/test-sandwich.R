# ============================================================================
# test-sandwich.R --- Tests for the Sandwich Variance Estimator Module
#
# Organisation:
#   A. compute_J_cluster --- Hand-computed known-answer tests
#   B. compute_J_cluster --- Edge cases
#   C. compute_der       --- Ratio computation and extreme warnings
#   D. H_obs block structure --- .block_diag verification
#   E. .build_param_labels --- Label construction
#   F. Input validation  --- Error messages
#   G. .check_pd / .ridge_regularize --- Robustness helpers
#   H. print.hbb_sandwich --- Print method robustness
#   I. DER ~ 1 under uniform weights --- Validation principle
#   J. Integration test  --- Full pipeline (requires CmdStan)
# ============================================================================


# ============================================================================
# A. compute_J_cluster --- Hand-computed known-answer tests
# ============================================================================

test_that("compute_J_cluster: hand-computed 2-stratum example", {
    # Setup: N=6, D=2, 2 strata, 2 PSUs each
    scores <- matrix(c(
        1,  2,   # obs 1, stratum 1 PSU 1
        3,  0,   # obs 2, stratum 1 PSU 1
        2,  1,   # obs 3, stratum 1 PSU 2
        0,  4,   # obs 4, stratum 2 PSU 3
        1, -1,   # obs 5, stratum 2 PSU 4
        1,  1    # obs 6, stratum 2 PSU 4
    ), nrow = 6, ncol = 2, byrow = TRUE)

    hbb_data <- list(
        stratum_idx = c(1L, 1L, 1L, 2L, 2L, 2L),
        psu_idx     = c(1L, 1L, 2L, 3L, 4L, 4L),
        w_tilde     = rep(1, 6),
        N           = 6L
    )

    J <- compute_J_cluster(scores, hbb_data)

    # Manual computation:
    # Stratum 1 (C_h = 2):
    #   PSU 1: s = (4, 2), PSU 2: s = (2, 1), s_bar = (3, 1.5)
    #   deltas: (1, 0.5), (-1, -0.5)
    #   J_1 = 2 * [outer((1,0.5)) + outer((-1,-0.5))]
    #       = 4 * outer((1,0.5)) = matrix(c(4, 2, 2, 1), 2, 2)
    #
    # Stratum 2 (C_h = 2):
    #   PSU 3: s = (0, 4), PSU 4: s = (2, 0), s_bar = (1, 2)
    #   deltas: (-1, 2), (1, -2)
    #   J_2 = 4 * outer((1,-2)) = matrix(c(4, -8, -8, 16), 2, 2)
    #
    # J = J_1 + J_2 = matrix(c(8, -6, -6, 17), 2, 2)
    J_expected <- matrix(c(8, -6, -6, 17), nrow = 2, ncol = 2)
    expect_equal(J, J_expected, tolerance = 1e-12)
})


test_that("compute_J_cluster: with non-uniform weights", {
    scores <- matrix(c(
        1,  0,
        0,  1,
        1,  1,
        2,  0,
        0,  2,
        1,  1
    ), nrow = 6, ncol = 2, byrow = TRUE)

    hbb_data <- list(
        stratum_idx = c(1L, 1L, 1L, 2L, 2L, 2L),
        psu_idx     = c(1L, 1L, 2L, 3L, 4L, 4L),
        w_tilde     = c(2, 1, 3, 1, 2, 2),
        N           = 6L
    )

    J <- compute_J_cluster(scores, hbb_data)

    # Manual:
    # Stratum 1: PSU1 = (2,1), PSU2 = (3,3), s_bar = (2.5,2)
    #   J_1 = 4 * outer((0.5,1)) = matrix(c(1, 2, 2, 4), 2, 2)
    # Stratum 2: PSU3 = (2,0), PSU4 = (2,6), s_bar = (2,3)
    #   J_2 = 4 * outer((0,3)) = matrix(c(0, 0, 0, 36), 2, 2)
    # J = matrix(c(1, 2, 2, 40), 2, 2)
    J_expected <- matrix(c(1, 2, 2, 40), nrow = 2, ncol = 2)
    expect_equal(J, J_expected, tolerance = 1e-12)
})


test_that("compute_J_cluster: D=1 case works (scalar)", {
    scores <- matrix(c(1, 3, 2, 4), nrow = 4, ncol = 1)
    hbb_data <- list(
        stratum_idx = c(1L, 1L, 1L, 1L),
        psu_idx     = c(1L, 1L, 2L, 2L),
        w_tilde     = rep(1, 4),
        N           = 4L
    )

    J <- compute_J_cluster(scores, hbb_data)

    # PSU 1 total: 4, PSU 2 total: 6, s_bar = 5
    # FPC = 2, J = 2 * (1 + 1) = 4
    expect_equal(J, matrix(4, 1, 1), tolerance = 1e-12)
})


test_that("compute_J_cluster: 3 PSUs per stratum", {
    scores <- matrix(c(1, 4, 7), nrow = 3, ncol = 1)
    hbb_data <- list(
        stratum_idx = c(1L, 1L, 1L),
        psu_idx     = c(1L, 2L, 3L),
        w_tilde     = rep(1, 3),
        N           = 3L
    )

    J <- compute_J_cluster(scores, hbb_data)

    # PSU totals: 1, 4, 7. Mean = 4. FPC = 3/2
    # J = 1.5 * (9 + 0 + 9) = 27
    expect_equal(J, matrix(27, 1, 1), tolerance = 1e-12)
})


test_that("compute_J_cluster: unequal PSUs per stratum", {
    scores <- matrix(c(
        1, 0,
        0, 1,
        1, 1,
        2, 0,
        0, 2
    ), nrow = 5, ncol = 2, byrow = TRUE)

    hbb_data <- list(
        stratum_idx = c(1L, 1L, 1L, 2L, 2L),
        psu_idx     = c(1L, 2L, 3L, 4L, 5L),
        w_tilde     = rep(1, 5),
        N           = 5L
    )

    J <- compute_J_cluster(scores, hbb_data)

    # Stratum 1 (C=3, FPC=3/2):
    #   PSU totals: (1,0), (0,1), (1,1). Mean: (2/3, 2/3)
    #   J_1 = matrix(c(1, -0.5, -0.5, 1), 2, 2)
    J_1 <- matrix(c(1, -0.5, -0.5, 1), 2, 2)

    # Stratum 2 (C=2, FPC=2):
    #   PSU totals: (2,0), (0,2). Mean: (1, 1)
    #   J_2 = matrix(c(4, -4, -4, 4), 2, 2)
    J_2 <- matrix(c(4, -4, -4, 4), 2, 2)

    J_expected <- J_1 + J_2
    expect_equal(J, J_expected, tolerance = 1e-12)
})


# ============================================================================
# B. compute_J_cluster --- Edge cases
# ============================================================================

test_that("compute_J_cluster: singleton strata are skipped with warning", {
    scores <- matrix(c(
        1, 2,
        3, 4,
        5, 6
    ), nrow = 3, ncol = 2, byrow = TRUE)

    hbb_data <- list(
        stratum_idx = c(1L, 1L, 2L),
        psu_idx     = c(1L, 2L, 3L),
        w_tilde     = rep(1, 3),
        N           = 3L
    )

    # Stratum 2 has only 1 PSU -> singleton, should warn
    expect_warning(
        J <- compute_J_cluster(scores, hbb_data),
        "singleton"
    )

    # Only stratum 1 contributes:
    #   PSU1=(1,2), PSU2=(3,4), s_bar=(2,3)
    #   J = 2 * [outer((-1,-1)) + outer((1,1))] = 4 * outer((1,1))
    J_expected <- matrix(c(4, 4, 4, 4), nrow = 2, ncol = 2)
    expect_equal(J, J_expected, tolerance = 1e-12)
})


test_that("compute_J_cluster: all-singleton strata gives zero matrix", {
    scores <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, byrow = TRUE)
    hbb_data <- list(
        stratum_idx = c(1L, 2L),
        psu_idx     = c(1L, 2L),
        w_tilde     = rep(1, 2),
        N           = 2L
    )

    suppressWarnings(
        J <- compute_J_cluster(scores, hbb_data)
    )
    expect_equal(J, matrix(0, 2, 2))
})


test_that("compute_J_cluster: zero scores produce zero J", {
    scores <- matrix(0, 4, 2)
    hbb_data <- list(
        stratum_idx = c(1L, 1L, 1L, 1L),
        psu_idx     = c(1L, 1L, 2L, 2L),
        w_tilde     = rep(1, 4),
        N           = 4L
    )

    J <- compute_J_cluster(scores, hbb_data)
    expect_equal(J, matrix(0, 2, 2), tolerance = 1e-15,
                 ignore_attr = TRUE)
})


test_that("compute_J_cluster: uniform scores with equal PSU sizes give zero J", {
    N <- 8L
    D <- 3L
    scores <- matrix(rep(c(1, 2, 3), each = N), nrow = N, ncol = D)

    hbb_data <- list(
        stratum_idx = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
        psu_idx     = c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L),
        w_tilde     = rep(1, N),
        N           = N
    )

    J <- compute_J_cluster(scores, hbb_data)
    expect_equal(max(abs(J)), 0, tolerance = 1e-12)
})


test_that("compute_J_cluster: weights scale PSU totals correctly", {
    N <- 4L
    scores <- matrix(c(1, 1, 1, 1), N, 1)

    base_data <- list(
        stratum_idx = c(1L, 1L, 1L, 1L),
        psu_idx     = c(1L, 1L, 2L, 2L),
        N           = N
    )

    # Uniform weights: PSU totals are equal -> J = 0
    data_unif <- base_data
    data_unif$w_tilde <- rep(1, N)
    J_unif <- compute_J_cluster(scores, data_unif)
    expect_equal(as.numeric(J_unif), 0, tolerance = 1e-12)

    # Unequal weights: PSU 1 gets 2, PSU 2 gets 0.5
    # PSU 1 total: 4, PSU 2 total: 1, s_bar = 2.5
    # J = 2 * (1.5^2 + 1.5^2) = 9
    data_wt <- base_data
    data_wt$w_tilde <- c(2, 2, 0.5, 0.5)
    J_wt <- compute_J_cluster(scores, data_wt)
    expect_equal(as.numeric(J_wt), 9, tolerance = 1e-12)
})


test_that("compute_J_cluster: weight doubling scales J by 4", {
    scores <- matrix(c(1, 0, 0, 1, 1, 1, -1, 0), nrow = 4, ncol = 2,
                     byrow = TRUE)

    base <- list(
        stratum_idx = c(1L, 1L, 1L, 1L),
        psu_idx     = c(1L, 1L, 2L, 2L),
        N           = 4L
    )

    base1 <- base; base1$w_tilde <- rep(1, 4)
    base2 <- base; base2$w_tilde <- rep(2, 4)

    J_1 <- compute_J_cluster(scores, base1)
    J_2 <- compute_J_cluster(scores, base2)

    # Doubling weights: PSU totals scale 2x, outer products 4x
    expect_equal(J_2, 4 * J_1, tolerance = 1e-12)
})


test_that("compute_J_cluster: symmetry of output", {
    set.seed(42)
    N <- 50
    D <- 5
    scores <- matrix(rnorm(N * D), nrow = N, ncol = D)

    hbb_data <- list(
        stratum_idx = rep(1:5, each = 10),
        psu_idx     = rep(1:10, each = 5),
        w_tilde     = runif(N, 0.5, 2),
        N           = N
    )

    J <- compute_J_cluster(scores, hbb_data)
    expect_equal(J, t(J), tolerance = 1e-14)
})


test_that("compute_J_cluster: PSD property", {
    set.seed(123)
    N <- 100
    D <- 4
    scores <- matrix(rnorm(N * D), nrow = N, ncol = D)

    hbb_data <- list(
        stratum_idx = rep(1:10, each = 10),
        psu_idx     = rep(1:20, each = 5),
        w_tilde     = runif(N, 0.5, 3),
        N           = N
    )

    J <- compute_J_cluster(scores, hbb_data)
    eigenvalues <- eigen(J, symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(eigenvalues >= -1e-10))
})


test_that("compute_J_cluster: large number of strata", {
    set.seed(99)
    H <- 50
    N <- H * 4  # 2 PSUs per stratum, 2 obs per PSU
    D <- 3
    scores <- matrix(rnorm(N * D), nrow = N, ncol = D)

    stratum_idx <- rep(seq_len(H), each = 4)
    psu_idx     <- rep(rep(1:2, each = 2), H)
    psu_idx     <- psu_idx + (stratum_idx - 1L) * 2L

    hbb_data <- list(
        stratum_idx = as.integer(stratum_idx),
        psu_idx     = as.integer(psu_idx),
        w_tilde     = runif(N, 0.5, 2),
        N           = as.integer(N)
    )

    J <- compute_J_cluster(scores, hbb_data)

    expect_equal(dim(J), c(D, D))
    expect_equal(J, t(J), tolerance = 1e-13)

    eig <- eigen(J, symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(eig >= -1e-10))
})


test_that("compute_J_cluster: default to uniform weights when NULL", {
    scores <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2, byrow = TRUE)

    hbb_data_w <- list(
        stratum_idx = c(1L, 1L),
        psu_idx     = c(1L, 2L),
        w_tilde     = c(1, 1),
        N           = 2L
    )

    hbb_data_null <- list(
        stratum_idx = c(1L, 1L),
        psu_idx     = c(1L, 2L),
        w_tilde     = NULL,
        N           = 2L
    )

    J_with_w <- compute_J_cluster(scores, hbb_data_w)
    J_null_w <- compute_J_cluster(scores, hbb_data_null)

    expect_equal(J_with_w, J_null_w, tolerance = 1e-14)
})


# ============================================================================
# C. compute_der --- Ratio computation and extreme warnings
# ============================================================================

test_that("compute_der: correct DER from known diagonal matrices", {
    mock_sandwich <- structure(
        list(
            V_sand       = diag(c(4, 9, 16)),
            H_obs_inv    = diag(c(1, 3, 4)),
            param_labels = c("alpha_1", "beta_1", "log_kappa")
        ),
        class = "hbb_sandwich"
    )

    der <- compute_der(mock_sandwich)
    expect_named(der, c("alpha_1", "beta_1", "log_kappa"))
    expect_equal(as.numeric(der), c(4, 3, 4), tolerance = 1e-12)
})


test_that("compute_der: DER = 1 when V_sand equals H_obs_inv", {
    D <- 5L
    V_sand <- diag(rep(2.5, D))
    H_obs_inv <- diag(rep(2.5, D))
    param_labels <- paste0("p", seq_len(D))

    mock_sandwich <- structure(
        list(V_sand = V_sand, H_obs_inv = H_obs_inv,
             param_labels = param_labels),
        class = "hbb_sandwich"
    )

    der <- compute_der(mock_sandwich)
    expect_equal(as.numeric(der), rep(1, D), tolerance = 1e-12)
})


test_that("compute_der: handles off-diagonal V_sand correctly", {
    V_sand <- matrix(c(4, 1, 1, 9), 2, 2)
    H_obs_inv <- diag(c(2, 3))

    mock_sandwich <- structure(
        list(V_sand = V_sand, H_obs_inv = H_obs_inv,
             param_labels = c("a1", "lk")),
        class = "hbb_sandwich"
    )

    der <- compute_der(mock_sandwich)
    expect_equal(as.numeric(der), c(4 / 2, 9 / 3), tolerance = 1e-12)
})


test_that("compute_der: errors on non-hbb_sandwich input", {
    expect_error(compute_der(list(V_sand = diag(3))), "hbb_sandwich")
    expect_error(compute_der("not a sandwich"), "hbb_sandwich")
})


test_that("compute_der: DER is exactly 1 when J_cluster equals H_obs", {
    D <- 5L
    set.seed(99)
    M <- matrix(rnorm(D * D), D, D)
    H_obs <- crossprod(M) + D * diag(D)
    H_obs_inv <- solve(H_obs)
    V_sand <- H_obs_inv  # because J = H means V = H^{-1}

    mock_sandwich <- structure(
        list(V_sand = V_sand, H_obs_inv = H_obs_inv,
             param_labels = paste0("p", seq_len(D))),
        class = "hbb_sandwich"
    )

    der <- compute_der(mock_sandwich)
    expect_equal(as.numeric(der), rep(1, D), tolerance = 1e-10)
})


test_that("compute_der: DER > 1 when J_cluster > H_obs (inflation)", {
    D <- 3L
    set.seed(77)
    M <- matrix(rnorm(D * D), D, D)
    H_obs <- crossprod(M) + D * diag(D)
    H_obs_inv <- solve(H_obs)
    V_sand <- H_obs_inv %*% (4 * H_obs) %*% H_obs_inv
    # V_sand = 4 * H_inv, so DER = 4

    mock_sandwich <- structure(
        list(V_sand = V_sand, H_obs_inv = H_obs_inv,
             param_labels = paste0("p", seq_len(D))),
        class = "hbb_sandwich"
    )

    der <- compute_der(mock_sandwich)
    expect_equal(as.numeric(der), rep(4, D), tolerance = 1e-10)
})


# ============================================================================
# D. H_obs block structure --- .block_diag
# ============================================================================

test_that(".block_diag: dimensions correct", {
    A <- matrix(1:4, 2, 2)
    B <- matrix(1:9, 3, 3)
    result <- hurdlebb:::.block_diag(A, B)
    expect_equal(dim(result), c(5, 5))
})


test_that(".block_diag: 1x1 blocks", {
    A <- matrix(5)
    B <- matrix(3)
    result <- hurdlebb:::.block_diag(A, B)
    expect_equal(result, matrix(c(5, 0, 0, 3), 2, 2))
})


test_that(".block_diag: preserves values and zero off-diagonal", {
    A <- matrix(c(1, 2, 3, 4), 2, 2)
    B <- matrix(c(5, 6, 7, 8, 9, 10, 11, 12, 13), 3, 3)
    result <- hurdlebb:::.block_diag(A, B)

    expect_equal(result[1:2, 1:2], A)
    expect_equal(result[3:5, 3:5], B)
    expect_equal(result[1:2, 3:5], matrix(0, 2, 3))
    expect_equal(result[3:5, 1:2], matrix(0, 3, 2))
})


test_that(".block_diag: preserves positive definiteness", {
    set.seed(42)
    A <- crossprod(matrix(rnorm(9), 3, 3)) + 3 * diag(3)
    B <- crossprod(matrix(rnorm(4), 2, 2)) + 2 * diag(2)

    result <- hurdlebb:::.block_diag(A, B)
    eig <- eigen(result, symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(eig > 0))
})


test_that(".block_diag: eigenvalues are union of block eigenvalues", {
    set.seed(123)
    A <- crossprod(matrix(rnorm(4), 2, 2)) + 2 * diag(2)
    B <- crossprod(matrix(rnorm(9), 3, 3)) + 3 * diag(3)

    result <- hurdlebb:::.block_diag(A, B)

    eig_A <- sort(eigen(A, symmetric = TRUE, only.values = TRUE)$values)
    eig_B <- sort(eigen(B, symmetric = TRUE, only.values = TRUE)$values)
    eig_R <- sort(eigen(result, symmetric = TRUE, only.values = TRUE)$values)

    expect_equal(eig_R, sort(c(eig_A, eig_B)), tolerance = 1e-10)
})


test_that("H_obs off-diagonal blocks are zero (block-diagonal structure)", {
    P <- 3
    D <- 2 * P + 1

    H_ext <- crossprod(matrix(runif(P * P), P, P))
    H_int <- crossprod(matrix(runif((P + 1)^2), P + 1, P + 1))

    H_obs <- hurdlebb:::.block_diag(H_ext, H_int)

    expect_equal(H_obs[1:P, (P + 1):D], matrix(0, P, P + 1))
    expect_equal(H_obs[(P + 1):D, 1:P], matrix(0, P + 1, P))
    expect_equal(H_obs[1:P, 1:P], H_ext)
    expect_equal(H_obs[(P + 1):D, (P + 1):D], H_int)
})


# ============================================================================
# E. .build_param_labels --- Label construction
# ============================================================================

test_that(".build_param_labels: standard case with column names", {
    hd <- list(
        P = 3L,
        X = matrix(0, nrow = 5, ncol = 3,
                    dimnames = list(NULL, c("(Intercept)", "poverty", "urban")))
    )

    labels <- hurdlebb:::.build_param_labels(hd)
    expect_equal(labels, c(
        "alpha_(Intercept)", "alpha_poverty", "alpha_urban",
        "beta_(Intercept)", "beta_poverty", "beta_urban",
        "log_kappa"
    ))
    expect_length(labels, 2 * 3 + 1)
})


test_that(".build_param_labels: no column names fallback", {
    hd <- list(P = 2L, X = matrix(0, nrow = 3, ncol = 2))
    labels <- hurdlebb:::.build_param_labels(hd)
    expect_equal(labels, c("alpha_x1", "alpha_x2",
                            "beta_x1", "beta_x2",
                            "log_kappa"))
})


test_that(".build_param_labels: D = 2P + 1 for various P", {
    for (P in 1:5) {
        hd <- list(
            P = P,
            X = matrix(0, nrow = 2, ncol = P,
                        dimnames = list(NULL, paste0("v", seq_len(P))))
        )
        labels <- hurdlebb:::.build_param_labels(hd)
        expect_length(labels, 2 * P + 1)
        expect_equal(labels[length(labels)], "log_kappa")
    }
})


# ============================================================================
# F. Input validation --- Error messages
# ============================================================================

test_that("sandwich_variance: errors on non-hbb_fit input", {
    expect_error(sandwich_variance(list(model_type = "weighted")), "hbb_fit")
    expect_error(sandwich_variance(42), "hbb_fit")
})


test_that("sandwich_variance: errors on unweighted (base) model", {
    mock_fit <- structure(
        list(model_type = "base", hbb_data = list(N = 10, P = 2)),
        class = "hbb_fit"
    )
    expect_error(sandwich_variance(mock_fit), "weighted")
})


test_that("sandwich_variance: errors on SVC (unweighted) model", {
    mock_fit <- structure(
        list(model_type = "svc", hbb_data = list(N = 10, P = 2)),
        class = "hbb_fit"
    )
    expect_error(sandwich_variance(mock_fit), "weighted")
})


test_that("sandwich_variance: errors when fit$fit is NULL", {
    mock_fit <- structure(
        list(
            model_type = "weighted",
            fit        = NULL,
            hbb_data   = structure(
                list(stratum_idx = 1:10, psu_idx = 1:10),
                class = "hbb_data"
            )
        ),
        class = "hbb_fit"
    )
    expect_error(sandwich_variance(mock_fit), "NULL")
})


test_that("sandwich_variance: errors when survey design missing", {
    mock_fit <- structure(
        list(
            model_type = "weighted",
            fit        = "placeholder",
            hbb_data   = list(stratum_idx = NULL, psu_idx = NULL)
        ),
        class = "hbb_fit"
    )
    expect_error(sandwich_variance(mock_fit), "stratum_idx")
})


test_that("compute_score_matrix: errors on non-hbb_fit", {
    expect_error(compute_score_matrix(list()), "hbb_fit")
})


test_that("compute_score_matrix: errors on non-weighted model", {
    fake_fit <- structure(list(model_type = "base"), class = "hbb_fit")
    expect_error(compute_score_matrix(fake_fit), "weighted")
})


test_that("compute_H_obs: errors on non-hbb_fit", {
    expect_error(compute_H_obs(list(a = 1)), "hbb_fit")
})


test_that("compute_J_cluster: errors on dimension mismatch", {
    scores <- matrix(1, 4, 2)
    bad_data <- list(
        stratum_idx = c(1L, 1L, 1L),
        psu_idx     = c(1L, 2L, 2L),
        w_tilde     = rep(1, 3),
        N           = 3L
    )
    expect_error(compute_J_cluster(scores, bad_data), "mismatch|Dimension")
})


test_that("compute_J_cluster: errors on missing stratum_idx", {
    scores <- matrix(1:4, nrow = 2, ncol = 2)
    hbb_data <- list(
        stratum_idx = NULL,
        psu_idx     = c(1L, 2L),
        w_tilde     = c(1, 1),
        N           = 2L
    )
    expect_error(compute_J_cluster(scores, hbb_data), "stratum_idx")
})


test_that("compute_J_cluster: errors on missing psu_idx", {
    scores <- matrix(1:4, nrow = 2, ncol = 2)
    hbb_data <- list(
        stratum_idx = c(1L, 1L),
        psu_idx     = NULL,
        w_tilde     = c(1, 1),
        N           = 2L
    )
    expect_error(compute_J_cluster(scores, hbb_data), "psu_idx")
})


test_that("compute_J_cluster: errors on w_tilde length mismatch", {
    scores <- matrix(1:6, nrow = 3, ncol = 2)
    hbb_data <- list(
        stratum_idx = c(1L, 1L, 1L),
        psu_idx     = c(1L, 1L, 2L),
        w_tilde     = c(1, 1),  # length 2 != 3
        N           = 3L
    )
    expect_error(compute_J_cluster(scores, hbb_data), "w_tilde")
})


# ============================================================================
# G. .check_pd / .ridge_regularize --- Robustness helpers
# ============================================================================

test_that(".check_pd correctly identifies PD matrix", {
    mat <- diag(5)
    result <- hurdlebb:::.check_pd(mat, label = "test")
    expect_true(result$is_pd)
    expect_equal(result$min_eig, 1, tolerance = 1e-14)
    expect_equal(result$max_eig, 1, tolerance = 1e-14)
    expect_equal(result$cond, 1, tolerance = 1e-10)
})


test_that(".check_pd correctly identifies non-PD matrix", {
    mat <- matrix(c(1, 2, 2, 1), 2, 2)  # eigenvalues: 3, -1
    result <- hurdlebb:::.check_pd(mat, label = "test")
    expect_false(result$is_pd)
    expect_true(result$min_eig < 0)
})


test_that(".check_pd handles zero matrix", {
    mat <- matrix(0, 3, 3)
    result <- hurdlebb:::.check_pd(mat, label = "test")
    expect_false(result$is_pd)
    expect_equal(result$min_eig, 0, tolerance = 1e-14)
})


test_that(".ridge_regularize makes non-PD matrix PD", {
    mat <- matrix(c(1, 2, 2, 1), 2, 2)
    result <- hurdlebb:::.ridge_regularize(mat)
    eig <- eigen(result, symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(eig > 0))
})


test_that(".ridge_regularize preserves symmetry", {
    mat <- matrix(c(1, 2, 2, 1), 2, 2)
    result <- hurdlebb:::.ridge_regularize(mat)
    expect_equal(result, t(result), tolerance = 1e-14)
})


test_that(".ridge_regularize with already-PD matrix still works", {
    mat <- diag(3)
    result <- hurdlebb:::.ridge_regularize(mat)
    eig <- eigen(result, symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(eig > 0))
})


# ============================================================================
# H. print.hbb_sandwich --- Print method robustness
# ============================================================================

test_that("print.hbb_sandwich: produces expected output sections", {
    D <- 3L
    V_sand <- diag(c(0.01, 0.02, 0.03))
    H_obs_inv <- diag(c(0.005, 0.01, 0.015))
    param_labels <- c("alpha_intercept", "beta_intercept", "log_kappa")

    mock_sandwich <- structure(
        list(
            V_sand       = V_sand,
            H_obs        = diag(1 / c(0.005, 0.01, 0.015)),
            H_obs_inv    = H_obs_inv,
            J_cluster    = diag(D),
            Sigma_MCMC   = diag(D),
            DER          = setNames(c(2, 2, 2), param_labels),
            param_labels = param_labels,
            D            = D,
            N            = 100L,
            model_type   = "weighted",
            survey_info  = list(n_strata = 10L, n_psu = 50L,
                                df = 40L, n_singleton = 0L),
            matrix_diagnostics = list(
                H_obs = list(min_eig = 1, max_eig = 200, cond = 200,
                             is_pd = TRUE),
                J_cluster = list(min_eig = 0.5, max_eig = 1.5,
                                 is_pd = TRUE),
                V_sand = list(min_eig = 0.005, max_eig = 0.03,
                              cond = 6, is_pd = TRUE),
                Sigma_MCMC = list(min_eig = 0.5, max_eig = 1.5,
                                  cond = 3, is_pd = TRUE)
            ),
            nearPD_applied = FALSE,
            call = quote(sandwich_variance(fit))
        ),
        class = "hbb_sandwich"
    )

    output <- capture.output(print(mock_sandwich))
    output_text <- paste(output, collapse = "\n")

    expect_true(grepl("Sandwich Variance", output_text))
    expect_true(grepl("Design Effect", output_text))
    expect_true(grepl("weighted", output_text))
    expect_true(grepl("alpha_intercept", output_text))
    expect_true(grepl("DER", output_text))
})


test_that("print.hbb_sandwich: returns x invisibly", {
    mock_sandwich <- structure(
        list(
            V_sand       = diag(2),
            H_obs        = diag(2),
            H_obs_inv    = diag(2),
            J_cluster    = diag(2),
            Sigma_MCMC   = diag(2),
            DER          = c(a = 1, b = 1),
            param_labels = c("a", "b"),
            D            = 2L,
            N            = 50L,
            model_type   = "weighted",
            survey_info  = list(n_strata = 1L, n_psu = 2L,
                                df = 1L, n_singleton = 0L),
            matrix_diagnostics = list(
                H_obs = list(min_eig = 1, max_eig = 1, cond = 1, is_pd = TRUE),
                V_sand = list(min_eig = 1, max_eig = 1, cond = 1, is_pd = TRUE),
                J_cluster = list(min_eig = 1, max_eig = 1, is_pd = TRUE),
                Sigma_MCMC = list(min_eig = 1, max_eig = 1, cond = 1, is_pd = TRUE)
            ),
            nearPD_applied = FALSE,
            call = quote(sandwich_variance(fit))
        ),
        class = "hbb_sandwich"
    )

    result <- withVisible(print(mock_sandwich))
    expect_false(result$visible)
    expect_s3_class(result$value, "hbb_sandwich")
})


test_that("print.hbb_sandwich: does not error on minimal object", {
    sand <- structure(
        list(
            V_sand       = NULL,
            DER          = NULL,
            param_labels = NULL,
            model_type   = "weighted",
            N            = 100L,
            D            = 5L,
            survey_info  = NULL,
            matrix_diagnostics = NULL,
            nearPD_applied = FALSE
        ),
        class = "hbb_sandwich"
    )
    expect_output(print(sand), "Sandwich Variance")
})


test_that("print.hbb_sandwich: does not error on corrupted object", {
    sand <- structure(
        list(V_sand = "not_a_matrix"),
        class = "hbb_sandwich"
    )
    expect_output(print(sand))
})


test_that("print.hbb_sandwich: shows nearPD note when applied", {
    D <- 2L
    mock_sandwich <- structure(
        list(
            V_sand       = diag(D),
            H_obs_inv    = diag(D),
            Sigma_MCMC   = diag(D),
            DER          = c(a = 1, b = 1),
            param_labels = c("a", "b"),
            D            = D,
            N            = 50L,
            model_type   = "weighted",
            survey_info  = list(n_strata = 1L, n_psu = 2L,
                                df = 1L, n_singleton = 0L),
            matrix_diagnostics = list(
                V_sand = list(min_eig = 1, max_eig = 1, cond = 1, is_pd = TRUE)
            ),
            nearPD_applied = TRUE
        ),
        class = "hbb_sandwich"
    )

    output <- capture.output(print(mock_sandwich))
    expect_true(any(grepl("nearPD", output)))
})


test_that("print.hbb_sandwich: shows singleton warning", {
    D <- 1L
    mock_sandwich <- structure(
        list(
            V_sand       = diag(D),
            H_obs_inv    = diag(D),
            Sigma_MCMC   = diag(D),
            DER          = c(p1 = 1),
            param_labels = "p1",
            D            = D,
            N            = 10L,
            model_type   = "weighted",
            survey_info  = list(n_strata = 5L, n_psu = 8L,
                                df = 3L, n_singleton = 2L),
            matrix_diagnostics = NULL,
            nearPD_applied = FALSE
        ),
        class = "hbb_sandwich"
    )

    output <- capture.output(print(mock_sandwich))
    expect_true(any(grepl("WARNING|singleton", output)))
})


# ============================================================================
# I. DER ~ 1 under uniform weights --- Validation principle
# ============================================================================

test_that("DER is near 1 with IID scores and uniform weights", {
    # Under the information identity H = sum s_i s_i', and if
    # scores are IID across PSUs, J_cluster should approximate H_obs,
    # giving DER near 1.
    set.seed(2024)
    N <- 120L
    D <- 3L
    n_strata <- 4L
    n_psu_per <- 3L

    stratum_idx <- rep(seq_len(n_strata), each = N / n_strata)
    psu_idx <- rep(seq_len(n_strata * n_psu_per),
                   each = N / (n_strata * n_psu_per))
    w_tilde <- rep(1, N)

    # IID scores (no within-PSU correlation)
    scores <- matrix(rnorm(N * D), N, D)

    mock_hbb_data <- list(
        stratum_idx = stratum_idx,
        psu_idx     = psu_idx,
        w_tilde     = w_tilde,
        N           = N
    )

    # H_obs via information identity
    H_obs <- crossprod(scores)
    H_obs_inv <- solve(H_obs)

    # J_cluster from our function
    J <- compute_J_cluster(scores, mock_hbb_data)

    # V_sand
    V_sand <- H_obs_inv %*% J %*% H_obs_inv

    # DER
    DER <- diag(V_sand) / diag(H_obs_inv)

    # With IID scores, DER should be between 0.5 and 2.0
    expect_true(mean(DER) > 0.5,
                label = sprintf("mean DER = %.3f should be > 0.5", mean(DER)))
    expect_true(mean(DER) < 2.0,
                label = sprintf("mean DER = %.3f should be < 2.0", mean(DER)))
})


# ============================================================================
# J. Integration test --- Full pipeline (requires CmdStan)
# ============================================================================

test_that("sandwich_variance: full integration with hbb() fit", {
    skip_if_no_cmdstan()
    skip_if_no_hbb_models()
    skip_on_cran_or_ci()

    data(nsece_synth_small, package = "hurdlebb")

    fit <- hbb(
        y | trials(n_trial) ~ poverty + urban,
        data      = nsece_synth_small,
        weights   = "weight",
        stratum   = "stratum",
        psu       = "psu",
        chains    = 2L,
        iter_warmup   = 200L,
        iter_sampling = 200L,
        seed      = 12345L,
        refresh   = 0L
    )

    expect_true(is.hbb_fit(fit))
    expect_equal(fit$model_type, "weighted")

    # Compute sandwich variance
    sand <- sandwich_variance(fit)

    # Check return class and structure
    expect_s3_class(sand, "hbb_sandwich")
    expect_true(!is.null(sand$V_sand))
    expect_true(!is.null(sand$H_obs))
    expect_true(!is.null(sand$H_obs_inv))
    expect_true(!is.null(sand$J_cluster))
    expect_true(!is.null(sand$Sigma_MCMC))
    expect_true(!is.null(sand$DER))
    expect_true(!is.null(sand$param_labels))

    # Dimensions: D = 2P + 1 = 2*3 + 1 = 7 (intercept + poverty + urban)
    P <- fit$hbb_data$P
    D <- 2L * P + 1L
    expect_equal(sand$D, D)
    expect_equal(dim(sand$V_sand), c(D, D))
    expect_equal(dim(sand$H_obs), c(D, D))
    expect_equal(dim(sand$J_cluster), c(D, D))
    expect_equal(dim(sand$Sigma_MCMC), c(D, D))
    expect_length(sand$DER, D)
    expect_length(sand$param_labels, D)

    # V_sand should be symmetric and positive (semi-)definite
    expect_true(isSymmetric(sand$V_sand, tol = 1e-8))
    eig_V <- eigen(sand$V_sand, symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(eig_V >= -1e-10))

    # DER should be positive
    expect_true(all(sand$DER > 0))

    # H_obs should be block-diagonal: off-diagonal blocks are zero
    off_diag_block <- sand$H_obs[1:P, (P + 1):D]
    expect_equal(max(abs(off_diag_block)), 0, tolerance = 1e-15)

    # compute_der should match stored DER
    der_computed <- compute_der(sand)
    expect_equal(der_computed, sand$DER, tolerance = 1e-12)

    # Print method should not error
    expect_output(print(sand), "Sandwich Variance")
})


test_that("sandwich_variance: score matrix has near-zero column means", {
    skip_if_no_cmdstan()
    skip_if_no_hbb_models()
    skip_on_cran_or_ci()

    data(nsece_synth_small, package = "hurdlebb")

    fit <- hbb(
        y | trials(n_trial) ~ poverty + urban,
        data      = nsece_synth_small,
        weights   = "weight",
        stratum   = "stratum",
        psu       = "psu",
        chains    = 2L,
        iter_warmup   = 200L,
        iter_sampling = 200L,
        seed      = 12345L,
        refresh   = 0L
    )

    S <- compute_score_matrix(fit)

    N <- fit$hbb_data$N
    P <- fit$hbb_data$P
    D <- 2L * P + 1L

    expect_equal(dim(S), c(N, D))

    # Score column means should be near zero at the posterior mean
    col_means <- colMeans(S)
    for (d in seq_len(D)) {
        expect_true(
            abs(col_means[d]) < 5,
            label = sprintf("Score mean for column %d = %.4f (should be near 0)",
                            d, col_means[d])
        )
    }
})


# ============================================================================
# K. Additional Coverage Tests --- Internal Helpers & Error Paths
# ============================================================================


# ---- K1: .extract_posterior_means (lines 1048-1089) ----

test_that(".extract_posterior_means: errors when draws extraction fails", {
    # Mock cmdstan_fit whose draws() always fails
    mock_cmdstan <- list(
        draws = function(variables = NULL, format = "draws_matrix") {
            stop("variable not found in Stan output")
        }
    )
    expect_error(
        hurdlebb:::.extract_posterior_means(
            mock_cmdstan, "score_ext", n_rows = 10, n_cols = 3,
            label = "score_ext"
        ),
        "score_ext"
    )
})


test_that(".extract_posterior_means: errors when draws are NULL or empty", {
    mock_cmdstan <- list(
        draws = function(variables = NULL, format = "draws_matrix") {
            matrix(numeric(0), nrow = 0, ncol = 0)
        }
    )
    expect_error(
        hurdlebb:::.extract_posterior_means(
            mock_cmdstan, "score_ext", n_rows = 10, n_cols = 3,
            label = "score_ext"
        ),
        "No draws|score_ext"
    )
})


test_that(".extract_posterior_means: errors on dimension mismatch", {
    # Returns 20 elements but we expect 30 (10 x 3)
    mock_cmdstan <- list(
        draws = function(variables = NULL, format = "draws_matrix") {
            matrix(rnorm(100), nrow = 5, ncol = 20)  # 20 columns != 30 expected
        }
    )
    expect_error(
        hurdlebb:::.extract_posterior_means(
            mock_cmdstan, "score_ext", n_rows = 10, n_cols = 3,
            label = "score_ext"
        ),
        "mismatch|Dimension|expected"
    )
})


test_that(".extract_posterior_means: returns correct matrix for valid input", {
    # 50 draws x 6 variables (N=3, P=2, so 3*2=6)
    set.seed(123)
    data_mat <- matrix(rnorm(50 * 6), nrow = 50, ncol = 6)
    mock_cmdstan <- list(
        draws = function(variables = NULL, format = "draws_matrix") data_mat
    )
    result <- hurdlebb:::.extract_posterior_means(
        mock_cmdstan, "test_var", n_rows = 3, n_cols = 2, label = "test_var"
    )
    expect_equal(dim(result), c(3, 2))
    # Values should be column means of data_mat, reshaped
    expected_means <- colMeans(data_mat)
    expect_equal(as.numeric(t(result)), expected_means, tolerance = 1e-10)
})


# ---- K2: .compute_sigma_mcmc (lines 1274-1299) ----

test_that(".compute_sigma_mcmc: returns D x D covariance matrix for valid input", {
    set.seed(456)
    P <- 2L
    D <- 2L * P + 1L
    # param_names: alpha[1], alpha[2], beta[1], beta[2], log_kappa
    draws_mat <- matrix(rnorm(100 * D), nrow = 100, ncol = D)
    colnames(draws_mat) <- c("alpha[1]", "alpha[2]",
                              "beta[1]", "beta[2]", "log_kappa")

    mock_cmdstan <- list(
        draws = function(variables = NULL, format = "draws_matrix") {
            if (!is.null(variables)) {
                idx <- match(variables, colnames(draws_mat))
                return(draws_mat[, idx, drop = FALSE])
            }
            draws_mat
        }
    )

    result <- hurdlebb:::.compute_sigma_mcmc(mock_cmdstan, P)
    expect_equal(dim(result), c(D, D))
    expect_true(isSymmetric(result, tol = 1e-10))
    # Diagonal should be positive (variances)
    expect_true(all(diag(result) > 0))
})


test_that(".compute_sigma_mcmc: returns NA matrix when draws fail", {
    P <- 2L
    D <- 2L * P + 1L

    mock_cmdstan <- list(
        draws = function(variables = NULL, format = "draws_matrix") {
            stop("cannot extract draws")
        }
    )

    expect_warning(
        result <- hurdlebb:::.compute_sigma_mcmc(mock_cmdstan, P),
        "Failed"
    )
    expect_equal(dim(result), c(D, D))
    expect_true(all(is.na(result)))
})


# ---- K3: .check_pd additional paths ----

test_that(".check_pd: handles eigenvalue computation failure", {
    # Create a "matrix" that causes eigen() to fail
    # We can mock this by passing something weird through tryCatch path
    # Actually, non-square matrix would cause eigen to fail
    mat <- matrix(c(1, NA, NA, 1), 2, 2)
    result <- hurdlebb:::.check_pd(mat, label = "test_na_matrix")
    # Eigen of a matrix with NA succeeds (returns NA eigenvalues)
    # so check the NA handling path
    expect_false(result$is_pd)
})


test_that(".check_pd: near-zero eigenvalue (PSD but not PD)", {
    # Matrix that is PSD but not PD: one eigenvalue = 0
    mat <- matrix(c(1, 1, 1, 1), 2, 2)  # eigenvalues: 2, 0
    result <- hurdlebb:::.check_pd(mat, label = "psd_matrix")
    expect_false(result$is_pd)
    expect_equal(result$min_eig, 0, tolerance = 1e-14)
    expect_equal(result$max_eig, 2, tolerance = 1e-14)
})


test_that(".check_pd: large well-conditioned matrix", {
    D <- 10L
    mat <- diag(seq(1, 10))
    result <- hurdlebb:::.check_pd(mat, label = "large_diag")
    expect_true(result$is_pd)
    expect_equal(result$min_eig, 1, tolerance = 1e-14)
    expect_equal(result$max_eig, 10, tolerance = 1e-14)
    expect_equal(result$cond, 10, tolerance = 1e-10)
})


# ---- K4: .ridge_regularize additional paths ----

test_that(".ridge_regularize: severely non-PD matrix", {
    # Matrix with very negative eigenvalue
    mat <- diag(c(-10, 1, 1))
    result <- hurdlebb:::.ridge_regularize(mat)
    eig <- eigen(result, symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(eig > 0))
})


test_that(".ridge_regularize: respects eps parameter", {
    mat <- matrix(c(1, 0.99, 0.99, 1), 2, 2)  # nearly singular
    result1 <- hurdlebb:::.ridge_regularize(mat, eps = 1e-6)
    result2 <- hurdlebb:::.ridge_regularize(mat, eps = 1e-2)
    # Larger eps should add more ridge
    eig1 <- min(eigen(result1, symmetric = TRUE, only.values = TRUE)$values)
    eig2 <- min(eigen(result2, symmetric = TRUE, only.values = TRUE)$values)
    expect_true(eig2 >= eig1)
})


# ---- K5: compute_der --- warnings (lines 797-820) ----

test_that("compute_der: warns on unusually low DER (< 0.5)", {
    D <- 3L
    param_labels <- c("alpha_int", "beta_int", "log_kappa")
    mock_sand <- structure(
        list(
            V_sand       = diag(c(0.001, 0.02, 0.03)),
            H_obs_inv    = diag(c(0.01, 0.02, 0.03)),  # DER[1] = 0.1 < 0.5
            param_labels = param_labels,
            D            = D
        ),
        class = "hbb_sandwich"
    )
    expect_warning(
        result <- compute_der(mock_sand),
        "low DER"
    )
    expect_equal(result[1], c(alpha_int = 0.1), tolerance = 1e-10)
})


test_that("compute_der: warns on unusually high DER (> 20)", {
    D <- 3L
    param_labels <- c("alpha_int", "beta_int", "log_kappa")
    mock_sand <- structure(
        list(
            V_sand       = diag(c(0.5, 0.02, 0.03)),
            H_obs_inv    = diag(c(0.01, 0.02, 0.03)),  # DER[1] = 50 > 20
            param_labels = param_labels,
            D            = D
        ),
        class = "hbb_sandwich"
    )
    expect_warning(
        result <- compute_der(mock_sand),
        "high DER"
    )
    expect_equal(result[1], c(alpha_int = 50), tolerance = 1e-10)
})


test_that("compute_der: errors on non-hbb_sandwich input", {
    expect_error(compute_der(list(V_sand = diag(2))), "hbb_sandwich")
    expect_error(compute_der(42), "hbb_sandwich")
})


test_that("compute_der: no warning when DER is in normal range", {
    D <- 3L
    param_labels <- c("alpha_int", "beta_int", "log_kappa")
    mock_sand <- structure(
        list(
            V_sand       = diag(c(0.02, 0.04, 0.06)),
            H_obs_inv    = diag(c(0.01, 0.02, 0.03)),  # DER = 2 everywhere
            param_labels = param_labels,
            D            = D
        ),
        class = "hbb_sandwich"
    )
    expect_no_warning(
        result <- compute_der(mock_sand)
    )
    expect_equal(unname(result), c(2, 2, 2), tolerance = 1e-10)
})


# ---- K6: .validate_sandwich_fit --- additional paths ----

test_that(".validate_sandwich_fit: errors when hbb_data is NULL", {
    mock_fit <- structure(
        list(
            model_type = "weighted",
            fit        = "present",
            hbb_data   = NULL
        ),
        class = "hbb_fit"
    )
    expect_error(
        hurdlebb:::.validate_sandwich_fit(mock_fit),
        "hbb_data"
    )
})


test_that(".validate_sandwich_fit: errors when psu_idx is NULL", {
    mock_fit <- structure(
        list(
            model_type = "weighted",
            fit        = "present",
            hbb_data   = list(stratum_idx = 1:10, psu_idx = NULL)
        ),
        class = "hbb_fit"
    )
    expect_error(
        hurdlebb:::.validate_sandwich_fit(mock_fit),
        "psu_idx"
    )
})


test_that(".validate_sandwich_fit: passes on valid weighted fit", {
    mock_fit <- structure(
        list(
            model_type = "weighted",
            fit        = "present",
            hbb_data   = list(stratum_idx = 1:10, psu_idx = 1:10)
        ),
        class = "hbb_fit"
    )
    expect_true(hurdlebb:::.validate_sandwich_fit(mock_fit))
})


test_that(".validate_sandwich_fit: accepts svc_weighted", {
    mock_fit <- structure(
        list(
            model_type = "svc_weighted",
            fit        = "present",
            hbb_data   = list(stratum_idx = 1:10, psu_idx = 1:10)
        ),
        class = "hbb_fit"
    )
    expect_true(hurdlebb:::.validate_sandwich_fit(mock_fit))
})


# ---- K7: print.hbb_sandwich additional coverage ----

test_that("print.hbb_sandwich: shows PSD and NO labels in matrix diagnostics", {
    D <- 2L
    mock_sandwich <- structure(
        list(
            V_sand       = diag(D),
            H_obs        = diag(D),
            H_obs_inv    = diag(D),
            J_cluster    = diag(D),
            Sigma_MCMC   = diag(D),
            DER          = c(a = 1.5, b = 2.5),
            param_labels = c("a", "b"),
            D            = D,
            N            = 50L,
            model_type   = "weighted",
            survey_info  = list(n_strata = 5L, n_psu = 20L,
                                df = 15L, n_singleton = 0L),
            matrix_diagnostics = list(
                V_sand = list(min_eig = -1e-11, max_eig = 1, cond = 1e11,
                              is_pd = FALSE),   # PSD (near zero neg)
                H_obs = list(min_eig = -5, max_eig = 10, cond = 2,
                             is_pd = FALSE),     # NO (definitely not PD)
                J_cluster = list(min_eig = 1, max_eig = 1, is_pd = TRUE),
                Sigma_MCMC = list(min_eig = 1, max_eig = 1, cond = 1,
                                  is_pd = TRUE)
            ),
            nearPD_applied = FALSE,
            call = quote(sandwich_variance(fit))
        ),
        class = "hbb_sandwich"
    )
    output <- capture.output(print(mock_sandwich))
    output_text <- paste(output, collapse = "\n")
    expect_true(grepl("PSD", output_text))
    expect_true(grepl("NO", output_text))
})


test_that("print.hbb_sandwich: shows flag markers for DER values", {
    D <- 3L
    param_labels <- c("alpha_int", "beta_int", "log_kappa")
    mock_sandwich <- structure(
        list(
            V_sand       = diag(c(0.001, 0.08, 0.5)),
            H_obs        = diag(c(100, 50, 2)),
            H_obs_inv    = diag(c(0.01, 0.01, 0.01)),
            J_cluster    = diag(D),
            Sigma_MCMC   = diag(D) * 0.1,
            DER          = c(alpha_int = 0.1, beta_int = 8.0,
                             log_kappa = 50.0),
            param_labels = param_labels,
            D            = D,
            N            = 100L,
            model_type   = "weighted",
            survey_info  = list(n_strata = 10L, n_psu = 50L,
                                df = 40L, n_singleton = 0L),
            matrix_diagnostics = list(
                V_sand = list(min_eig = 0.001, max_eig = 0.5, cond = 500,
                              is_pd = TRUE),
                H_obs = list(min_eig = 2, max_eig = 100, cond = 50,
                             is_pd = TRUE),
                J_cluster = list(min_eig = 1, max_eig = 1, is_pd = TRUE),
                Sigma_MCMC = list(min_eig = 0.1, max_eig = 0.1, cond = 1,
                                  is_pd = TRUE)
            ),
            nearPD_applied = FALSE,
            call = quote(sandwich_variance(fit))
        ),
        class = "hbb_sandwich"
    )
    output <- capture.output(print(mock_sandwich))
    output_text <- paste(output, collapse = "\n")
    # [LOW], [*], [HIGH] flags
    expect_true(grepl("\\[LOW\\]", output_text))
    expect_true(grepl("\\[\\*\\]", output_text))
    expect_true(grepl("\\[HIGH\\]", output_text))
})


# ---- K8: compute_score_matrix error paths ----

test_that("compute_score_matrix: errors on non-weighted model types", {
    fake_fit <- structure(list(model_type = "svc"), class = "hbb_fit")
    expect_error(compute_score_matrix(fake_fit), "weighted")

    fake_fit2 <- structure(list(model_type = "base"), class = "hbb_fit")
    expect_error(compute_score_matrix(fake_fit2), "weighted")
})


# ---- K9: compute_H_obs error path ----

test_that("compute_H_obs: errors on missing weights", {
    mock_fit <- structure(
        list(
            model_type = "weighted",
            hbb_data   = list(
                N = 5L, P = 2L,
                X = matrix(1, 5, 2),
                w_tilde = NULL   # no weights
            )
        ),
        class = "hbb_fit"
    )
    scores <- matrix(rnorm(5 * 5), nrow = 5, ncol = 5)
    expect_error(
        compute_H_obs(mock_fit, scores = scores),
        "weights|w_tilde"
    )
})


# ---- K10: .block_diag ----

test_that(".block_diag: assembles two blocks correctly", {
    A <- matrix(c(1, 2, 3, 4), 2, 2)
    B <- matrix(c(5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16), 4, 3)
    # B is 4x3 so it should fail since .block_diag expects square matrices?
    # Let me use square matrices
    B <- matrix(c(5, 6, 7, 8, 9, 10, 11, 12, 13), 3, 3)
    result <- hurdlebb:::.block_diag(A, B)
    expect_equal(dim(result), c(5, 5))
    # Top-left block = A
    expect_equal(result[1:2, 1:2], A)
    # Bottom-right block = B
    expect_equal(result[3:5, 3:5], B)
    # Off-diagonal blocks = 0
    expect_true(all(result[1:2, 3:5] == 0))
    expect_true(all(result[3:5, 1:2] == 0))
})
