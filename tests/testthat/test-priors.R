# ============================================================================
# test-priors.R --- Tests for hbb_prior(), default_prior(),
#                   and print.hbb_prior()
#
# Covers:
#   Section 1: default_prior()
#   Section 2: hbb_prior() custom specifications
#   Section 3: Error paths
#   Section 4: print.hbb_prior()
# ============================================================================


# ============================================================================
# Section 1: default_prior()
# ============================================================================

test_that("default_prior: returns 'hbb_prior' class", {
    p <- default_prior()

    expect_s3_class(p, "hbb_prior")
    expect_true(inherits(p, "hbb_prior"))
})


test_that("default_prior: has all expected fields", {
    p <- default_prior()

    expected_fields <- c("alpha", "beta", "log_kappa", "gamma", "tau",
                         "lkj_eta")
    for (field in expected_fields) {
        expect_true(field %in% names(p),
                    info = paste("Missing field:", field))
    }
})


test_that("default_prior: default values match specification", {
    p <- default_prior()

    # alpha ~ Normal(0, 2)
    expect_equal(p$alpha$dist, "normal")
    expect_equal(p$alpha$mean, 0)
    expect_equal(p$alpha$sd, 2)

    # beta ~ Normal(0, 2)
    expect_equal(p$beta$dist, "normal")
    expect_equal(p$beta$mean, 0)
    expect_equal(p$beta$sd, 2)

    # log_kappa ~ Normal(2, 1.5)
    expect_equal(p$log_kappa$dist, "normal")
    expect_equal(p$log_kappa$mean, 2)
    expect_equal(p$log_kappa$sd, 1.5)

    # gamma ~ Normal(0, 1)
    expect_equal(p$gamma$dist, "normal")
    expect_equal(p$gamma$mean, 0)
    expect_equal(p$gamma$sd, 1)

    # tau ~ Normal(0, 1)
    expect_equal(p$tau$dist, "normal")
    expect_equal(p$tau$mean, 0)
    expect_equal(p$tau$sd, 1)

    # lkj_eta = 2
    expect_equal(p$lkj_eta, 2)
})


# ============================================================================
# Section 2: hbb_prior() custom specifications
# ============================================================================

test_that("hbb_prior: custom alpha prior", {
    p <- hbb_prior(alpha = list(dist = "normal", mean = 1, sd = 5))

    expect_s3_class(p, "hbb_prior")
    expect_equal(p$alpha$dist, "normal")
    expect_equal(p$alpha$mean, 1)
    expect_equal(p$alpha$sd, 5)
    # Other fields should retain defaults
    expect_equal(p$beta$dist, "normal")
    expect_equal(p$beta$mean, 0)
    expect_equal(p$beta$sd, 2)
})


test_that("hbb_prior: custom log_kappa prior", {
    p <- hbb_prior(log_kappa = list(dist = "normal", mean = 3, sd = 2))

    expect_equal(p$log_kappa$dist, "normal")
    expect_equal(p$log_kappa$mean, 3)
    expect_equal(p$log_kappa$sd, 2)
    # alpha should retain default
    expect_equal(p$alpha$mean, 0)
    expect_equal(p$alpha$sd, 2)
})


test_that("hbb_prior: custom lkj_eta", {
    p <- hbb_prior(lkj_eta = 5)

    expect_equal(p$lkj_eta, 5)
    # Rest should be defaults
    expect_equal(p$alpha$mean, 0)
})


test_that("hbb_prior: partial specification fills rest from defaults", {
    # Specify only alpha and lkj_eta; the rest should come from defaults
    p <- hbb_prior(
        alpha   = list(dist = "normal", mean = -1, sd = 3),
        lkj_eta = 10
    )

    expect_s3_class(p, "hbb_prior")

    # Custom values
    expect_equal(p$alpha$mean, -1)
    expect_equal(p$alpha$sd, 3)
    expect_equal(p$lkj_eta, 10)

    # Default values (unchanged)
    expect_equal(p$beta$dist, "normal")
    expect_equal(p$beta$mean, 0)
    expect_equal(p$beta$sd, 2)
    expect_equal(p$log_kappa$mean, 2)
    expect_equal(p$log_kappa$sd, 1.5)
    expect_equal(p$gamma$mean, 0)
    expect_equal(p$gamma$sd, 1)
    expect_equal(p$tau$mean, 0)
    expect_equal(p$tau$sd, 1)
})


test_that("hbb_prior: all parameters custom", {
    p <- hbb_prior(
        alpha     = list(dist = "normal", mean = 1, sd = 10),
        beta      = list(dist = "normal", mean = -1, sd = 5),
        log_kappa = list(dist = "normal", mean = 0, sd = 3),
        gamma     = list(dist = "normal", mean = 0.5, sd = 2),
        tau       = list(dist = "normal", mean = 0, sd = 0.5),
        lkj_eta   = 1
    )

    expect_s3_class(p, "hbb_prior")
    expect_equal(p$alpha$sd, 10)
    expect_equal(p$beta$mean, -1)
    expect_equal(p$log_kappa$sd, 3)
    expect_equal(p$gamma$mean, 0.5)
    expect_equal(p$tau$sd, 0.5)
    expect_equal(p$lkj_eta, 1)
})


test_that("hbb_prior: returns identical to default_prior when called with no args", {
    p_default <- default_prior()
    p_empty   <- hbb_prior()

    # Compare all fields
    expect_equal(p_default$alpha, p_empty$alpha)
    expect_equal(p_default$beta, p_empty$beta)
    expect_equal(p_default$log_kappa, p_empty$log_kappa)
    expect_equal(p_default$gamma, p_empty$gamma)
    expect_equal(p_default$tau, p_empty$tau)
    expect_equal(p_default$lkj_eta, p_empty$lkj_eta)
})


# ============================================================================
# Section 3: Error paths
# ============================================================================

test_that("hbb_prior: errors on negative sd", {
    expect_error(
        hbb_prior(alpha = list(dist = "normal", mean = 0, sd = -1)),
        regexp = "sd|negative|positive|greater"
    )
})


test_that("hbb_prior: errors on zero sd", {
    expect_error(
        hbb_prior(beta = list(dist = "normal", mean = 0, sd = 0)),
        regexp = "sd|zero|positive|greater"
    )
})


test_that("hbb_prior: errors on non-numeric mean", {
    expect_error(
        hbb_prior(alpha = list(dist = "normal", mean = "zero", sd = 2)),
        regexp = "mean|numeric|character"
    )
})


test_that("hbb_prior: errors on negative lkj_eta", {
    expect_error(
        hbb_prior(lkj_eta = -1),
        regexp = "lkj_eta|negative|positive|greater"
    )
})


test_that("hbb_prior: errors on zero lkj_eta", {
    expect_error(
        hbb_prior(lkj_eta = 0),
        regexp = "lkj_eta|zero|positive|greater"
    )
})


test_that("hbb_prior: errors on invalid dist name", {
    expect_error(
        hbb_prior(alpha = list(dist = "cauchy", mean = 0, sd = 2)),
        regexp = "dist|distribution|normal|supported|invalid"
    )
})


test_that("hbb_prior: errors when alpha is not a list", {
    expect_error(
        hbb_prior(alpha = c(0, 2)),
        regexp = "list|alpha"
    )
})


test_that("hbb_prior: errors when alpha is a numeric scalar", {
    expect_error(
        hbb_prior(alpha = 5),
        regexp = "list|alpha"
    )
})


test_that("hbb_prior: errors when beta is a character string", {
    expect_error(
        hbb_prior(beta = "normal(0,2)"),
        regexp = "list|beta"
    )
})


test_that("hbb_prior: errors on missing sd in list", {
    expect_error(
        hbb_prior(alpha = list(dist = "normal", mean = 0)),
        regexp = "sd|missing|required"
    )
})


test_that("hbb_prior: errors on missing mean in list", {
    expect_error(
        hbb_prior(alpha = list(dist = "normal", sd = 2)),
        regexp = "mean|missing|required"
    )
})


test_that("hbb_prior: errors on missing dist in list", {
    expect_error(
        hbb_prior(alpha = list(mean = 0, sd = 2)),
        regexp = "dist|missing|required"
    )
})


test_that("hbb_prior: errors on non-numeric lkj_eta", {
    expect_error(
        hbb_prior(lkj_eta = "two"),
        regexp = "lkj_eta|numeric|character"
    )
})


test_that("hbb_prior: errors on Inf sd", {
    expect_error(
        hbb_prior(alpha = list(dist = "normal", mean = 0, sd = Inf)),
        regexp = "sd|finite|Inf"
    )
})


test_that("hbb_prior: errors on NA mean", {
    expect_error(
        hbb_prior(alpha = list(dist = "normal", mean = NA, sd = 2)),
        regexp = "mean|NA|missing"
    )
})


# ============================================================================
# Section 4: print.hbb_prior()
# ============================================================================

test_that("print.hbb_prior: prints without error", {
    p <- default_prior()
    expect_output(print(p))
})


test_that("print.hbb_prior: output contains distribution names", {
    p <- default_prior()

    out <- capture.output(print(p))
    combined <- paste(out, collapse = "\n")

    # Should mention "normal" (the distribution for all default priors)
    expect_true(grepl("normal", combined, ignore.case = TRUE),
                info = "'normal' not found in print output")

    # Should mention key parameter names
    expect_true(grepl("alpha", combined, ignore.case = TRUE),
                info = "'alpha' not found in print output")
    expect_true(grepl("beta", combined, ignore.case = TRUE),
                info = "'beta' not found in print output")
    expect_true(grepl("kappa", combined, ignore.case = TRUE),
                info = "'kappa' not found in print output")
    expect_true(grepl("lkj", combined, ignore.case = TRUE),
                info = "'lkj' not found in print output")
})


test_that("print.hbb_prior: custom prior shows custom values", {
    p <- hbb_prior(alpha = list(dist = "normal", mean = 5, sd = 10))

    out <- capture.output(print(p))
    combined <- paste(out, collapse = "\n")

    # The custom mean=5 and sd=10 should appear
    expect_true(grepl("5", combined),
                info = "Custom mean not found in print output")
    expect_true(grepl("10", combined),
                info = "Custom sd not found in print output")
})


test_that("print.hbb_prior: returns the object invisibly", {
    p <- default_prior()

    result <- withr::with_output_sink(tempfile(), print(p))
    # print methods should return the object invisibly
    expect_s3_class(result, "hbb_prior")
})
