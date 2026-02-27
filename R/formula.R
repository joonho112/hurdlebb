# ============================================================================
# formula.R --- Formula Parser for Hurdle Beta-Binomial Models
#
# Parses brms-style formulas into a structured S3 object (class "hbb_formula")
# that downstream fitting functions use to construct model matrices.
#
# Supported formula syntax:
#   y | trials(n) ~ x1 + x2                          # fixed effects only
#   y | trials(n) ~ x1 + (1 | group)                 # random intercept
#   y | trials(n) ~ x1 + (x1 + x2 | group)           # SVC (random slopes)
#   y | trials(n) ~ x1 + (x1 | group) + state_level(v1 + v2)  # + policy mods
#
# Key design:
#   - Fixed effects enter BOTH hurdle margins (alpha and beta share X)
#   - Intercept is ALWAYS implicitly included (not in $fixed)
#   - state_level() terms become the V matrix (policy moderators)
#   - Random slopes must be a subset of fixed effects
#   - Dot expansion (y ~ .) is NOT supported
#
# Contents:
#   1. hbb_formula            --- Main parser (exported)
#   2. validate_hbb_formula   --- Validate formula against data (exported)
#   3. print.hbb_formula      --- Print method for hbb_formula objects
#   4. .parse_lhs             --- Parse left-hand side (response + trials)
#   5. .parse_rhs             --- Parse right-hand side (orchestrator)
#   6. .parse_state_level     --- Extract state_level() terms
#   7. .parse_random_effects  --- Extract (... | group) terms
#   8. .parse_fixed_effects   --- Extract remaining fixed effect terms
#   9. .parse_bar_terms       --- Parse a single (... | group) expression
# ============================================================================


# ============================================================================
# 1. hbb_formula — Main parser (exported)
# ============================================================================

#' Parse a Hurdle Beta-Binomial Model Formula
#'
#' Converts a brms-style formula into a structured specification for a
#' two-part hurdle Beta-Binomial model. The fixed effects specified in
#' the formula enter **both** model margins (extensive and intensive)
#' with separate coefficient vectors \eqn{\alpha} and \eqn{\beta}.
#'
#' @param formula A two-sided formula of the form
#'   `y | trials(n) ~ predictors`. See **Details** for the full syntax.
#'
#' @return An S3 object of class `"hbb_formula"` with components:
#' \describe{
#'   \item{`response`}{Character. The response variable name.}
#'   \item{`trials`}{Character. The trials variable name.}
#'   \item{`fixed`}{Character vector of fixed effect term names,
#'     **excluding** the intercept (which is always implicit).}
#'   \item{`group`}{Character or `NULL`. The grouping variable for
#'     random effects.}
#'   \item{`random`}{Character vector of random effect terms, or `NULL`
#'     if no random effects. Contains `"1"` for intercept-only random
#'     effects, or variable names for state-varying coefficients.}
#'   \item{`policy`}{Character vector of policy moderator term names,
#'     or `NULL` if none specified.}
#'   \item{`formula`}{The original formula object.}
#'   \item{`svc`}{Logical. `TRUE` if the model includes random slopes
#'     (state-varying coefficients), not just a random intercept.}
#' }
#'
#' @details
#' ## Formula syntax
#'
#' The left-hand side **must** specify both the response variable and
#' the trials variable using the `| trials()` syntax:
#' ```
#' y | trials(n_trial) ~ ...
#' ```
#'
#' The right-hand side supports four types of terms:
#'
#' 1. **Fixed effects**: Standard R formula terms (e.g., `x1 + x2`).
#'    An intercept is always included automatically. Use `y | trials(n) ~ 1`
#'    for an intercept-only model.
#'
#' 2. **Random intercept**: `(1 | group)` specifies a random intercept
#'    grouped by `group`.
#'
#' 3. **State-varying coefficients (SVC)**: `(x1 + x2 | group)` specifies
#'    random slopes for `x1` and `x2` (plus the implicit random intercept).
#'    All random slope variables must also appear as fixed effects.
#'
#' 4. **Policy moderators**: `state_level(v1 + v2)` specifies cross-level
#'    interaction terms that predict the random effects. These require
#'    that random effects are also specified.
#'
#' ## Model mapping
#'
#' The parsed formula maps to the model as follows:
#' - Fixed effects \eqn{\to} the \eqn{X} matrix (with auto-prepended
#'   intercept column)
#' - The grouping variable \eqn{\to} the state index vector
#' - Policy moderators \eqn{\to} the \eqn{V} matrix (with auto-prepended
#'   intercept column)
#' - Both hurdle margins share the same \eqn{X} design matrix
#'
#' @examples
#' # Intercept-only model
#' f1 <- hbb_formula(y | trials(n_trial) ~ 1)
#' f1
#'
#' # Fixed effects only
#' f2 <- hbb_formula(y | trials(n_trial) ~ poverty + urban)
#' f2
#'
#' # Random intercept
#' f3 <- hbb_formula(y | trials(n_trial) ~ poverty + (1 | state_id))
#' f3
#'
#' # SVC with policy moderators
#' f4 <- hbb_formula(
#'   y | trials(n_trial) ~ poverty + urban +
#'     (poverty + urban | state_id) +
#'     state_level(mr_pctile + tiered_reim)
#' )
#' f4
#'
#' @seealso [validate_hbb_formula()] to check a formula against a dataset.
#' @family formula
#' @export
hbb_formula <- function(formula) {

    # -- Input validation -----------------------------------------------------
    if (!inherits(formula, "formula")) {
        cli_abort(c(
            "Expected a formula object.",
            "x" = "Got {.cls {class(formula)}}.",
            "i" = "Use brms-style syntax: {.code y | trials(n) ~ x1 + x2}."
        ))
    }
    if (length(formula) != 3L) {
        cli_abort(c(
            "Formula must be two-sided.",
            "x" = "A right-hand side only formula was supplied.",
            "i" = "Use: {.code y | trials(n) ~ predictors}."
        ))
    }

    # -- Parse left-hand side -------------------------------------------------
    lhs <- .parse_lhs(formula)

    # -- Parse right-hand side ------------------------------------------------
    rhs <- .parse_rhs(formula)

    # -- Cross-component validation -------------------------------------------

    # Random slopes must be a subset of fixed effects
    if (!is.null(rhs$random) && !identical(rhs$random, "1")) {
        bad <- setdiff(rhs$random, rhs$fixed)
        if (length(bad) > 0L) {
            cli_abort(c(
                "{qty(length(bad))}Random slope variable{?s} must also appear
                 as fixed effect{?s}.",
                "x" = "Not in fixed effects: {.val {bad}}.",
                "i" = "Add {.code {paste(bad, collapse = ' + ')}} to the
                       right-hand side of the formula."
            ))
        }
    }

    # Policy moderators require random effects
    if (!is.null(rhs$policy) && is.null(rhs$group)) {
        cli_abort(c(
            "{.code state_level()} requires random effects.",
            "x" = "Policy moderators were specified but no
                   {.code (... | group)} term was found.",
            "i" = "Add a random effects term, e.g.,
                   {.code (1 | state_id)} or
                   {.code (x1 | state_id)}."
        ))
    }

    # -- Determine SVC flag ---------------------------------------------------
    svc <- !is.null(rhs$random) && !identical(rhs$random, "1")

    # -- Construct and return S3 object ---------------------------------------
    out <- structure(
        list(
            response = lhs$response,
            trials   = lhs$trials,
            fixed    = rhs$fixed,
            group    = rhs$group,
            random   = rhs$random,
            policy   = rhs$policy,
            formula  = formula,
            svc      = svc
        ),
        class = "hbb_formula"
    )

    out
}


# ============================================================================
# 2. validate_hbb_formula — Validate formula against data (exported)
# ============================================================================

#' Validate an hbb_formula Against a Dataset
#'
#' Checks that all variables referenced in the formula exist in the
#' supplied data frame(s) and that they have appropriate types and
#' ranges. This is typically called internally by the fitting function,
#' but can be used directly for early error checking.
#'
#' @param object An object of class `"hbb_formula"`, as returned by
#'   [hbb_formula()].
#' @param data A data frame containing provider-level variables
#'   (response, trials, fixed effects, grouping variable).
#' @param state_data An optional data frame containing state-level
#'   policy variables. Required if `object$policy` is not `NULL`.
#'   Must contain a column matching `object$group` for merging.
#'
#' @return Invisibly returns `TRUE` if all checks pass. Throws an
#'   informative error otherwise.
#'
#' @details
#' The following checks are performed:
#' \enumerate{
#'   \item `object` is an `hbb_formula`.
#'   \item `data` is a data frame with at least one row.
#'   \item The response variable exists and contains non-negative integers.
#'   \item The trials variable exists and contains positive integers.
#'   \item Response does not exceed trials for any observation.
#'   \item All fixed effect variables exist and are numeric.
#'   \item The grouping variable (if any) exists and is coercible to
#'     integer.
#'   \item If `state_data` is supplied: all policy variables exist in
#'     `state_data`, the grouping variable exists in `state_data`, and
#'     policy variables are numeric.
#' }
#'
#' @examples
#' f <- hbb_formula(y | trials(n_trial) ~ poverty + urban)
#' data(nsece_synth_small)
#' validate_hbb_formula(f, nsece_synth_small)
#'
#' @seealso [hbb_formula()]
#' @family formula
#' @export
validate_hbb_formula <- function(object, data, state_data = NULL) {

    # -- Check object class ---------------------------------------------------
    if (!inherits(object, "hbb_formula")) {
        cli_abort(c(
            "Expected an {.cls hbb_formula} object.",
            "x" = "Got {.cls {class(object)}}.",
            "i" = "Use {.fun hbb_formula} to create a formula specification."
        ))
    }

    # -- Check data -----------------------------------------------------------
    assert_data_frame(data, min.rows = 1L, .var.name = "data")

    # -- Response variable ----------------------------------------------------
    if (!object$response %in% names(data)) {
        cli_abort(c(
            "Response variable {.val {object$response}} not found in {.arg data}.",
            "i" = "Available columns: {.val {head(names(data), 20)}}."
        ))
    }
    resp <- data[[object$response]]
    if (!is.numeric(resp)) {
        cli_abort(c(
            "Response variable {.val {object$response}} must be numeric.",
            "x" = "Got {.cls {class(resp)}}."
        ))
    }
    if (any(is.nan(resp))) {
        cli_abort("Response variable {.val {object$response}} contains NaN values.")
    }
    if (any(is.infinite(resp))) {
        cli_abort("Response variable {.val {object$response}} contains infinite values.")
    }
    if (anyNA(resp)) {
        n_na <- sum(is.na(resp))
        cli_warn(c(
            "Response variable {.val {object$response}} contains {n_na}
             missing value{?s}.",
            "i" = "Missing values will be excluded during fitting."
        ))
    }
    if (any(resp < 0, na.rm = TRUE)) {
        cli_abort("Response variable {.val {object$response}} contains negative values.")
    }
    if (any(resp != floor(resp), na.rm = TRUE)) {
        cli_abort("Response variable {.val {object$response}} must contain integers (whole numbers).")
    }

    # -- Trials variable ------------------------------------------------------
    if (!object$trials %in% names(data)) {
        cli_abort(c(
            "Trials variable {.val {object$trials}} not found in {.arg data}.",
            "i" = "Available columns: {.val {head(names(data), 20)}}."
        ))
    }
    trials <- data[[object$trials]]
    if (!is.numeric(trials)) {
        cli_abort(c(
            "Trials variable {.val {object$trials}} must be numeric.",
            "x" = "Got {.cls {class(trials)}}."
        ))
    }
    if (any(is.nan(trials))) {
        cli_abort("Trials variable {.val {object$trials}} contains NaN values.")
    }
    if (any(is.infinite(trials))) {
        cli_abort("Trials variable {.val {object$trials}} contains infinite values.")
    }
    if (anyNA(trials)) {
        n_na <- sum(is.na(trials))
        cli_warn(c(
            "Trials variable {.val {object$trials}} contains {n_na}
             missing value{?s}.",
            "i" = "Missing values will be excluded during fitting."
        ))
    }
    if (any(trials < 1, na.rm = TRUE)) {
        cli_abort("Trials variable {.val {object$trials}} must be >= 1 for all observations.")
    }
    if (any(trials != floor(trials), na.rm = TRUE)) {
        cli_abort("Trials variable {.val {object$trials}} must contain integers (whole numbers).")
    }

    # -- Response <= Trials ---------------------------------------------------
    violations <- sum(resp > trials, na.rm = TRUE)
    if (violations > 0L) {
        cli_abort(c(
            "Response exceeds trials for {violations} observation{?s}.",
            "x" = "{.val {object$response}} > {.val {object$trials}} is not allowed.",
            "i" = "The response must satisfy 0 <= y <= n for all observations."
        ))
    }

    # -- Fixed effect variables -----------------------------------------------
    if (length(object$fixed) > 0L) {
        missing_fe <- setdiff(object$fixed, names(data))
        if (length(missing_fe) > 0L) {
            cli_abort(c(
                "{qty(length(missing_fe))}Fixed effect variable{?s} not found in {.arg data}.",
                "x" = "Missing: {.val {missing_fe}}.",
                "i" = "Available columns: {.val {head(names(data), 20)}}."
            ))
        }
        for (v in object$fixed) {
            if (!is.numeric(data[[v]])) {
                cli_abort(c(
                    "Fixed effect variable {.val {v}} must be numeric.",
                    "x" = "Got {.cls {class(data[[v]])}}."
                ))
            }
        }
    }

    # -- Grouping variable ----------------------------------------------------
    if (!is.null(object$group)) {
        if (!object$group %in% names(data)) {
            cli_abort(c(
                "Grouping variable {.val {object$group}} not found in {.arg data}.",
                "i" = "Available columns: {.val {head(names(data), 20)}}."
            ))
        }
        grp <- data[[object$group]]
        if (!is.numeric(grp) && !is.integer(grp) && !is.factor(grp) &&
            !is.character(grp)) {
            cli_abort(c(
                "Grouping variable {.val {object$group}} must be numeric, integer,
                 factor, or character.",
                "x" = "Got {.cls {class(grp)}}."
            ))
        }
    }

    # -- Policy moderators (state_data) ---------------------------------------
    if (!is.null(object$policy)) {
        if (is.null(state_data)) {
            cli_abort(c(
                "Policy moderators require {.arg state_data}.",
                "x" = "{.code state_level()} terms were specified but
                       {.arg state_data} was not provided.",
                "i" = "Supply a data frame of state-level policy variables."
            ))
        }
        assert_data_frame(state_data, min.rows = 1L,
                          .var.name = "state_data")

        # Grouping variable must exist in state_data for merging
        if (!object$group %in% names(state_data)) {
            cli_abort(c(
                "Grouping variable {.val {object$group}} not found in
                 {.arg state_data}.",
                "i" = "The grouping variable must appear in both {.arg data}
                       and {.arg state_data} for merging."
            ))
        }

        # State coverage: all states in data must appear in state_data
        data_states <- unique(data[[object$group]])
        sd_states <- unique(state_data[[object$group]])
        missing_states <- setdiff(data_states, sd_states)
        if (length(missing_states) > 0L) {
            cli_abort(c(
                "{qty(length(missing_states))}{length(missing_states)}
                 group{?s} in {.arg data} not found in {.arg state_data}.",
                "x" = "Missing: {.val {head(missing_states, 10)}}.",
                "i" = "{.arg state_data} must contain a row for every unique
                       value of {.val {object$group}} in {.arg data}."
            ))
        }

        # No duplicate states in state_data
        if (anyDuplicated(state_data[[object$group]])) {
            cli_abort(c(
                "{.arg state_data} contains duplicate entries for
                 {.val {object$group}}.",
                "i" = "Each group should appear exactly once in
                       {.arg state_data}."
            ))
        }

        # Policy variables must exist and be numeric
        missing_pol <- setdiff(object$policy, names(state_data))
        if (length(missing_pol) > 0L) {
            cli_abort(c(
                "{qty(length(missing_pol))}Policy variable{?s} not found in {.arg state_data}.",
                "x" = "Missing: {.val {missing_pol}}.",
                "i" = "Available columns: {.val {head(names(state_data), 20)}}."
            ))
        }
        for (v in object$policy) {
            if (!is.numeric(state_data[[v]])) {
                cli_abort(c(
                    "Policy variable {.val {v}} must be numeric.",
                    "x" = "Got {.cls {class(state_data[[v]])}}."
                ))
            }
        }
    }

    invisible(TRUE)
}


# ============================================================================
# 3. print.hbb_formula — Print method
# ============================================================================

#' @export
print.hbb_formula <- function(x, ...) {

    cat("Hurdle Beta-Binomial Formula\n")
    cat("----------------------------\n")
    cat("  Response :", x$response, "\n")
    cat("  Trials   :", x$trials, "\n")

    if (length(x$fixed) == 0L) {
        cat("  Fixed    : (intercept only)\n")
    } else {
        cat("  Fixed    : intercept +",
            paste(x$fixed, collapse = " + "), "\n")
    }

    if (is.null(x$group)) {
        cat("  Random   : none\n")
    } else if (identical(x$random, "1")) {
        cat("  Random   : (1 |", x$group, ")   [random intercept]\n")
    } else {
        cat("  Random   : (1 +",
            paste(x$random, collapse = " + "), "|",
            x$group, ")   [SVC]\n")
    }

    if (is.null(x$policy)) {
        cat("  Policy   : none\n")
    } else {
        cat("  Policy   : state_level(",
            paste(x$policy, collapse = " + "), ")\n")
    }

    cat("  SVC      :", x$svc, "\n")

    invisible(x)
}


# ============================================================================
# 4. .parse_lhs — Parse left-hand side
# ============================================================================

#' Parse the left-hand side of an hbb formula
#'
#' Extracts the response variable name and the trials variable name
#' from the LHS. Expected form: `y | trials(n)`.
#'
#' @param formula A two-sided formula.
#' @return A named list with `response` (character) and `trials` (character).
#' @noRd
.parse_lhs <- function(formula) {

    lhs_expr <- formula[[2L]]
    lhs_text <- deparse(lhs_expr, width.cutoff = 500L)
    # Collapse multi-line deparse output
    lhs_text <- paste(lhs_text, collapse = " ")
    # Remove whitespace for consistent matching
    lhs_clean <- gsub("\\s+", "", lhs_text)

    # Match: response | trials(variable)
    # The response can be any valid R name (letters, digits, dots, underscores)
    # The trials variable similarly
    pattern <- "^([a-zA-Z._][a-zA-Z0-9._]*)\\|trials\\(([a-zA-Z._][a-zA-Z0-9._]*)\\)$"
    m <- regmatches(lhs_clean, regexec(pattern, lhs_clean))[[1L]]

    if (length(m) != 3L) {
        cli_abort(c(
            "Invalid left-hand side.",
            "x" = "Got: {.code {lhs_text}}.",
            "i" = "Expected: {.code y | trials(n)}.",
            "i" = "The LHS must specify a response variable and a trials
                   variable using the {.code | trials()} syntax."
        ))
    }

    list(
        response = m[2L],
        trials   = m[3L]
    )
}


# ============================================================================
# 5. .parse_rhs — Parse right-hand side (orchestrator)
# ============================================================================

#' Parse the right-hand side of an hbb formula
#'
#' Orchestrates extraction of state_level(), random effects, and fixed
#' effects from the RHS of the formula.
#'
#' @param formula A two-sided formula.
#' @return A named list with `fixed`, `group`, `random`, `policy`.
#' @noRd
.parse_rhs <- function(formula) {

    # Get the RHS as a character string
    rhs_expr <- formula[[3L]]
    rhs_text <- deparse(rhs_expr, width.cutoff = 500L)
    rhs_text <- paste(rhs_text, collapse = " ")

    # -- Reject dot expansion -------------------------------------------------
    # Check for standalone "." (not inside a variable name)
    rhs_clean_check <- gsub("\\s+", " ", trimws(rhs_text))
    # Split into tokens by + and -, check if any token is exactly "."
    tokens_check <- trimws(strsplit(rhs_clean_check, "[+]")[[1L]])
    if ("." %in% tokens_check) {
        cli_abort(c(
            "Dot expansion ({.code ~ .}) is not supported.",
            "i" = "List predictor variables explicitly:
                   {.code y | trials(n) ~ x1 + x2}."
        ))
    }

    # -- Reject explicit intercept removal (- 1 or + 0) ----------------------
    # The intercept is always included in hbb models.
    # Check for standalone "0" as a term and "- 1" independently to avoid
    # operator-precedence bugs when RE terms containing "|" are present.
    tokens_for_zero <- trimws(strsplit(rhs_clean_check, "[+]")[[1L]])
    has_zero <- "0" %in% tokens_for_zero
    has_minus_one <- any(grepl("^-\\s*1$", tokens_for_zero)) ||
        grepl("-\\s*1\\s*$", rhs_clean_check) ||
        grepl("-\\s*1\\s*\\+", rhs_clean_check)

    if (has_minus_one || has_zero) {
        cli_abort(c(
            "Removing the intercept is not supported.",
            "x" = "Found {.code {if (has_minus_one) '- 1' else '0'}} in
                   the formula.",
            "i" = "The intercept is always automatically included in
                   hurdle Beta-Binomial models."
        ))
    }

    # -- Step 1: Extract state_level() terms ----------------------------------
    sl_result <- .parse_state_level(rhs_text)
    rhs_text <- sl_result$rhs_remaining

    # -- Step 2: Extract random effects (... | group) -------------------------
    re_result <- .parse_random_effects(rhs_text)
    rhs_text <- re_result$rhs_remaining

    # -- Step 3: Parse remaining fixed effects --------------------------------
    fixed <- .parse_fixed_effects(rhs_text)

    list(
        fixed  = fixed,
        group  = re_result$group,
        random = re_result$random,
        policy = sl_result$policy
    )
}


# ============================================================================
# 6. .parse_state_level — Extract state_level() terms
# ============================================================================

#' Extract state_level() terms from the RHS string
#'
#' Finds `state_level(v1 + v2 + ...)`, parses the variable names,
#' and removes the term from the RHS string.
#'
#' @param rhs_text Character. The RHS of the formula as a string.
#' @return A list with `policy` (character vector or NULL) and
#'   `rhs_remaining` (character string with state_level removed).
#' @noRd
.parse_state_level <- function(rhs_text) {

    # Pattern to match state_level(...)
    # Use a greedy match for the contents inside parentheses
    pattern <- "state_level\\(([^)]+)\\)"
    m <- regmatches(rhs_text, regexec(pattern, rhs_text))[[1L]]

    if (length(m) == 0L) {
        return(list(policy = NULL, rhs_remaining = rhs_text))
    }

    # Check for multiple state_level() calls
    all_matches <- gregexpr(pattern, rhs_text)[[1L]]
    if (length(all_matches) > 1L) {
        cli_abort(c(
            "Multiple {.code state_level()} terms are not allowed.",
            "i" = "Combine all policy variables into a single
                   {.code state_level()} call:
                   {.code state_level(v1 + v2 + v3)}."
        ))
    }

    # Parse the inner content: "v1 + v2 + v3"
    inner <- m[2L]
    inner_clean <- gsub("\\s+", "", inner)
    policy_vars <- strsplit(inner_clean, "\\+")[[1L]]

    # Remove empty strings from leading/trailing +
    policy_vars <- policy_vars[nzchar(policy_vars)]

    if (length(policy_vars) == 0L) {
        cli_abort(c(
            "{.code state_level()} must contain at least one variable.",
            "x" = "Got empty {.code state_level()}.",
            "i" = "Example: {.code state_level(mr_pctile + tiered_reim)}."
        ))
    }

    # Validate variable names
    bad_names <- policy_vars[!grepl("^[a-zA-Z._][a-zA-Z0-9._]*$", policy_vars)]
    if (length(bad_names) > 0L) {
        cli_abort(c(
            "{qty(length(bad_names))}Invalid variable name{?s} in {.code state_level()}.",
            "x" = "Invalid: {.val {bad_names}}.",
            "i" = "Variable names must be valid R identifiers."
        ))
    }

    # Check for duplicates
    if (anyDuplicated(policy_vars)) {
        dups <- policy_vars[duplicated(policy_vars)]
        cli_abort(c(
            "{qty(length(unique(dups)))}Duplicate variable{?s} in {.code state_level()}.",
            "x" = "Duplicated: {.val {unique(dups)}}."
        ))
    }

    # Reject "1" (intercept) in state_level — it is always implicit
    if ("1" %in% policy_vars) {
        cli_abort(c(
            "Do not include {.code 1} in {.code state_level()}.",
            "i" = "The intercept is always automatically included in the
                   V matrix."
        ))
    }

    # Remove state_level(...) from rhs_text
    rhs_remaining <- sub(pattern, "", rhs_text)
    # Clean up orphaned + signs: e.g., "x1 +  + x2" or leading/trailing +
    rhs_remaining <- gsub("\\+\\s*\\+", "+", rhs_remaining)
    rhs_remaining <- gsub("^\\s*\\+", "", rhs_remaining)
    rhs_remaining <- gsub("\\+\\s*$", "", rhs_remaining)
    rhs_remaining <- trimws(rhs_remaining)

    list(
        policy        = policy_vars,
        rhs_remaining = rhs_remaining
    )
}


# ============================================================================
# 7. .parse_random_effects — Extract (... | group) terms
# ============================================================================

#' Extract random effects terms from the RHS string
#'
#' Finds `(terms | group)` expressions, extracts the random effect
#' variables and grouping factor, then removes the term from the
#' RHS string.
#'
#' @param rhs_text Character. The RHS of the formula as a string
#'   (with state_level already removed).
#' @return A list with `group` (character or NULL), `random` (character
#'   vector or NULL), and `rhs_remaining` (character string).
#' @noRd
.parse_random_effects <- function(rhs_text) {

    # Pattern: match outermost parentheses containing a pipe
    # This finds (stuff | group) allowing nested content
    # Use a non-greedy match inside parens, then the bar, then the group
    pattern <- "\\(([^)]+)\\)"

    # Find all parenthesized expressions
    all_matches <- gregexpr(pattern, rhs_text, perl = TRUE)[[1L]]
    if (all_matches[1L] == -1L) {
        return(list(group = NULL, random = NULL, rhs_remaining = rhs_text))
    }

    # Extract each match and filter for those containing "|"
    match_texts <- regmatches(rhs_text, gregexpr(pattern, rhs_text, perl = TRUE))[[1L]]
    bar_matches <- match_texts[grepl("\\|", match_texts)]

    if (length(bar_matches) == 0L) {
        return(list(group = NULL, random = NULL, rhs_remaining = rhs_text))
    }

    if (length(bar_matches) > 1L) {
        cli_abort(c(
            "Multiple random effects terms are not allowed.",
            "x" = "Found {length(bar_matches)} {.code (... | group)} terms.",
            "i" = "Specify all random slopes in a single term:
                   {.code (x1 + x2 | state_id)}."
        ))
    }

    # Parse the single random effects term
    parsed <- .parse_bar_terms(bar_matches[1L])

    # Remove the random effects term from rhs_text
    # Escape the match for use in sub()
    escaped <- gsub("([\\(\\)\\|\\+\\.])", "\\\\\\1", bar_matches[1L])
    rhs_remaining <- sub(escaped, "", rhs_text)
    # Clean up orphaned + signs
    rhs_remaining <- gsub("\\+\\s*\\+", "+", rhs_remaining)
    rhs_remaining <- gsub("^\\s*\\+", "", rhs_remaining)
    rhs_remaining <- gsub("\\+\\s*$", "", rhs_remaining)
    rhs_remaining <- trimws(rhs_remaining)

    list(
        group         = parsed$group,
        random        = parsed$random,
        rhs_remaining = rhs_remaining
    )
}


# ============================================================================
# 8. .parse_fixed_effects — Extract remaining fixed effect terms
# ============================================================================

#' Parse fixed effect terms from the remaining RHS string
#'
#' After state_level() and random effects have been removed, this
#' function extracts the fixed effect variable names from the
#' remaining string.
#'
#' @param rhs_text Character. The RHS string with state_level and
#'   random effects already removed.
#' @return Character vector of fixed effect term names (excluding
#'   the intercept). May be `character(0)` for intercept-only models.
#' @noRd
.parse_fixed_effects <- function(rhs_text) {

    rhs_text <- trimws(rhs_text)

    # Handle empty RHS (after removing random effects / state_level)
    if (nchar(rhs_text) == 0L || rhs_text == "1") {
        return(character(0L))
    }

    # Split by +
    terms <- trimws(strsplit(rhs_text, "\\+")[[1L]])
    terms <- terms[nzchar(terms)]

    # Remove the intercept marker "1" if present (it's always implicit)
    terms <- terms[terms != "1"]

    if (length(terms) == 0L) {
        return(character(0L))
    }

    # Reject interaction terms (we don't support them at the formula level)
    has_interaction <- grepl("[*:]", terms)
    if (any(has_interaction)) {
        cli_abort(c(
            "Interaction terms are not supported in the formula.",
            "x" = "Found: {.val {terms[has_interaction]}}.",
            "i" = "Create interaction variables manually and include them
                   as regular predictors."
        ))
    }

    # Reject transformation functions (e.g., log(x), poly(x, 2))
    has_func <- grepl("\\(", terms)
    if (any(has_func)) {
        cli_abort(c(
            "Function calls in formula terms are not supported
             (except {.code state_level()}).",
            "x" = "Found: {.val {terms[has_func]}}.",
            "i" = "Apply transformations to your data before specifying the
                   formula."
        ))
    }

    # Validate variable names
    bad_names <- terms[!grepl("^[a-zA-Z._][a-zA-Z0-9._]*$", terms)]
    if (length(bad_names) > 0L) {
        cli_abort(c(
            "{qty(length(bad_names))}Invalid variable name{?s} in fixed effects.",
            "x" = "Invalid: {.val {bad_names}}.",
            "i" = "Variable names must be valid R identifiers."
        ))
    }

    # Check for duplicates
    if (anyDuplicated(terms)) {
        dups <- terms[duplicated(terms)]
        cli_abort(c(
            "{qty(length(unique(dups)))}Duplicate fixed effect variable{?s}.",
            "x" = "Duplicated: {.val {unique(dups)}}.",
            "i" = "Each variable should appear only once in the formula."
        ))
    }

    terms
}


# ============================================================================
# 9. .parse_bar_terms — Parse a single (... | group) expression
# ============================================================================

#' Parse a single random effects expression
#'
#' Takes a string like `"(x1 + x2 | group)"` and extracts the random
#' effect terms and grouping variable.
#'
#' @param term_text Character. A parenthesized bar expression, e.g.,
#'   `"(1 | state_id)"` or `"(x1 + x2 | state_id)"`.
#' @return A list with `group` (character) and `random` (character vector:
#'   `"1"` for intercept-only, or variable names for SVC).
#' @noRd
.parse_bar_terms <- function(term_text) {

    # Strip outer parentheses
    inner <- gsub("^\\(|\\)$", "", trimws(term_text))
    inner <- trimws(inner)

    # Split on "|"
    parts <- strsplit(inner, "\\|")[[1L]]

    if (length(parts) != 2L) {
        cli_abort(c(
            "Invalid random effects specification.",
            "x" = "Got: {.code {term_text}}.",
            "i" = "Expected: {.code (terms | group)}, e.g.,
                   {.code (1 | state_id)} or
                   {.code (x1 + x2 | state_id)}."
        ))
    }

    lhs <- trimws(parts[1L])
    rhs <- trimws(parts[2L])

    # -- Parse grouping variable (RHS of |) -----------------------------------
    group <- gsub("\\s+", "", rhs)

    if (!grepl("^[a-zA-Z._][a-zA-Z0-9._]*$", group)) {
        cli_abort(c(
            "Invalid grouping variable name.",
            "x" = "Got: {.val {group}}.",
            "i" = "The grouping variable must be a valid R identifier,
                   e.g., {.code state_id}."
        ))
    }

    # -- Parse random effect terms (LHS of |) ---------------------------------
    re_terms <- trimws(strsplit(lhs, "\\+")[[1L]])
    re_terms <- re_terms[nzchar(re_terms)]

    if (length(re_terms) == 0L) {
        cli_abort(c(
            "Empty left-hand side in random effects specification.",
            "x" = "Got: {.code {term_text}}.",
            "i" = "Specify at least {.code (1 | group)} for a random
                   intercept."
        ))
    }

    # Check if this is intercept-only: (1 | group)
    if (length(re_terms) == 1L && re_terms == "1") {
        return(list(group = group, random = "1"))
    }

    # For SVC: remove "1" if present (intercept is always implicit in RE)
    re_terms <- re_terms[re_terms != "1"]

    if (length(re_terms) == 0L) {
        # Edge case: user wrote (1 | group) but we already handled that above
        return(list(group = group, random = "1"))
    }

    # Validate variable names
    bad_names <- re_terms[!grepl("^[a-zA-Z._][a-zA-Z0-9._]*$", re_terms)]
    if (length(bad_names) > 0L) {
        cli_abort(c(
            "{qty(length(bad_names))}Invalid variable name{?s} in random effects.",
            "x" = "Invalid: {.val {bad_names}}.",
            "i" = "Variable names must be valid R identifiers."
        ))
    }

    # Check for duplicates
    if (anyDuplicated(re_terms)) {
        dups <- re_terms[duplicated(re_terms)]
        cli_abort(c(
            "{qty(length(unique(dups)))}Duplicate variable{?s} in random effects.",
            "x" = "Duplicated: {.val {unique(dups)}}."
        ))
    }

    list(group = group, random = re_terms)
}
