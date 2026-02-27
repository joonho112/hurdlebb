# Validate an hbb_formula Against a Dataset

Checks that all variables referenced in the formula exist in the
supplied data frame(s) and that they have appropriate types and ranges.
This is typically called internally by the fitting function, but can be
used directly for early error checking.

## Usage

``` r
validate_hbb_formula(object, data, state_data = NULL)
```

## Arguments

- object:

  An object of class `"hbb_formula"`, as returned by
  [`hbb_formula()`](https://joonho112.github.io/hurdlebb/reference/hbb_formula.md).

- data:

  A data frame containing provider-level variables (response, trials,
  fixed effects, grouping variable).

- state_data:

  An optional data frame containing state-level policy variables.
  Required if `object$policy` is not `NULL`. Must contain a column
  matching `object$group` for merging.

## Value

Invisibly returns `TRUE` if all checks pass. Throws an informative error
otherwise.

## Details

The following checks are performed:

1.  `object` is an `hbb_formula`.

2.  `data` is a data frame with at least one row.

3.  The response variable exists and contains non-negative integers.

4.  The trials variable exists and contains positive integers.

5.  Response does not exceed trials for any observation.

6.  All fixed effect variables exist and are numeric.

7.  The grouping variable (if any) exists and is coercible to integer.

8.  If `state_data` is supplied: all policy variables exist in
    `state_data`, the grouping variable exists in `state_data`, and
    policy variables are numeric.

## See also

[`hbb_formula()`](https://joonho112.github.io/hurdlebb/reference/hbb_formula.md)

Other formula:
[`hbb_formula()`](https://joonho112.github.io/hurdlebb/reference/hbb_formula.md)

## Examples

``` r
f <- hbb_formula(y | trials(n_trial) ~ poverty + urban)
data(nsece_synth_small)
validate_hbb_formula(f, nsece_synth_small)
```
