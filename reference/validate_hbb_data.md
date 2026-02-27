# Validate an hbb_data Object

Performs comprehensive checks on a Stan data list to ensure internal
consistency before passing to CmdStan. Verifies field existence,
dimensional consistency, value ranges, and cross-field coherence.

## Usage

``` r
validate_hbb_data(data_list)
```

## Arguments

- data_list:

  An object of class `"hbb_data"`, as returned by
  [`prepare_stan_data()`](https://joonho112.github.io/hurdlebb/reference/prepare_stan_data.md).

## Value

Invisibly returns `TRUE` if all checks pass. Throws an informative error
otherwise.

## See also

[`prepare_stan_data()`](https://joonho112.github.io/hurdlebb/reference/prepare_stan_data.md)

Other data:
[`prepare_stan_data()`](https://joonho112.github.io/hurdlebb/reference/prepare_stan_data.md)

## Examples

``` r
data(nsece_synth_small, package = "hurdlebb")
f <- hbb_formula(y | trials(n_trial) ~ poverty + urban)
d <- prepare_stan_data(f, nsece_synth_small)
validate_hbb_data(d)
```
