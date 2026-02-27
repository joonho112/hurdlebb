# Test if an Object is an hbb_fit

Returns `TRUE` if `x` inherits from class `"hbb_fit"`, `FALSE`
otherwise. This is the recommended way to check whether an object was
produced by
[`hbb()`](https://joonho112.github.io/hurdlebb/reference/hbb.md).

## Usage

``` r
is.hbb_fit(x)
```

## Arguments

- x:

  Any R object.

## Value

Logical scalar.

## See also

[`hbb()`](https://joonho112.github.io/hurdlebb/reference/hbb.md) for
fitting models.

Other fitting:
[`hbb()`](https://joonho112.github.io/hurdlebb/reference/hbb.md),
[`print.hbb_fit()`](https://joonho112.github.io/hurdlebb/reference/print.hbb_fit.md)

## Examples

``` r
is.hbb_fit(1)          # FALSE
#> [1] FALSE
is.hbb_fit(list())     # FALSE
#> [1] FALSE

if (FALSE) { # \dontrun{
fit <- hbb(y | trials(n_trial) ~ poverty, data = my_data)
is.hbb_fit(fit)        # TRUE
} # }
```
