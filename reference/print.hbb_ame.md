# Print Method for hbb_ame Objects

Displays a structured summary of the AME decomposition, including the
decomposition table with sign patterns, reversal probabilities, and
optional Wald comparison.

## Usage

``` r
# S3 method for class 'hbb_ame'
print(x, digits = 4, ...)
```

## Arguments

- x:

  An object of class `"hbb_ame"` returned by
  [`ame`](https://joonho112.github.io/hurdlebb/reference/ame.md).

- digits:

  Integer; number of significant digits for numeric output. Default is
  4.

- ...:

  Additional arguments (currently unused).

## Value

The object `x`, invisibly.

## See also

[`ame`](https://joonho112.github.io/hurdlebb/reference/ame.md),
[`ame_decomposition`](https://joonho112.github.io/hurdlebb/reference/ame_decomposition.md)
