# Build the AME decomposition table

Build the AME decomposition table

## Usage

``` r
.ame_build_decomp(ext_summ, int_summ, total_summ, labels, idx)
```

## Arguments

- ext_summ:

  Data frame: extensive summary from .ame_summarize.

- int_summ:

  Data frame: intensive summary from .ame_summarize.

- total_summ:

  Data frame: total summary from .ame_summarize.

- labels:

  Character vector of length P.

- idx:

  Integer vector: indices of non-intercept covariates.

## Value

Data frame with 13 columns.
