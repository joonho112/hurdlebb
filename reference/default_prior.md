# Default Prior Specification for Hurdle Beta-Binomial Models

Returns the package default priors, which are weakly informative and
match the Stan model code:

## Usage

``` r
default_prior()
```

## Value

An S3 object of class `"hbb_prior"`. See
[`hbb_prior()`](https://joonho112.github.io/hurdlebb/reference/hbb_prior.md)
for the full structure.

## Details

|               |                  |
|---------------|------------------|
| **Parameter** | **Prior**        |
| `alpha`       | Normal(0, 2)     |
| `beta`        | Normal(0, 2)     |
| `log_kappa`   | Normal(2, 1.5)   |
| `gamma`       | Normal(0, 1)     |
| `tau`         | HalfNormal(0, 1) |
| `L_Omega`     | LKJ(2)           |

## See also

[`hbb_prior()`](https://joonho112.github.io/hurdlebb/reference/hbb_prior.md)
for custom prior specification.

Other priors:
[`hbb_prior()`](https://joonho112.github.io/hurdlebb/reference/hbb_prior.md)

## Examples

``` r
dp <- default_prior()
dp
#> Hurdle Beta-Binomial Prior Specification
#> ----------------------------------------
#>   alpha     ~ Normal(0, 2) 
#>   beta      ~ Normal(0, 2) 
#>   log_kappa ~ Normal(2, 1.5) 
#>   gamma     ~ Normal(0, 1) 
#>   tau       ~ HalfNormal(0, 1) 
#>   L_Omega   ~ LKJ(2) 

# Inspect a specific component
dp$alpha
#> $dist
#> [1] "normal"
#> 
#> $mean
#> [1] 0
#> 
#> $sd
#> [1] 2
#> 
dp$lkj_eta
#> [1] 2
```
