# Small Synthetic NSECE 2019 Dataset

A stratified subsample of
[`nsece_synth`](https://joonho112.github.io/hurdlebb/reference/nsece_synth.md)
with approximately 500 rows. Designed for quick demonstrations, unit
tests, and package examples where computational speed matters more than
statistical power.

## Usage

``` r
nsece_synth_small
```

## Format

A data frame with approximately 500 rows and the same 13 columns as
[`nsece_synth`](https://joonho112.github.io/hurdlebb/reference/nsece_synth.md):

- provider_id:

  Integer. Unique provider identifier (re-indexed 1 to N for this
  subset).

- state_id:

  Integer. State identifier (1 to 51, including DC).

- y:

  Integer. Infant/toddler (IT) enrollment count.

- n_trial:

  Integer. Total enrollment for children ages 0–5.

- z:

  Integer. Participation indicator: 1 if `y > 0`, 0 otherwise.

- it_share:

  Numeric. IT enrollment share `y / n_trial`.

- poverty:

  Numeric. Community poverty rate.

- urban:

  Numeric. Community urbanization rate.

- black:

  Numeric. Community percent Black population.

- hispanic:

  Numeric. Community percent Hispanic population.

- weight:

  Numeric. Survey sampling weight.

- stratum:

  Integer. Sampling stratum identifier.

- psu:

  Integer. Primary sampling unit (PSU) identifier.

## Source

Stratified subsample of
[`nsece_synth`](https://joonho112.github.io/hurdlebb/reference/nsece_synth.md).
See `data-raw/generate_synthetic.R` for the generation script.

## Details

The subsample is drawn via stratified sampling on `state_id` to preserve
the state size distribution of the full dataset. All 51 state
identifiers are represented. The structural zero rate and covariate
distributions approximately match the full dataset, though with more
sampling variability due to the smaller size.

This dataset is **not suitable for model fitting**. With only \\\sim
10\\ observations per state on average, state-varying coefficient models
will not converge reliably. Use
[`nsece_synth`](https://joonho112.github.io/hurdlebb/reference/nsece_synth.md)
for any serious modeling.

## See also

[`nsece_synth`](https://joonho112.github.io/hurdlebb/reference/nsece_synth.md)
for the full dataset,
[`nsece_state_policy`](https://joonho112.github.io/hurdlebb/reference/nsece_state_policy.md)
for state-level policy variables.

## Examples

``` r
data(nsece_synth_small)

# Quick overview
cat("N =", nrow(nsece_synth_small), "\n")
#> N = 504 
cat("States:", length(unique(nsece_synth_small$state_id)), "\n")
#> States: 51 
cat("Zero rate:", round(1 - mean(nsece_synth_small$z), 3), "\n")
#> Zero rate: 0.379 

# Useful for quick tests and demonstrations
table(nsece_synth_small$z)
#> 
#>   0   1 
#> 191 313 
```
