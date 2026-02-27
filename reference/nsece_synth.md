# Synthetic NSECE 2019 Center-Based Provider Data

A synthetic dataset mimicking the structure and statistical properties
of the 2019 National Survey of Early Care and Education (NSECE)
center-based provider data. Generated via Gaussian copula to preserve
marginal distributions and rank-correlation structure of the original
restricted-use data.

## Usage

``` r
nsece_synth
```

## Format

A data frame with approximately 6,785 rows and 13 columns:

- provider_id:

  Integer. Unique provider identifier (1 to N).

- state_id:

  Integer. State identifier (1 to 51, including DC).

- y:

  Integer. Infant/toddler (IT) enrollment count. Zero for centers not
  serving IT (`z = 0`).

- n_trial:

  Integer. Total enrollment for children ages 0–5. Always \\\ge 1\\.

- z:

  Integer. Participation indicator: 1 if the center serves any IT
  children, 0 otherwise. Equal to `I(y > 0)`.

- it_share:

  Numeric. IT enrollment share `y / n_trial`. Zero when `z = 0`.

- poverty:

  Numeric. Community poverty rate (percent below poverty line). Range
  approximately 2–54.

- urban:

  Numeric. Community urbanization rate (percent urban population). Range
  0–100, heavily right-skewed.

- black:

  Numeric. Community percent Black population. Range 0–97.

- hispanic:

  Numeric. Community percent Hispanic population. Range 0–98.

- weight:

  Numeric. Survey sampling weight. Positive values; larger weights
  indicate the provider represents more units in the population.

- stratum:

  Integer. Sampling stratum identifier.

- psu:

  Integer. Primary sampling unit (PSU) identifier, nested within strata.

## Source

Generated from the NSECE 2019 restricted-use data via Gaussian copula
(Miratrix, 2025) with hurdle Beta-Binomial outcomes calibrated to the
empirical moments of the original data. See
`data-raw/generate_synthetic.R` for the full generation script.

## Details

The data-generating process uses a two-part hurdle Beta-Binomial model:

- Part 1 (Extensive margin):

  Whether a center serves IT children at all, modeled as \\z_i \sim
  \textrm{Bernoulli}(q_i)\\ where \\\textrm{logit}(q_i) = X_i \alpha +
  \delta\_{1,s_i}\\.

- Part 2 (Intensive margin):

  Among servers, the IT enrollment count follows a zero-truncated
  Beta-Binomial: \\y_i \mid z_i = 1 \sim \textrm{ZT-BetaBin}(n_i, \mu_i,
  \kappa)\\ where \\\textrm{logit}(\mu_i) = X_i \beta +
  \delta\_{2,s_i}\\.

The structural zero rate is approximately 35 percent, and the mean IT
share among servers is approximately 48 percent. These match the
empirical moments of the original NSECE 2019 data.

## References

National Survey of Early Care and Education (NSECE), 2019. U. S.
Department of Health and Human Services.

## See also

[`nsece_synth_small`](https://joonho112.github.io/hurdlebb/reference/nsece_synth_small.md)
for a smaller version suitable for quick examples and tests,
[`nsece_state_policy`](https://joonho112.github.io/hurdlebb/reference/nsece_state_policy.md)
for state-level policy variables.

## Examples

``` r
data(nsece_synth)
head(nsece_synth)
#>   provider_id state_id  y n_trial z   it_share  poverty    urban      black
#> 1           1       33  1      24 1 0.04166667 29.70972 99.95720  4.7455889
#> 2           2       37 30      58 1 0.51724138 19.30702 99.96461  4.5828925
#> 3           3       26  0     101 0 0.00000000 12.71033 99.92934 18.5843027
#> 4           4        5  0      19 0 0.00000000 18.54633 15.12154  0.8891587
#> 5           5        5 31     125 1 0.24800000 10.63834 99.94319  1.5584889
#> 6           6       10  0     172 0 0.00000000 19.32121 99.96716 32.6650238
#>    hispanic    weight stratum psu
#> 1 60.887394  3.276547       3 138
#> 2 42.464795  6.939392      27 402
#> 3  6.466195  4.953660      20 367
#> 4  9.839974 16.319319       1  42
#> 5  9.304892  5.363512       1  45
#> 6 20.409599  2.164282       5 202

# Basic summary
cat("N =", nrow(nsece_synth), "\n")
#> N = 6785 
cat("Zero rate:", round(1 - mean(nsece_synth$z), 3), "\n")
#> Zero rate: 0.354 
cat("IT share (servers):",
    round(mean(nsece_synth$it_share[nsece_synth$z == 1]), 3), "\n")
#> IT share (servers): 0.496 

# State distribution
cat("States:", length(unique(nsece_synth$state_id)), "\n")
#> States: 51 
```
