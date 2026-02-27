# State-Level Childcare Policy Variables

State-level policy indicators for all 51 jurisdictions (50 states + DC).
Used as cross-level moderators in the hurdle Beta-Binomial model with
state-varying coefficients (Module D).

## Usage

``` r
nsece_state_policy
```

## Format

A data frame with 51 rows and 5 columns:

- state_id:

  Integer. State identifier matching `nsece_synth$state_id` (1 to 51).

- state_name:

  Character. Generic state label (`"State_01"` through `"State_51"`).
  Names are anonymized to prevent identification of actual states.

- mr_pctile:

  Numeric. CCDF market rate percentile (standardized, mean \\\approx
  0\\, SD \\\approx 1\\). Higher values indicate more generous
  reimbursement rates.

- tiered_reim:

  Integer. 1 if the state uses tiered reimbursement
  (quality-differentiated CCDF rates), 0 otherwise. Approximately 84
  percent of states have tiered reimbursement.

- it_addon:

  Integer. 1 if the state provides an IT add-on payment (supplemental
  rate for infant/toddler care), 0 otherwise. Approximately 24 percent
  of states have IT add-ons.

## Source

Synthetic values calibrated to the distribution of actual U.S. childcare
subsidy policies. See the CCDF Policies Database maintained by the Urban
Institute for original policy data. Generated in
`data-raw/generate_synthetic.R`.

## Details

These policy variables serve as cross-level moderators in the
state-varying coefficient (SVC) model. They enter the model as
predictors of the state random effects: \$\$\delta_s = \Gamma v_s +
\eta_s\$\$ where \\v_s = (1, \texttt{mr\\pctile}\_s,
\texttt{tiered\\reim}\_s, \texttt{it\\addon}\_s)'\\ and \\\Gamma\\
captures the cross-level interaction effects.

The `mr_pctile` variable is standardized (centered and scaled) to
improve MCMC sampling and interpretability. The binary indicators
`tiered_reim` and `it_addon` are left uncentered.

## References

National Survey of Early Care and Education (NSECE), 2019. U.S.
Department of Health and Human Services.

## See also

[`nsece_synth`](https://joonho112.github.io/hurdlebb/reference/nsece_synth.md)
for the provider-level data that can be merged with this table via
`state_id`.

## Examples

``` r
data(nsece_state_policy)
head(nsece_state_policy)
#>   state_id state_name  mr_pctile tiered_reim it_addon
#> 1        1   State_01  0.2151828           1        0
#> 2        2   State_02 -1.2477813           0        0
#> 3        3   State_03 -1.3731782           1        0
#> 4        4   State_04  1.6363478           1        1
#> 5        5   State_05  0.5913735           0        0
#> 6        6   State_06 -1.2895803           1        1

# Policy summary
cat("Tiered reimbursement:",
    round(mean(nsece_state_policy$tiered_reim), 2), "\n")
#> Tiered reimbursement: 0.84 
cat("IT add-on:",
    round(mean(nsece_state_policy$it_addon), 2), "\n")
#> IT add-on: 0.24 

# Merge with provider data
data(nsece_synth)
merged <- merge(nsece_synth, nsece_state_policy, by = "state_id")
cat("Merged rows:", nrow(merged), "\n")
#> Merged rows: 6785 
```
