
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mPower: A Real Data-based Power Analysis Tool for Microbiome Study Design

This function estimates the power for microbiome study design based on
various experimental designs and parameters. For details of the proposed
method, please refer to the following paper:

Lu Yang, Jun Chen. mPower: A Real Data-based Power Analysis Tool for
Microbiome Study Design

## Get the overview of mPower: ![This is a caption for the figure](workflows.png)

## Use our Shiny App without coding skills: <https://microbiomestat.shinyapps.io/mPower/>

## Installation

You can install mPower as follows:

``` r
# install.packages("devtools")
devtools::install_github("chloelulu/mPower")
```

## Example

1.  Assume you are designing a case-control study (i.e., cancer vs
    healthy individuals), you are not

``` r
library(mPower)
## Estimate taxa-level power for case control study for 50 samples, 10% differential taxa
data(feature.dat)
## Exclude feature exist in less than 2 samples
feature.dat <- feature.dat[rowSums(feature.dat != 0) > 2, ]
## Estimate the parameters
model.paras <- EstPara(ref.otu.tab = feature.dat)
#> Iteration 1: Log-likelihood value: -305650.388324419
#> Iteration 2: Log-likelihood value: -304380.969463483
#> Iteration 3: Log-likelihood value: -303439.452508666
#> Iteration 4: Log-likelihood value: -302892.343764286
#> Iteration 5: Log-likelihood value: -302636.990397356
#> Iteration 6: Log-likelihood value: -302552.4224922
#> Iteration 7: Log-likelihood value: -302537.930631092
#> Iteration 8: Log-likelihood value: -302537.268905297
#> Iteration 9: Log-likelihood value: -302537.266697766
#> Iteration 10: Log-likelihood value: -302537.266697725
## Estimate taxa-level power for case control study
res1 <- mPower(feature.dat = feature.dat, model.paras = model.paras,
               test = 'Community', design = 'CaseControl',
               nSams = 50, grp.ratio = 0.5,
               iters = 10, alpha = 0.05, distance = 'BC',
               diff.otu.pct = 0.1, diff.otu.direct = 'balanced',diff.otu.mode = 'random',
               covariate.eff.min = 0, covariate.eff.maxs = c(1, 2, 3),
               confounder = 'no', depth.mu = 10000, depth.theta = 5)
#> ----Calculating in progressing----
#>  nSam= 50  & covariate.eff.max= 1 
#> mPower is calculating: 1 of 10 iterations.
#> mPower is calculating: 2 of 10 iterations.
#> mPower is calculating: 3 of 10 iterations.
#> mPower is calculating: 4 of 10 iterations.
#> mPower is calculating: 5 of 10 iterations.
#> mPower is calculating: 6 of 10 iterations.
#> mPower is calculating: 7 of 10 iterations.
#> mPower is calculating: 8 of 10 iterations.
#> mPower is calculating: 9 of 10 iterations.
#> mPower is calculating: 10 of 10 iterations.
#> 
#>  nSam= 50  & covariate.eff.max= 2 
#> mPower is calculating: 1 of 10 iterations.
#> mPower is calculating: 2 of 10 iterations.
#> mPower is calculating: 3 of 10 iterations.
#> mPower is calculating: 4 of 10 iterations.
#> mPower is calculating: 5 of 10 iterations.
#> mPower is calculating: 6 of 10 iterations.
#> mPower is calculating: 7 of 10 iterations.
#> mPower is calculating: 8 of 10 iterations.
#> mPower is calculating: 9 of 10 iterations.
#> mPower is calculating: 10 of 10 iterations.
#> 
#>  nSam= 50  & covariate.eff.max= 3 
#> mPower is calculating: 1 of 10 iterations.
#> mPower is calculating: 2 of 10 iterations.
#> mPower is calculating: 3 of 10 iterations.
#> mPower is calculating: 4 of 10 iterations.
#> mPower is calculating: 5 of 10 iterations.
#> mPower is calculating: 6 of 10 iterations.
#> mPower is calculating: 7 of 10 iterations.
#> mPower is calculating: 8 of 10 iterations.
#> mPower is calculating: 9 of 10 iterations.
#> mPower is calculating: 10 of 10 iterations.
res1$power
#>   nSam covariate.eff.max value median        SD      ymax      ymin        SE
#> 1   50                 1   0.2      0 0.4216370 0.4613333 0.0000000 0.1333333
#> 2   50                 2   0.6      1 0.5163978 0.9200667 0.2799333 0.1632993
#> 3   50                 3   0.9      1 0.3162278 1.0960000 0.7040000 0.1000000
#>   ErrRate variable
#> 1       0     beta
#> 2       0     beta
#> 3       0     beta
## Estimate community-level power for case control study
res2 <- mPower(feature.dat = feature.dat, model.paras = model.paras,
               test = 'Taxa', design = 'CaseControl',
               nSams = c(20, 40, 80), grp.ratio = 0.5,
               iters = 100, alpha = 0.05, distance = 'BC',
               diff.otu.pct = 0.1, diff.otu.direct = 'balanced',diff.otu.mode = 'random',
               covariate.eff.min = 0, covariate.eff.maxs = 2,
               prev.filter = 0.1, max.abund.filter = 0.002,
               confounder = 'yes', depth.mu = 100000, depth.theta = 5)
#> ----Calculating in progressing----
#>  nSam= 20  & covariate.eff.max= 2 
#> mPower is calculating: 1 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  196 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 2 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  193 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 3 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  195 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 4 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  190 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 5 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  184 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 6 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  203 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 7 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  195 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 8 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  186 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 9 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  199 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 10 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  178 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 11 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  192 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 12 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  173 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 13 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  191 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 14 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  189 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 15 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  200 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 16 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  198 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 17 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  196 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 18 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  211 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 19 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  206 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 20 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  194 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 21 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  201 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 22 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  194 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 23 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  196 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 24 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  183 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 25 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  183 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 26 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  178 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 27 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  200 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 28 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  187 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 29 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  185 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 30 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  192 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 31 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  201 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 32 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  177 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 33 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  191 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 34 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  187 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 35 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  173 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 36 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  203 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 37 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  182 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 38 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  184 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 39 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  203 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 40 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  200 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 41 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  193 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 42 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  205 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 43 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  191 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 44 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  195 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 45 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  204 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 46 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  157 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 47 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  172 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 48 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  185 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 49 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  174 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 50 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  167 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 51 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  207 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 52 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  181 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 53 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  209 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 54 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  172 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 55 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  196 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 56 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  217 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 57 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  188 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 58 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  171 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 59 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  170 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 60 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  195 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 61 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  167 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 62 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  196 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 63 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  191 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 64 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  194 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 65 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  191 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 66 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  193 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 67 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  190 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 68 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  189 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 69 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  203 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 70 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  192 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 71 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  196 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 72 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  190 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 73 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  177 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 74 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  186 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 75 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  203 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 76 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  169 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 77 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  172 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 78 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  180 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 79 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  200 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 80 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  203 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 81 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  187 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 82 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  193 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 83 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  177 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 84 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  178 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 85 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  192 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 86 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  179 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 87 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  200 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 88 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  181 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 89 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  181 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 90 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  194 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 91 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  185 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 92 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  188 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 93 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  192 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 94 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  191 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 95 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  163 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 96 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  191 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 97 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  180 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 98 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  197 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 99 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  159 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> mPower is calculating: 100 of 100 iterations.
#> Differential abundance analysis is testing:  20  samples and  194 out of 382  features!
#> Warning in linda(meta.dat = meta.dat, feature.dat = feature.dat, formula = formula, : Some features have less than 3 nonzero values!
#>                      They have virtually no statistical power. You may consider filtering them in the analysis!
#> 
#>  nSam= 40  & covariate.eff.max= 2 
#> mPower is calculating: 1 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  232 out of 382  features!
#> mPower is calculating: 2 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  207 out of 382  features!
#> mPower is calculating: 3 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  221 out of 382  features!
#> mPower is calculating: 4 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  238 out of 382  features!
#> mPower is calculating: 5 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  212 out of 382  features!
#> mPower is calculating: 6 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  233 out of 382  features!
#> mPower is calculating: 7 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  226 out of 382  features!
#> mPower is calculating: 8 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  230 out of 382  features!
#> mPower is calculating: 9 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  222 out of 382  features!
#> mPower is calculating: 10 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  217 out of 382  features!
#> mPower is calculating: 11 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  230 out of 382  features!
#> mPower is calculating: 12 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  228 out of 382  features!
#> mPower is calculating: 13 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  218 out of 382  features!
#> mPower is calculating: 14 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  223 out of 382  features!
#> mPower is calculating: 15 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  248 out of 382  features!
#> mPower is calculating: 16 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  218 out of 382  features!
#> mPower is calculating: 17 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  234 out of 382  features!
#> mPower is calculating: 18 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  217 out of 382  features!
#> mPower is calculating: 19 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  237 out of 382  features!
#> mPower is calculating: 20 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  217 out of 382  features!
#> mPower is calculating: 21 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  230 out of 382  features!
#> mPower is calculating: 22 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  223 out of 382  features!
#> mPower is calculating: 23 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  215 out of 382  features!
#> mPower is calculating: 24 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  211 out of 382  features!
#> mPower is calculating: 25 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  228 out of 382  features!
#> mPower is calculating: 26 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  227 out of 382  features!
#> mPower is calculating: 27 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  234 out of 382  features!
#> mPower is calculating: 28 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  223 out of 382  features!
#> mPower is calculating: 29 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  241 out of 382  features!
#> mPower is calculating: 30 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  211 out of 382  features!
#> mPower is calculating: 31 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  222 out of 382  features!
#> mPower is calculating: 32 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  224 out of 382  features!
#> mPower is calculating: 33 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  218 out of 382  features!
#> mPower is calculating: 34 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  209 out of 382  features!
#> mPower is calculating: 35 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  221 out of 382  features!
#> mPower is calculating: 36 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  225 out of 382  features!
#> mPower is calculating: 37 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  225 out of 382  features!
#> mPower is calculating: 38 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  229 out of 382  features!
#> mPower is calculating: 39 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  217 out of 382  features!
#> mPower is calculating: 40 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  221 out of 382  features!
#> mPower is calculating: 41 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  189 out of 382  features!
#> mPower is calculating: 42 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  196 out of 382  features!
#> mPower is calculating: 43 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  216 out of 382  features!
#> mPower is calculating: 44 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  228 out of 382  features!
#> mPower is calculating: 45 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  234 out of 382  features!
#> mPower is calculating: 46 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  231 out of 382  features!
#> mPower is calculating: 47 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  220 out of 382  features!
#> mPower is calculating: 48 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  211 out of 382  features!
#> mPower is calculating: 49 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  223 out of 382  features!
#> mPower is calculating: 50 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  225 out of 382  features!
#> mPower is calculating: 51 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  236 out of 382  features!
#> mPower is calculating: 52 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  214 out of 382  features!
#> mPower is calculating: 53 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  213 out of 382  features!
#> mPower is calculating: 54 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  231 out of 382  features!
#> mPower is calculating: 55 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  231 out of 382  features!
#> mPower is calculating: 56 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  234 out of 382  features!
#> mPower is calculating: 57 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  233 out of 382  features!
#> mPower is calculating: 58 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  237 out of 382  features!
#> mPower is calculating: 59 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  236 out of 382  features!
#> mPower is calculating: 60 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  223 out of 382  features!
#> mPower is calculating: 61 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  227 out of 382  features!
#> mPower is calculating: 62 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  234 out of 382  features!
#> mPower is calculating: 63 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  228 out of 382  features!
#> mPower is calculating: 64 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  229 out of 382  features!
#> mPower is calculating: 65 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  237 out of 382  features!
#> mPower is calculating: 66 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  213 out of 382  features!
#> mPower is calculating: 67 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  220 out of 382  features!
#> mPower is calculating: 68 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  227 out of 382  features!
#> mPower is calculating: 69 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  203 out of 382  features!
#> mPower is calculating: 70 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  205 out of 382  features!
#> mPower is calculating: 71 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  237 out of 382  features!
#> mPower is calculating: 72 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  220 out of 382  features!
#> mPower is calculating: 73 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  233 out of 382  features!
#> mPower is calculating: 74 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  201 out of 382  features!
#> mPower is calculating: 75 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  241 out of 382  features!
#> mPower is calculating: 76 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  200 out of 382  features!
#> mPower is calculating: 77 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  194 out of 382  features!
#> mPower is calculating: 78 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  225 out of 382  features!
#> mPower is calculating: 79 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  215 out of 382  features!
#> mPower is calculating: 80 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  229 out of 382  features!
#> mPower is calculating: 81 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  224 out of 382  features!
#> mPower is calculating: 82 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  213 out of 382  features!
#> mPower is calculating: 83 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  203 out of 382  features!
#> mPower is calculating: 84 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  213 out of 382  features!
#> mPower is calculating: 85 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  205 out of 382  features!
#> mPower is calculating: 86 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  228 out of 382  features!
#> mPower is calculating: 87 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  206 out of 382  features!
#> mPower is calculating: 88 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  242 out of 382  features!
#> mPower is calculating: 89 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  223 out of 382  features!
#> mPower is calculating: 90 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  219 out of 382  features!
#> mPower is calculating: 91 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  214 out of 382  features!
#> mPower is calculating: 92 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  229 out of 382  features!
#> mPower is calculating: 93 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  231 out of 382  features!
#> mPower is calculating: 94 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  227 out of 382  features!
#> mPower is calculating: 95 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  208 out of 382  features!
#> mPower is calculating: 96 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  228 out of 382  features!
#> mPower is calculating: 97 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  230 out of 382  features!
#> mPower is calculating: 98 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  225 out of 382  features!
#> mPower is calculating: 99 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  233 out of 382  features!
#> mPower is calculating: 100 of 100 iterations.
#> Differential abundance analysis is testing:  40  samples and  222 out of 382  features!
#> 
#>  nSam= 80  & covariate.eff.max= 2 
#> mPower is calculating: 1 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  231 out of 382  features!
#> mPower is calculating: 2 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  241 out of 382  features!
#> mPower is calculating: 3 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  250 out of 382  features!
#> mPower is calculating: 4 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  257 out of 382  features!
#> mPower is calculating: 5 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  247 out of 382  features!
#> mPower is calculating: 6 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  245 out of 382  features!
#> mPower is calculating: 7 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  239 out of 382  features!
#> mPower is calculating: 8 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  252 out of 382  features!
#> mPower is calculating: 9 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  249 out of 382  features!
#> mPower is calculating: 10 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  249 out of 382  features!
#> mPower is calculating: 11 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  255 out of 382  features!
#> mPower is calculating: 12 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  258 out of 382  features!
#> mPower is calculating: 13 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  247 out of 382  features!
#> mPower is calculating: 14 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  265 out of 382  features!
#> mPower is calculating: 15 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  254 out of 382  features!
#> mPower is calculating: 16 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  251 out of 382  features!
#> mPower is calculating: 17 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  267 out of 382  features!
#> mPower is calculating: 18 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  242 out of 382  features!
#> mPower is calculating: 19 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  252 out of 382  features!
#> mPower is calculating: 20 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  243 out of 382  features!
#> mPower is calculating: 21 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  246 out of 382  features!
#> mPower is calculating: 22 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  252 out of 382  features!
#> mPower is calculating: 23 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  252 out of 382  features!
#> mPower is calculating: 24 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  251 out of 382  features!
#> mPower is calculating: 25 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  255 out of 382  features!
#> mPower is calculating: 26 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  242 out of 382  features!
#> mPower is calculating: 27 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  246 out of 382  features!
#> mPower is calculating: 28 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  252 out of 382  features!
#> mPower is calculating: 29 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  240 out of 382  features!
#> mPower is calculating: 30 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  253 out of 382  features!
#> mPower is calculating: 31 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  244 out of 382  features!
#> mPower is calculating: 32 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  246 out of 382  features!
#> mPower is calculating: 33 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  251 out of 382  features!
#> mPower is calculating: 34 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  261 out of 382  features!
#> mPower is calculating: 35 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  233 out of 382  features!
#> mPower is calculating: 36 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  261 out of 382  features!
#> mPower is calculating: 37 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  240 out of 382  features!
#> mPower is calculating: 38 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  245 out of 382  features!
#> mPower is calculating: 39 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  249 out of 382  features!
#> mPower is calculating: 40 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  262 out of 382  features!
#> mPower is calculating: 41 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  243 out of 382  features!
#> mPower is calculating: 42 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  237 out of 382  features!
#> mPower is calculating: 43 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  235 out of 382  features!
#> mPower is calculating: 44 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  237 out of 382  features!
#> mPower is calculating: 45 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  253 out of 382  features!
#> mPower is calculating: 46 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  263 out of 382  features!
#> mPower is calculating: 47 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  250 out of 382  features!
#> mPower is calculating: 48 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  262 out of 382  features!
#> mPower is calculating: 49 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  248 out of 382  features!
#> mPower is calculating: 50 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  255 out of 382  features!
#> mPower is calculating: 51 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  251 out of 382  features!
#> mPower is calculating: 52 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  243 out of 382  features!
#> mPower is calculating: 53 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  254 out of 382  features!
#> mPower is calculating: 54 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  246 out of 382  features!
#> mPower is calculating: 55 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  250 out of 382  features!
#> mPower is calculating: 56 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  247 out of 382  features!
#> mPower is calculating: 57 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  250 out of 382  features!
#> mPower is calculating: 58 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  251 out of 382  features!
#> mPower is calculating: 59 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  241 out of 382  features!
#> mPower is calculating: 60 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  257 out of 382  features!
#> mPower is calculating: 61 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  248 out of 382  features!
#> mPower is calculating: 62 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  252 out of 382  features!
#> mPower is calculating: 63 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  262 out of 382  features!
#> mPower is calculating: 64 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  240 out of 382  features!
#> mPower is calculating: 65 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  255 out of 382  features!
#> mPower is calculating: 66 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  237 out of 382  features!
#> mPower is calculating: 67 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  253 out of 382  features!
#> mPower is calculating: 68 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  256 out of 382  features!
#> mPower is calculating: 69 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  256 out of 382  features!
#> mPower is calculating: 70 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  255 out of 382  features!
#> mPower is calculating: 71 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  243 out of 382  features!
#> mPower is calculating: 72 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  249 out of 382  features!
#> mPower is calculating: 73 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  254 out of 382  features!
#> mPower is calculating: 74 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  246 out of 382  features!
#> mPower is calculating: 75 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  249 out of 382  features!
#> mPower is calculating: 76 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  242 out of 382  features!
#> mPower is calculating: 77 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  239 out of 382  features!
#> mPower is calculating: 78 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  245 out of 382  features!
#> mPower is calculating: 79 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  264 out of 382  features!
#> mPower is calculating: 80 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  254 out of 382  features!
#> mPower is calculating: 81 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  244 out of 382  features!
#> mPower is calculating: 82 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  250 out of 382  features!
#> mPower is calculating: 83 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  236 out of 382  features!
#> mPower is calculating: 84 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  251 out of 382  features!
#> mPower is calculating: 85 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  244 out of 382  features!
#> mPower is calculating: 86 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  244 out of 382  features!
#> mPower is calculating: 87 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  240 out of 382  features!
#> mPower is calculating: 88 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  244 out of 382  features!
#> mPower is calculating: 89 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  250 out of 382  features!
#> mPower is calculating: 90 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  250 out of 382  features!
#> mPower is calculating: 91 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  248 out of 382  features!
#> mPower is calculating: 92 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  257 out of 382  features!
#> mPower is calculating: 93 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  253 out of 382  features!
#> mPower is calculating: 94 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  242 out of 382  features!
#> mPower is calculating: 95 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  240 out of 382  features!
#> mPower is calculating: 96 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  244 out of 382  features!
#> mPower is calculating: 97 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  253 out of 382  features!
#> mPower is calculating: 98 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  261 out of 382  features!
#> mPower is calculating: 99 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  255 out of 382  features!
#> mPower is calculating: 100 of 100 iterations.
#> Differential abundance analysis is testing:  80  samples and  236 out of 382  features!
res2$pOCR
#>   nSam covariate.eff.max value variable
#> 1   20                 2  0.30    LinDA
#> 2   40                 2  0.70    LinDA
#> 3   80                 2  0.89    LinDA
res2$aTPR
#>   variable nSam covariate.eff.max      value    median         SD       ymax
#> 1    LinDA   20                 2 0.01699349 0.0000000 0.02907630 0.02269245
#> 2    LinDA   40                 2 0.05685735 0.0400000 0.05732261 0.06809258
#> 3    LinDA   80                 2 0.10708756 0.1071429 0.07801662 0.12237882
#>         ymin          SE ErrRate
#> 1 0.01129454 0.002907630       0
#> 2 0.04562212 0.005732261       0
#> 3 0.09179630 0.007801662       0
```
