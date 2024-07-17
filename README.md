
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mPower: A Real Data-based Power Analysis Tool for Microbiome Study Design

This function estimates the power for microbiome study design based on
various experimental designs and parameters. For details of the proposed
method, please refer to the following paper:

Lu Yang, Jun Chen. mPower: A Real Data-based Power Analysis Tool for
Microbiome Study Design

![This is a caption for the figure](workflows.png)

## Usage

### 1. Local package installation

You can install mPower as follows:

``` r
# install.packages("devtools")
devtools::install_github("chloelulu/mPower")
```

### 2. Shiny R App

You can also use our Shiny App without coding skills, please access the
Shiny App at: <https://microbiomestat.shinyapps.io/mPower/>

## Example: Case-Control Gut Microbiome Study

Assume you are designing a case-control gut microbiome study (e.g.,
cancer vs healthy individuals) and want to determine the sample size
needed to achieve 80% community-level and taxa-level power. Assume a log
2 fold change of 2 based on the literature.

### Estimate Taxa-Level Power for Case-Control Study

For a study with 50 samples and 10% differential taxa:

``` r
library(mPower)
data(feature.dat)
feature.dat <- feature.dat[rowSums(feature.dat != 0) > 2, ]
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
```

``` r
res1 <- mPower(feature.dat = feature.dat, model.paras = model.paras,
               test = 'Taxa', design = 'CaseControl',
               nSams = 50, grp.ratio = 0.5,
               iters = 10, alpha = 0.05, distance = 'BC',
               diff.otu.pct = 0.1, diff.otu.direct = 'balanced',diff.otu.mode = 'random',
               covariate.eff.min = 0, covariate.eff.maxs = c(1, 2),
               confounder = 'no', depth.mu = 10000, depth.theta = 5, verbose = F)
```

``` r
knitr::kable(res1$power)
```

``` r
res1$plot
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />
