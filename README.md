mPower: A Real Data-based Power Analysis Tool for Microbiome Study
Design
================
2024-07-17

This package estimates the power for microbiome study design based on
various experimental designs and parameters.

![](workflows.png)

## 1. Usage

### 1.1 Package installation

You can install mPower as follows:

``` r
# install.packages("devtools")
devtools::install_github("chloelulu/mPower")
```

### 1.2 Shiny R App

You can also use our Shiny App without coding skills, please access the
Shiny App at: <https://microbiomestat.shinyapps.io/mPower/>

## 2. Examples

Assume you are designing a case-control gut microbiome study (e.g.,
cancer vs healthy individuals) and want to determine the sample size
needed to achieve 80% community-level and taxa-level power. Assume a log
2 fold change of 2 based on the literature.

``` r
library(mPower)
data(feature.dat)
```

Preprocessing: exclude features present in fewer than 2 samples.

``` r
feature.dat <- feature.dat[rowSums(feature.dat != 0) > 2, ]
```

Estimate the parameters

``` r
model.paras <- EstPara(ref.otu.tab = feature.dat)
```

#### 2.1 Estimate Community-Level Power for Case-Control Study

``` r
res1 <- mPower(feature.dat = feature.dat, model.paras = model.paras,
               test = 'Community', design = 'CaseControl',
               nSams = 50, grp.ratio = 0.5,
               iters = 500, alpha = 0.05, distance = 'BC',
               diff.otu.pct = 0.1, diff.otu.direct = 'balanced',diff.otu.mode = 'random',
               covariate.eff.min = 0, covariate.eff.maxs = c(1, 2),
               confounder = 'no', depth.mu = 10000, depth.theta = 5, verbose = F)
```

##### 2.1.1 Output1: Community-level power table

“power” column indicates community-level power: the probability of
rejecting the null hypothesis when the null hypothesis is false.

| Sample size | max log2 fold change | power |        SD |      ymax |      ymin |
|:------------|:---------------------|------:|----------:|----------:|----------:|
| 50          | 1                    | 0.386 | 0.4873181 | 0.4287153 | 0.3432847 |
| 50          | 2                    | 0.748 | 0.4345961 | 0.7860940 | 0.7099060 |

##### 2.1.2 Output2: $R^2$ (variance explained) and community-level power curve

$R^2$ is generated from PERMANOVA, denotes the proportion of the total
variation in the response data that is explained by the explanatory
variables. ![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

#### 2.2 Estimate Taxa-Level Power for Case-Control Study

``` r
res2 <- mPower(feature.dat = feature.dat, model.paras = model.paras,
               test = 'Taxa', design = 'CaseControl',
               nSams = c(20, 80), grp.ratio = 0.5,
               iters = 50, alpha = 0.05, distance = 'BC',
               diff.otu.pct = 0.1, diff.otu.direct = 'balanced',diff.otu.mode = 'random',
               covariate.eff.min = 0, covariate.eff.maxs = 2,
               prev.filter = 0.1, max.abund.filter = 0.002,
               confounder = 'yes', depth.mu = 100000, depth.theta = 5, verbose = F)
```

##### 2.2.1 Output1: Taxa-level power table

“pOCR” indicates the probability of making at least one correct
rejection, meaning at least one taxon has a FDR-adjusted pvalue less
than alpha(0.05).

| Sample size | max log2 fold change | pOCR |
|:------------|:---------------------|-----:|
| 20          | 2                    | 0.40 |
| 80          | 2                    | 0.86 |

##### 2.2.2 Output2: Taxa-level power table

“aTPR” column indicates the average true positive rate. SD(Standard
deviation) for each setting, i.e., sample size = 20 & max log2 fold
change =2. ymax and ymin represents the upper and lower bound of the 95%
confidence interval for the aTPR.

| Sample size | max log2 fold change |      aTPR |        SD |      ymax |      ymin |
|:------------|:---------------------|----------:|----------:|----------:|----------:|
| 20          | 2                    | 0.0265493 | 0.0395909 | 0.0375233 | 0.0155752 |
| 80          | 2                    | 0.0916704 | 0.0632646 | 0.1092064 | 0.0741344 |

##### 2.2.3 Output3: aTPR power curve(left) and pOCR power curve(right)

![](README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

### References

(Yang and Chen 2022, 2023, n.d.)

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-yang2022comprehensive" class="csl-entry">

Yang, Lu, and Jun Chen. 2022. “A Comprehensive Evaluation of Microbial
Differential Abundance Analysis Methods: Current Status and Potential
Solutions.” *Microbiome* 10 (1): 130.

</div>

<div id="ref-yang2023benchmarking" class="csl-entry">

———. 2023. “Benchmarking Differential Abundance Analysis Methods for
Correlated Microbiome Sequencing Data.” *Briefings in Bioinformatics* 24
(1): bbac607.

</div>

<div id="ref-yang2024mpower" class="csl-entry">

———. n.d. “mPower: A Real Data-Based Power Analysis Tool for Microbiome
Study Design.”

</div>

</div>
