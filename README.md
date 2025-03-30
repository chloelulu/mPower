mPower: A Real Data-based Power Analysis Tool for Microbiome Study
Design
================

This package estimates the power for microbiome study design based on
various experimental designs and parameters.

![](workflows.png)

# 1. Usage

## 1.1 Package installation

You can install mPower as follows:

``` r
# install.packages("devtools")
devtools::install_github("chloelulu/mPower")
```

## 1.2 Shiny R App

You can also use our Shiny App without coding skills, please access the
Shiny App at: <https://microbiomestat.shinyapps.io/mPower/>

# 2. Examples

You need to have a reference feature table, it should be a collection of
microbiome sequencing samples from a study population at a specific
sampling site. It should be large enough to capture the main
compositional variation in the population of interest.

``` r
library(mPower)
```

    ## Warning: replacing previous import 'MASS::select' by 'dplyr::select' when
    ## loading 'mPower'

    ## Warning: replacing previous import 'dplyr::lag' by 'stats::lag' when loading
    ## 'mPower'

    ## Warning: replacing previous import 'dplyr::filter' by 'stats::filter' when
    ## loading 'mPower'

``` r
data(feature.dat)
```

Preprocessing: exclude features present in fewer than 2 samples.

``` r
feature.dat <- feature.dat[rowSums(feature.dat != 0) > 2, ]
```

Estimate the parameters: We will obtain posterior samples of the
underlying composition based on an empirical Bayes model, and generate
the absolute abundance of the reference feature table.

``` r
model.paras <- EstPara(ref.otu.tab = feature.dat)
```

### 2.1 Estimate Community-Level Power for Case-Control Study

Assume you are designing a case-control gut microbiome study (e.g.,
cancer vs healthy individuals) and want to determine the sample size
needed to achieve 80% community-level and taxa-level power.

- Set the sample size: if you are not sure how many sample size could be
  enough to achieve desirable power, i.e., 90% power, you can set a
  range of sample sizes, such as (20, 60, 100).

- Set the iterations: 500 at least for community-level power estimate.

- Set alpha: alpha for community-level means the probability of
  rejecting the null hypothesis when it is true. Default 0.05 is chosen.

- Set distance: Bray-Curtis or Jaccard distance can be chosen. Here we
  choose “BC” considers both presence and species abundance.

- Set the effect sizes:

  - 1)  the percentage of differential taxa: can be estimated based on
        the association p-value distribution.
  - 2)  “mPower” allows the effect size (log2 fold change, LFC, between
        two groups, or in response to 1 S.D. change of a continuous
        variable) to come from a probabilistic distribution. In current
        implementation, it assumes a uniform distribution on the
        interval \[min LFC , max LFC\]. The max LFC can be easily
        estimated from real data while min LFC is set to 0 by default.
        When min LFC = max LFC, a fixed LFC is set for all differential
        taxa.

- Set the differential setting to set the direction of differential
  taxa, either “unbalanced” for creating strong compositional effect or
  “balanced” for moderate compositional effect. Here we choose
  “balanced”.

- Set the differential taxa from “rare” or “abundant” or “random”,
  indicating the direction of change for these differential taxa.

``` r
res1 <- mPower(feature.dat = feature.dat, model.paras = model.paras,
               test = 'Community', design = 'CaseControl',
               nSams = c(20, 60, 100), grp.ratio = 0.5,
               iters = 500, alpha = 0.05, distance = 'BC',
               diff.otu.pct = 0.1, 
               covariate.eff.min = 0, covariate.eff.maxs = 2,
               diff.otu.direct = 'balanced',diff.otu.mode = 'random',
               confounder = 'no', depth.mu = 10000, depth.sd = 4000, verbose = F)
```

#### 2.1.1 Output1: Community-level power table

“power”(community-level power): the probability of rejecting the null
hypothesis when the null hypothesis is false. An pOCR at least 90% will
ensure a high likelihood of making some discoveries. Thus in this
example, around 100 samples should be well-powered.

``` r
knitr::kable(res1$power, format = "markdown")
```

| Sample size | max log2 fold change | power |
|:------------|:---------------------|------:|
| 20          | 2                    | 0.428 |
| 60          | 2                    | 0.854 |
| 100         | 2                    | 0.930 |

#### 2.1.2 Output2: $R^2$ (variance explained) and community-level power curve

$R^2$: the proportion of the total variation in the response data that
is explained by the explanatory variables.

``` r
res1$plot
```

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

### 2.2 Estimate Taxa-Level Power for Case-Control Study

``` r
res2 <- mPower(feature.dat = feature.dat, model.paras = model.paras,
               test = 'Taxa', design = 'CaseControl',
               nSams = c(20, 60, 100), grp.ratio = 0.5,
               iters = 50, alpha = 0.05, distance = 'BC',
               diff.otu.pct = 0.1, diff.otu.direct = 'balanced',diff.otu.mode = 'random',
               covariate.eff.min = 0, covariate.eff.maxs = 2,
               prev.filter = 0.1, max.abund.filter = 0.002,
               confounder = 'yes', depth.mu = 100000, depth.sd = 4000, verbose = F)
```

    ## Differential abundance analysis is testing:  20  samples and  145 out of 204  features!
    ## Differential abundance analysis is testing:  20  samples and  174 out of 237  features!
    ## Differential abundance analysis is testing:  20  samples and  169 out of 225  features!
    ## Differential abundance analysis is testing:  20  samples and  168 out of 249  features!
    ## Differential abundance analysis is testing:  20  samples and  192 out of 258  features!
    ## Differential abundance analysis is testing:  20  samples and  159 out of 203  features!
    ## Differential abundance analysis is testing:  20  samples and  177 out of 216  features!
    ## Differential abundance analysis is testing:  20  samples and  166 out of 230  features!
    ## Differential abundance analysis is testing:  20  samples and  185 out of 248  features!
    ## Differential abundance analysis is testing:  20  samples and  183 out of 246  features!
    ## Differential abundance analysis is testing:  20  samples and  171 out of 239  features!
    ## Differential abundance analysis is testing:  20  samples and  159 out of 236  features!
    ## Differential abundance analysis is testing:  20  samples and  170 out of 225  features!
    ## Differential abundance analysis is testing:  20  samples and  171 out of 234  features!
    ## Differential abundance analysis is testing:  20  samples and  159 out of 227  features!
    ## Differential abundance analysis is testing:  20  samples and  171 out of 235  features!
    ## Differential abundance analysis is testing:  20  samples and  162 out of 211  features!
    ## Differential abundance analysis is testing:  20  samples and  170 out of 232  features!
    ## Differential abundance analysis is testing:  20  samples and  185 out of 259  features!
    ## Differential abundance analysis is testing:  20  samples and  173 out of 230  features!
    ## Differential abundance analysis is testing:  20  samples and  168 out of 239  features!
    ## Differential abundance analysis is testing:  20  samples and  174 out of 251  features!
    ## Differential abundance analysis is testing:  20  samples and  195 out of 247  features!
    ## Differential abundance analysis is testing:  20  samples and  172 out of 233  features!
    ## Differential abundance analysis is testing:  20  samples and  156 out of 225  features!
    ## Differential abundance analysis is testing:  20  samples and  152 out of 219  features!
    ## Differential abundance analysis is testing:  20  samples and  178 out of 244  features!
    ## Differential abundance analysis is testing:  20  samples and  168 out of 230  features!
    ## Differential abundance analysis is testing:  20  samples and  166 out of 231  features!
    ## Differential abundance analysis is testing:  20  samples and  169 out of 237  features!
    ## Differential abundance analysis is testing:  20  samples and  205 out of 264  features!
    ## Differential abundance analysis is testing:  20  samples and  164 out of 227  features!
    ## Differential abundance analysis is testing:  20  samples and  195 out of 260  features!
    ## Differential abundance analysis is testing:  20  samples and  164 out of 211  features!
    ## Differential abundance analysis is testing:  20  samples and  192 out of 260  features!
    ## Differential abundance analysis is testing:  20  samples and  184 out of 248  features!
    ## Differential abundance analysis is testing:  20  samples and  172 out of 234  features!
    ## Differential abundance analysis is testing:  20  samples and  172 out of 245  features!
    ## Differential abundance analysis is testing:  20  samples and  171 out of 244  features!
    ## Differential abundance analysis is testing:  20  samples and  189 out of 255  features!
    ## Differential abundance analysis is testing:  20  samples and  197 out of 266  features!
    ## Differential abundance analysis is testing:  20  samples and  187 out of 263  features!
    ## Differential abundance analysis is testing:  20  samples and  181 out of 267  features!
    ## Differential abundance analysis is testing:  20  samples and  196 out of 242  features!
    ## Differential abundance analysis is testing:  20  samples and  191 out of 265  features!
    ## Differential abundance analysis is testing:  20  samples and  174 out of 250  features!
    ## Differential abundance analysis is testing:  20  samples and  175 out of 238  features!
    ## Differential abundance analysis is testing:  20  samples and  173 out of 240  features!
    ## Differential abundance analysis is testing:  20  samples and  177 out of 224  features!
    ## Differential abundance analysis is testing:  20  samples and  158 out of 233  features!
    ## Differential abundance analysis is testing:  60  samples and  236 out of 355  features!
    ## Differential abundance analysis is testing:  60  samples and  259 out of 357  features!
    ## Differential abundance analysis is testing:  60  samples and  256 out of 356  features!
    ## Differential abundance analysis is testing:  60  samples and  248 out of 357  features!
    ## Differential abundance analysis is testing:  60  samples and  224 out of 341  features!
    ## Differential abundance analysis is testing:  60  samples and  253 out of 360  features!
    ## Differential abundance analysis is testing:  60  samples and  220 out of 346  features!
    ## Differential abundance analysis is testing:  60  samples and  249 out of 370  features!
    ## Differential abundance analysis is testing:  60  samples and  237 out of 356  features!
    ## Differential abundance analysis is testing:  60  samples and  236 out of 346  features!
    ## Differential abundance analysis is testing:  60  samples and  230 out of 353  features!
    ## Differential abundance analysis is testing:  60  samples and  257 out of 359  features!
    ## Differential abundance analysis is testing:  60  samples and  238 out of 349  features!
    ## Differential abundance analysis is testing:  60  samples and  229 out of 345  features!
    ## Differential abundance analysis is testing:  60  samples and  243 out of 353  features!
    ## Differential abundance analysis is testing:  60  samples and  245 out of 357  features!
    ## Differential abundance analysis is testing:  60  samples and  264 out of 363  features!
    ## Differential abundance analysis is testing:  60  samples and  223 out of 347  features!
    ## Differential abundance analysis is testing:  60  samples and  249 out of 364  features!
    ## Differential abundance analysis is testing:  60  samples and  240 out of 353  features!
    ## Differential abundance analysis is testing:  60  samples and  246 out of 362  features!
    ## Differential abundance analysis is testing:  60  samples and  247 out of 354  features!
    ## Differential abundance analysis is testing:  60  samples and  251 out of 364  features!
    ## Differential abundance analysis is testing:  60  samples and  256 out of 370  features!
    ## Differential abundance analysis is testing:  60  samples and  239 out of 357  features!
    ## Differential abundance analysis is testing:  60  samples and  245 out of 361  features!
    ## Differential abundance analysis is testing:  60  samples and  250 out of 354  features!
    ## Differential abundance analysis is testing:  60  samples and  249 out of 371  features!
    ## Differential abundance analysis is testing:  60  samples and  239 out of 357  features!
    ## Differential abundance analysis is testing:  60  samples and  257 out of 361  features!
    ## Differential abundance analysis is testing:  60  samples and  228 out of 349  features!
    ## Differential abundance analysis is testing:  60  samples and  225 out of 348  features!
    ## Differential abundance analysis is testing:  60  samples and  248 out of 359  features!
    ## Differential abundance analysis is testing:  60  samples and  259 out of 357  features!
    ## Differential abundance analysis is testing:  60  samples and  235 out of 358  features!
    ## Differential abundance analysis is testing:  60  samples and  248 out of 360  features!
    ## Differential abundance analysis is testing:  60  samples and  238 out of 344  features!
    ## Differential abundance analysis is testing:  60  samples and  263 out of 364  features!
    ## Differential abundance analysis is testing:  60  samples and  232 out of 341  features!
    ## Differential abundance analysis is testing:  60  samples and  244 out of 367  features!
    ## Differential abundance analysis is testing:  60  samples and  244 out of 350  features!
    ## Differential abundance analysis is testing:  60  samples and  237 out of 359  features!
    ## Differential abundance analysis is testing:  60  samples and  244 out of 361  features!
    ## Differential abundance analysis is testing:  60  samples and  237 out of 364  features!
    ## Differential abundance analysis is testing:  60  samples and  250 out of 355  features!
    ## Differential abundance analysis is testing:  60  samples and  248 out of 355  features!
    ## Differential abundance analysis is testing:  60  samples and  232 out of 339  features!
    ## Differential abundance analysis is testing:  60  samples and  242 out of 365  features!
    ## Differential abundance analysis is testing:  60  samples and  247 out of 360  features!
    ## Differential abundance analysis is testing:  60  samples and  245 out of 362  features!
    ## Differential abundance analysis is testing:  100  samples and  253 out of 375  features!
    ## Differential abundance analysis is testing:  100  samples and  264 out of 379  features!
    ## Differential abundance analysis is testing:  100  samples and  274 out of 382  features!
    ## Differential abundance analysis is testing:  100  samples and  261 out of 376  features!
    ## Differential abundance analysis is testing:  100  samples and  261 out of 372  features!
    ## Differential abundance analysis is testing:  100  samples and  266 out of 380  features!
    ## Differential abundance analysis is testing:  100  samples and  257 out of 376  features!
    ## Differential abundance analysis is testing:  100  samples and  287 out of 380  features!
    ## Differential abundance analysis is testing:  100  samples and  268 out of 379  features!
    ## Differential abundance analysis is testing:  100  samples and  260 out of 371  features!
    ## Differential abundance analysis is testing:  100  samples and  258 out of 381  features!
    ## Differential abundance analysis is testing:  100  samples and  270 out of 381  features!
    ## Differential abundance analysis is testing:  100  samples and  258 out of 376  features!
    ## Differential abundance analysis is testing:  100  samples and  269 out of 378  features!
    ## Differential abundance analysis is testing:  100  samples and  266 out of 373  features!
    ## Differential abundance analysis is testing:  100  samples and  259 out of 377  features!
    ## Differential abundance analysis is testing:  100  samples and  276 out of 381  features!
    ## Differential abundance analysis is testing:  100  samples and  268 out of 377  features!
    ## Differential abundance analysis is testing:  100  samples and  253 out of 376  features!
    ## Differential abundance analysis is testing:  100  samples and  283 out of 380  features!
    ## Differential abundance analysis is testing:  100  samples and  272 out of 379  features!
    ## Differential abundance analysis is testing:  100  samples and  258 out of 370  features!
    ## Differential abundance analysis is testing:  100  samples and  269 out of 378  features!
    ## Differential abundance analysis is testing:  100  samples and  269 out of 381  features!
    ## Differential abundance analysis is testing:  100  samples and  255 out of 377  features!
    ## Differential abundance analysis is testing:  100  samples and  255 out of 378  features!
    ## Differential abundance analysis is testing:  100  samples and  279 out of 378  features!
    ## Differential abundance analysis is testing:  100  samples and  258 out of 378  features!
    ## Differential abundance analysis is testing:  100  samples and  266 out of 377  features!
    ## Differential abundance analysis is testing:  100  samples and  268 out of 374  features!
    ## Differential abundance analysis is testing:  100  samples and  259 out of 376  features!
    ## Differential abundance analysis is testing:  100  samples and  275 out of 382  features!
    ## Differential abundance analysis is testing:  100  samples and  281 out of 380  features!
    ## Differential abundance analysis is testing:  100  samples and  270 out of 377  features!
    ## Differential abundance analysis is testing:  100  samples and  254 out of 376  features!
    ## Differential abundance analysis is testing:  100  samples and  257 out of 377  features!
    ## Differential abundance analysis is testing:  100  samples and  254 out of 379  features!
    ## Differential abundance analysis is testing:  100  samples and  261 out of 373  features!
    ## Differential abundance analysis is testing:  100  samples and  266 out of 378  features!
    ## Differential abundance analysis is testing:  100  samples and  263 out of 378  features!
    ## Differential abundance analysis is testing:  100  samples and  266 out of 380  features!
    ## Differential abundance analysis is testing:  100  samples and  258 out of 375  features!
    ## Differential abundance analysis is testing:  100  samples and  255 out of 377  features!
    ## Differential abundance analysis is testing:  100  samples and  262 out of 382  features!
    ## Differential abundance analysis is testing:  100  samples and  243 out of 374  features!
    ## Differential abundance analysis is testing:  100  samples and  250 out of 374  features!
    ## Differential abundance analysis is testing:  100  samples and  268 out of 378  features!
    ## Differential abundance analysis is testing:  100  samples and  268 out of 376  features!
    ## Differential abundance analysis is testing:  100  samples and  262 out of 378  features!
    ## Differential abundance analysis is testing:  100  samples and  263 out of 379  features!

#### 2.2.1 Output1: Taxa-level power table - aTPR

“aTPR”(average true positive rate): represents the average proportion of
truly differential taxa that are correctly identified as such. SD, ymax
and ymin represents standard deviation, the upper and lower bound of the
95% confidence interval for the aTPR, respectively.

``` r
knitr::kable(res2$aTPR, format = "markdown")
```

| Sample size | max log2 fold change |      aTPR |
|:------------|:---------------------|----------:|
| 20          | 2                    | 0.0113019 |
| 60          | 2                    | 0.0790285 |
| 100         | 2                    | 0.1400519 |

#### 2.2.2 Output1: Taxa-level power table - pOCR

“pOCR”(probability of making at least one correct rejection): akin to
the conventional understanding of power, differs in that the specific
taxa rejected need not be consistent.

``` r
knitr::kable(res2$pOCR, format = "markdown")
```

| Sample size | max log2 fold change | pOCR |
|:------------|:---------------------|-----:|
| 20          | 2                    | 0.16 |
| 60          | 2                    | 0.88 |
| 100         | 2                    | 1.00 |

#### 2.2.3 Output3: aTPR power curve(left) and pOCR power curve(right)

``` r
res2$plot
```

![](README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

## References

(Yang and Chen 2022, 2023, n.d.)

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

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
