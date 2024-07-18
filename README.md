
# mPower: A Real Data-based Power Analysis Tool for Microbiome Study Design

This package estimates the power for microbiome study design based on
various experimental designs and parameters. For details of the proposed
method, please refer to the following paper: Lu Yang, Jun Chen. mPower:
A Real Data-based Power Analysis Tool for Microbiome Study Design

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
library(reactable) # show table in Markdown
library(mPower)
data(feature.dat)
```

Exclude features present in fewer than 2 samples

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

##### 2.1.1 Output1: Community-level power table (“power” column indicates community-level power)

<div class="reactable html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-34857f0edcbfbb2c0cba" style="width:auto;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-34857f0edcbfbb2c0cba">{"x":{"tag":{"name":"Reactable","attribs":{"data":{"Sample size":["50","50"],"max log2 fold change":["1","2"],"power":[0.378,0.732],"SD":[0.485373232006768,0.443361182645205],"ymax":[0.420544839616056,0.770862321124791],"ymin":[0.335455160383944,0.693137678875208]},"columns":[{"id":"Sample size","name":"Sample size","type":"factor"},{"id":"max log2 fold change","name":"max log2 fold change","type":"factor"},{"id":"power","name":"power","type":"numeric"},{"id":"SD","name":"SD","type":"numeric"},{"id":"ymax","name":"ymax","type":"numeric"},{"id":"ymin","name":"ymin","type":"numeric"}],"resizable":true,"dataKey":"50c82439d041f8f057a878bdcd02b1f7"},"children":[]},"class":"reactR_markup"},"evals":[],"jsHooks":[]}</script>

##### 2.1.2 Output2: R2 (variance explained) and community-level power curve

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

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

##### 2.2.1 Output1: Taxa-level power table (“pOCR” column indicates the probability of making at least one correct rejection)

<div class="reactable html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-595501ac9c6ffe9b0f1c" style="width:auto;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-595501ac9c6ffe9b0f1c">{"x":{"tag":{"name":"Reactable","attribs":{"data":{"Sample size":["20","80"],"max log2 fold change":["2","2"],"pOCR":[0.24,1]},"columns":[{"id":"Sample size","name":"Sample size","type":"factor"},{"id":"max log2 fold change","name":"max log2 fold change","type":"factor"},{"id":"pOCR","name":"pOCR","type":"numeric"}],"resizable":true,"dataKey":"a8e841316291e83d954e6838c8a712a1"},"children":[]},"class":"reactR_markup"},"evals":[],"jsHooks":[]}</script>

##### 2.2.2 Output2: Taxa-level power table (“aTPR” column indicates the average true positive rate)

<div class="reactable html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-e0b9afa38c2b5d42ced5" style="width:auto;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-e0b9afa38c2b5d42ced5">{"x":{"tag":{"name":"Reactable","attribs":{"data":{"Sample size":["20","80"],"max log2 fold change":["2","2"],"aTPR":[0.0183058194454075,0.115363484340006],"SD":[0.0381776159825601,0.07665312933019],"ymax":[0.0288881146965228,0.136610647779036],"ymin":[0.00772352419429232,0.0941163209009754]},"columns":[{"id":"Sample size","name":"Sample size","type":"factor"},{"id":"max log2 fold change","name":"max log2 fold change","type":"factor"},{"id":"aTPR","name":"aTPR","type":"numeric"},{"id":"SD","name":"SD","type":"numeric"},{"id":"ymax","name":"ymax","type":"numeric"},{"id":"ymin","name":"ymin","type":"numeric"}],"resizable":true,"dataKey":"b641aebe47037b9226d3d4746ac65025"},"children":[]},"class":"reactR_markup"},"evals":[],"jsHooks":[]}</script>

##### 2.2.3 Output3: aTPR power curve(left) and pOCR power curve(right)

![](README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->
