% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mPower.R
\name{mPower}
\alias{mPower}
\title{A Real Data-based Power Analysis Tool for Microbiome Study Design}
\usage{
mPower(
  feature.dat,
  model.paras,
  test = c("Taxa", "Community"),
  design = c("CaseControl", "CrossSectional", "MatchedPair"),
  nSams = 50,
  grp.ratio = 0.5,
  iters = 50,
  alpha = 0.05,
  distance = c("BC", "Jaccard"),
  diff.otu.pct = 0.1,
  diff.otu.direct = c("unbalanced", "balanced"),
  diff.otu.mode = c("random", "abundant", "rare"),
  covariate.eff.min = 0,
  covariate.eff.maxs = c(1, 2),
  prev.filter = 0.1,
  max.abund.filter = 0.002,
  confounder = c("no", "yes"),
  conf.cov.cor = 0.6,
  conf.diff.otu.pct = 0.1,
  conf.nondiff.otu.pct = 0.1,
  confounder.eff.min = 0,
  confounder.eff.max = 1,
  depth.mu = 10000,
  depth.sd = 4000,
  verbose = TRUE
)
}
\arguments{
\item{feature.dat}{data frame where taxa are in rows and samples are in columns.}

\item{model.paras}{list of model parameters output from the \code{EstPara} function.}

\item{test}{should be "Taxa" or "Community", indicating estimate taxa-level power or community-level power.}

\item{design}{the experimental design should be "CaseControl" or "CrossSectional" or "MatchedPair".}

\item{nSams}{an integer or a vector of integers, indicating the number of samples or subjects for power estimate.}

\item{grp.ratio}{a value between 0 and 1, only need to be supplied for "CaseControl" design. i.e., group ratio=0.6 (60/(60+40)=0.6) for 100 samples means comparing group1: 60 vs group2: 40.}

\item{iters}{the number of iterations for simulating microbiome data for power estimation. Community-level test: at least 500; Taxon-level test: at least 50.}

\item{alpha}{a value between 0 and 1 indicating the significance level for detecting differential taxa. Alpha in community-level test: this is the probability of rejecting the null hypothesis when it is true. The significance level is typically set at 0.05, meaning that there is a 5 percent chance we would reject the null hypothesis when it's true. Alpha in taxon-level test: False discovery rate(FDR), the expected proportion of non-differential taxa in the discovery. A typical FDR level can be 5 percent or 10 percent.}

\item{distance}{the type of distance metric used in community-level power estimation: "BC" or "Jaccard". Should be ignored for Taxa-level power estimation. Jaccard distance assesses community similarity by species presence or absence; Bray-Curtis (BC) distance considers both presence and species abundance.}

\item{diff.otu.pct}{a numeric value between 0 and 1 indicating the percentage of differential taxa to be estimated. If 0, global null setting is simulated. The default is 0.1.}

\item{diff.otu.direct}{should be "balanced" or "unbalanced". "balanced": the direction of change for these differential taxa is random, "unbalanced": direction of change is the same. The default is "balanced".}

\item{diff.otu.mode}{should be "rare" or "abundant" or "random". "abundant": differential taxa come from the top quartile of the abundance distribution, "rare": differential OTU come from the bottom quartile of the abundance distribution, and "random": random set. The default is "abundant".}

\item{covariate.eff.min}{a number indicating the minimal log2 fold change associated with the condition ("CaseControl" design), or the expected minimal change in the condition for a one-unit change in the taxa ("CrossSectional" design), or the minimal log2 fold change associated with the condition ("MatchedPair" design).}

\item{covariate.eff.maxs}{a number or a vector indicating the lower and upper range for max log2 fold change. The definition of max log2 fold change: In "CaseControl" design (or "CrossSectional" with binary covariate), it represents the log2 fold change (LFC) between the conditions; In "CrossSectional" desgin with continuous covariate, it represents the expected LFC in response to one-unit change of the condition.}

\item{prev.filter}{a value between 0 and 1. The prevalence cutoff, defined as the percentage of non-zero values, determines the threshold below which taxa will be filtered out.}

\item{max.abund.filter}{a value between 0 and 1. The maximum relative abundance threshold determines the point below which taxa will be excluded from the analysis.}

\item{confounder}{the presence of confounders: "yes" or "no". If you have confounder, power is mostly affected by the correlation between the condition and the confounder, other parameters have much smaller effects and can be set as the default.}

\item{conf.cov.cor}{a value between 0 and 1, indicating the correlation between the condition and the confounder. A large value indicates stronger correlation. Condition is generated based on: condition= (correlation^2/(1-correlation^2))^0.5 x confounder + N(0,1). Binary condition can be generated by dichotomizing condition using group ratio.}

\item{conf.diff.otu.pct}{percentage of taxa associated with the confounder and the condition. A numeric value less than percentage of differential taxa associated with the condition. The percentage of taxa affected by the confounder and the condition.}

\item{conf.nondiff.otu.pct}{percentage of taxa only associated with the confounder. The percentage of taxa affected by the confounder but not the condition.}

\item{confounder.eff.min}{the minimum log2 fold change of taxa associated with the confounder.}

\item{confounder.eff.max}{the maximum log2 fold change of taxa associated with the confounder.}

\item{depth.mu}{the mean sequencing depth for their experiment could be reasonably determined. The default is 10,000.}

\item{depth.sd}{the standard deviation of the sequencing depth. The default is 4,000.}

\item{verbose}{Logical. If TRUE, prints progress messages. Default is TRUE.}
}
\value{
A list with the elements
\item{call}{the call}
\item{plot}{the power curves}
\item{pOCR}{the probability of making at least one correct rejection. Only returns for Taxa-level power estimate.}
\item{aTPR}{the average true positive rate. Only returns for Taxa-level power estimate.}
\item{power}{Community-level power.}
}
\description{
This function estimates the power for microbiome study design based on various experimental designs and parameters.
}
\examples{
data(feature.dat)
## Exclude feature exist in less than 2 samples
feature.dat <- feature.dat[rowSums(feature.dat != 0) > 2, ]
## Estimate the parameters
model.paras <- EstPara(ref.otu.tab = feature.dat)
## Estimate Community-level power for a case-control study
res1 <- mPower(feature.dat = feature.dat, model.paras = model.paras,
              test = 'Community', design = 'CaseControl',
              nSams = 50, grp.ratio = 0.5,
              iters = 500, alpha = 0.05, distance = 'BC',
              diff.otu.pct = 0.1, diff.otu.direct = 'balanced',diff.otu.mode = 'random',
              covariate.eff.min = 0, covariate.eff.maxs = c(1, 2, 3),
              confounder = 'no', depth.mu = 10000, depth.sd = 4000)
## Estimate Taxa-level power for a case-control study
res2 <- mPower(feature.dat = feature.dat, model.paras = model.paras,
              test = 'Taxa', design = 'CaseControl',
              nSams = c(20, 40, 80), grp.ratio = 0.5,
              iters = 100, alpha = 0.05, distance = 'BC',
              diff.otu.pct = 0.1, diff.otu.direct = 'balanced',diff.otu.mode = 'random',
              covariate.eff.min = 0, covariate.eff.maxs = 2,
              prev.filter = 0.1, max.abund.filter = 0.002,
              confounder = 'yes', depth.mu = 100000, depth.sd = 4000)
}
\references{
Lu Yang, Jun Chen. mPower: a Real Data-based Power Analysis Tool for Microbiome Study Design.
}
\author{
Lu Yang \email{luyang1005@gmail.com} & Jun Chen \email{Chen.Jun2@mayo.edu}
}
