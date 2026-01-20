
utils::globalVariables(c("otu_table", "sample_data", "lmer", "mclapply", "%do%", "foreach", "mlv", "value", "nSam", "variable", "ymin", "ymax","Sample size","max log2 fold change"))

rowMeans <- function(data){
  return(apply(data, 1, function(x) mean(x)))
}

rowSums <- function(data){
  return(apply(data, 1, function(x) sum(x)))
}

rowMaxs <- function(data){
  return(apply(data, 1, function(x) x[which.max(x)]))
}

colSums <- function(data){
  return(apply(data, 2, function(x) sum(x)))
}

colMeans <- function(data){
  return(apply(data, 2, function(x) mean(x)))
}


generate_vector <- function(lower_bound, upper_bound, n) {
  if(n < 2) {
    stop("n must be at least 2 to create a range.")
  }
  vector <- round(seq(lower_bound, upper_bound, length.out = n), 2)
  return(vector)
}


winsor.fun <- function(Y, quan, feature.dat.type) {
  if (feature.dat.type == 'count') {
    N <- colSums(Y)
    P <- t(t(Y) / N)
    cut <- apply(P, 1, quantile, quan)
    Cut <- matrix(rep(cut, ncol(Y)), nrow(Y))
    ind <- P > Cut
    P[ind] <- Cut[ind]
    Y <- round(t(t(P) * N))
  }

  if (feature.dat.type == 'proportion') {
    cut <- apply(Y, 1, quantile, quan)
    Cut <- matrix(rep(cut, ncol(Y)), nrow(Y))
    ind <- Y > Cut
    Y[ind] <- Cut[ind]
  }
  return(Y)
}

linda <- function(feature.dat, meta.dat, phyloseq.obj = NULL, formula, feature.dat.type = c('count', 'proportion'),
                  prev.filter = 0, mean.abund.filter = 0, max.abund.filter = 0,
                  is.winsor = TRUE, outlier.pct = 0.03,
                  adaptive = TRUE, zero.handling = c('pseudo-count', 'imputation'),
                  pseudo.cnt = 0.5, corr.cut = 0.1,
                  p.adj.method = 'BH', alpha = 0.05,
                  n.cores = 1, verbose = TRUE) {

  feature.dat.type <- match.arg(feature.dat.type)

  if (!is.null(phyloseq.obj)) {
    feature.dat.type <- 'count'

    feature.dat <- otu_table(phyloseq.obj)@.Data
    meta.dat <- data.frame(sample_data(phyloseq.obj))
  }

  if(any(is.na(feature.dat))) {
    stop('The feature table contains NAs! Please remove!\n')
  }
  allvars <- all.vars(as.formula(formula))
  Z <- as.data.frame(meta.dat[, allvars])
  m0 <- nrow(feature.dat)
  ###############################################################################
  # Filter sample
  keep.sam <- which(rowSums(is.na(Z)) == 0)
  Y <- feature.dat[, keep.sam]
  Z <- as.data.frame(Z[keep.sam, ])
  names(Z) <- allvars

  # Filter features
  temp <- t(t(Y) / colSums(Y))

  keep.tax <- rowMeans(temp != 0) >= prev.filter & rowMeans(temp) >= mean.abund.filter & rowMaxs(temp) >= max.abund.filter
  names(keep.tax) <- rownames(Y)
  rm(temp)
  # if (verbose)  cat(sum(!keep.tax), ' features are filtered!\n')
  Y <- Y[keep.tax, ]

  n <- ncol(Y)
  m <- nrow(Y)

  ## some samples may have zero total counts after screening taxa
  if(any(colSums(Y) == 0)) {
    ind <- which(colSums(Y) > 0)
    Y <- Y[, ind]
    Z <- as.data.frame(Z[ind, ])
    names(Z) <- allvars
    keep.sam <- keep.sam[ind]
    n <- ncol(Y)
  }

  if (verbose) cat('Differential abundance analysis is testing: ', n, ' samples and ', m,'out of', m0, ' features!\n' )

  if (sum(rowSums(Y != 0) <= 2) != 0) {
    warning('Some features have less than 3 nonzero values!
						They have virtually no statistical power. You may consider filtering them in the analysis!\n')
  }

  ###############################################################################
  ## scaling numerical variables
  ind <- sapply(1 : ncol(Z), function(i) is.numeric(Z[, i]))
  Z[, ind] <- scale(Z[, ind])

  ## winsorization
  if(is.winsor) {
    Y <- winsor.fun(Y, 1 - outlier.pct, feature.dat.type)
  }

  ##
  if(grepl('\\(', formula)) {
    random.effect <- TRUE
  } else {
    random.effect <- FALSE
  }

  if(is.null(rownames(feature.dat))) {
    taxa.name <- (1 : nrow(feature.dat))[keep.tax]
  } else {
    taxa.name <- rownames(feature.dat)[keep.tax]
  }
  if(is.null(rownames(meta.dat))) {
    samp.name <- (1 : nrow(meta.dat))[keep.sam]
  } else {
    samp.name <- rownames(meta.dat)[keep.sam]
  }

  ## handling zeros
  if (feature.dat.type == 'count') {
    if(any(Y == 0)) {
      N <- colSums(Y)
      if(adaptive) {
        logN <- log(N)
        if(random.effect) {
          tmp <- lmer(as.formula(paste0('logN', formula)), Z)
        } else {
          tmp <- lm(as.formula(paste0('logN', formula)), Z)
        }
        corr.pval <- coef(summary(tmp))[-1, "Pr(>|t|)"]
        if(any(corr.pval <= corr.cut)) {
          # if (verbose) cat('Imputation approach is used.\n')
          zero.handling <- "Imputation"
        } else {
          # if (verbose) cat('Pseudo-count approach is used.\n')
          zero.handling <- "Pseudo-count"
        }
      }
      if(zero.handling == 'imputation') {
        N.mat <- matrix(rep(N, m), nrow = m, byrow = TRUE)
        N.mat[Y > 0] <- 0
        tmp <- N[max.col(N.mat)]
        Y <- Y + N.mat / tmp
      } else {
        Y <- Y + pseudo.cnt
      }
    }
  }

  if (feature.dat.type == 'proportion') {
    if(any(Y == 0)) {
      # Half minimum approach

      Y <- t(apply(Y, 1, function (x) {
        x[x == 0] <- 0.5 * min(x[x != 0])
        return(x)
      }))
      colnames(Y) <- samp.name
      rownames(Y) <- taxa.name

    }
  }

  ## CLR transformation
  logY <- log2(Y)
  W <- t(logY) - colMeans(logY)

  ## linear regression

  #	oldw <- getOption('warn')
  #	options(warn = -1)
  if(!random.effect) {
    # if (verbose)  cat("Fit linear models ...\n")
    suppressMessages(fit <- lm(as.formula(paste0('W', formula)), Z))
    res <- do.call(rbind, coef(summary(fit)))
    df <- rep(n - ncol(model.matrix(fit)), m)
  } else {
    # if (verbose)  cat("Fit linear mixed effects models ...\n")
    fun <- function(i) {
      w <- W[, i]
      suppressMessages(fit <- lmer(as.formula(paste0('w', formula)), Z))
      coef(summary(fit))
    }
    if(n.cores > 1) {
      res <- mclapply(c(1 : m), function(i) fun(i), mc.cores = n.cores)
    } else {
      suppressMessages(res <- foreach(i = 1 : m) %do% fun(i))
    }
    res <- do.call(rbind, res)
  }
  #	options(warn = oldw)

  res.intc <- res[which(rownames(res) == '(Intercept)'), ]
  rownames(res.intc) <- NULL
  baseMean <- 2 ^ res.intc[, 1]
  baseMean <- baseMean / sum(baseMean) * 1e6

  output.fun <- function(x) {
    res.voi <- res[which(rownames(res) == x), ]
    rownames(res.voi) <- NULL

    if(random.effect) {
      df <- res.voi[, 3]
    }

    log2FoldChange <- res.voi[, 1]
    lfcSE <- res.voi[, 2]
    #		oldw <- getOption('warn')
    #		options(warn = -1)
    suppressMessages(bias <- mlv(sqrt(n) * log2FoldChange,
                                 method = 'meanshift', kernel = 'gaussian') / sqrt(n))
    #		options(warn = oldw)
    log2FoldChange <- log2FoldChange - bias
    stat <- log2FoldChange / lfcSE

    pvalue <- 2 * pt(-abs(stat), df)
    padj <- p.adjust(pvalue, method = p.adj.method)
    reject <- padj <= alpha
    output <- cbind.data.frame(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, reject, df)
    rownames(output) <- taxa.name
    return(list(bias = bias, output = output))
  }

  variables <- unique(rownames(res))[-1]
  variables.n <- length(variables)
  bias <- rep(NA, variables.n)
  output <- list()
  for(i in 1 : variables.n) {
    tmp <- output.fun(variables[i])
    output[[i]] <- tmp[[2]]
    bias[i] <- tmp[[1]]
  }
  names(output) <- variables

  rownames(Y) <- taxa.name
  colnames(Y) <- samp.name
  rownames(Z) <- samp.name
  # if (verbose)  cat("Completed.\n")
  return(list(variables = variables, bias = bias, output = output, feature.dat.use = Y, meta.dat.use = Z))
}



perform_DAA <- function(feature.dat, meta.dat, grp.name, adj.name, design, prev.filter, max.abund.filter,
                        method = c('LinDA'), verbose=T){
  if(method =='LinDA'){
    if(design %in% c('CaseControl','CrossSectional')){
      if (is.null(adj.name)) {
        formula <- paste0('~', grp.name)
      } else {
        formula <- paste0('~', paste(grp.name, adj.name, sep = '+'))
      }

      obj <- linda(meta.dat = meta.dat, feature.dat = feature.dat, formula=formula, feature.dat.type = "count",prev.filter=prev.filter,max.abund.filter=max.abund.filter, verbose = verbose)
      p.adj.fdr <- obj$output[[grep(grp.name,names(obj$output))]][,'padj'] # Default is BH
      names(p.adj.fdr) <- rownames(obj$output[[grep(grp.name,names(obj$output))]])
      P <- obj$output[[grep(grp.name,names(obj$output))]][,c('pvalue','padj')]
      colnames(P)[grep('padj',colnames(P))] <- 'p.adj.fdr'
    }

    if(design =='ReplicateSampling'){
      if (is.null(adj.name)) {
        formula <- paste0('~', grp.name,'+(1|SubjectID)')
      } else {
        formula <- paste0('~', paste(grp.name, adj.name, sep = '+'), '+(1|SubjectID)')
      }

      obj <- linda(meta.dat = meta.dat, feature.dat = feature.dat, prev.filter=prev.filter,max.abund.filter=max.abund.filter,
                   formula=formula, feature.dat.type = "count", verbose = verbose)
      p.adj.fdr <- obj$output[[grep(grp.name,names(obj$output))]][,'padj'] # Default is BH
      names(p.adj.fdr) <- rownames(obj$output[[grep(grp.name,names(obj$output))]])
      P <- obj$output[[grep(grp.name,names(obj$output))]][,c('pvalue','padj')]
      colnames(P)[2] <- 'p.adj.fdr'
    }

    if(design =='MatchedPair'){
      formula <- paste0('~',grp.name, '+ SubjectID')
      obj <- linda(meta.dat = meta.dat, feature.dat = feature.dat, formula=formula, feature.dat.type = "count",prev.filter=prev.filter,max.abund.filter=max.abund.filter, verbose = verbose)
      p.adj.fdr <- obj$output[[grep(grp.name,names(obj$output))]][,'padj'] # Default is BH
      names(p.adj.fdr) <- rownames(obj$output[[grep(grp.name,names(obj$output))]])
      P <- obj$output[[grep(grp.name,names(obj$output))]][,c('pvalue','padj')]
      colnames(P)[2] <- 'p.adj.fdr'
    }

    if(design =='longitudinal'){
      formula <- '~X*time+Z+(time|SubjectID)'
      obj <- linda(meta.dat = meta.dat, feature.dat = feature.dat,prev.filter=prev.filter,max.abund.filter=max.abund.filter,
                   formula=formula, feature.dat.type = "count", verbose = verbose)

      p.adj.fdr.X <- obj$output[['X1']][,'padj'] # Default is BH
      names(p.adj.fdr.X) <- rownames(obj$output[['X1']])

      p.adj.fdr.time <- obj$output[['time']][,'padj'] # Default is BH
      names(p.adj.fdr.time) <- rownames(obj$output[['time']])

      p.adj.fdr.interaction <- obj$output[['X1:time']][,'padj'] # Default is BH
      names(p.adj.fdr.interaction) <- rownames(obj$output[['X1:time']])

      P.X <- obj$output$`X1`[, c('pvalue', 'padj')]
      colnames(P.X) <- c('pvalue.X','p.adj.fdr.X')
      P.time <- obj$output$time[, c('pvalue', 'padj')]
      colnames(P.time) <- c('pvalue.T','p.adj.fdr.T')
      P.interaction <- obj$output$`X1:time`[, c('pvalue', 'padj')]
      colnames(P.interaction) <- c('pvalue.XT','p.adj.fdr.XT')

      P <- merge(P.X, P.time, by =0, all = T) %>% column_to_rownames('Row.names') %>% merge(P.interaction, by = 0, all=T)%>% column_to_rownames('Row.names')
    }
  }
  return(P)
}


rdirichlet.m <- function (alpha) {
  Gam <- matrix(rgamma(length(alpha), shape = alpha), nrow(alpha), ncol(alpha))
  t(t(Gam) / colSums(Gam))
}



#' @title estimate the true proportion of ref.otu.tab
#' @param ref.otu.tab a matrix, the reference OTU count table (row - taxa, column - samples), serving as the template for synthetic sample generation.
#' @return A list with the elements
#' \item{mu}{true proportion}
#' \item{ref.otu.tab}{absolute abundance table}
#' @rdname EstPara
#' @import dirmult
#' @import stats
#' @export

EstPara <- function (ref.otu.tab) {

  if (is.null(rownames(ref.otu.tab))) {
    rownames(ref.otu.tab) <- paste0('OTU', 1 : nrow(ref.otu.tab))
  } # otu * sample
  samplenames <- colnames(ref.otu.tab)
  taxnames <- rownames(ref.otu.tab)

  dirmult.paras <- dirmult(t(ref.otu.tab))

  gamma <- dirmult.paras$gamma
  names(gamma) <- names(dirmult.paras$pi)

  # Add pseduo count(each OTU add gamma estimated from dirmult)
  ref.otu.tab <- sapply(1:ncol(ref.otu.tab), function (i) gamma + ref.otu.tab[,i]) # C_ij otu * sample

  # back to dirchlet, calculate the true proportion
  ref.otu.tab.p <- rdirichlet.m(ref.otu.tab) # P_ij nOTU*nSam
  colnames(ref.otu.tab.p) <- samplenames
  rownames(ref.otu.tab.p) <- taxnames

  # order OTUs by mean OTU proportion, for later selection
  ord <- order(rowMeans(ref.otu.tab.p), decreasing = TRUE)
  ref.otu.tab.p <-  ref.otu.tab.p[ord,]

  # apply size factor
  Si <- exp(rnorm(ncol(ref.otu.tab.p)))
  ref.otu.tab0 <- t(t(ref.otu.tab.p)*Si)
  colnames(ref.otu.tab0) <- colnames(ref.otu.tab.p)
  return(list(mu = ref.otu.tab.p, ref.otu.tab = ref.otu.tab0))
}

EstPara_abs <- function (feature.dat, gamma) {

  if (is.null(rownames(feature.dat))) {
    rownames(feature.dat) <- paste0('OTU', 1 : nrow(feature.dat))
  } # otu * sample
  samplenames <- colnames(feature.dat)
  taxnames <- rownames(feature.dat)

  # Add pseduo count(each OTU add gamma estimated from dirmult)
  ref.otu.tab <- sapply(1:ncol(feature.dat), function (i) gamma + feature.dat[,i]) # C_ij otu * sample

  # back to dirchlet, calculate the true proportion
  ref.otu.tab.p <- rdirichlet.m(ref.otu.tab) # P_ij nOTU*nSam
  colnames(ref.otu.tab.p) <- samplenames
  rownames(ref.otu.tab.p) <- taxnames

  # order OTUs by mean OTU proportion, for later selection
  ord <- order(rowMeans(ref.otu.tab.p), decreasing = TRUE)
  ref.otu.tab.p <-  ref.otu.tab.p[ord,]

  # apply size factor
  Si <- exp(rnorm(ncol(ref.otu.tab.p)))
  ref.otu.tab0 <- t(t(ref.otu.tab.p)*Si)
  colnames(ref.otu.tab0) <- colnames(ref.otu.tab.p)

  return(ref.otu.tab0)
}


SimulateMSeqU <- function (ref.otu.tab, model.paras,
                           nSam, nOTU, diff.otu.pct = 0.1,
                           diff.otu.direct = c("balanced", "unbalanced"),
                           diff.otu.mode = 1,
                           covariate.type = c("binary", "continuous"),
                           grp.ratio = 0.5,
                           covariate.eff.min = 0, covariate.eff.max = 1,
                           confounder = c("no", "yes"),
                           confounder.type = c("none", "binary", "continuous", "both"),
                           conf.cov.cor = 0.6, conf.diff.otu.pct = 0.05, conf.nondiff.otu.pct = 0.1,
                           confounder.eff.min = 0, confounder.eff.max = 1,  # for uniform distribution range
                           error.sd = 0,
                           depth.mu = 10000, depth.sd = 4000, depth.conf.factor = 0) {

  diff.otu.direct <- match.arg(diff.otu.direct)
  covariate.type <- match.arg(covariate.type)
  confounder <- match.arg(confounder)
  confounder.type <- match.arg(confounder.type)

  ref.otu.tab0 <- ref.otu.tab # raw otu table
  # ref.otu.tab <- EstPara_abs(feature.dat = ref.otu.tab, gamma = gamma) # absolute abundance table
  ref.otu.tab <- model.paras$ref.otu.tab #EstPara(feature.dat = ref.otu.tab)
  sample.names <- colnames(ref.otu.tab)
  ref.otu.tab <- ref.otu.tab[(1:(nOTU)), ]
  idx.sample <- sample(sample.names, nSam, replace = F)
  idx.nonsample <- colnames(ref.otu.tab)[!(colnames(ref.otu.tab) %in% idx.sample)]
  ref.otu.tab <- ref.otu.tab[, idx.sample]
  ## filter with prevalence for avoiding rare taxa =0 issue
  idx.nozero <- rownames(ref.otu.tab0)[rowSums(ref.otu.tab0[,idx.sample]>0) > 1]
  ref.otu.tab <- ref.otu.tab[idx.nozero,]
  nOTU <- nrow(ref.otu.tab) # update on nOTU as 0 OTUs were filtered out, only apply to power estimation, not for other simulation
  idx.otu <- rownames(ref.otu.tab)
  if(confounder =='no'){
    confounder.type <- "continuous"
    confounder.eff.mean <- 0
    confounder.eff.sd <- 0
    Z <- cbind(rnorm(nSam))
  }

  if(confounder =='yes'){
    confounder.type <- "continuous"
    if (confounder.type == "continuous")
      Z <- cbind(rnorm(nSam))
    if (confounder.type == "binary")
      Z <- cbind(c(rep(0, nSam%/%2), rep(1, nSam - nSam%/%2)))
    if (confounder.type == "both")
      Z <- cbind(rnorm(nSam), c(rep(0, nSam%/%2), rep(1, nSam - nSam%/%2)))
  }

  rho <- sqrt(conf.cov.cor^2/(1 - conf.cov.cor^2))
  X <- rho * scale(scale(Z) %*% rep(1, ncol(Z))) + rnorm(nSam)
  if(covariate.type == "binary") X <- cbind(ifelse(X <= quantile(X, grp.ratio), 0, 1))
  rownames(X) <- colnames(ref.otu.tab)

  if(diff.otu.direct == "balanced") eta.diff <- sample(c(runif(floor(nOTU/2), min = covariate.eff.min, max = covariate.eff.max), runif(nOTU - floor(nOTU/2), min = -covariate.eff.max, max = -covariate.eff.min))) %*% t(scale(X))
  if(diff.otu.direct == "unbalanced") eta.diff <- sample(runif(nOTU, min = covariate.eff.min, max = covariate.eff.max)) %*% t(scale(X))

  # only balanced for confounder eff
  eta.conf <- sample(c(runif(floor(nOTU/2), min = confounder.eff.min, max = confounder.eff.max),
                       runif(nOTU - floor(nOTU/2), min = -confounder.eff.max, max = -confounder.eff.min))) %*% t(scale(scale(Z) %*% rep(1, ncol(Z))))

  otu.ord <- 1:(nOTU)
  diff.otu.ind <- NULL
  diff.otu.num <- round(diff.otu.pct * nOTU)
  if(diff.otu.mode == 1) {
    diff.ind <- sample(otu.ord, diff.otu.num)
    diff.otu.ind <- c(diff.otu.ind, diff.ind)
  }


  if(diff.otu.mode < 0.5){
    diff.ind <- sample(otu.ord[1:round(length(otu.ord)* diff.otu.mode)], diff.otu.num)
    diff.otu.ind <- c(diff.otu.ind, diff.ind)
  }

  if(diff.otu.mode > 0.5) {
    diff.ind <- sample(otu.ord[round(diff.otu.mode * length(otu.ord)):length(otu.ord)], diff.otu.num)
    diff.otu.ind <- c(diff.otu.ind, diff.ind)
  }

  if(length(diff.otu.ind) >= round(nOTU * conf.diff.otu.pct)) {
    conf.otu.ind1 <- sample(diff.otu.ind, round(nOTU * conf.diff.otu.pct))
  }else {
    conf.otu.ind1 <- diff.otu.ind
  }

  conf.otu.ind <- c(conf.otu.ind1, sample(setdiff(1:(nOTU),diff.otu.ind), round(conf.nondiff.otu.pct * nOTU)))
  eta.diff[setdiff(1:(nOTU), diff.otu.ind), ] <- 0
  eta.conf[setdiff(1:(nOTU), conf.otu.ind), ] <- 0
  eta.error <- matrix(rnorm(nOTU * nSam, 0, error.sd), nOTU, nSam)
  eta.exp <- exp(t(eta.diff + eta.conf + eta.error))
  eta.exp <- eta.exp * t(ref.otu.tab)
  ref.otu.tab.prop <- eta.exp/rowSums(eta.exp)
  ref.otu.tab.prop <- t(ref.otu.tab.prop)
  depth.theta <- (depth.mu^2) / (depth.sd^2 - depth.mu)
  nSeq <- rnegbin(nSam, mu = depth.mu * exp(scale(X) * depth.conf.factor), theta = depth.theta)
  otu.tab.sim <- sapply(1:ncol(ref.otu.tab.prop), function(i) rmultinom(1, nSeq[i], ref.otu.tab.prop[, i]))
  colnames(otu.tab.sim) <- rownames(eta.exp)
  rownames(otu.tab.sim) <- rownames(ref.otu.tab)
  diff.otu.ind = (1:nOTU) %in% diff.otu.ind
  conf.otu.ind = (1:nOTU) %in% conf.otu.ind
  return(list(otu.tab.sim = otu.tab.sim, covariate = X, confounder = Z,
              diff.otu.ind = diff.otu.ind, otu.names = idx.otu, conf.otu.ind = conf.otu.ind))
}



SimulateMSeqCU <- function (ref.otu.tab, model.paras,#gamma,
                            nSubject = 40, nOTU = 50, nTime = 2, error.sd = 1,
                            MgX.min = 0, MgX.max = 1,
                            X.diff.otu.pct = 0.1, grp.ratio = 1,
                            diff.otu.mode = 1,
                            balanced.X = TRUE, SbT = 0, T.diff.otu.pct = 0,
                            MgT.min = 0, MgT.max = 1,
                            balanced.T = TRUE, XT.diff.otu.pct = 0,
                            MgXT.min = 0, MgXT.max = 0,
                            balanced.XT = TRUE, conf.cov.cor = 0.6, confounder = c("none","X", "T"),
                            MgZ.min = 0, MgZ.max = 0,
                            Z.diff.otu.pct = 0.05, Z.nondiff.otu.pct = 0.1,
                            depth.mu = 10000, depth.sd = 4000, depth.conf.factor = 0) {

  confounder <- match.arg(confounder)
  if (confounder == "X") {
    if (Z.diff.otu.pct > X.diff.otu.pct)
      stop("Z.diff.otu.pct should not be larger than X.diff.otu.pct!\n")
  }
  if (confounder == "T") {
    if (Z.diff.otu.pct > T.diff.otu.pct)
      stop("Z.diff.otu.pct should not be larger than T.diff.otu.pct!\n")
  }

  ref.otu.tab <- model.paras$ref.otu.tab #EstPara_abs(feature.dat = ref.otu.tab, gamma = gamma) # absolute abundance table
  sample.names <- colnames(ref.otu.tab)
  ref.otu.tab <- ref.otu.tab[(1:(nOTU)), ]

  idx.sample <- sample(sample.names, nSubject, replace = T)
  OTU <- model.paras$ref.otu.tab[, idx.sample, drop = F]
  otu.names <- names(sort(rowMeans(OTU), decreasing = T))[1:nOTU]
  otu_tab <- OTU[otu.names, , drop = F]
  ## filter with prevalence for avoiding rare taxa =0 issue
  idx.nozero <- rowSums(ref.otu.tab[,idx.sample, drop = F]>0) > 1
  otu_tab <- otu_tab[idx.nozero,]
  nOTU <- nrow(otu_tab) # update on nOTU as 0 OTUs were filtered out, only apply to power estimation, not for other simulation
  idx.otu <- rownames(otu_tab)

  error.mean <- 0
  otu.tab <- otu_tab[, rep(1:nSubject, each = nTime)]
  colnames(otu.tab) <- paste0(paste0("rep", 1:nTime), "_", rep(paste0("Subject", 1:nSubject), each = nTime))
  nSam <- ncol(otu.tab)
  SubjectID <- as.numeric(gsub(".*Subject", "", colnames(otu.tab)))
  time <- as.numeric(as.factor(gsub("_.*", "", colnames(otu.tab)))) - 1
  X <- as.vector(ifelse(SubjectID <= quantile(unique(SubjectID), grp.ratio/(1 + grp.ratio)), 0, 1))
  rho <- sqrt(conf.cov.cor^2/(1 - conf.cov.cor^2))
  if (confounder == "X") {
    Z <- as.vector(rho * scale(X) + rnorm(nSam))
  }else if(confounder=="T"){
    Z <- as.vector(rho * scale(time) + rnorm(nSam))
  }

  if(confounder %in% c("X","T")){
    meta <- cbind.data.frame(X, Z, SubjectID, time)
  }else{
    meta <- cbind.data.frame(X, SubjectID, time)
  }

  rownames(meta) <- colnames(otu.tab)
  error <- replicate(nSam, rnorm(nOTU, error.mean, error.sd))
  otu.ord <- 1:nOTU
  X.diff.otu.ind <- T.diff.otu.ind <- XT.diff.otu.ind <- NULL
  X.diff.otu.num <- round(X.diff.otu.pct * nOTU)
  T.diff.otu.num <- round(T.diff.otu.pct * nOTU)

  if(diff.otu.mode <0.5){
    X.diff.otu.ind.tmp <- sample(otu.ord[1:round(length(otu.ord) * diff.otu.mode)], X.diff.otu.num)
    T.diff.otu.ind.tmp <- sample(otu.ord[1:round(length(otu.ord) * diff.otu.mode)], T.diff.otu.num)
  }

  if(diff.otu.mode >0.5){
    X.diff.otu.ind.tmp <- sample(otu.ord[round(length(otu.ord) * diff.otu.mode):length(otu.ord)], X.diff.otu.num)
    T.diff.otu.ind.tmp <- sample(otu.ord[round(length(otu.ord) * length(otu.ord)):length(otu.ord)], T.diff.otu.num)
  }

  if(diff.otu.mode ==1){
    X.diff.otu.ind.tmp <- sample(otu.ord, X.diff.otu.num)
    T.diff.otu.ind.tmp <- sample(otu.ord, T.diff.otu.num)
  }

  X.diff.otu.ind <- c(X.diff.otu.ind, X.diff.otu.ind.tmp)
  T.diff.otu.ind <- c(T.diff.otu.ind, T.diff.otu.ind.tmp)

  if(confounder %in% c("X","T")){
    Z.diff.otu.num <- round(Z.diff.otu.pct * nOTU)
    Z.nondiff.otu.num <- round(Z.nondiff.otu.pct * nOTU)
  }

  if (confounder == "T") {
    Z.diff.otu.ind <- c(sample(T.diff.otu.ind, Z.diff.otu.num),
                        sample(setdiff(otu.ord, T.diff.otu.ind), Z.nondiff.otu.num))
  }else if(confounder == "X") {
    Z.diff.otu.ind <- c(sample(X.diff.otu.ind, Z.diff.otu.num),
                        sample(setdiff(otu.ord, X.diff.otu.ind), Z.nondiff.otu.num))
  }

  if (balanced.X) {
    coef.X <- sample(c(runif(floor(nOTU/2), min = MgX.min, max = MgX.max), runif(nOTU - floor(nOTU/2), min = -MgX.max, max = -MgX.min)))
  }else {
    coef.X <- runif(nOTU, min = 0, max = MgX.max)
  }

  coef.X[setdiff(otu.ord, X.diff.otu.ind)] <- 0
  eta.X <- coef.X %*% t(scale(X))

  if (balanced.T) {
    coef.T <- sample(c(runif(floor(nOTU/2), min = MgT.min, max = MgT.max), runif(nOTU - floor(nOTU/2), min = -MgT.max, max = -MgT.min)))
  }else {
    coef.T <- runif(nOTU, min = 0, max = MgT.max)
  }

  coef.T[setdiff(otu.ord, T.diff.otu.ind)] <- 0
  coef.T.Subject <- replicate(nSubject, rnorm(nOTU, coef.T, SbT))
  coef.T <- matrix(data = apply(coef.T.Subject, 2, function(x) rep(x, nTime)), ncol = ncol(coef.T.Subject) * nTime)
  eta.T <- t(t(coef.T) * as.vector(scale(time)))

  if (balanced.XT) {
    coef.XT <- sample(c(runif(floor(nOTU/2), min = MgXT.min, max = MgXT.max), runif(nOTU - floor(nOTU/2), min = -MgXT.max, max = -MgXT.min)))
  }else {
    coef.XT <- runif(nOTU, min = 0, max = MgXT.max)
  }

  coef.XT[setdiff(otu.ord, XT.diff.otu.ind)] <- 0
  coef.XT.Subject <- replicate(nSubject, rnorm(nOTU, coef.XT, 0))
  coef.XT <- matrix(data = apply(coef.XT.Subject, 2, function(x) rep(x, nTime)), ncol = ncol(coef.XT.Subject) * nTime)
  eta.XT <- t(t(coef.XT) * as.vector(scale(time) * scale(X)))

  if(confounder %in% c("X","T")){
    coef.Z <- sample(c(runif(floor(nOTU/2), min = MgZ.min, max = MgZ.max), runif(nOTU - floor(nOTU/2), min = -MgZ.max, max = -MgZ.min)))
    coef.Z[setdiff(otu.ord, Z.diff.otu.ind)] <- 0
    eta.Z <- coef.Z %*% t(scale(Z))
    eta.exp <- exp(t(eta.X + eta.Z + eta.T + eta.XT + error))
  }else{
    eta.exp <- exp(t(eta.X + eta.T + eta.XT + error))
  }

  eta.exp <- eta.exp * t(otu.tab)
  otu.tab.prop <- t(eta.exp/rowSums(eta.exp))
  # nSeq <- rnorm(ncol(otu.tab.prop), mean = depth.mu * exp(scale(X) * depth.conf.factor), sd = depth.theta)
  depth.theta <- (depth.mu^2) / (depth.sd^2 - depth.mu)
  nSeq <- rnegbin(ncol(otu.tab.prop), mu = depth.mu * exp(scale(X) * depth.conf.factor), theta = depth.theta)
  otu.tab.sim <- sapply(1:ncol(otu.tab.prop), function(i) rmultinom(1, nSeq[i], otu.tab.prop[, i]))
  colnames(otu.tab.sim) <- rownames(eta.exp)
  rownames(otu.tab.sim) <- colnames(eta.exp)
  otu.names <- colnames(eta.exp)
  X.diff.otu.ind = (1:nOTU) %in% X.diff.otu.ind
  T.diff.otu.ind = (1:nOTU) %in% T.diff.otu.ind
  XT.diff.otu.ind = (1:nOTU) %in% XT.diff.otu.ind
  if(confounder %in% c("X","T")){
    Z.diff.otu.ind = (1:nOTU) %in% Z.diff.otu.ind
  }else{
    Z.diff.otu.ind = NULL
  }
  return(list(otu.tab.sim = otu.tab.sim, meta = meta, otu.names = otu.names,
              X.diff.otu.ind = X.diff.otu.ind, T.diff.otu.ind = T.diff.otu.ind,
              XT.diff.otu.ind = XT.diff.otu.ind, Z.diff.otu.ind = Z.diff.otu.ind))
}


data_summary <- function(data, formula){
  ## formula: value ~ variables to be aggregated
  error.info <- aggregate(as.formula(formula), data, function(x) mean(is.na(x)))
  m <- aggregate(as.formula(formula), data, function(x) mean(x[!is.na(x)]))
  med <- aggregate(as.formula(formula), data, function(x) median(x[!is.na(x)]))
  se <- aggregate(as.formula(formula), data, function(x) {
    ind <- !is.na(x)
    sd(x[ind]) / sqrt(length(x[ind]))})
  sd <- aggregate(as.formula(formula), data, function(x) {
    ind <- !is.na(x)
    sd(x[ind])})
  ymin <- m[, ncol(m)] - 1.96 * se[, ncol(m)]
  ymin <- ifelse(ymin > 0, ymin, 0)
  ymax <- m[, ncol(m)] + 1.96 * se[, ncol(m)]
  sum <- cbind(m, SD = sd[, ncol(sd)], ymax=ymax, ymin=ymin)
  return(sum)
}

cal.tpr.fdr <- function(p.adj.fdr, truth, cutoff){
  TP <- sum(p.adj.fdr <= cutoff & truth, na.rm = T)
  FP <- sum(p.adj.fdr <= cutoff & !truth, na.rm = T)
  FN <- sum(p.adj.fdr > cutoff & truth, na.rm = T)
  TPR <- ifelse((TP + FN)==0, 0, TP / (TP + FN))
  FDR <- ifelse((FP + TP)==0, 0, FP / (FP + TP))
  res <- c(TPR, FDR)
  names(res) <- c('TPR','FDR')
  return(res)
}

list_to_df <- function(list.df){
  if(length(list.df)>1){
    df <- as.data.frame(list.df[[1]])
    df$grp <- names(list.df)[1]
    for(i in 2:length(list.df)){
      df.tmp <- as.data.frame(list.df[[i]])
      df.tmp$grp <- names(list.df)[i]
      df <- rbind(df, df.tmp)
    }
  }else{
    df <- as.data.frame(list.df[[1]])
    df$grp <- names(list.df)[1]
  }
  return(df)
}

list_vector2df <- function(ListOfVector){
  df <- NULL
  for(i in 1:length(ListOfVector)){df <- cbind(df, ListOfVector[[i]])}
  colnames(df) <- names(ListOfVector)
  return(as.data.frame(df))
}


generate_plot3 <- function(data, effect.size, ylab, error.bar =T, R2=F, matched.pair=F){
  if(matched.pair==T){
    xlab = 'Subject number'
  }else{
    xlab = 'Sample size'
  }
  if(length(unique(data$nSam))==1 & length(unique(data[,effect.size]))>1){
    if(matched.pair==T){
      data$nSam <- paste0('Subject number=', data$nSam)
    }else{
      data$nSam <- paste0('Sample size=', data$nSam)
    }

    if(R2==T){
      # cat('--- category 1---\n')
      plt <- ggplot(data, aes(x=!!as.name(effect.size), y = value, fill = nSam)) +
        geom_boxplot(size = 1.3, outlier.shape = NA)+
        geom_jitter(shape = 21, size = 4, width = 0.2) +
        facet_grid(~nSam) +
        scale_fill_brewer(palette = 'Accent') +
        scale_color_brewer(palette = 'Accent') +
        theme_classic() +
        theme(axis.title = element_text(color = 'black', size =24),
              axis.text = element_text(color = 'black', size = 24),
              legend.title = element_text(color = 'black', size =24),
              legend.text = element_text(color = 'black', size =24),
              strip.text = element_text(color = 'black', size =24),
              legend.position = 'none') +
        labs(x = 'max log2 fold change', y =ylab, fill = '', color = '')
    }else{
      # cat('--- category 2---\n')
      plt <- ggplot(data, aes(x=!!as.name(effect.size), y = value)) +
        geom_line(aes(group = variable, color = variable), size = 1.3)+
        geom_point(shape = 21, aes(fill = variable), size = 4) +
        facet_grid(~nSam) +
        scale_fill_brewer(palette = 'Set1') +
        scale_color_brewer(palette = 'Set1') +
        theme_classic() +
        theme(axis.title = element_text(color = 'black', size =24),
              axis.text = element_text(color = 'black', size = 24),
              legend.title = element_text(color = 'black', size =24),
              legend.text = element_text(color = 'black', size =24),
              strip.text = element_text(color = 'black', size =24),
              legend.position = 'none') +
        labs(x = 'max log2 fold change', y = ylab, fill = '', color = '')
      if(error.bar ==T){
        plt <- plt + geom_errorbar(mapping=aes(ymin=ymin, ymax=ymax,color = variable), width=0.07, size=0.5)
      }
    }
  }else if(length(unique(data$nSam))>1 & length(unique(data[,effect.size]))==1){
    data[,effect.size] <- paste0('max log2 fold change=', data[,effect.size])

    if(R2==T){
      # cat('--- category 3---\n')
      plt <- ggplot(data, aes(x=nSam, y = value, fill = nSam)) +
        geom_boxplot(size = 1.3, outlier.shape = NA)+
        geom_jitter(shape = 21, size = 4, width = 0.2) +
        facet_grid(as.formula(paste0('~',effect.size))) +
        scale_fill_brewer(palette = 'Accent') +
        scale_color_brewer(palette = 'Accent') +
        theme_classic() +
        theme(axis.title = element_text(color = 'black', size =24),
              axis.text = element_text(color = 'black', size = 24),
              legend.title = element_text(color = 'black', size =24),
              legend.text = element_text(color = 'black', size =24),
              strip.text = element_text(color = 'black', size =24),
              legend.position = 'none') +
        labs(x =xlab, y =ylab, fill = '', color = '')
    }else{
      # cat('--- category 4---\n')
      plt <- ggplot(data, aes(x=nSam, y = value)) +
        geom_line(aes(group = variable, color = variable), size = 1.3)+
        geom_point(shape = 21, aes(fill = variable), size = 4) +
        facet_grid(as.formula(paste0('~',effect.size))) +
        scale_fill_brewer(palette = 'Set1') +
        scale_color_brewer(palette = 'Set1') +
        theme_classic() +
        theme(axis.title = element_text(color = 'black', size =24),
              axis.text = element_text(color = 'black', size =24),
              legend.title = element_text(color = 'black', size =24),
              legend.text = element_text(color = 'black', size =24),
              strip.text = element_text(color = 'black', size =24),
              legend.position = 'none') +
        labs(x = xlab, y = ylab, fill = '', color = '')
      if(error.bar ==T){
        plt <- plt + geom_errorbar(mapping=aes(ymin=ymin, ymax=ymax,color = variable),width=0.07, size=0.5)
      }
    }

  }else if(length(unique(data$nSam))==1 & length(unique(data$effect.size))==1){
    data[,effect.size] <- paste0('max log2 fold change=', data[,effect.size])
    if(matched.pair==T){
      data$nSam <- paste0('Subject number=', data$nSam)
    }else{
      data$nSam <- paste0('Sample size=', data$nSam)
    }
    if(R2==T){
      # cat('--- category 5---\n')
      plt <- ggplot(data, aes(x=nSam, y = value, fill = nSam)) +
        geom_boxplot(size = 1.3, outlier.shape = NA)+
        geom_jitter(shape = 21, size = 4, width = 0.2) +
        facet_grid(as.formula(paste0('~',effect.size))) +
        scale_fill_brewer(palette = 'Accent') +
        scale_color_brewer(palette = 'Accent') +
        theme_classic() +
        theme(axis.title = element_text(color = 'black', size =24),
              axis.text = element_text(color = 'black', size = 24),
              legend.title = element_text(color = 'black', size =24),
              legend.text = element_text(color = 'black', size =24),
              strip.text = element_text(color = 'black', size =24),
              legend.position = 'none') +
        labs(x = xlab, y =ylab, fill = '', color = '')
    }else{
      # cat('--- category 6---\n')
      plt <- ggplot(data, aes(x=nSam, y = value, fill = nSam)) +
        geom_bar(position="dodge",stat = 'identity', aes(fill = variable)) +
        facet_grid(as.formula(paste0('~',effect.size))) +
        scale_fill_brewer(palette = 'Accent') +
        scale_color_brewer(palette = 'Accent') +
        theme_classic() +
        theme(axis.title = element_text(color = 'black', size =24),
              axis.text = element_text(color = 'black', size =24),
              legend.title = element_text(color = 'black', size =24),
              legend.text = element_text(color = 'black', size =24),
              strip.text = element_text(color = 'black', size =24),
              legend.position = 'none') +
        labs(x = xlab, y = ylab, fill = '', color = '')
      if(error.bar ==T){
        plt <- plt + geom_errorbar(mapping=aes(ymin=ymin, ymax=ymax,color = variable),width=0.07, size=0.5,position=position_dodge(.9))
      }
    }
  }else{
    # cat('--- category 7---\n')
    plt <- ggplot(data, aes(x=nSam, y = value)) +
      geom_line(aes(group = variable, color = variable), size = 1.3)+
      geom_point(shape = 21, aes(fill = variable), size = 4) +
      facet_grid(as.formula(paste0('~',effect.size))) +
      scale_fill_brewer(palette = 'Set1') +
      scale_color_brewer(palette = 'Set1') +
      theme_classic() +
      theme(axis.title = element_text(color = 'black', size =24),
            axis.text = element_text(color = 'black', size =24),
            legend.title = element_text(color = 'black', size =24),
            legend.text = element_text(color = 'black', size =24),
            strip.text = element_text(color = 'black', size =24),
            legend.position = 'none') +
      labs(x = xlab, y = ylab, fill = '', color = '')
    if(error.bar ==T){
      plt <- plt + geom_errorbar(mapping=aes(ymin=ymin, ymax=ymax,color = variable), width=0.07, size=0.5)
    }
  }
  return(plt)
}

