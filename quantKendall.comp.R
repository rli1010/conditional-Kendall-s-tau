

#################################################################################
################# conditional Kendall's tau for complete data ###################
#################################################################################
K.est.comp <- function(data, formula1, formula2, tauseq = seq(0.01, 0.99, 0.01), ...) {
  mf1 <- model.frame(formula1, data)    ## prepare for extracting y
  mf2 <- model.frame(formula2, data)
  Y1 <- model.response(mf1)
  Y2 <- model.response(mf2)

  fit1 <- quantreg::rq(formula1, tau = tauseq, data, ...)     ## quantile reg for y1
  fit2 <- quantreg::rq(formula2, tau = tauseq, data, ...)
  qmat1 <- cbind(-Inf, fit1$fitted.values, Inf)             ## quantile(0) = -Inf; quantile(1) = Inf
  qmat2 <- cbind(-Inf, fit2$fitted.values, Inf)

  lt <- length(tauseq) + 2            ## length of tauseq + 2
  n.da <- nrow(data)
  Y1mat <- matrix(rep(Y1, lt), n.da)
  Y2mat <- matrix(rep(Y2, lt), n.da)
  compare1.mat <- (Y1mat <= qmat1)
  compare2.mat <- (Y2mat <= qmat2)

  Cmat <- matrix(0, lt, lt)
  Cmat = .Fortran("cmatfortran", Cmat = matrix(as.double(Cmat), lt, lt),
                  mat1 = as.double(compare1.mat),
                  mat2 = as.double(compare2.mat),
                  m = as.integer(n.da),
                  n = as.integer(lt))[[1]]

  dCmat <- Cmat[-1, -1] - Cmat[-1, -lt] - Cmat[-lt, -1] + Cmat[-lt, -lt]
  C.est <- (Cmat[-1, -1] + Cmat[-1, -lt] + Cmat[-lt, -1] + Cmat[-lt, -lt])/4
  K.est <- 4 * sum(C.est * dCmat) - 1
  result <- list(Est = K.est)
  result
}


fisher <- function(x) {0.5*log((1+x)/(1-x))} # fisher transformation

quantKendall.comp <- function(data, formula1, formula2, tauseq = seq(0.01, 0.99, 0.01), n.bt = 200, ...){
  if (!requireNamespace("quantreg", quietly = TRUE)) {install.packages("quantreg")}  ## install "quantreg"
  call <- match.call()
  n.da <- nrow(data)
  result <- K.est.comp(data, formula1, formula2, tauseq, ...)
  K.bt <- rep(NA, n.bt)
  for(i in 1:n.bt){
    idx <- sample(1:n.da, replace=TRUE)
    K.bt[i] <- K.est.comp(data[idx, ], formula1, formula2, tauseq, ...)$Est
  }
  ztest <- fisher(result$Est)/sd(fisher(K.bt))
  result$call <- call
  result$BTsample <- K.bt
  result$p.value <- 2*(1 - pnorm(abs(ztest)))
  result
}
