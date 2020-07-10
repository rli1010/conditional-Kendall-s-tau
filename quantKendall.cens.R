

quantKendall.cens <- function(data, formula1, formula2, upper=NULL, tauseq = seq(0.01, 0.99, 0.01), n.bt = 200, ...){
  if(!requireNamespace("quantreg", quietly = T)) {install.packages("quantreg")}  ## install "quantreg"
  if(!requireNamespace("survival", quietly = T)) {install.packages("survival")}  ## install "quantreg"
  call <- match.call()
  formula1 <- paste0("survival::", format(formula1))  ## save formula as string with "survival::"
  formula2 <- paste0("survival::", format(formula2))
  n.da <- nrow(data)
  result <- K.est.cens(data, formula1, formula2, upper, tauseq, ...)
  K.bt <- rep(NA, n.bt)
  for(i in 1:n.bt){
    idx <- sample(1:n.da,replace=TRUE)
    K.bt[i] <- K.est.cens(data[idx, ], formula1, formula2, upper, tauseq, Boot=TRUE, ...)$Est
  }
  ztest <- fisher(result$Est)/sd(fisher(K.bt))
  result$BTsample <- K.bt
  result$p.value <- 2*(1 - pnorm(abs(ztest)))
  class(result) <- "quantKendall"
  result
}

fisher <- function(x) {0.5*log((1+x)/(1-x))} # fisher transformation

#################################################################################
################# conditional Kendall's tau for censored data ###################
#################################################################################
K.est.cens <- function(data, formula1, formula2, upper, tauseq = seq(0.01, 0.99, 0.01), Boot=FALSE, ...) {
  mf1 <- model.frame(formula1, data)    ## prepare for extracting y
  mf2 <- model.frame(formula2, data)
  Y1 <- model.response(mf1)             ## for survival data, 2 columns (T, delta)
  Y2 <- model.response(mf2)
  Zmat1 <- model.matrix(as.formula(formula1), data)  ## model.matrix doesn't work with string formula
  Zmat2 <- model.matrix(as.formula(formula2), data)
  n.b.tau <- dim(Zmat1)[2]

  f1 <- quantreg::crq(formula1, grid = tauseq, data = data, method = "PengHuang", ...)     ## quantile reg for y1
  f1$sol <- f1$sol[, complete.cases(t(f1$sol))]  ## coefs may be NaN, remove columns with NaN
  f2 <- quantreg::crq(formula2, grid = tauseq, data = data, method = "PengHuang", ...)
  f2$sol <- f2$sol[, complete.cases(t(f2$sol))]
  if(is.null(upper)){
    upper <- c(max(f1$sol[1, ]), max(f2$sol[1, ]))
  } else {
    if(length(upper) == 1) upper <- rep(upper, 2)
    if(upper[1] > max(f1$sol[1, ])) {
      upper[1] <- max(f1$sol[1, ])
      if(!Boot){warning(paste0("Too large upper1, ", upper[1], " is used instead"))}
    }
    if(upper[2] > max(f2$sol[1, ])) {
      upper[2] <- max(f2$sol[1, ])
      if(!Boot){warning(paste0("Too large upper2, ", upper[2], " is used instead"))}
    }
  }
  taus <- min(which(f1$sol[1, ] >= upper[1]))    ## to select taus within (0, upper1)
  fit1 <- list(tauseq = f1$sol[1, 1:taus], b.tau = f1$sol[2:(n.b.tau + 1), 1:taus])
  pred.q1 = cbind(-Inf, monotone.respect(fit1, Zmat1))
  taus <- min(which(f2$sol[1, ] >= upper[2]))
  fit2 <- list(tauseq = f2$sol[1, 1:taus], b.tau = f2$sol[2:(n.b.tau + 1), 1:taus])
  pred.q2 = cbind(-Inf, monotone.respect(fit2, Zmat2))

  obj1 = compare.func.new(Y1[,1], pred.q1, Y1[,2], fit1$tauseq)
  obj2 = compare.func.new(Y2[,1], pred.q2, Y2[,2], fit2$tauseq)
  dtime.1 = obj1$dtime
  dtime.2 = obj2$dtime

  l = as.integer(length(dtime.1))
  m = as.integer(length(dtime.2))
  n = length(Y1)
  B.mat = matrix(0, l, m)
  eta1.dt.seq = obj1$eta.d
  eta2.dt.seq = obj2$eta.d

  atrisk.1 = obj1$atrisk
  atrisk.2 = obj2$atrisk
  M1.dt.mat = t(t(obj1$dindi) - t(atrisk.1) * eta1.dt.seq)
  M2.dt.mat = t(t(obj2$dindi) - t(atrisk.2) * eta2.dt.seq)

  #calculate B matrix per eqn (11)
  B.mat = .Fortran("pc", B = matrix(as.double(B.mat), l, m), l = l, m = m, n = n,
    atrisk1 = matrix(as.integer(obj1$atrisk), n, l),
    atrisk2 = matrix(as.integer(obj2$atrisk), n, m),
    d1 = as.double(eta1.dt.seq),
    d2 = as.double(eta2.dt.seq),
    M1 = matrix(as.double(M1.dt.mat), n, l),
    M2 = matrix(as.double(M2.dt.mat), n, m))[[1]]

  B.mat = matrix(B.mat, l, m)
  B.mat[is.na(B.mat)] = 0

  #calculate the L matrix based on eqn (9)
  Obj.mat = matrix(0, l, m)

  Obj.mat[1, ] = 1 + cumsum(B.mat[1, ])
  Obj.mat[, 1] = 1 + cumsum(B.mat[, 1])

  if ((dim(B.mat)[1] > 1) & (dim(B.mat)[2] > 1)) {
    Obj.mat = matrix(.Fortran("objsum", obj = as.double(Obj.mat), l = as.integer(length(dtime.1)),
                              m = as.integer(length(dtime.2)), Bmat = as.double(B.mat))[[1]],
                     length(dtime.1), length(dtime.2))
  }

  F1 <- 1 - dtime.1
  F2 <- 1 - dtime.2
  F12.mat <- F1 %*% t(F2)

  F.mat <- F12.mat * Obj.mat
  F.mat <- rbind(c(1, F2), cbind(F1, F.mat))
  rownames(F.mat) <- paste0("t=", c(0, dtime.1))
  colnames(F.mat) <- paste0("t=", c(0, dtime.2))

  #### Calculate C-index using bivariate survival function ####
  Dim1 <- dim(F.mat)[1]
  Dim2 <- dim(F.mat)[2]
  M <- sum(F.mat[-Dim1, -Dim2] * (F.mat[-1, -1] - F.mat[-Dim1, -1] - F.mat[-1, -Dim2] + F.mat[-Dim1, -Dim2]))
  N <- sum((F.mat[-Dim1, -1] - F.mat[-Dim1, -Dim2]) * (F.mat[-1, -Dim2] - F.mat[-Dim1, -Dim2]))
  return(list(Est = (M - N) / (M + N), fit1 = f1, fit2 = f2, upper = upper,
              upper.model = c(max(f1$sol[1, ]), max(f2$sol[1, ]))))
}

compare.func.new = function(Y, qmat, eta, tauseq) {
  n = length(Y)
  m = dim(qmat)[2]
  Ymat = round(matrix(rep(Y, m - 1), n), 7)
  Yeval = (Ymat <= qmat[, -1]) - (Ymat <= qmat[, -m])
  dindi = Yeval * eta # n*(m-1)
  jump = !(apply(dindi, 2, sum) == 0) # check whether observed for specific quantiles
  dtime = tauseq[jump]
  atrisk = (Ymat > qmat[, -m])

  eta.d = apply(dindi, 2, mean) / apply(atrisk, 2, mean)
  surv.obj <- 1 - apply(dindi[, jump], 2, sum) / apply(atrisk[, jump], 2, sum)
  surv.margin <- cumprod(surv.obj)
  return(list(dindi=dindi[,jump], atrisk=atrisk[,jump], eta.d=eta.d[jump], m=length(dtime),
              dtime=dtime, surv.margin=surv.margin))
}

monotone.respect = function(fit, Zmat) {    ## to make the predicted quantiles monotone
  pred.q = cbind(as.matrix(Zmat) %*% fit$b.tau)
  m = length(fit$tauseq)
  si = floor(m / 2)
  respect.seq = respect.sum = rep(FALSE, m)
  pred.q.curr = pred.q[, si]
  respect.seq[si] = TRUE

  n0 <- dim(Zmat)[1]
  for (indx in (si + 1):m) {
    respect = (pred.q.curr <= pred.q[, indx])
    respect.sum[indx] = sum(respect)
    if (sum(respect) == n0) {
      pred.q.curr = pred.q[, indx]
      respect.seq[indx] = TRUE
    }
  }

  pred.q.curr = pred.q[, si]
  for (indx in (si - 1):1) {
    respect = (pred.q.curr >= pred.q[, indx])
    respect.sum[indx] = sum(respect)
    if (sum(respect) == n0) {
      pred.q.curr = pred.q[, indx]
      respect.seq[indx] = TRUE
    }
  }
  min.respect = min((1:m)[respect.seq])
  max.respect = max((1:m)[respect.seq])
  respect.seq[1:min.respect] = TRUE
  respect.seq[max.respect:m] = TRUE

  tauseq.respect = fit$tauseq[respect.seq]
  if (is.null(dim(fit$b.tau))) {
    b.new1 = approxfun(tauseq.respect, fit$b.tau[respect.seq], "linear", rule = 2)
  } else{
    b.new1 = approxfun(tauseq.respect, fit$b.tau[1, respect.seq], "linear", rule = 2)
  }
  b.tau.new = b.new1(fit$tauseq)
  Zcol.n <- dim(Zmat)[2]
  if (Zcol.n > 1) {
    for (i in 2:Zcol.n) {
      b.newi = approxfun(tauseq.respect, fit$b.tau[i, respect.seq], "linear", rule = 2)
      b.tau.new = rbind(b.tau.new, b.newi(fit$tauseq))
    }
  }
  pred.respect = as.matrix(Zmat) %*% b.tau.new
  return(pred.respect)
}
