pscore <- function(A,X,weights=NULL,stabilized=FALSE){
  if (!is.null(weights)){
    w0 = (1-A)*weights
    w1 = A*weights
    ips = cbind(w0,w1)
    return(ips)
  }
  if (is.null(X)){
    w0 = (1-A)/mean(1-A)
    w1 = A/mean(A)
  } else{
    score = predict(glm(A~X, family='binomial'),type='response')
    w0 = (1-A)/(1-score)
    w1 = A/score
    if (stabilized){
      w0 = w0/mean(w0)*(1-A)
      w1 = w1/mean(w1)*A
    }
  }
  ips = cbind(w0,w1)
  return(ips)
}

matchy <- function(yvec,xvec,newx){
  ivec = sapply(newx, function(x) max(which(xvec<=x)))
  return(yvec[ivec])
}

geteta <- function(x, beta){
  drop(x %*% beta)
}

semi_markov <- function(A,T,dT,R,dR,X=NULL,weights=NULL,a,conf.int=NULL){
  w = pscore(A,X,weights)
  dat = list(A=A,T=T,dT=dT,R=R,dR=dR,X=X,weight=weights)
  if (length(a)==1) a = rep(a,3)
  tseq = sort(unique(c(T[dT==1],R[dR==1])))
  K = length(tseq)
  haz1 = haz2 = haz3 = haz1t = haz2t = haz3t = 0
  F1 = F2 = F3 = F1t = F2t = F3t = F3t.n = 0
  G1.1 = G1.2 = G1.3 = G2.1 = G2.2 = G2.3 = G3 = 0
  for (k in 1:K){
    t0 = tseq[k]
    Y1 = sum(w[,a[1]+1]*(T>=t0)*(R>=t0))
    dhaz1 = sum(w[,a[1]+1]*(T==t0)*dT*(1-dR))/Y1
    if (Y1==0) dhaz1 = 0
    Y2 = sum(w[,a[2]+1]*(T>=t0)*(R>=t0))
    dhaz2 = sum(w[,a[2]+1]*(R==t0)*dR)/Y2
    if (Y2==0) dhaz2 = 0
    Y3 = sum(w[,a[3]+1]*(T>=t0)*(R<=t0)*dR)
    dhaz3 = sum(w[,a[3]+1]*(T==t0)*dT*dR)/Y3
    if (Y3==0) dhaz3 = 0
    F1t = F1t + exp(-haz1t-haz2t)*dhaz1
    F2t = F2t + exp(-haz1t-haz2t)*dhaz2
    F3t.n = F3t.n + exp(-haz1t-haz2t+haz3t)*dhaz2
    F3t = F2t - F3t.n*exp(-haz3t)
    if (!is.null(conf.int)){
    Y1w = sum(w[,a[1]+1]^2*(T>=t0)*(R>=t0))
    Y2w = sum(w[,a[2]+1]^2*(T>=t0)*(R>=t0))
    Y3w = sum(w[,a[3]+1]^2*(T>=t0)*(R<=t0)*dR)
    dG1.1 = Y1w*dhaz1/Y1^2
    dG1.2 = (F2t-F3t)*exp(haz3t)*Y1w*dhaz1/Y1^2
    dG1.3 = (F2t-F3t)^2*exp(haz3t)^2*Y1w*dhaz1/Y1^2
    dG2.1 = Y2w*dhaz2/Y2^2
    dG2.2 = (1-F1t-F3t)*exp(haz3t)*Y2w*dhaz2/Y2^2
    dG2.3 = (1-F1t-F3t)^2*exp(haz3t)^2*Y2w*dhaz2/Y2^2
    dG3 = (F2t-F3t)^2*exp(haz3t)^2*Y3w*dhaz3/Y3^2
    if (Y1==0) dG1.1 = dG1.2 = dG1.3 = 0
    if (Y2==0) dG2.1 = dG2.2 = dG2.3 = 0
    if (Y3==0) dG3 = 0
    G1.1 = append(G1.1, dG1.1)
    G1.2 = append(G1.2, dG1.2)
    G1.3 = append(G1.3, dG1.3)
    G2.1 = append(G2.1, dG2.1)
    G2.2 = append(G2.2, dG2.2)
    G2.3 = append(G2.3, dG2.3)
    G3 = append(G3, dG3)
    }
    haz1t = haz1t + dhaz1
    haz2t = haz2t + dhaz2
    haz3t = haz3t + dhaz3
    haz1 = append(haz1, haz1t)
    haz2 = append(haz2, haz2t)
    haz3 = append(haz3, haz3t)
    F1 = append(F1, F1t)
    F2 = append(F2, F2t)
    F3 = append(F3, F3t)
  }
  if (!is.null(conf.int)){
  G1 = (1-F1-F3)^2*cumsum(G1.1)-2*(1-F1-F3)*exp(-haz3)*cumsum(G1.2)+exp(-haz3)^2*cumsum(G1.3)
  G2 = (1-F1-F3)^2*cumsum(G2.1)-2*(1-F1-F3)*exp(-haz3)*cumsum(G2.2)+exp(-haz3)^2*cumsum(G2.3)
  G3 = exp(-haz3)^2*cumsum(G3)
  G1[is.na(G1)] = 0
  G2[is.na(G2)] = 0
  se = sqrt(G1+G2+G3)
  }
  if (tseq[1]==0){
    tseq = tseq[-1]
    F1 = F1[-1]
    F2 = F2[-1]
    F3 = F3[-1]
    if (!is.null(conf.int)) se = se[-1]
  }
  if (!is.null(conf.int)){
    return(list(dat=dat,time=c(0,tseq),F1=F1,F2=F2,F3=F3,a=a,se=se))
  }
  return(list(dat=dat,time=c(0,tseq),F1=F1,F2=F2,F3=F3,a=a))
}

semi_semimarkov <- function(A,T,dT,R,dR,X=NULL,weights=NULL,a,conf.int=NULL){
  w = pscore(A,X,weights)
  dat = list(A=A,T=T,dT=dT,R=R,dR=dR,X=X,weight=weights)
  if (length(a)==1) a = rep(a,3)
  haz1 = haz2 = haz3 = haz1t = haz2t = haz3t = M3 = 0
  F1 = F2 = F3 = F1t = F2t = F3t = 0
  G1 = G2 = G3 = 0
  M1 = M2 = NULL
  TR = sort(unique((T-R)[dR==1]))
  K = length(TR)
  for (k in 1:K){
    t0 = TR[k]
    Y3 = sum(w[,a[3]+1]*(T-R>=t0)*dR)
    Y3w = sum(w[,a[3]+1]^2*(T-R>=t0)*dR)
    dhaz3 = sum(w[,a[3]+1]*(T-R==t0)*dT*dR)/Y3
    dM3 = Y3w*dhaz3/Y3^2
    if (Y3==0) dhaz3 = dM3 = 0
    haz3 = append(haz3, haz3t + dhaz3)
    haz3t = haz3t + dhaz3
    M3 = append(M3, dM3)
  }
  if (TR[1]==0){
    haz3 = haz3[-1]
    M3 = M3[-1]
  } else {
    TR = c(0,TR)
  }
  M3 = cumsum(M3)
  tseq = sort(unique(c(T[dT==1],R[dR==1])))
  K = length(tseq)
  for (k in 1:K){
    t0 = tseq[k]
    Y1 = sum(w[,a[1]+1]*(T>=t0)*(R>=t0))
    dhaz1 = sum(w[,a[1]+1]*(T==t0)*dT*(1-dR))/Y1
    if (Y1==0) dhaz1 = 0
    Y2 = sum(w[,a[2]+1]*(T>=t0)*(R>=t0))
    dhaz2 = sum(w[,a[2]+1]*(R==t0)*dR)/Y2
    if (Y2==0) dhaz2 = 0
    hazu = matchy(haz3,TR,t0-tseq[1:k])
    M = matchy(M3,TR,t0-tseq[1:k])
    if (!is.null(conf.int)){
    Y1w = sum(w[,a[1]+1]^2*(T>=t0)*(R>=t0))
    Y2w = sum(w[,a[2]+1]^2*(T>=t0)*(R>=t0))
    if (Y1!=0){
      M1 = append(M1, Y1w*dhaz1/Y1^2)
    } else {
      M1 = append(M1, 0)
    }
    if (Y2!=0){
      M2 = append(M2, Y2w*dhaz2/Y2^2)
    } else {
      M2 = append(M2, 0)
    }
    G1 = append(G1, sum((1-F1t-F2t+sum(F2*exp(-hazu))-
                           cumsum(F2*exp(-hazu)))^2*M1))
    G2 = append(G2, sum((1-F1t-F2t-exp(-haz1-haz2)*exp(-hazu)+
                           sum(F2*exp(-hazu))-cumsum(F2*exp(-hazu)))^2*M2))
    G3 = append(G3, sum((F2*exp(-hazu)*(2*cumsum(F2*exp(-hazu))-F2*exp(-hazu)))*M))
    }
    F1t = F1t + exp(-haz1t-haz2t)*dhaz1
    F2t = F2t + exp(-haz1t-haz2t)*dhaz2
    F3t = sum(F2*(1-exp(-hazu)))
    F1 = append(F1, exp(-haz1t-haz2t)*dhaz1)
    F2 = append(F2, exp(-haz1t-haz2t)*dhaz2)
    haz1t = haz1t + dhaz1
    haz2t = haz2t + dhaz2
    haz1 = append(haz1, haz1t)
    haz2 = append(haz2, haz2t)
    F3 = append(F3, F3t)
  }
  F1 = cumsum(F1)
  F2 = cumsum(F2)
  if (!is.null(conf.int)) se = sqrt(G1+G2+G3)
  if (tseq[1]==0){
    tseq = tseq[-1]
    F1 = F1[-1]
    F2 = F2[-1]
    F3 = F3[-1]
    if (!is.null(conf.int)) se = se[-1]
  }
  if (!is.null(conf.int)){
  return(list(dat=dat,time=c(0,tseq),F1=F1,F2=F2,F3=F3,a=a,se=se))
  }
  return(list(dat=dat,time=c(0,tseq),F1=F1,F2=F2,F3=F3,a=a))
}

semi_sensitivity <- function(A,T,dT,R,dR,X=NULL,weights=NULL,a,sens=0){
  w = pscore(A,X,weights)
  if (length(a)==1) a = rep(a,3)
  haz1 = haz2 = haz3 = haz1t = haz2t = haz3t = haz3mt = haz3m = 0
  F1 = F2 = F3 = F1t = F2t = F3t = 0
  TR = sort(unique((T-R)[dR==1]))
  K = length(TR)
  for (k in 1:K){
    t0 = TR[k]
    Y3 = sum(w[,a[3]+1]*(T-R>=t0)*dR)
    dhaz3 = sum(w[,a[3]+1]*(T-R==t0)*dT*dR)/Y3
    if (Y3==0) dhaz3 = 0
    haz3 = append(haz3, haz3t + dhaz3)
    haz3t = haz3t + dhaz3
  }
  if (TR[1]==0){
    haz3 = haz3[-1]
  } else {
    TR = c(0,TR)
  }
  tseq = sort(unique(c(T[dT==1],R[dR==1])))
  K = length(tseq)
  for (k in 1:K){
    t0 = tseq[k]
    Y1 = sum(w[,a[1]+1]*(T>=t0)*(R>=t0))
    dhaz1 = sum(w[,a[1]+1]*(T==t0)*dT*(1-dR))/Y1
    if (Y1==0) dhaz1=0
    Y2 = sum(w[,a[2]+1]*(T>=t0)*(R>=t0))
    dhaz2 = sum(w[,a[2]+1]*(R==t0)*dR)/Y2
    if (Y2==0) dhaz2=0
    Y3 = sum(w[,a[3]+1]*(T>=t0)*(R<=t0)*dR)
    dhaz3 = sum(w[,a[3]+1]*(T==t0)*dT*dR)/Y3
    if (Y3==0) dhaz3=0
    hazu = matchy(haz3,TR,t0-tseq[1:k])
    hazu = 1 - exp((1-sens)*(-haz3mt+haz3m) - sens*hazu)
    haz3mt = haz3mt + dhaz3
    haz3m = append(haz3m, haz3mt)
    F1t = F1t + exp(-haz1t-haz2t)*dhaz1
    F2t = F2t + exp(-haz1t-haz2t)*dhaz2
    F3t = sum(F2*hazu)
    F1 = append(F1, exp(-haz1t-haz2t)*dhaz1)
    F2 = append(F2, exp(-haz1t-haz2t)*dhaz2)
    haz1t = haz1t + dhaz1
    haz2t = haz2t + dhaz2
    haz1 = append(haz1, haz1t)
    haz2 = append(haz2, haz2t)
    F3 = append(F3, F3t)
  }
  F1 = cumsum(F1)
  F2 = cumsum(F2)
  if (tseq[1]==0){
    tseq = tseq[-1]
    F1 = F1[-1]
    F2 = F2[-1]
    F3 = F3[-1]
  }
  return(list(time=c(0,tseq),F1=F1,F2=F2,F3=F3,a=a,sens=sens))
}

treatmentpolicy <- function(A,T,dT,R=NULL,dR=NULL,X=NULL,weights=NULL,a,conf.int=NULL){
  w = pscore(A,X,weights)
  tseq = sort(unique(T[dT==1]))
  if (!is.null(dR)) tseq = sort(unique(c(T[dT==1],R[dR==1])))
  K = length(tseq)
  haz = hazt = 0
  F = Ft = G = 0
  for (k in 1:K){
    t0 = tseq[k]
    Y = sum(w[,a+1]*(T>=t0))
    dhaz = sum(w[,a+1]*(T==t0)*dT)/Y
    if (Y==0) dhaz = 0
    Ft = Ft + exp(-hazt)*dhaz
    if (!is.null(conf.int)){
      Yw = sum(w[,a+1]^2*(T>=t0))
      dG1 = Yw*dhaz/Y^2
      if (Y==0) dG1 = 0
      G = append(G, dG1)
    }
    hazt = hazt + dhaz
    haz = append(haz, hazt)
    F = append(F, Ft)
  }
  if (!is.null(conf.int)){
    G[is.na(G)] = 0
    se = sqrt(G)
  }
  if (tseq[1]==0){
    tseq = tseq[-1]
    F = F[-1]
    if (!is.null(conf.int)) se = se[-1]
  }
  if (!is.null(conf.int)){
    return(list(time=c(0,tseq),F=F,a=a,se=se))
  }
  return(list(time=c(0,tseq),F=F,a=a))
}

huang.nx <- function(A,T,dT,R,dR,subset=NULL,a){
  if (!is.null(subset)){
    A = A[subset]
    T = T[subset]
    R = R[subset]
    dT = dT[subset]
    dR = dR[subset]
  }
  if (length(a)==1) a = rep(a,2)
  if (length(a)>2) a = a[1:2]
  a1 = a[1]
  a2 = a[2]
  tseq = sort(unique(c(T[dT==1],R[dR==1])))
  K = length(tseq)
  haz = NULL
  for (k in 1:K){
    t0 = tseq[k]
    Y0 = sum(T>=t0 & R>t0 & A==a2)
    Y1 = sum(T>=t0 & R<=t0 & dR==1 & A==a2)
    YT0 = sum(T>=t0 & R>t0 & A==a1)
    YT1 = sum(T>=t0 & R<=t0 & dR==1 & A==a1)
    N0 = sum(T==t0 & R>t0 & A==a1 & dT==1)
    N1 = sum(T==t0 & R<=t0 & A==a1 & dT==1)
    haz0 = Y0/(Y0+Y1)*N0/YT0
    haz1 = Y1/(Y0+Y1)*N1/YT1
    if (YT0*(Y1+Y0)==0) haz0 = 0
    if (YT1*(Y1+Y0)==0) haz1 = 0
    haz = append(haz, haz1+haz0)
  }
  F = c(0, 1 - exp(-cumsum(haz)))
  Time = c(0, tseq)
  return(list(Time=Time,F=F,a=a))
}

huang <- function(A,T,dT,R,dR,X,a){
  X.un = unique(X)
  P = sapply(X.un,function(u) mean(X==u))
  tseq = sort(unique(c(0,T[dT==1],R[dR==1])))
  F = rep(0,length(tseq))
  for (k in 1:length(X.un)){
    x = X.un[k]
    fit = huang.nx(A,T,dT,R,dR,subset=(X==x),a)
    t = fit$Time
    Fx = fit$F
    F = F + matchy(Fx,t,tseq)*P[k]
  }
  return(list(Time=tseq,F=F,a=a))
}

testpath <- function(A,T,dT,R,dR,X=NULL,weights=NULL){
  w = pscore(A,X,weights)
  tseq = sort(unique(c(T[dT==1],R[dR==1])))
  TR = sort(unique((T-R)[dR==1]))
  K = length(tseq)
  test0 = test1 = test2 = test3 = test3s = test4 = test5 = 0
  var0 = var1 = var2 = var3 = var3s = var4 = var5 = 0
  for (k in 1:K){
    t0 = tseq[k]
    Y1.1 = sum(w[,2]*(T>=t0)*(R>=t0))
    Y1.0 = sum(w[,1]*(T>=t0)*(R>=t0))
    if (Y1.1*Y1.0!=0){
      N1.1 = sum(w[,2]*(T==t0)*dT*(1-dR))
      N2.1 = sum(w[,2]*(R==t0)*dR)
      N1.0 = sum(w[,1]*(T==t0)*dT*(1-dR))
      N2.0 = sum(w[,1]*(R==t0)*dR)
      Y1.1w = sum(w[,2]^2*(T>=t0)*(R>=t0))
      Y1.0w = sum(w[,1]^2*(T>=t0)*(R>=t0))
      dtest1 = (Y1.0*N1.1-Y1.1*N1.0)/(Y1.1+Y1.0)
      dtest2 = (Y1.0*N2.1-Y1.1*N2.0)/(Y1.1+Y1.0)
      dvar1 = (Y1.0^2*Y1.1w+Y1.1^2*Y1.0w)*(N1.1+N1.0)/(Y1.1+Y1.0)^3
      dvar2 = (Y1.0^2*Y1.1w+Y1.1^2*Y1.0w)*(N2.1+N2.0)/(Y1.1+Y1.0)^3
      test1 = test1 + dtest1
      test2 = test2 + dtest2
      var1 = var1 + dvar1
      var2 = var2 + dvar2
    }
    Y0.1 = sum(w[,2]*(T>=t0))
    Y0.0 = sum(w[,1]*(T>=t0))
    if (Y0.1*Y0.0!=0){
      N0.1 = sum(w[,2]*(T==t0)*dT)
      N0.0 = sum(w[,1]*(T==t0)*dT)
      N4.1 = sum(w[,2]*(T==t0)*dR*dT)
      N4.0 = sum(w[,1]*(T==t0)*dR*dT)
      N5.1 = sum(w[,2]*(T==t0)*(1-dR)*dT)
      N5.0 = sum(w[,1]*(T==t0)*(1-dR)*dT)
      Y0.1w = sum(w[,2]^2*(T>=t0))
      Y0.0w = sum(w[,1]^2*(T>=t0))
      dtest0 = (Y0.0*N0.1-Y0.1*N0.0)/(Y0.1+Y0.0)
      dtest4 = (Y0.0*N4.1-Y0.1*N4.0)/(Y0.1+Y0.0)
      dtest5 = (Y0.0*N5.1-Y0.1*N5.0)/(Y0.1+Y0.0)
      dvar0 = (Y0.0^2*Y0.1w+Y0.1^2*Y0.0w)*(N0.1+N0.0)/(Y0.1+Y0.0)^3
      dvar4 = (Y0.0^2*Y0.1w+Y0.1^2*Y0.0w)*(N4.1+N4.0)/(Y0.1+Y0.0)^3
      dvar5 = (Y0.0^2*Y0.1w+Y0.1^2*Y0.0w)*(N5.1+N5.0)/(Y0.1+Y0.0)^3
      test0 = test0 + dtest0
      test4 = test4 + dtest4
      test5 = test5 + dtest5
      var0 = var0 + dvar0
      var4 = var4 + dvar4
      var5 = var5 + dvar5
    }
    Y3.1 = sum(w[,2]*(T>=t0)*(R<=t0)*dR)
    Y3.0 = sum(w[,1]*(T>=t0)*(R<=t0)*dR)
    if (Y3.1*Y3.0!=0){
      N3.1 = sum(w[,2]*(T==t0)*dT*dR)
      N3.0 = sum(w[,1]*(T==t0)*dT*dR)
      Y3.1w = sum(w[,2]^2*(T>=t0)*(R<=t0)*dR)
      Y3.0w = sum(w[,1]^2*(T>=t0)*(R<=t0)*dR)
      dtest3 = (Y3.0*N3.1-Y3.1*N3.0)/(Y3.1+Y3.0)
      dvar3 = (Y3.0^2*Y3.1w+Y3.1^2*Y3.0w)*(N3.1+N3.0)/(Y3.1+Y3.0)^3
      test3 = test3 + dtest3
      var3 = var3 + dvar3
    }
  }
  K = length(TR)
  for (k in 1:K){
    t0 = TR[k]
    Y3.1 = sum(w[,2]*(T-R>=t0)*dR)
    Y3.0 = sum(w[,1]*(T-R>=t0)*dR)
    if (Y3.1*Y3.0!=0){
      N3.1 = sum(w[,2]*(T-R==t0)*dR*dT)
      N3.0 = sum(w[,1]*(T-R==t0)*dR*dT)
      Y3.1w = sum(w[,2]^2*(T-R>=t0)*dR)
      Y3.0w = sum(w[,1]^2*(T-R>=t0)*dR)
      dtest3 = (Y3.0*N3.1-Y3.1*N3.0)/(Y3.1+Y3.0)
      dvar3 = (Y3.0^2*Y3.1w+Y3.1^2*Y3.0w)*(N3.1+N3.0)/(Y3.1+Y3.0)^3
      test3s = test3s + dtest3
      var3s = var3s + dvar3
    }
  }
    p1 = 2*pnorm(-abs(test1)/sqrt(var1))
    p2 = 2*pnorm(-abs(test2)/sqrt(var2))
    p3 = 2*pnorm(-abs(test3)/sqrt(var3))
    p3s = 2*pnorm(-abs(test3s)/sqrt(var3s))
    p4 = 2*pnorm(-abs(test4)/sqrt(var4))
    p5 = 2*pnorm(-abs(test5)/sqrt(var5))
    p0 = 2*pnorm(-abs(test0)/sqrt(var0))
    if (var1==0) p1=1
    if (var2==0) p2=1
    if (var3==0) p3=1
    if (var3s==0) p3s=1
    if (var4==0) p4=1
    if (var5==0) p5=1
    if (var0==0) p0=1
    return(list(p1=p1,p2=p2,p3=c(p3,p3s),p0=p0,p23=p4,p01=p5))
}

cox_est <- function(A,T,R,dT,dR,X=NULL,max_iter=50){
  t <- T[dR == 1]
  r <- R[dR == 1]
  a <- A[dR == 1]
  dtr <- dT[dR == 1]
  X_all = cbind(R, A, X)
  X <- X_all[dR == 1,]
  o <- order(t)
  to <- t[o]
  ro <- r[o]
  ao <- a[o]
  dto <- dtr[o]
  Xo <- X[o,]
  beta_initial <- rep(0,ncol(X))
  L0 <- 0
  for(i in 1:max_iter){
    eta = geteta(Xo, beta_initial)
    exp_eta = exp(eta)
    K = length(to)
    L1 <- rep(0, ncol(X))
    L2 <- matrix(0, nrow = ncol(X), ncol = ncol(X))
    l1 <- H1 <- matrix(NA, nrow = K, ncol = ncol(X))
    H2 <- rep(NA, length = K)
    l2 <- H3 <- matrix(NA, nrow = K*ncol(X),ncol = ncol(X) )
    for(k in 1:K){
      t0 <- to[k]
      H1[k,] <- apply((exp_eta*Xo)*(to >= t0)*(ro <= t0),2,sum)
      H2[k] <- sum((exp_eta)[(to >= t0) & (ro <= t0)])
      l1[k,] <- dto[k]*(Xo[k,] - H1[k,]/H2[k])
      L1 = L1 + l1[k,]
      H3[((k-1)*ncol(Xo)+1):(k*ncol(Xo)),] <- t(Xo[to >= t0 & ro <= t0,]) %*%
        (exp_eta*Xo)[to >= t0 & ro <= t0,]
      l2[((k-1)*ncol(Xo)+1):(k*ncol(Xo)),] <- H3[((k-1)*ncol(Xo)+1):(k*ncol(Xo)),]/H2[k] -
        H1[k,]%*%t(H1[k,])/((H2[k])^2)
      L2 = L2 + dto[k]*l2[((k-1)*ncol(Xo)+1):(k*ncol(Xo)),]
    }
    L <- sum(dto*(eta - log(H2)))/K
    beta <- beta_initial + solve(L2,L1)
    if(abs(L-L0)/(0.1 + abs(L)) < 1e-4){
      Inform = solve(L2)
      fit <- NULL
      fit$se <- sqrt(diag(Inform))
      p.value = 2*pnorm( - abs(beta/fit$se))
      fit$beta <- beta
      fit$p.value <- p.value
      return(fit)
      break
    }
    L0 <- L
    beta_initial = (beta + beta_initial)/2
    if(i == max_iter){
      warning("Algorithm doesn't Converge!")
      Inform = solve(L2)
      fit <- NULL
      fit$se <- sqrt(diag(Inform))
      p.value = 2*pnorm( - abs(beta/fit$se))
      fit$beta <- beta
      fit$p.value <- p.value
      return(fit)
    }
  }
}

semicr <- function(A, T, dT=rep(1,length(T)), R, dR=rep(1,length(R)), X=NULL, weights=NULL, a,
                   asm=c('markov','semimarkov'), sens=NULL, conf.int=NULL, nboot=0, subset=NULL){
  if (length(asm)>1) asm = asm[1]
  if (!asm %in% c('markov','semimarkov')) {
    return(cat('Please specify "asm" = "markov" or "semimarkov"!'))
  }
  if (!is.null(X)) X = as.matrix(X)
  if (!is.null(subset)){
    A = A[subset]
    T = T[subset]
    R = R[subset]
    dT = dT[subset]
    dR = dR[subset]
    if (!is.null(X)) X = X[subset,]
  }
  if (!is.null(sens)) fit=semi_sensitivity(A,T,dT,R,dR,X,weights,a,sens)
  if (asm=='markov') fit=semi_markov(A,T,dT,R,dR,X,weights,a,conf.int)
  if (asm=='semimarkov') fit=semi_semimarkov(A,T,dT,R,dR,X,weights,a,conf.int)
  Time = fit$time
  if (!is.null(conf.int)){
    if (nboot<=1){
      ci_l = fit$F1+fit$F3+qnorm((1-conf.int)/2)*fit$se
      ci_u = fit$F1+fit$F3-qnorm((1-conf.int)/2)*fit$se
      ci_l[is.na(ci_l)] = 0
      ci_u[is.na(ci_u)] = 0
    } else {
    boot.matrix = NULL
    for (b in 1:nboot){
      fit.b = semicr(A,T,dT,R,dR,X,weights,a,asm,sens=sens,subset=sample(length(A),replace=TRUE))
      tm = fit.b$time
      F = fit.b$F1+fit.b$F3
      boot.matrix = rbind(boot.matrix, matchy(c(0,F),c(0,tm),Time))
    }
    ci = apply(boot.matrix,2,function(x) quantile(x,c((1-conf.int)/2,1-(1-conf.int)/2)))
    ci[is.na(ci)] = 0
    ci_l = as.numeric(ci[1,])
    ci_u = as.numeric(ci[2,])
    }
    return(c(fit,list(ci_l=ci_l,ci_u=ci_u)))
  }
  return(fit)
}

semicr.causal <- function(A, T, dT, R, dR, X=NULL, weights=NULL, a1=c(1,0,0), a0=c(0,0,0), 
                          asm='markov', sens=NULL, conf.int=NULL, nboot=0,
                          xlab='Time', ylim=c(0,1), legend=c('A = 1','A = 0'), cex=0.8, ...){
  fit1 = semicr(A,T,dT,R,dR,X,weights,a1,asm,sens,conf.int,nboot)
  fit0 = semicr(A,T,dT,R,dR,X,weights,a0,asm,sens,conf.int,nboot)
  plot(fit1$time,fit1$F1,type='s',lwd=2,ylim=c(0,1),col='brown',
       ylab='Cumulative incidence',xlab=xlab,main='Cumulative incidences of each state')
  points(fit1$time,1-fit1$F2,type='s',lwd=2,lty=5,col='brown')
  points(fit1$time,1-fit1$F3,type='s',lwd=2,lty=4,col='brown')
  points(fit0$time,fit0$F1,type='s',lwd=2,col='darkcyan')
  points(fit0$time,1-fit0$F2,type='s',lwd=2,lty=5,col='darkcyan')
  points(fit0$time,1-fit0$F3,type='s',lwd=2,lty=4,col='darkcyan')
  legend('left',cex=cex,legend=legend,col=c('brown','darkcyan'),
         lwd=c(2,2),...)
  plot(fit1$time,fit1$F1+fit1$F3,type='s',lwd=2,ylim=ylim,col='brown',
       ylab='Cumulative incidence',xlab=xlab,main='Incidence of the primary event')
  points(fit0$time,fit0$F1+fit0$F3,type='s',lwd=2,col='darkcyan')
  if (is.null(sens) & !is.null(conf.int)){
    points(fit1$time,fit1$ci_u,type='s',lwd=1,lty=2,col='brown')
    points(fit1$time,fit1$ci_l,type='s',lwd=1,lty=2,col='brown')
    points(fit0$time,fit0$ci_u,type='s',lwd=1,lty=2,col='darkcyan')
    points(fit0$time,fit0$ci_l,type='s',lwd=1,lty=2,col='darkcyan')
  }
  legend('topleft',cex=cex,legend=legend,col=c('brown','darkcyan'),
         lwd=c(2,2),...)
}

semicr.sensitivity <- function(A, T, dT, R, dR, X=NULL, a1=1, weights=NULL, a0=0, sens=seq(0,1,0.2),
                               xlab='Time', ylim=c(-1,1), cex=0.8, ...){
  plot(NULL,NULL,xlim=c(0,max(c(T[dT==1],R[dR==1]))),ylim=ylim,xlab=xlab,ylab='Separable effect',
       main='Sensitivity Analysis')
  mtext('0: markov, 1: semimarkov')
  abline(h=0,lty=3)
  m = length(sens)
  for (s in 1:m){
    fit1 = semi_sensitivity(A,T,dT,R,dR,X,weights,a1,sens[s])
    fit0 = semi_sensitivity(A,T,dT,R,dR,X,weights,a0,sens[s])
    tm = fit1$time
    y = fit1$F1+fit1$F3-fit0$F1-fit0$F3
    points(tm,y,type='s',col=s+1,lwd=1.5)
  }
  legend('topleft',cex=cex,legend=sens,col=(1:m)+1,lwd=rep(1.5,m),...)
}

#Time = ((R+T)-abs(R-T))/2
#cstatus = dT + 2*dR
#cstatus[cstatus>2] = 1