# Confidence interval coverage

source('semicompeting.R')
expit <- function(x) 1/(1+exp(-x))
xseq = seq(0,8,length=100)
generatedata <- function(N,setting,Z1=1,Z2=1,Z3=1){
  X1 = rbinom(N,1,0.5)/2+0.5
  X2 = rbinom(N,1,0.5)/2+0.5
  X = cbind(X1,X2)
  A = rbinom(N,1,expit(0.4*X1+0.8*X2-0.6))
  if (setting==1){ #markov&semimarkov
    TorR0 = rexp(N,0.25*X1)
    TorR1 = rexp(N,0.25*X1+0.15*Z1+0.1*Z2)
    T0d = rbinom(N,1,15/25)
    T1d = rbinom(N,1,0.15*(X1+Z1)/(0.25*X1+0.15*Z1+0.1*Z2))
    # If Tad=1, T happens; If Tad=0, R happens.
    # Next generate T-R when Tad=0.
    TwR0 = rexp(N,0.2*X2)
    TwR1 = rexp(N,0.2*(X2+Z3))
  }
  if (setting==2){ #markov
    TorR0 = rweibull(N,2,sqrt(2/(0.06*X1)))
    TorR1 = rweibull(N,2,sqrt(2/(0.06*X1+0.04*Z1+0.02*Z2)))
    T0d = rbinom(N,1,2/3)
    T1d = rbinom(N,1,(2*X1+2*Z1)/(3*X1+2*Z1+Z2))
    # If Tad=1, T happens; If Tad=0, R happens.
    # Next generate T-R when Tad=0.
    TwR0 = rexp(N,0.05*X2*1/2)
    TwR1 = rexp(N,0.05*(X2+Z3)/2)
    TwR0 = sqrt(TwR0+TorR0^2) - TorR0
    TwR1 = sqrt(TwR1+TorR1^2) - TorR1
  }
  if (setting==3){ #semimarkov
    TorR0 = rweibull(N,2,sqrt(2/(0.06*X1)))
    TorR1 = rweibull(N,2,sqrt(2/(0.06*X1+0.04*Z1+0.02*Z2)))
    T0d = rbinom(N,1,2/3)
    T1d = rbinom(N,1,2*(X1+Z1)/(3*X1+2*Z1+Z2))
    # If Tad=1, T happens; If Tad=0, R happens.
    # Next generate T-R when Tad=0.
    TwR0 = rweibull(N,2,sqrt(2/(0.1*X2)))
    TwR1 = rweibull(N,2,sqrt(2/(0.1*X2+0.1*Z3)))
  }
  R1 = TorR1 + T1d*99
  R0 = TorR0 + T0d*99
  T1 = TorR1 + (1-T1d)*TwR1
  T0 = TorR0 + (1-T0d)*TwR0
  cif1 = sapply(xseq, function(x) mean(T1<=x))
  cif0 = sapply(xseq, function(x) mean(T0<=x))
  T = A*T1 + (1-A)*T0
  R = A*R1 + (1-A)*R0
  C = runif(N,6,10)
  dR = as.numeric(R<=C)
  dT = as.numeric(T<=C)
  Rt = dR*R + (1-dR)*C
  Tt = dT*T + (1-dT)*C
  return(list(A=A,T=Tt,dT=dT,R=Rt,dR=dR,X=X,cif1=cif1,cif0=cif0))
}

cif <- function(setting,a){
  xseq = seq(0,8,length=100)
  N = 10000
  Z1 = a[1]
  Z2 = a[2]
  Z3 = a[3]
  X1 = rbinom(N,1,0.5)/2+0.5
  X2 = rbinom(N,1,0.5)/2+0.5
  if (setting==1){ #markov&semimarkov
    TorR1 = rexp(N,0.25*X1+0.15*Z1+0.1*Z2)
    T1d = rbinom(N,1,0.15*(X1+Z1)/(0.25*X1+0.15*Z1+0.1*Z2))
    # If Tad=1, T happens; If Tad=0, R happens.
    # Next generate T-R when Tad=0.
    TwR1 = rexp(N,0.2*(X2+Z3))
  }
  if (setting==2){ #markov
    TorR1 = rweibull(N,2,sqrt(2/(0.06*X1+0.04*Z1+0.02*Z2)))
    T1d = rbinom(N,1,(2*X1+2*Z1)/(3*X1+2*Z1+Z2))
    # If Tad=1, T happens; If Tad=0, R happens.
    # Next generate T-R when Tad=0.
    TwR1 = rexp(N,0.05*(X2+Z3)/2)
    TwR1 = sqrt(TwR1+TorR1^2) - TorR1
  }
  if (setting==3){ #semimarkov
    TorR1 = rweibull(N,2,sqrt(2/(0.06*X1+0.04*Z1+0.02*Z2)))
    T1d = rbinom(N,1,2*(X1+Z1)/(3*X1+2*Z1+Z2))
    # If Tad=1, T happens; If Tad=0, R happens.
    # Next generate T-R when Tad=0.
    TwR1 = rweibull(N,2,sqrt(2/(0.1*X2+0.1*Z3)))
  }
  R1 = TorR1 + T1d*99
  T1 = TorR1 + (1-T1d)*TwR1
  cif.a = sapply(xseq, function(x) mean(T1<=x))
  return(cif.a)
}


B = 10000
cvtab = NULL
a = c(1,0,1)
Z1 = a[1]
Z2 = a[2]
Z3 = a[3]
xseq = 1:8

cat('Coverage (markov/semimarkov), Width (markov/semimarkov)\n')

for (setting in 1:3){
  for (N in c(100,500)){
    cat('Setting',setting,'\n')
    cat('N =',N,'\n')
    set.seed(2024)
    count.mar = count.smm = rep(0,length(xseq))
    width.mar = width.smm = rep(0,length(xseq))
    #yt = meaninc(xseq,setting,a)
    yt = matchy(cif(setting,a),seq(0,8,length=100),xseq)
    for (i in 1:B){
      dat = generatedata(N,setting,Z1,Z2,Z3)
      fit = semicr(dat$A,dat$T,dat$dT,dat$R,dat$dR,dat$X,a=a,asm='markov',conf.int=.95)
      x = fit$time
      index = sapply(xseq, function(v) max(which(x<=v)))
      upper = fit$ci_u[index]
      lower = fit$ci_l[index]
      upper[is.na(upper)] = 1
      lower[is.na(lower)] = 0
      count.mar = count.mar + (upper>=yt & lower<=yt)
      width.mar = width.mar + (upper-lower)
      fit = semicr(dat$A,dat$T,dat$dT,dat$R,dat$dR,dat$X,a=a,asm='semimarkov',conf.int=.95)
      x = fit$time
      index = sapply(xseq, function(v) max(which(x<=v)))
      upper = fit$ci_u[index]
      lower = fit$ci_l[index]
      upper[is.na(upper)] = 1
      lower[is.na(lower)] = 0
      count.smm = count.smm + (upper>=yt & lower<=yt)
      width.smm = width.smm + (upper-lower)
    }
    cat(a,'\n')
    cv.mar = count.mar/B
    cv.smm = count.smm/B
    wd.mar = round(width.mar/B,4)
    wd.smm = round(width.smm/B,4)
    cvtab = rbind(cvtab, rbind(cv.mar,wd.mar,cv.smm,wd.smm))
    print(cvtab)
  }
}

library(xtable)
xtable(cvtab, digits=4)
