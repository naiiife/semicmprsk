# Testing separable effects

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
tietab = NULL
for (N in c(100,500)){
  cat('Power of dhaz1, dhaz2, dhaz3(markov), dhaz3(semimarkov)\n')
  cat('N =',N,'\n')
  for (setting in 1:3){
    cat('Setting',setting,'\n')
    for (s in 1:6){
      p1 = p2 = p3 = p4 = 0
      Z1 = 1-as.numeric(s==2|s==6)
      Z2 = 1-as.numeric(s==3|s>=5)
      Z3 = 1-as.numeric(s==4|s>=5)
      set.seed(2024)
      for (b in 1:B){
        dat = generatedata(N,setting,Z1,Z2,Z3)
        #weights = dat$A/(dat$X/3+0.25) + (1-dat$A)/(1-dat$X/3-0.25)
        fit = testpath(dat$A,dat$T,dat$dT,dat$R,dat$dR,dat$X)
        p1 = p1 + (fit$p1<0.05)
        p2 = p2 + (fit$p2<0.05)
        p3 = p3 + (fit$p3[1]<0.05)
        p4 = p4 + (fit$p3[2]<0.05)
      }
      tie = c(p1,p2,p3,p4)/B
      print(tie)
      tietab = rbind(tietab, tie)
    }
  }
}

tietab = cbind(tietab[1:18,],tietab[19:36,])
library(xtable)
xtable(tietab, digits=4)
