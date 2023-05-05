generatedata <- function(N,setting,h1,h0){
  X = rbinom(N,1,0.5)/2 + 0.5
  A = rbinom(N,1,X/3+0.25)
  if (setting==1){ #markov&semimarkov
    TorR1 = rexp(N,(h1[1]+h1[2])*.75)
    TorR0 = rexp(N,(h0[1]+h0[2])*.75)
    T1d = rbinom(N,1,h1[1]/(h1[1]+h1[2]))
    T0d = rbinom(N,1,h0[1]/(h0[1]+h0[2]))
    # If Tad=1, T happens; If Tad=0, R happens.
    # Next generate T-R when Tad=0.
    TwR1 = rexp(N,h1[3]*.75)
    TwR0 = rexp(N,h0[3]*.75)
  }
  if (setting==2){ #markov
    TorR1 = rweibull(N,2,sqrt(2/((h1[1]+h1[2])*.75)))
    TorR0 = rweibull(N,2,sqrt(2/((h0[1]+h0[2])*.75)))
    T1d = rbinom(N,1,h1[1]/(h1[1]+h1[2]))
    T0d = rbinom(N,1,h0[1]/(h0[1]+h0[2]))
    # If Tad=1, T happens; If Tad=0, R happens.
    # Next generate T-R when Tad=0.
    TwR1 = rexp(N,h1[3]*.75/2)
    TwR0 = rexp(N,h0[3]*.75/2)
    TwR1 = sqrt(TwR1+TorR1^2) - TorR1
    TwR0 = sqrt(TwR0+TorR0^2) - TorR0
  }
  if (setting==3){ #semimarkov
    TorR1 = rweibull(N,2,sqrt(2/((h1[1]+h1[2])*.75)))
    TorR0 = rweibull(N,2,sqrt(2/((h0[1]+h0[2])*.75)))
    T1d = rbinom(N,1,h1[1]/(h1[1]+h1[2]))
    T0d = rbinom(N,1,h0[1]/(h0[1]+h0[2]))
    # If Tad=1, T happens; If Tad=0, R happens.
    # Next generate T-R when Tad=0.
    TwR1 = rweibull(N,2,sqrt(2/(h1[3]*.75)))
    TwR0 = rweibull(N,2,sqrt(2/(h0[3]*.75)))
  }
  R1 = TorR1 + T1d*99
  R0 = TorR0 + T0d*99
  T1 = TorR1 + (1-T1d)*TwR1
  T0 = TorR0 + (1-T0d)*TwR0
  T = A*T1 + (1-A)*T0
  R = A*R1 + (1-A)*R0
  C = runif(N,6,10)
  dR = as.numeric(R<=C)
  dT = as.numeric(T<=C)
  Rt = dR*R + (1-dR)*C
  Tt = dT*T + (1-dT)*C
  return(list(A=A,T=Tt,dT=dT,R=Rt,dR=dR,X=X))
}


# Testing separable effects

set.seed(2023)
B = 10000
for (N in c(100,500,2000)){
cat('Power of dhaz1, dhaz2, dhaz3(memeryless), dhaz3(semimarkov)\n')
cat('N =',N,'\n')
for (setting in 1:3){
cat('Setting',setting,'\n')
  if (setting==1){
    h1 = rbind(c(0.2,0.2,0.2),c(0.2,0.2,0.2),c(0.2,0.2,0.2),c(0.2,0.2,0.2),c(0.2,0.2,0.2))
    h0 = rbind(c(0.5,0.5,0.5),c(0.2,0.5,0.5),c(0.5,0.2,0.5),c(0.5,0.5,0.2),c(0.2,0.2,0.2))
  }
  if (setting==2){
    h1 = rbind(c(0.1,0.1,0.1),c(0.1,0.1,0.1),c(0.1,0.1,0.1),c(0.1,0.1,0.1),c(0.1,0.1,0.1))
    h0 = rbind(c(0.2,0.2,0.2),c(0.1,0.2,0.2),c(0.2,0.1,0.2),c(0.2,0.2,0.1),c(0.1,0.1,0.1))
  }
  if (setting==3){
    h1 = rbind(c(0.1,0.1,0.1),c(0.1,0.1,0.1),c(0.1,0.1,0.1),c(0.1,0.1,0.1),c(0.1,0.1,0.1))
    h0 = rbind(c(0.2,0.2,0.2),c(0.1,0.2,0.2),c(0.2,0.1,0.2),c(0.2,0.2,0.1),c(0.1,0.1,0.1))
  }
for (s in 1:5){
  p1 = p2 = p3 = p4 = 0
  for (b in 1:B){
    dat = generatedata(N,setting,h1[s,],h0[s,])
    #weights = dat$A/(dat$X/3+0.25) + (1-dat$A)/(1-dat$X/3-0.25)
    fit = testpath(dat$A,dat$T,dat$dT,dat$R,dat$dR,dat$X)
    p1 = p1 + (fit$p1<0.05)
    p2 = p2 + (fit$p2<0.05)
    p3 = p3 + (fit$p3[1]<0.05)
    p4 = p4 + (fit$p3[2]<0.05)
  }
  tie = c(p1,p2,p3,p4)/B
  print(tie)
}
}
}


# Confidence interval coverage

B = 10000
a = 1
xseq = 1:8
for (N in c(100,500,2000)){
cat('Coverage (markov/semimarkov), Width (markov/semimarkov)\n')
cat('N =',N,'\n')
for (setting in 1:3){
cat('Setting',setting,'\n')
  if (setting==1){
    h1 = c(0.2,0.2,0.2)
    h0 = c(0.5,0.5,0.5)
  }
  if (setting==2){
    h1 = c(0.1,0.1,0.1)
    h0 = c(0.2,0.2,0.2)
  }
  if (setting==3){
    h1 = c(0.1,0.1,0.1)
    h0 = c(0.2,0.2,0.2)
  }
set.seed(2023)
count.mem = count.mar = rep(0,length(xseq))
width.mem = width.mar = rep(0,length(xseq))
if (a==1) {
  h1a = h1[1]
  h2a = h1[2]
  h3a = h1[3]
} else {
  h1a = h0[1]
  h2a = h0[2]
  h3a = h0[3]
}
yt = meaninc(xseq,setting,h1a,h2a,h3a)
for (i in 1:B){
  dat = generatedata(N,setting,h1,h0)
  fit = semicr(dat$A,dat$T,dat$dT,dat$R,dat$dR,dat$X,a=a,asm='markov',conf.int=.95)
  x = fit$time
  index = sapply(xseq, function(v) max(which(x<=v)))
  upper = fit$ci_u[index]
  lower = fit$ci_l[index]
  count.mem = count.mem + (upper>=yt & lower<=yt)
  width.mem = width.mem + (upper-lower)
  fit = semicr(dat$A,dat$T,dat$dT,dat$R,dat$dR,dat$X,a=a,asm='semimarkov',conf.int=.95)
  x = fit$time
  index = sapply(xseq, function(v) max(which(x<=v)))
  upper = fit$ci_u[index]
  lower = fit$ci_l[index]
  count.mar = count.mar + (upper>=yt & lower<=yt)
  width.mar = width.mar + (upper-lower)
}
cat(a,'\n')
print(count.mem/B)
print(count.mar/B)
print(width.mem/B)
print(width.mar/B)
}
}

# Test semi-semimarkov

set.seed(2023)
B = 500
for (N in c(100,500,2000)){
  cat('Power of testing semi-semimarkov\n')
  cat('N =',N,'\n')
  for (setting in 1:3){
    cat('Setting',setting,'\n')
    if (setting==1){
      h1 = rbind(c(0.2,0.2,0.2),c(0.2,0.2,0.2),c(0.2,0.2,0.2),c(0.2,0.2,0.2),c(0.2,0.2,0.2))
      h0 = rbind(c(0.5,0.5,0.5),c(0.2,0.5,0.5),c(0.5,0.2,0.5),c(0.5,0.5,0.2),c(0.2,0.2,0.2))
    }
    if (setting==2){
      h1 = rbind(c(0.1,0.1,0.1),c(0.1,0.1,0.1),c(0.1,0.1,0.1),c(0.1,0.1,0.1),c(0.1,0.1,0.1))
      h0 = rbind(c(0.2,0.2,0.2),c(0.1,0.2,0.2),c(0.2,0.1,0.2),c(0.2,0.2,0.1),c(0.1,0.1,0.1))
    }
    if (setting==3){
      h1 = rbind(c(0.1,0.1,0.1),c(0.1,0.1,0.1),c(0.1,0.1,0.1),c(0.1,0.1,0.1),c(0.1,0.1,0.1))
      h0 = rbind(c(0.2,0.2,0.2),c(0.1,0.2,0.2),c(0.2,0.1,0.2),c(0.2,0.2,0.1),c(0.1,0.1,0.1))
    }
    for (s in 1:5){
    tie = 0
    for (b in 1:B){
      dat = generatedata(N,setting,h1[s,],h0[s,])
      #weights = dat$A/(dat$X/3+0.25) + (1-dat$A)/(1-dat$X/3-0.25)
      TR = dat$T - dat$R
      W = pscore(dat$A,dat$X)
      W = as.numeric(W[,1]+W[,2])
      fit = coxph(Surv(TR,dat$dT)~dat$R+dat$A,weights=W,subset=(dat$dR==1))
      tie = tie + (summary(fit)$coefficients[1,6]<0.05)
    }
    tie = tie/B
    print(tie)
    }
  }
}
