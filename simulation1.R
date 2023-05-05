## Estimation

generatedata <- function(N,setting,Z1=1,Z2=1,Z3=1){
  X = rbinom(N,1,0.5)/2 + 0.5
  A = rbinom(N,1,X/3+0.25)
  if (setting==1){ #markov&semimarkov
    TorR0 = rexp(N,0.25*X)
    TorR1 = rexp(N,0.25*X+0.15*Z1+0.1*Z2)
    T0d = rbinom(N,1,15/25)
    T1d = rbinom(N,1,0.15*(X+Z1)/(0.25*X+0.15*Z1+0.1*Z2))
    # If Tad=1, T happens; If Tad=0, R happens.
    # Next generate T-R when Tad=0.
    TwR0 = rexp(N,0.2)
    TwR1 = rexp(N,0.2*(1+Z3))
  }
  if (setting==2){ #markov
    TorR0 = rweibull(N,2,sqrt(2/(0.06*X)))
    TorR1 = rweibull(N,2,sqrt(2/(0.06*X+0.04*Z1+0.02*Z2)))
    T0d = rbinom(N,1,2/3)
    T1d = rbinom(N,1,(2*X+2*Z1)/(3*X+2*Z1+Z2))
    # If Tad=1, T happens; If Tad=0, R happens.
    # Next generate T-R when Tad=0.
    TwR0 = rexp(N,0.05*1/2)
    TwR1 = rexp(N,0.05*(1+Z3)/2)
    TwR0 = sqrt(TwR0+TorR0^2) - TorR0
    TwR1 = sqrt(TwR1+TorR1^2) - TorR1
  }
  if (setting==3){ #semimarkov
    TorR0 = rweibull(N,2,sqrt(2/(0.06*X)))
    TorR1 = rweibull(N,2,sqrt(2/(0.06*X+0.04*Z1+0.02*Z2)))
    T0d = rbinom(N,1,2/3)
    T1d = rbinom(N,1,2*(X+Z1)/(3*X+2*Z1+Z2))
    # If Tad=1, T happens; If Tad=0, R happens.
    # Next generate T-R when Tad=0.
    TwR0 = rweibull(N,2,sqrt(2/(0.1*1)))
    TwR1 = rweibull(N,2,sqrt(2/(0.1*1+0.1*Z3)))
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

Erf <- function(x) 2*pnorm(x,0,1/sqrt(2))-1
set1n <- function(t,x,a){
  ax = (a+c(x,x,1))*c(0.75,0.5,1)
  if (ax[1]+ax[2]==ax[3]) {
    y = 1 - exp(-0.2*t*(ax[1]+ax[2])) - 
      exp(-0.2*t*(ax[1]+ax[2]))*t*0.2*(ax[2])
    return(y)
  }
  y1 = 1 - exp(-0.2*t*(ax[1]+ax[2]))
  y2 = (ax[2])/(ax[1]+ax[2]-ax[3])
  y3 = exp(-0.2*t*(ax[3])) - exp(-0.2*t*(ax[1]+ax[2]))
  return(y1-y2*y3)
}
set2n <- function(t,x,a){
  ax = (a+c(x,x,1))*c(0.4,0.2,0.5)
  if (ax[1]+ax[2]==ax[3]) {
    y = 1 - exp(-0.05*t^2*(ax[1]+ax[2])) * (1+0.05*t^2*(ax[2]))
    return(y)
  }
  y1 = 1 - exp(-0.05*t^2*(ax[1]+ax[2]))
  y2 = (ax[2])*exp(-0.05*t^2*(ax[1]+ax[2]))/(ax[1]+ax[2]-ax[3])
  y3 = -1+exp(0.05*t^2*(ax[1]+ax[2]-ax[3]))
  return(y1-y2*y3)
}
set3n <- function(t,x,a){
  ax = (a+c(x,x,1))*c(0.4,0.2,1)
  y1 = 1 - exp(-0.05*t^2*(ax[1]+ax[2]))
  y2 = 0.1*(ax[2])*exp(-0.05*t^2*(ax[1]+ax[2]+2*ax[3]))/sqrt(sum(ax)^3)
  y3 = 10*(exp(0.1*t^2*(ax[3]))-exp(0.05*t^2*sum(ax)))*sqrt(sum(ax))
  y4 = exp(0.05*t^2*(ax[3])^2/sum(ax)+0.05*t^2*sum(ax))*t*sqrt(5*pi)*(ax[3])
  y5 = Erf(t*sqrt(0.05)*(ax[1]+ax[2])/sqrt(sum(ax)))
  y6 = Erf(t*sqrt(0.05)*(ax[3])/sqrt(sum(ax)))
  return(y1-y2*(-y3+y4*(y5+y6)))
}
meaninc <- function(t,setting,a){
  if (setting==1) y = (set1n(t,1,a)+set1n(t,0.5,a))/2
  if (setting==2) y = (set2n(t,1,a)+set2n(t,0.5,a))/2
  if (setting==3) y = (set3n(t,1,a)+set3n(t,0.5,a))/2
  return(y)
}

# Draw cumulative incidence
B = 100
a = c(1,0,0)
for (N in c(100,500,2000)){
  par(mfcol=c(2,3))
  for (setting in 1:3){
    cat('N =',N,' Setting',setting,' a =',a,'\n')
    for (asm in c('markov','semimarkov')){
      set.seed(2024)
      plot(NULL,NULL,xlim=c(0,8),ylim=c(0,1),xlab='Time',ylab='Cumulative incidence',
           main=paste('Setting', setting))
      mtext(asm, cex=0.8)
      for (i in 1:B){
        dat = generatedata(N,setting)
        fit = semicr(dat$A,dat$T,dat$dT,dat$R,dat$dR,dat$X,a=a,asm=asm)
        xt = fit$time
        y = fit$F1+fit$F3
        x = c(xt[xt<8],8)
        y = c(y[xt<8],max(y[xt<8]))
        points(x,y,type='l',col='grey')
      }
      # Asymptotic CI
      fit = semicr(dat$A,dat$T,dat$dT,dat$R,dat$dR,dat$X,a=a,asm=asm,conf.int=.95,nboot=0)
      xt = fit$time
      y = fit$F1+fit$F3
      x = c(xt[xt<8],8)
      y = c(y[xt<8],max(y[xt<8]))
      points(x,y,type='l',col='darkcyan')
      ciu = c(fit$ci_u[xt<8],max(fit$ci_u[xt<8]))
      cil = c(fit$ci_l[xt<8],max(fit$ci_l[xt<8]))
      points(x,ciu,type='l',lwd=1,lty=5,col='darkcyan')
      points(x,cil,type='l',lwd=1,lty=5,col='darkcyan')
      # Bootstrap CI
      fit = semicr(dat$A,dat$T,dat$dT,dat$R,dat$dR,dat$X,a=a,asm=asm,conf.int=.95,nboot=200)
      xt = fit$time
      y = fit$F1+fit$F3
      x = c(xt[xt<8],8)
      y = c(y[xt<8],max(y[xt<8]))
      ciu = c(fit$ci_u[xt<8],max(fit$ci_u[xt<8]))
      cil = c(fit$ci_l[xt<8],max(fit$ci_l[xt<8]))
      points(x,ciu,type='l',lwd=1,lty=4,col='brown')
      points(x,cil,type='l',lwd=1,lty=4,col='brown')
      x = seq(0,8,length=100)
      y = meaninc(x,setting,a)
      points(x,y,type='l')
    }
  }
}

# Draw bias
B = 1000
xseq = seq(0,8,length=100)
for (N in c(100,500,2000)){
  for (alist in 1:3){
    if (alist==1) a = c(0,0,0)
    if (alist==2) a = c(1,0,0)
    if (alist==3) a = c(1,0,1)
    par(mfcol=c(1,3))
    for (setting in 1:3){
      cat('N =',N,' Setting',setting,' a =',a,'\n')
      set.seed(2024)
      plot(NULL,NULL,xlim=c(0,8),ylim=c(-.04,.07),xlab='Time',
           ylab='Bias',main=paste('Setting', setting))
      abline(h=0,lty=2)
      if (alist==1) mtext('a = (0,0,0)', cex=0.8)
      if (alist==2) mtext('a = (1,0,0)', cex=0.8)
      if (alist==3) mtext('a = (1,0,1)', cex=0.8)
      yseq = meaninc(xseq,setting,a)
      bias1 = bias2 = bias3 = bias4 = rep(0,100)
      rmse1 = rmse2 = rmse3 = rmse4 = rep(0,100)
      for (i in 1:B){
        dat = generatedata(N,setting)
        fit = semicr(dat$A,dat$T,dat$dT,dat$R,dat$dR,dat$X,a=a,asm='markov')
        bs = matchy(fit$F1+fit$F3,fit$time,xseq) - yseq
        bias1 = bias1 + bs
        rmse1 = rmse1 + bs^2
        fit = semicr(dat$A,dat$T,dat$dT,dat$R,dat$dR,dat$X,a=a,asm='semimarkov')
        bs = matchy(fit$F1+fit$F3,fit$time,xseq) - yseq
        bias2 = bias2 + bs
        rmse2 = rmse2 + bs^2
        fit = huang(dat$A,dat$T,dat$dT,dat$R,dat$dR,dat$X,a=a[1:2])
        bs = matchy(fit$F,fit$Time,xseq) - yseq
        bias3 = bias3 + bs
        rmse3 = rmse3 + bs^2
      }
      ibias1 = cumsum(bias1/B)*8/99
      ibias2 = cumsum(bias2/B)*8/99
      ibias3 = cumsum(bias3/B)*8/99
      points(xseq,bias1/B,type='l',lwd=1.2,col='blue')
      points(xseq,bias2/B,type='l',lwd=1.2,col='darkcyan')
      points(xseq,bias3/B,type='l',lwd=1.2,col='brown')
      points(xseq,ibias1,type='l',lwd=1.2,lty=5,col='blue')
      points(xseq,ibias2,type='l',lwd=1.2,lty=5,col='darkcyan')
      points(xseq,ibias3,type='l',lwd=1.2,lty=5,col='brown')
      legend('topleft',cex=1,lwd=c(1.2,1.2,1.2),
             col=c('blue','darkcyan','brown'),
             legend=c('Proposed (ma.)','Proposed (sm.)','Huang'))
    }
  }
}

