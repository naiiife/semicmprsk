
## Simulation (dismissible treatment decomposition violated)

generatedata <- function(N,setting,h1,h0){
  X = rbinom(N,1,0.5)/2 + 0.5
  A = rbinom(N,1,X/3+0.25)
  if (setting==1){ #markov&semimarkov
    TorR1 = rexp(N,(h1[1]+h1[2])*X)
    TorR0 = rexp(N,(h0[1]+h0[2])*X)
    T1d = rbinom(N,1,h1[1]/(h1[1]+h1[2]))
    T0d = rbinom(N,1,h0[1]/(h0[1]+h0[2]))
    # If Tad=1, T happens; If Tad=0, R happens.
    # Next generate T-R when Tad=0.
    TwR1 = rexp(N,h1[3])
    TwR0 = rexp(N,h0[3])
  }
  if (setting==2){ #markov
    TorR1 = rweibull(N,2,sqrt(2/((h1[1]+h1[2])*X)))
    TorR0 = rweibull(N,2,sqrt(2/((h0[1]+h0[2])*X)))
    T1d = rbinom(N,1,h1[1]/(h1[1]+h1[2]))
    T0d = rbinom(N,1,h0[1]/(h0[1]+h0[2]))
    # If Tad=1, T happens; If Tad=0, R happens.
    # Next generate T-R when Tad=0.
    TwR1 = rexp(N,h1[3]/2)
    TwR0 = rexp(N,h0[3]/2)
    TwR1 = sqrt(TwR1+TorR1^2) - TorR1
    TwR0 = sqrt(TwR0+TorR0^2) - TorR0
  }
  if (setting==3){ #semimarkov
    TorR1 = rweibull(N,2,sqrt(2/((h1[1]+h1[2])*X)))
    TorR0 = rweibull(N,2,sqrt(2/((h0[1]+h0[2])*X)))
    T1d = rbinom(N,1,h1[1]/(h1[1]+h1[2]))
    T0d = rbinom(N,1,h0[1]/(h0[1]+h0[2]))
    # If Tad=1, T happens; If Tad=0, R happens.
    # Next generate T-R when Tad=0.
    TwR1 = rweibull(N,2,sqrt(2/(h1[3])))
    TwR0 = rweibull(N,2,sqrt(2/(h0[3])))
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

set1n <- function(t,x,a,b,c){
  a = a*x
  b = b*x
  y1 = 1-exp(-(a+b)*t)
  y2 = b/(a+b-c)*(exp(-c*t)-exp(-(a+b)*t))
  if (a+b==c) y2 = b*t^2*exp(-c*t)/2
  return(y1-y2)
}
set2n <- function(t,x,a,b,c){
  a = a*x
  b = b*x
  y1 = 1-exp(-(a+b)*t^2/2)
  y2 = b/(a+b-c)*(exp(-c*t^2/2)-exp(-(a+b)*t^2/2))
  if (a+b==c) y2 = b*t^2*exp(-c*t^2/2)/2
  return(y1-y2)
}
set3n <- function(t,x,a,b,c){
  a = a*x
  b = b*x
  y1 = 1-exp(-(a+b)*t^2/2)
  y2 = b*exp(c^2*t^2/(a+b+c)/2-c*t^2/2)/(a+b+c)
  y3 = exp(-c^2*t^2/(a+b+c)/2)-exp(-(a+b)^2*t^2/(a+b+c)/2)
  y4 = sqrt(2*pi/(a+b+c))*c*t*
    (pnorm(t,c*t/(a+b+c),1/sqrt(a+b+c))-pnorm(0,c*t/(a+b+c),1/sqrt(a+b+c)))
  return(y1-y2*(y3+y4))
}
meaninc <- function(t,setting,a,b,c){
  if (setting==1) y = set1n(t,1,a,b,c)/2+set1n(t,0.5,a,b,c)/2
  if (setting==2) y = set2n(t,1,a,b,c)/2+set2n(t,0.5,a,b,c)/2
  if (setting==3) y = set3n(t,1,a,b,c)/2+set3n(t,0.5,a,b,c)/2
  return(y)
}


## Draw cumulative incidence

B = 100
a = c(1,0,1)
for (N in c(100,500,2000)){
  par(mfcol=c(2,3))
  for (setting in 1:3){
    if (setting==1){
      h0 = c(0.1,0.2,0.25)
      h1 = c(0.2,0.3,0.35)
    }
    if (setting==2){
      h0 = c(0.05,0.15,0.15)
      h1 = c(0.15,0.25,0.25)
    }
    if (setting==3){
      h0 = c(0.05,0.15,0.15)
      h1 = c(0.15,0.25,0.25)
    }
    cat('N =',N,' Setting',setting,'\n')
    for (asm in c('markov','semimarkov')){
    set.seed(2023)
    plot(NULL,NULL,xlim=c(0,8),ylim=c(0,1),xlab='Time',ylab='Cumulative incidence',
         main=paste('Setting', setting))
    mtext(asm, cex=0.8)
    for (i in 1:B){
      dat = generatedata(N,setting,h1,h0)
      fit = semicr(dat$A,dat$T,dat$dT,dat$R,dat$dR,dat$X,a=a,asm=asm)
      x = fit$time
      y = fit$F1+fit$F3
      x[x>8] = NA
      points(x,y,type='l',col='grey')
    }
    fit = semicr(dat$A,dat$T,dat$dT,dat$R,dat$dR,dat$X,a=a,asm=asm,conf.int=.95,nboot=0)
    x = fit$time
    x[x>8] = NA
    points(x,fit$F1+fit$F3,type='l',col='darkcyan')
    points(x,fit$ci_u,type='l',lty=2,col='darkcyan')
    points(x,fit$ci_l,type='l',lty=2,col='darkcyan')
    x = seq(0,8,length=100)
    h1a = h1[1]
    h2a = h1[2]
    h3a = h1[3]
    if (a[1]==0)  h1a = h0[1]
    if (a[2]==0)  h2a = h0[2]
    if (a[3]==0)  h3a = h0[3]
    y = meaninc(x,setting,h1a,h2a,h3a)
    points(x,y,type='l')
    }
  }
}


## Draw bias

B = 1000
xseq = seq(0,8,length=100)
for (N in c(100,500,2000)){
  for (alist in 1:3){
    if (alist==1) a = c(0,0,0)
    if (alist==2) a = c(1,0,0)
    if (alist==3) a = c(1,0,1)
    par(mfcol=c(1,3))
    for (setting in 1:3){
      if (setting==1){
        h0 = c(0.15,0.15,0.4)
        h1 = h0*2
      }
      if (setting==2){
        h0 = c(0.1,0.1,0.25)
        h1 = h0*2
      }
      if (setting==3){
        h0 = c(0.1,0.15,0.4)
        h1 = h0*2
      }
    cat('N =',N,' Setting',setting,' a =',a,'\n')
    set.seed(2024)
    plot(NULL,NULL,xlim=c(0,8),ylim=c(-.04,.07),xlab='Time',ylab='Bias',
             main=paste('Setting', setting))
    abline(h=0,lty=2)
        if (alist==1) mtext('a = (0,0,0)', cex=0.8)
        if (alist==2) mtext('a = (1,0,0)', cex=0.8)
        if (alist==3) mtext('a = (1,0,1)', cex=0.8)
        h1a = h1[1]
        h2a = h1[2]
        h3a = h1[3]
        if (a[1]==0) h1a = h0[1]
        if (a[2]==0) h2a = h0[2]
        if (a[3]==0) h3a = h0[3]
        yseq = meaninc(xseq,setting,h1a,h2a,h3a)
        bias1 = bias2 = bias3 = rep(0,100)
        rmse1 = rmse2 = rmse3 = rep(0,100)
        for (i in 1:B){
          dat = generatedata(N,setting,h1,h0)
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
        legend('topleft',cex=1,lwd=c(1.2,1.2,1.2),col=c('blue','darkcyan','brown'),
               legend=c('Proposed (ma.)','Proposed (sm.)','Huang'))
      }
  }
}