# Testing separable effects

B = 10000
for (N in c(100,500,2000)){
  cat('Power of dhaz1, dhaz2, dhaz3(markov), dhaz3(semimarkov)\n')
  cat('N =',N,'\n')
  for (setting in 1:3){
    cat('Setting',setting,'\n')
    for (s in 1:5){
      p1 = p2 = p3 = p4 = 0
      Z1 = 1-as.numeric(s==2|s==5)
      Z2 = 1-as.numeric(s==3|s==5)
      Z3 = 1-as.numeric(s==4|s==5)
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
    }
  }
}


# Confidence interval coverage

B = 10000
a = c(0,0,0)
xseq = 1:8
for (N in c(100,500,2000)){
  cat('Coverage (markov/semimarkov), Width (markov/semimarkov)\n')
  cat('N =',N,'\n')
  for (setting in 1:3){
    cat('Setting',setting,'\n')
    set.seed(2024)
    count.mem = count.mar = rep(0,length(xseq))
    width.mem = width.mar = rep(0,length(xseq))
    yt = meaninc(xseq,setting,a)
    for (i in 1:B){
      dat = generatedata(N,setting)
      fit = semicr(dat$A,dat$T,dat$dT,dat$R,dat$dR,dat$X,a=a,asm='markov',conf.int=.95)
      x = fit$time
      index = sapply(xseq, function(v) max(which(x<=v)))
      upper = fit$ci_u[index]
      lower = fit$ci_l[index]
      if (is.na(upper)) upper = 0
      if (is.na(lower)) lower = 0
      count.mem = count.mem + (upper>=yt & lower<=yt)
      width.mem = width.mem + (upper-lower)
      fit = semicr(dat$A,dat$T,dat$dT,dat$R,dat$dR,dat$X,a=a,asm='semimarkov',conf.int=.95)
      x = fit$time
      index = sapply(xseq, function(v) max(which(x<=v)))
      upper = fit$ci_u[index]
      lower = fit$ci_l[index]
      if (is.na(upper)) upper = 0
      if (is.na(lower)) lower = 0
      count.mar = count.mar + (upper>=yt & lower<=yt)
      width.mar = width.mar + (upper-lower)
    }
    cat(a,'\n')
    print(count.mem/B)
    print(count.mar/B)
    print(round(width.mem/B,4))
    print(round(width.mar/B,4))
  }
}
