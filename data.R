setwd('D:/PAPER/Semi-competing risks/SPE_SemiCom0420/Codes')
source('semicompeting.R')
dat = read.csv('leukemiaPKU.csv')
dat = dat[dat$MRD==1&dat$CR==1,]
A = dat$TRANSPLANT
T = dat$OST
R = dat$RELAPSET
dT = dat$OS
dR = dat$RELAPSE
X = cbind(dat$AGE,dat$SEX,dat$ALL)

tp = testpath(A,T,dT,R,dR,X)

# Compare counterfactual cumulative incidences
semicr.causal(A,T,dT,R,dR,X,a1=c(1,0,0),a0=c(0,0,0), 
              asm='semimarkov', sens=NULL, conf.int=.95, nboot=0,
              legend=c('NRM hazard under Haplo-SCT','NRM hazard under MSDT'),
              cex=1)
text(x=200,y=0.8, paste0('p = ',round(tp$p1,4)))
semicr.causal(A,T,dT,R,dR,X,a1=c(1,1,0),a0=c(1,0,0), 
              asm='semimarkov', sens=NULL, conf.int=.95, nboot=0,
              legend=c('Relapse hazard under Haplo-SCT','Relapse hazard under MSDT'),
              cex=1)
text(x=200,y=0.8, paste0('p = ',round(tp$p2,4)))
semicr.causal(A,T,dT,R,dR,X,a1=c(1,1,1),a0=c(1,1,0), 
              asm='semimarkov', sens=NULL, conf.int=.95, nboot=0,
              legend=c('RRM hazard under Haplo-SCT','RRM hazard under MSDT'),
              cex=1)
text(x=200,y=0.8, paste0('p = ',round(tp$p3[2],4)))

# Sensitivity analysis
semicr.sensitivity(A,T,dT,R,dR,X,a1=c(1,0,0),a0=c(0,0,0),
                   ylim=c(-0.3,0.3),cex=1)
semicr.sensitivity(A,T,dT,R,dR,X,a1=c(1,1,0),a0=c(1,0,0),
                   ylim=c(-0.3,0.3),cex=1)
semicr.sensitivity(A,T,dT,R,dR,X,a1=c(1,1,1),a0=c(1,1,0),
                   ylim=c(-0.3,0.3),cex=1)

fit.km = survfit(Surv(T,dT)~A, weights=rowSums(pscore(A,X)))
tt=fit.km$time
tt0=tt[1:54]
tt1=tt[55:226]
tt=sort(unique(tt))
cif0=matchy(cif0,tt0,tt)
cif1=matchy(cif1,tt1,tt)
cif0=1-fit.km$surv[1:54]
cif1=1-fit.km$surv[55:226]
semicr.causal(A,T,dT,R,dR,X,a1=c(1,1,1),a0=c(0,0,0), 
              asm='markov', sens=NULL, conf.int=.95, nboot=0,
              legend=c('Hazards under Haplo-SCT','Hazards under MSDT'),
              cex=1)
points(tt[tt<2500],cif0[tt<2500],type='s',lwd=1.5)
points(tt[tt<2500],cif1[tt<2500],type='s',lwd=1.5)
semicr.causal(A,T,dT,R,dR,X,a1=c(1,1,1),a0=c(0,0,0), 
              asm='semimarkov', sens=NULL, conf.int=.95, nboot=0,
              legend=c('Hazards under Haplo-SCT','Hazards under MSDT'),
              cex=1)
points(tt[tt<2500],cif0[tt<2500],type='s',lwd=1.5)
points(tt[tt<2500],cif1[tt<2500],type='s',lwd=1.5)
semicr.sensitivity(A,T,dT,R,dR,X,a1=c(1,1,1),a0=c(0,0,0),
                   ylim=c(-0.3,0.3),cex=1)
points(tt[tt<2500],(cif1-cif0)[tt<2500],type='s',lwd=1.5)

# Partial isolation
semixcr.causal(A,T,dT,R,dR,X,a1=c(1,0,0),a0=c(0,0,0), 
              conf.int=.95, nboot=100,
              legend=c('NRM hazard under Haplo-SCT','NRM hazard under MSDT'),
              cex=1)
semixcr.causal(A,T,dT,R,dR,X,a1=c(1,1,0),a0=c(1,0,0), 
              conf.int=.95, nboot=100,
              legend=c('Relapse hazard under Haplo-SCT','Relapse hazard under MSDT'),
              cex=1)
semixcr.causal(A,T,dT,R,dR,X,a1=c(1,1,1),a0=c(1,1,0), 
              conf.int=.95, nboot=100,
              legend=c('RRM hazard under Haplo-SCT','RRM hazard under MSDT'),
              cex=1)
