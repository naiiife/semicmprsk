setwd('D:/PAPER/Semi-competing risks/SPE_SemiCom0127/Codes')
dat = read.csv('leukemiaPKU.csv')
dat = dat[dat$MRD==1,]
A = 1-dat$TRANSPLANT
T = dat$OST
R = dat$RELAPSET
dT = dat$OS
dR = dat$RELAPSE
X = cbind(dat$AGE,dat$CR,dat$ALL)

testpath(A,T,dT,R,dR,X)
# 0.3060, 0.0105, 0.0513~0.3908, 0.0037, 0.0315

semicr.causal(A,T,dT,R,dR,X,a1=c(1,0,0),a0=c(0,0,0), 
              asm='semimarkov', sens=NULL, conf.int=.95, nboot=0,
              legend=c('NRM hazard under MSDT','NRM hazard under Haplo-SCT'),
              cex=1)
text(x=200,y=0.8, 'p = 0.3060')
semicr.causal(A,T,dT,R,dR,X,a1=c(1,1,0),a0=c(1,0,0), 
              asm='semimarkov', sens=NULL, conf.int=.95, nboot=0,
              legend=c('Relapse hazard under MSDT','Relapse hazard under Haplo-SCT'),
              cex=1)
text(x=200,y=0.8, 'p = 0.0105')
semicr.causal(A,T,dT,R,dR,X,a1=c(1,1,1),a0=c(1,1,0), 
              asm='semimarkov', sens=NULL, conf.int=.95, nboot=0,
              legend=c('RRM hazard under MSDT','RRM hazard under Haplo-SCT'),
              cex=1)
text(x=200,y=0.8, 'p = 0.3908')

semicr.sensitivity(A,T,dT,R,dR,X,a1=c(1,0,0),a0=c(0,0,0),
                   ylim=c(-0.1,0.3),cex=1)
semicr.sensitivity(A,T,dT,R,dR,X,a1=c(1,1,0),a0=c(1,0,0),
                   ylim=c(-0.1,0.3),cex=1)
semicr.sensitivity(A,T,dT,R,dR,X,a1=c(1,1,1),a0=c(1,1,0),
                   ylim=c(-0.1,0.3),cex=1)
