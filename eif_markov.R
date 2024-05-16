## Partial isolation, with baseline covariates

setwd('D:/PAPER/Semi-competing risks/SPE_SemiCom0420/Codes')
dat = read.csv('leukemiaPKU.csv')
dim(dat)
dat$A = dat$TRANSPLANT
dat$RRMT = dat$OST
dat$RRM = dat$OS*dat$RELAPSE
dat = dat[dat$MRD==1&dat$CR==1,]
dim(dat)
X = cbind(dat$ALL,dat$AGE,dat$SEX)
n = length(dat$A)
attach(dat)

library(survival)

phfit_3 = function(Tr,Dr,Td,Dd,A,X,a){
  Tr = Tr[A==a]
  Td = Td[A==a]
  Dr = Dr[A==a]
  Dd = Dd[A==a]
  X = X[A==a,]
  tt = sort(unique(Td[Dd==1]))
  l = length(tt)
  k = length(Td)
  lam = rep(1/l, l)
  p = ncol(X)
  beta = beta0 = rep(0, p)
  while(TRUE){
    Xb = as.numeric(X%*%beta)
    Lam = sapply(1:k, function(i) sum((tt<=Td[i])*(tt>=Tr[i])*(Dr[i]==1)*lam))
    Lam = Lam * exp(Xb)
    dbeta = t(X)%*%(Dd*Dr-Lam)
    ddbeta = -t(X)%*%diag(Lam)%*%X
    beta = beta - solve(ddbeta) %*% dbeta
    Xb = as.numeric(X%*%beta)
    lam = sapply(tt, function(t) sum(Dd*(Td==t)*Dr)/sum((Td>=t)*Dr*(Tr<=t)*
                                                          exp(Xb)))
    lam[is.nan(lam)] = 0
    tol = max(abs(beta-beta0))
    if (tol<0.0001) break
    beta0 = beta
  }
  return(list(beta=beta,tt=tt,lam=lam))
}

matchy = function(x,y,newx){
  newy = sapply(newx, function(t) y[which(x==t)])
  newy = as.numeric(newy)
  newy[is.na(newy)] = 0
  return(newy)
}


Td = dat$OST
Dd = dat$OS
Tr = dat$RELAPSET
Dr = dat$RELAPSE
dat$X = X
A = dat$A
p = ncol(X)
tt = sort(unique(c(Td,Tr)))
l = length(tt)

## bootstrap
F000.s = F100.s = F110.s = F111.s = NULL
for (b in 1:200){
ss = sample(n,n,replace=TRUE)
Td = dat$OST[ss]
Dd = dat$OS[ss]
Tr = dat$RELAPSET[ss]
Dr = dat$RELAPSE[ss]
A = dat$A[ss]
X = dat$X[ss,]

# Estimating RRM hazard
fit_d = phfit_3(Tr,Dr,Td,Dd,A,X,a=1)
Xb = as.numeric(X%*%fit_d$beta)
lam_d = matchy(fit_d$tt, fit_d$lam, tt)
lam_ord1 = sapply(1:l, function(t) lam_d[t]*exp(Xb))
fit_d = phfit_3(Tr,Dr,Td,Dd,A,X,a=0)
Xb = as.numeric(X%*%fit_d$beta)
lam_d = matchy(fit_d$tt, fit_d$lam, tt)
lam_ord0 = sapply(1:l, function(t) lam_d[t]*exp(Xb))

# Estimating TRM hazard
fit.TRM.1 = coxph(Surv(TRMT,TRM)~ALL+AGE+SEX, dat[ss,], subset=(A==1))
summary(fit.TRM.1)
exb = as.numeric(exp(X%*%fit.TRM.1$coefficients))
lam.TRM.1 = sapply(tt, function(t) sum((TRMT==t)*TRM*A)/sum((TRMT>=t)*exb*A))
lam_od1 = t(sapply(1:n, function(i) lam.TRM.1*exb[i]))
fit.TRM.0 = coxph(Surv(TRMT,TRM)~ALL+AGE+SEX, dat[ss,], subset=(A==0))
summary(fit.TRM.0)
exb = as.numeric(exp(X%*%fit.TRM.0$coefficients))
lam.TRM.0 = sapply(tt, function(t) sum((TRMT==t)*TRM*(1-A))/sum((TRMT>=t)*exb*(1-A)))
lam_od0 = t(sapply(1:n, function(i) lam.TRM.0*exb[i]))

# Estimating REL hazard
fit.REL.1 = coxph(Surv(RELAPSET,RELAPSE)~ALL+AGE+SEX, dat[ss,], subset=(A==1))
summary(fit.REL.1)
exb = as.numeric(exp(X%*%fit.REL.1$coefficients))
lam.REL.1 = sapply(tt, function(t) sum((RELAPSET==t)*RELAPSE*A)/
                     sum((RELAPSET>=t)*exb*A))
lam_or1 = t(sapply(1:n, function(i) lam.REL.1*exb[i]))
fit.REL.0 = coxph(Surv(RELAPSET,RELAPSE)~ALL+AGE+SEX, dat[ss,], subset=(A==0))
summary(fit.TRM.0)
exb = as.numeric(exp(X%*%fit.REL.0$coefficients))
lam.REL.0 = sapply(tt, function(t) sum((RELAPSET==t)*RELAPSE*(1-A))/
                     sum((RELAPSET>=t)*exb*(1-A)))
lam_or0 = t(sapply(1:n, function(i) lam.REL.0*exb[i]))

# Censoring probability
fit.c = coxph(Surv(OST,1-OS)~A+ALL+AGE+SEX, dat[ss,])
summary(fit.c)
exb = as.numeric(exp(cbind(A,X)%*%fit.c$coefficients))
lam.c = sapply(tt, function(t) sum((OST==t)*(1-OS))/
                 sum((OST>=t)*exb))
lam.c = t(sapply(1:n, function(i) lam.c*exb[i]))
SC = t(exp(-apply(lam.c, 1, cumsum)))

# Propensity score
fit = glm(A~X, family='binomial')
ps = matrix(1,nrow=n,ncol=2)
pscore = predict(fit, type='response')
ps[,1] = (1 - pscore) #* mean((1-A)/(1-pscore))
ps[,2] = pscore #* mean(A/pscore)

# observable incidence
lam_od_A = A*lam_od1 + (1-A)*lam_od0
lam_or_A = A*lam_or1 + (1-A)*lam_or0
lam_ord_A = A*lam_ord1 + (1-A)*lam_ord0
Lam_o_A = t(apply(lam_od_A+lam_or_A, 1, cumsum))
Lam_or_A = t(apply(lam_ord_A, 1, cumsum))

dF_or_A = exp(-Lam_o_A)*lam_or_A
dF_od_A = exp(-Lam_o_A)*lam_od_A
F_or_A = t(apply(dF_or_A, 1, cumsum))
F_od_A = t(apply(dF_od_A, 1, cumsum))
dF_or__A = t(apply(dF_or_A*exp(Lam_or_A), 1, cumsum))*exp(-Lam_or_A)
dF_ord_A = dF_or__A*lam_ord_A
F_ord_A = t(apply(dF_ord_A, 1, cumsum))
F_A = colMeans(F_od_A+F_ord_A)

EIF_d = function(a){
  # counterfactual incidence
  lam_od_a = a['od']*lam_od1 + (1-a['od'])*lam_od0
  lam_or_a = a['or']*lam_or1 + (1-a['or'])*lam_or0
  lam_ord_a = a['rd']*lam_ord1 + (1-a['rd'])*lam_ord0
  Lam_o_a = t(apply(lam_od_a+lam_or_a, 1, cumsum))
  Lam_or_a = t(apply(lam_ord_a, 1, cumsum))
  
  dF_or_a = exp(-Lam_o_a)*lam_or_a
  dF_od_a = exp(-Lam_o_a)*lam_od_a
  F_od_a = t(apply(dF_od_a, 1, cumsum))
  F_or_a = t(apply(dF_or_a, 1, cumsum))
  dF_or__a = t(apply(dF_or_a*exp(Lam_or_a), 1, cumsum))*exp(-Lam_or_a)
  dF_ord_a = dF_or__a*lam_ord_a
  F_ord_a = t(apply(dF_ord_a, 1, cumsum))
  Fr = F_od_a + F_ord_a
  
  # od
  Y_o = sapply(tt, function(t) (Td>=t)*(Tr>=t))
  PY_o = (1-F_od_A-F_or_A) * SC
  dM_od = sapply(tt, function(t) (Td==t)*Dd*(Tr>=t)*(1-Dr)) - Y_o*lam_od_A
  dQ_od = dM_od * (A==a['od'])/ps[,a['od']+1] / PY_o
  dM_or = sapply(tt, function(t) (Tr==t)*Dr*(Td>=t)) - Y_o*lam_or_A
  dQ_or = dM_or * (A==a['or'])/ps[,a['or']+1] / PY_o
  dQ_od[PY_o==0] = 0
  dQ_or[PY_o==0] = 0
  G1_od = dQ_od - t(apply(dQ_od+dQ_or,1,cumsum))*lam_od_a
  G_od = t(apply(G1_od*exp(-Lam_o_a), 1, cumsum))
  
  # ord
  Y_or = sapply(tt, function(t) (Td>=t)*(Tr<=t)*Dr)
  PY_or = (F_or_A - F_ord_A) * SC
  dM_ord = sapply(tt, function(t) (Td==t)*Dd*Dr) - Y_or*lam_ord_A
  dQ_ord = dM_ord * (A==a['rd'])/ps[,a['rd']+1] / PY_or
  dQ_ord[PY_or==0] = 0
  G1_or = dQ_or - t(apply(dQ_od+dQ_or,1,cumsum))*lam_or_a
  G1_or = t(apply(G1_or*exp(Lam_or_a-Lam_o_a), 1, cumsum))
  G2_or = dQ_ord - t(apply(dQ_ord,1,cumsum))*lam_ord_a
  G2_or1 = t(apply(exp(Lam_or_a-Lam_o_a)*lam_or_a, 1, cumsum)) 
  G2_or2 = G2_or1 * G2_or
  G_ord = t(apply((G1_or*lam_ord_a+G2_or2)*exp(-Lam_or_a), 1, cumsum))
  
  EIF = G_od + G_ord
  se = sqrt(apply(EIF^2, 2, mean) / n)
  Feff = Fr + EIF
  return(list(Freg=Fr, EIF=EIF, Feff=Feff, se=se))
}

a000 = c(0,0,0)
a100 = c(1,0,0)
a110 = c(1,1,0)
a111 = c(1,1,1)
names(a000)=names(a100)=names(a110)=names(a111)=c('od','or','rd')
fit000 = EIF_d(a000)
F000 = colMeans(fit000$Feff)
#se000 = fit000$se
fit100 = EIF_d(a100)
F100 = colMeans(fit100$Feff)
#se100 = fit100$se
fit110 = EIF_d(a110)
F110 = colMeans(fit110$Feff)
#se110 = fit110$se
fit111 = EIF_d(a111)
F111 = colMeans(fit111$Feff)
#se111 = fit111$se
F000.s = rbind(F000.s,F000)
F100.s = rbind(F100.s,F100)
F110.s = rbind(F110.s,F110)
F111.s = rbind(F111.s,F111)
}

se000 = apply(F000.s,2,sd,na.rm=TRUE)
se100 = apply(F100.s,2,sd,na.rm=TRUE)
se110 = apply(F110.s,2,sd,na.rm=TRUE)
se111 = apply(F111.s,2,sd,na.rm=TRUE)

ind = (tt<=2500)
plot(tt[ind],F000[ind],type='s',ylim=c(0,1),lwd=2,col='darkcyan',xlab='Time',
     ylab='Cumulative incidence', main='Incidence of the primary event')
points(tt[ind],F000[ind]+1.96*se000[ind],col='darkcyan',type='s',lty=2)
points(tt[ind],F000[ind]-1.96*se000[ind],col='darkcyan',type='s',lty=2)
points(tt[ind],F100[ind],type='s',lwd=2,col='brown')
points(tt[ind],F100[ind]+1.96*se100[ind],col='brown',type='s',lty=2)
points(tt[ind],F100[ind]-1.96*se100[ind],col='brown',type='s',lty=2)
legend('topleft', lwd=c(2,2), col=c('brown','darkcyan'),lty=c(1,1),
       legend=c('NRM hazard under Haplo-SCT','NRM hazard under MSDT'))

plot(tt[ind],F100[ind],type='s',ylim=c(0,1),lwd=2,col='darkcyan',xlab='Time',
     ylab='Cumulative incidence', main='Incidence of the primary event')
points(tt[ind],F100[ind]+1.96*se000[ind],col='darkcyan',type='s',lty=2)
points(tt[ind],F100[ind]-1.96*se000[ind],col='darkcyan',type='s',lty=2)
points(tt[ind],F110[ind],type='s',lwd=2,col='brown')
points(tt[ind],F110[ind]+1.96*se100[ind],col='brown',type='s',lty=2)
points(tt[ind],F110[ind]-1.96*se100[ind],col='brown',type='s',lty=2)
legend('topleft', lwd=c(2,2), col=c('brown','darkcyan'),lty=c(1,1),
       legend=c('Relapse hazard under Haplo-SCT','Relapse hazard under MSDT'))

plot(tt[ind],F110[ind],type='s',ylim=c(0,1),lwd=2,col='darkcyan',xlab='Time',
     ylab='Cumulative incidence', main='Incidence of the primary event')
points(tt[ind],F110[ind]+1.96*se000[ind],col='darkcyan',type='s',lty=2)
points(tt[ind],F110[ind]-1.96*se000[ind],col='darkcyan',type='s',lty=2)
points(tt[ind],F111[ind],type='s',lwd=2,col='brown')
points(tt[ind],F111[ind]+1.96*se100[ind],col='brown',type='s',lty=2)
points(tt[ind],F111[ind]-1.96*se100[ind],col='brown',type='s',lty=2)
legend('topleft', lwd=c(2,2), col=c('brown','darkcyan'),lty=c(1,1),
       legend=c('RRM hazard under Haplo-SCT','RRM hazard under MSDT'))
