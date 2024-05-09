
library(ggplot2)
library(Rfast)
library(geostats)
library(AMISEpi) # https://github.com/rretkute/AMISEpi
library(matrixStats)


##  Simulate angles
set.seed(1)
sm.ang.sngl <- rvonmises(1000, pi, 5, rads = TRUE)
sm.ang.sngl<-sm.ang.sngl*180/pi

hist(sm.ang.sngl, 100, xlim=c(0, 360), main="",
     xlab = "Angle")

##  Importance sampling
param.lim<-data.frame(min=c(0, 0), max=c(360, 700))

ans<-data.frame(m=c(), k=c(), LL=c())

for(ii in 1:(10^5)){
  m<-runif(1, min=param.lim$min[1], max=param.lim$max[1])
  k<-runif(1, min=param.lim$min[1], max=700)
  dens<-geostats::vonMises(sm.ang.sngl, m, k)
  y=sum(log(dens))-log(1/param.lim$max[1])-log(1/param.lim$max[2])
  ans<-rbind(ans, data.frame(m=m, k=k, LL=y))
  if(ii %% 1000==0) cat(c(ii, ""))
}

ess<-data.frame(nn=c(), ESS=c())
for(n in 10^seq(2, log10(nrow(ans)), 0.1)){
  wl<- ans$LL[1:round(n)]
  CC<-logSumExp(wl)
  WW<-exp(wl-CC)
  WW<-WW/sum(WW)
  ess<-rbind(ess, data.frame(nn=round(n), ESS=(sum((WW)^2))^(-1)))
}

ggplot(ess, aes(x=nn, y=ESS))+
  geom_point(fill="lightblue", pch=21, size=2) +
  scale_x_continuous(trans='log10')+
  xlab('Number of samples')+ ylab('ESS')+
  theme_bw()


##  Adaptive Importance Sampling

ESS.R<-100
n.param<-2
NN<-1000

# Prior PDF value 
dprop0<-function(pp, param.lim){
  ans<-dunif(pp[1], min=param.lim$min[1], max=param.lim$max[1], log=FALSE)
  if(length(pp)>1){
    for(i in 2:length(pp)){
      ans<-ans*dunif(pp[i], min=param.lim$min[i], max=param.lim$max[i], log=FALSE)
    }
  }
  return(ans)
}

prior.dns<-log(1/param.lim$max[1])+log(1/param.lim$max[2])

#  Set up for mixture function as a multivariate Students distribution
proposal=mvtComp(df=3); mixture=mclustMix();
dprop <- proposal$d
rprop <- proposal$r


## initial iteration 
it=1
ans<-data.frame(m=c(), k=c(), LL=c(), it=c())

while(nrow(ans)<NN){
  m<-runif(1, min=0, max=360)
  k<-runif(1, min=0, max=700)
  dens<-geostats::vonMises(sm.ang.sngl, m, k)
  y=sum(log(dens))-log(1/param.lim$max[1])-log(1/param.lim$max[2])
  if(y>-Inf){
    ans<-rbind(ans, data.frame(m=m, k=k, LL=y, it=it))
  }
}

ess<-data.frame(nn=c(), ESS=c())
for(n in seq(NN, nrow(ans), NN)){
  wl<- ans$LL[1:round(n)]
  CC<-logSumExp(wl)
  WW<-exp(wl-CC)
  WW<-WW/sum(WW)
  ess<-rbind(ess, data.frame(nn=round(n), ESS=(sum((WW)^2))^(-1), it=1))
}

ess

param<-ans[, 1:2]
wl<- ans$LL
CC<-logSumExp(wl)
WW<-exp(wl-CC)
WW<-WW/sum(WW)


# Iterations 2+
stop<-0
while(stop==0){
  it<-it+1
  cat(c("\n Started iteration ", it,"\n"))
  J<-unique(sample(1:nrow(param), NN, prob=WW, replace=TRUE))
  if(length(J)<3) J<-c(J, sample(1:nrow(param), 2))
  xx<-as.matrix(param[J,1:n.param])
  clustMix <- mixture(xx)
  G <- clustMix$G
  cluster <- clustMix$cluster
  ### Components of the mixture
  ppt <- clustMix$alpha
  muHatt <- clustMix$muHat
  varHatt <- clustMix$SigmaHat
  
  param<-data.frame(m=c(), k=c())
  # Draw new parameters and calculate log likelihood
  while(nrow(param)<NN){
    compo <- sample(1:G,1,prob=ppt)
    x1 <- t(rprop(1,muHatt[compo,], varHatt[,,compo]))
    new.param<-as.numeric(x1)
    if(dprop0(new.param, param.lim)>0){
      m<-new.param[1]; k<-new.param[2]
      dens<-geostats::vonMises(sm.ang.sngl, m, k)
      y=sum(log(dens))-log(1/param.lim$max[1])-log(1/param.lim$max[2])
      if(y>-Inf){
        ans<-rbind(ans, data.frame(m=m, k=k, LL=y, it=it))
        param<-rbind(param, data.frame(m=m, k=k))
      }
    }
  }
  
  wl<- ans$LL[ans$it==it]
  CC<-logSumExp(wl)
  WW<-exp(wl-CC)
  WW<-WW/sum(WW)
  ess<-rbind(ess, data.frame(nn=nrow(ans), ESS=(sum((WW)^2))^(-1), it=it))
  lst.ess<-ess$ESS[ess$it==it]
  print(paste0("ESS=", lst.ess))
  if(min(lst.ess)>=ESS.R) stop<-1
  if(it>= 100) stop<-1
}

ggplot(ess, aes(x=nn, y=ESS))+
  geom_path(col="lightblue") +
  geom_point(fill="lightblue", pch=21, size=2) +
  geom_hline(yintercept=ESS.R, linetype="dashed", color = "black")+
  scale_x_continuous(trans='log10')+
  xlab('Number of samples')+ ylab('ESS')+
  theme_bw() 




