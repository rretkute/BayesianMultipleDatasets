
library(ggplot2)
library(Rfast)
library(geostats)
library(AMISEpi) # https://github.com/rretkute/AMISEpi
library(matrixStats)

##  Simulate angles

nn<-1000
prm.m.true<-c(30, 110, 240)
prm.k.true<-c(100, 20, 5)
sim.angl<-matrix(NA, ncol=nn, nrow=length(prm.m.true))
set.seed(1)
for(ii in 1:length(prm.m.true)){
  m<-prm.m.true[ii]*pi/180
  k<-prm.k.true[ii]
  x <- rvonmises(nn, m, k, rads = TRUE)
  sim.angl[ii,]<- x*180/pi
}

sim.angl.all<-rbind(data.frame(Angle=sim.angl[1,], Sets="Set 1"),
                    data.frame(Angle=sim.angl[2,], Sets="Set 2"),
                    data.frame(Angle=sim.angl[3,], Sets="Set 3"))

ggplot(sim.angl.all, aes(x=Angle))+
  geom_histogram(aes(fill=Sets), bins=100) +
  theme_bw()


## Adaptive Multiple Importance Sampling

# Number of data in each dataset
n.dat<-nrow(sim.angl)
# Number of parameter sets to sample each iteration
nn.parsets<-10^3
# Max number of iterations
NN<-rep(nn.parsets,100)
# Required effective sample size
ESS.R<-100
# Dimension of parameter vector
n.param<-2
# Range of parameters
param.lim<-data.frame(min=c(0, 0), max=c(360, 700))


#  Set up for mixture function as a multivariate Students distribution
proposal=mvtComp(df=3); mixture=mclustMix();
dprop <- proposal$d
rprop <- proposal$r


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

# Sampling from prior
rprop0<-function(param.lim){
  ans<-sapply(1:nrow(param.lim), function(a) runif(1, min=param.lim$min[a],
                                                   max=param.lim$max[a]))
  return(ans)
}

prior.dns<-log(1/param.lim$max[1])+log(1/param.lim$max[2])

Sigma <- list(NA)
Mean<-list(NA)
PP<-list(NA)
GG<-list(NA)
ess<-data.frame(it=c(), Sets=c(), ESS=c())

# Initial iteration
set.seed(1)
it<-1
cat(c("Started iteration ", it,"\n"))

param<-data.frame(m=c(), k=c())
DENS<-matrix(NA, nrow=n.dat, ncol=0)

for(ii in 1:NN[it]){
  new.param<-rprop0(param.lim)
  m<-new.param[1]; k<-new.param[2]
  param<-rbind(param, data.frame(m=m, k=k))
  dens<-geostats::vonMises(unlist(as.list(t(sim.angl*2*pi/360))),
                           m*2*pi/360, k)
  dens<-matrix(dens, ncol=nn, byrow=TRUE)
  DENS<-cbind(DENS, rowSums(log(dens)))
}

WW<-0*DENS
for(ii in 1:3){
  wl<- DENS[ii,]-prior.dns
  CC<-logSumExp(wl)
  ww<-exp(wl-CC)
  ww<-ww/sum(ww)
  WW[ii, ]<-ww
  ess<-rbind(ess, data.frame(it=it, Sets=paste("Set ", ii), ESS=(sum((ww)^2))^(-1)))
}
ess
lst.ess<-ess$ESS

it<-1
# Iterations 2+
stop<-0
while(stop==0){
  it<-it+1
  cat(c("\n Started iteration ", it,"\n"))
  wh.dt<-which(lst.ess<ESS.R)
  if(length(wh.dt)>1){
    J<-unique(sample(1:sum(NN[1:(it-1)]), 1000, prob=colMeans(WW[wh.dt,]), replace=TRUE))
  } else {
    J<-unique(sample(1:sum(NN[1:(it-1)]), 1000, prob=WW[wh.dt,], replace=TRUE))
  }
  xx<-as.matrix(param[J,1:n.param])
  clustMix <- mixture(xx)
  G <- clustMix$G
  cluster <- clustMix$cluster
  ### Components of the mixture
  ppt <- clustMix$alpha
  muHatt <- clustMix$muHat
  varHatt <- clustMix$SigmaHat
  GG[[it-1]]<-G
  G1<-0; G2<-G
  if(it>2) {
    G1<-sum(sapply(1:(it-2), function(a) GG[[a]]))
    G2<-sum(sapply(1:(it-1), function(a) GG[[a]]))
  }
  for(i in 1:G){
    Sigma[[i+G1]] <- varHatt[,,i]
    Mean[[i+G1]] <- muHatt[i,]
    PP[[i+G1]]<-ppt[i]
  }
  
  # Draw new parameters and calculate log likelihood
  while(nrow(param)<sum(NN[1:it])){
    compo <- sample(1:G,1,prob=ppt)
    x1 <- t(rprop(1,muHatt[compo,], varHatt[,,compo]))
    new.param<-as.numeric(x1)
    if(dprop0(new.param, param.lim)>0){
      m<-new.param[1]; k<-new.param[2]
      dens<-geostats::vonMises(unlist(as.list(t(sim.angl*2*pi/360))),
                               m*2*pi/360, k)
      dens<-matrix(dens, ncol=nn, byrow=TRUE)
      LL<-rowSums(log(dens))
      if(min(LL)>-Inf){
        DENS<-cbind(DENS, LL)
        param<-rbind(param, data.frame(m=m, k=k))
      }
    }
  }
  
  q <- (NN[1]/sum(NN[1:it]))*exp(prior.dns) +
    (sum(NN[2:it])/sum(NN[1:it]))* rowSums(as.matrix(sapply(1:G2, 
                                                            function(a) PP[[a]] * dprop(param[,1:n.param],mu= Mean[[a]], Sig=Sigma[[a]]))))
  
  #Calculate ESS and weighst
  WW<-0*DENS
  for(ii in 1:3){
    wl<- DENS[ii,]-log(q)
    CC<-logSumExp(wl)
    ww<-exp(wl-CC)
    ww<-ww/sum(ww)
    WW[ii, ]<-ww
    ess<-rbind(ess, data.frame(it=it, Sets=paste("Set ", ii), ESS=(sum((ww)^2))^(-1)))
  }
  lst.ess<-round(ess$ESS[ess$it==it])
  print(paste0("ESS1=", lst.ess[1], ", ESS2=", lst.ess[2], ", ESS3=", lst.ess[3]))
  if(min(lst.ess)>=ESS.R) stop<-1
  if(it>= length(NN)) stop<-1
}

ggplot(ess, aes(x=it*1000, y=ESS))+
  geom_path(aes(col=Sets))+
  geom_point(aes(fill=Sets), pch=21)+
  geom_hline(yintercept=ESS.R, linetype="dashed", color = "black")+
  scale_x_continuous(trans='log10')+
  xlab('Number of samples')+ ylab('ESS')+
  theme_bw()

