#########NDCUSUM###########
epsilon;sigma0;sigma1
Kernelf1<-function(x){
  dnorm(x,sd=1)
}

fhat<-function(x,samplex,wh){
  return(sum((Kernelf1((samplex-x)/wh)))/(length(samplex)*wh) )
}

N<-1000

samplex<-numeric(N)
for(i in 1:N){
  samplex[i]<-rnorm(n=1, mean =0,sd=sample(c(sigma0, sigma1), size = 1, prob = c(1-epsilon, epsilon), replace = T))
}

wh=0.5

fhat(0,samplex,wh)

(1-epsilon)*dnorm(0)+epsilon*(dnorm(0,sd=3))

npf<-function(x){
  return(fhat(x,samplex,wh))
}
#plot(npf,type='l')

mu1<-2
npl<-function(x){
  log((npf(x-mu1))/(npf(x)))
}

#######plot###
xxx<-seq(-10,10,0.01)
yyy<-numeric(length(xxx))
for(i in 1:length(xxx)){yyy[i]<-npl(xxx[i])}
library(ggplot2)



xx<-seq(-20,20,by=0.01)
yy<-lr(xxx)
npyy<-yyy
#dataxy<-as.data.frame(cbind(xxx,yyy))
#rootnp<-as.data.frame(cbind(root_n,root_p))

ltype=c('npllr'=1,'llr'=2)
ONE<-ggplot(data=NULL,aes(xx,npyy),axis.line = element_line(colour='white'))+
  geom_line(aes(linetype='npllr'))+
  geom_line(aes(x=xx,y=yy,linetype='llr'))+
  scale_linetype_manual(name=NULL,values = ltype,labels=c("real LLR","LLR of kernel density"))+
  xlab("x")+
  ylab("LLR")+
  coord_cartesian(xlim =c(-9,11), ylim = c(-4,4))+
  theme_bw()


ONE+theme(legend.position=c(0.25,0.8),legend.key.size = unit(35, "pt"))

xx<-seq(-20,20,by=0.01)
yy<-lr(xxx)
npyy<-yyy
#dataxy<-as.d
ltype=c('npllr'=1,'llr'=2)
ONE<-ggplot(data=NULL,aes(xx,npyy),axis.line = element_line(colour='white'))+
  geom_line(aes(linetype='npllr'))+
  geom_line(aes(x=xx,y=yy,linetype='llr'))+
  scale_linetype_manual(name=NULL,values = ltype,labels=c("LLR","RLLR"))+
  xlab("x")+
  ylab("LLR")+
  coord_cartesian(xlim =c(-9,11), ylim = c(-4,4))+
  theme_bw()


ONE+theme(legend.position="none")

##############NDLLR########


##########ARL-monte-carlo#######

arl.cusum<-function(h,mu)
{
  t<-1000;
  half.w<-h/(2*t+1)
  r.mat<-matrix(0,nrow = t, ncol = t)
  r.f<-rep(0,2*t-1)
  
  dist.cdf<-function(x) (1-epsilon)*pnorm(x, mean = mu, sd = sigma0)+epsilon*pnorm(x, mean = mu, sd = sigma1)
  
  for(i in seq(from = 1, to = 2*t-1))
  {
    y<-(i-t)*2*half.w+half.w
    lr.y<-function(x) rnpl(x)-y
    temp<-as.numeric(uniroot(f = lr.y, lower = -20, upper = 20)$root)
    r.f[i]<-dist.cdf(temp)
  }
  
  for(i in seq(1,t))
    r.mat[i,1]<-r.f[t+1-i]
  
  for(j in seq(2,t))
  {
    for(i in seq(1,t))
      r.mat[i,j]<-r.f[j-i+t]-r.f[j-i-1+t]
  }
  
  one.vec<-rep(1,t)
  unit.mat<-diag(t)
  arl<-solve(unit.mat-r.mat,one.vec)
  return(arl[1])
}

minh=0
maxh=10
h=(minh+maxh)/2
arls<-arl.cusum(h=h,mu=0)
while((abs(arls-1000))>0.5){
  if(arls<1000){minh=h
  h=(minh+maxh)/2
  arls<-arl.cusum(h=h,mu=0)}
  else{maxh=h
  h=(minh+maxh)/2
  arls<-arl.cusum(h=h,mu=0)
  }
}
print(c(mu1,h,arls))
#hvector<-c(hvector,h)
#arlvector<-c(arlvector,arls)

testshift<-seq(0.5,10,0.5)
arlvector.rlcusum.fd<-matrix(data=NA,nrow=length(testshift),ncol=2)
arlvector.rlcusum.fd[,1]<-testshift

for(i in 1:length(testshift)){
  arlvector.rlcusum.fd[i,2]<-arl.cusum(h=h,mu=testshift[i])
}
arlvector.rlcusum.fd

###############cont normal##################

###############################
#############1000##############
#mu1=1
mu0<-0; mu1<-1;sigma0<-1;sigma1<-3;epsilon1<-0.1;
N<-1000


samplex<-read.csv('C:\\Users\\www23\\Desktop\\Research\\Paper1\\data_contnormal_1000.csv')
samplex<-samplex[,2]

wh=0.5
Kernelf1<-function(x){
  dnorm(x,sd=1)
}

fhat<-function(x,samplex,wh){
  return(sum((Kernelf1((samplex-x)/wh)))/(length(samplex)*wh) )
}


npf0<-function(x){
  return(fhat(x,samplex,wh))
}
den<-density(samplex,bw=wh,kernel='gaussian')
npf1<- approxfun(den$x,den$y)

npf<-function(x){
  if(is.na(npf1(x))==1){
    return(npf0(x))
  }else{
    return(npf1(x))
  }
}

npl<-function(x){
  return(log(npf(x-mu1)/npf(x)))
}

xxx<-seq(-10,10,0.01)
yyy<-numeric(length(xxx))
for(i in 1:length(xxx)){yyy[i]<-npl(xxx[i])}

plot(xxx,yyy,type='l',ylim=c(-4,4),xlim=c(-10,10))
##################################################

lr<-function(x)
{
  log((1-epsilon)*dnorm(x, mean=mu1, sd=sigma0)+epsilon*dnorm(x, mean=mu1, sd=sigma1))-
    log((1-epsilon)*dnorm(x, mean=mu0, sd=sigma0)+epsilon*dnorm(x, mean=mu0, sd=sigma1))
}


xx<-seq(-10,10,0.01)
yy<-lr(xx)
lines(xx,yy,lty=2)

abline(h=0,v=0,lty=3)

##################################################

x1=-5.1441930;x2=-1.889062;x3=2.98091;x4=4.4069490;x5=5.3023680;x6=6.1977880
rnpl<-function(x){
  if(x<x1){return(npl(x)-(npl(x1)-npl(x2)))}
  else if(x<x2){return(npl(x2))}
  else if(x<x3){return(npl(x))}
  else if(x<x4){return(npl(x3))}
  else if(x<x5){return(npl(x)+npl(x3)-npl(x4))}
  else if(x<x6){return(npl(x5)+npl(x3)-npl(x4))}
  else {return(npl(x)+npl(x5)-npl(x6)+npl(x3)-npl(x4))}
}

xxx<-seq(-10,10,0.01)
yyy<-numeric(length(xxx))
for(i in 1:length(xxx)){yyy[i]<-rnpl(xxx[i])}

plot(xxx,yyy,type='l',ylim=c(-4,4))

###############################################

arl.cusum<-function(h,mu)
{
  t<-1000;
  half.w<-h/(2*t+1)
  r.mat<-matrix(0,nrow = t, ncol = t)
  r.f<-rep(0,2*t-1)
  
  dist.cdf<-function(x) (1-epsilon)*pnorm(x, mean = mu, sd = sigma0)+epsilon*pnorm(x, mean = mu, sd = sigma1)
  
  for(i in seq(from = 1, to = 2*t-1))
  {
    y<-(i-t)*2*half.w+half.w
    lr.y<-function(x) rnpl(x)-y
    temp<-as.numeric(uniroot(f = lr.y, lower = -20, upper = 20)$root)
    r.f[i]<-dist.cdf(temp)
  }
  
  for(i in seq(1,t))
    r.mat[i,1]<-r.f[t+1-i]
  
  for(j in seq(2,t))
  {
    for(i in seq(1,t))
      r.mat[i,j]<-r.f[j-i+t]-r.f[j-i-1+t]
  }
  
  one.vec<-rep(1,t)
  unit.mat<-diag(t)
  arl<-solve(unit.mat-r.mat,one.vec)
  return(arl[1])
}

#####################################################

NBoot=10000
rlv=numeric(NBoot)
h=4.985
S=0
num=0

for (i in 1:NBoot){
  S=0
  num=0
  while(S<=h){
    num=num+1
    znow<-samplex[ceiling(runif(1,min=0,max=1000))]
    S=max(0,S+rnpl(znow))
  }
  rlv[i]=num
}
print(c(mean(rlv),sd(rlv)))

arl.cusum(h=h,mu=0)





#########################################################
#########################################################
#########################################################

###############cont normal##################

###############################
#############1000##############
#mu1=2
mu0<-0; mu1<-2;sigma0<-1;sigma1<-3;epsilon1<-0.1;
N<-1000


samplex<-read.csv('C:\\Users\\www23\\Desktop\\Research\\Paper1\\data_contnormal_1000.csv')
samplex<-samplex[,2]

wh=0.5
Kernelf1<-function(x){
  dnorm(x,sd=1)
}

fhat<-function(x,samplex,wh){
  return(sum((Kernelf1((samplex-x)/wh)))/(length(samplex)*wh) )
}


npf0<-function(x){
  return(fhat(x,samplex,wh))
}
den<-density(samplex,bw=wh,kernel='gaussian')
npf1<- approxfun(den$x,den$y)

npf<-function(x){
  if(is.na(npf1(x))==1){
    return(npf0(x))
  }else{
    return(npf1(x))
  }
}

npl<-function(x){
  return(log(npf(x-mu1)/npf(x)))
}

xxx<-seq(-10,10,0.01)
yyy<-numeric(length(xxx))
for(i in 1:length(xxx)){yyy[i]<-npl(xxx[i])}

plot(xxx,yyy,type='l',ylim=c(-4,4),xlim=c(-10,10))
##################################################

lr<-function(x)
{
  log((1-epsilon)*dnorm(x, mean=mu1, sd=sigma0)+epsilon*dnorm(x, mean=mu1, sd=sigma1))-
    log((1-epsilon)*dnorm(x, mean=mu0, sd=sigma0)+epsilon*dnorm(x, mean=mu0, sd=sigma1))
}


xx<-seq(-10,10,0.01)
yy<-lr(xx)
lines(xx,yy,lty=2)

abline(h=0,v=0,lty=3)

##################################################

x1=-5.0345240;x2=-1.718155;x3=3.246219;x4=5.832987;x5=6.595752;x6=6.970732
rnpl<-function(x){
  if(x<x1){return(npl(x)-(npl(x1)-npl(x2)))}
  else if(x<x2){return(npl(x2))}
  else if(x<x3){return(npl(x))}
  else if(x<x4){return(npl(x3))}
  else if(x<x5){return(npl(x)+npl(x3)-npl(x4))}
  else if(x<x6){return(npl(x5)+npl(x3)-npl(x4))}
  else {return(npl(x)+npl(x5)-npl(x6)+npl(x3)-npl(x4))}
}

xxx<-seq(-10,10,0.01)
yyy<-numeric(length(xxx))
for(i in 1:length(xxx)){yyy[i]<-rnpl(xxx[i])}

plot(xxx,yyy,type='l',ylim=c(-4,4))

###############################################

arl.cusum<-function(h,mu)
{
  t<-1000;
  half.w<-h/(2*t+1)
  r.mat<-matrix(0,nrow = t, ncol = t)
  r.f<-rep(0,2*t-1)
  
  dist.cdf<-function(x) (1-epsilon)*pnorm(x, mean = mu, sd = sigma0)+epsilon*pnorm(x, mean = mu, sd = sigma1)
  
  for(i in seq(from = 1, to = 2*t-1))
  {
    y<-(i-t)*2*half.w+half.w
    lr.y<-function(x) rnpl(x)-y
    temp<-as.numeric(uniroot(f = lr.y, lower = -20, upper = 20)$root)
    r.f[i]<-dist.cdf(temp)
  }
  
  for(i in seq(1,t))
    r.mat[i,1]<-r.f[t+1-i]
  
  for(j in seq(2,t))
  {
    for(i in seq(1,t))
      r.mat[i,j]<-r.f[j-i+t]-r.f[j-i-1+t]
  }
  
  one.vec<-rep(1,t)
  unit.mat<-diag(t)
  arl<-solve(unit.mat-r.mat,one.vec)
  return(arl[1])
}

#####################################################

NBoot=100000
rlv=numeric(NBoot)
h=5.443
S=0
num=0

for (i in 1:NBoot){
  S=0
  num=0
  while(S<=h){
    num=num+1
    znow<-samplex[ceiling(runif(1,min=0,max=1000))]
    S=max(0,S+rnpl(znow))
  }
  rlv[i]=num
}
print(c(mean(rlv),sd(rlv)))

arl.cusum(h=h,mu=0)

#########################################################
#########################################################
#########################################################

###############t##################

N<-1000
v=3;tmu0=0;tmu1=1;tsigma=1/sqrt(3)

samplex<-read.csv('C:\\Users\\www23\\Desktop\\Research\\Paper1\\data_t_1000.csv')
samplex<-samplex[,2]






















