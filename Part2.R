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








