#########0##########
mu0<-0; mu1<-2;sigma0<-1;sigma1<-3;epsilon1<-0.1;
epsilon=epsilon1
hvector<-numeric(0);arlvector<-numeric(0)
#########cusum#########
####log-likelihood ratio function
lr<-function(x)
{
  log((1-epsilon)*dnorm(x, mean=mu1, sd=sigma0)+epsilon*dnorm(x, mean=mu1, sd=sigma1))-
    log((1-epsilon)*dnorm(x, mean=mu0, sd=sigma0)+epsilon*dnorm(x, mean=mu0, sd=sigma1))
}

####ln(F(x|mu1)/F(x|mu0))
fn.cdf.lr <- function(x)
{
  log((1-epsilon)*pnorm(x, mean=mu1, sd=sigma0)+epsilon*pnorm(x, mean=mu1, sd=sigma1))-
    log((1-epsilon)*pnorm(x, mean=mu0, sd=sigma0)+epsilon*pnorm(x, mean=mu0, sd=sigma1))
}

####ln((1-F(x|mu1))/(1-F(x|mu0)))
fn.surv.lr <- function(x)
{
  log((1-epsilon)*pnorm(x, mean=mu1, sd=sigma0, lower.tail=FALSE)+epsilon*pnorm(x, mean=mu1, sd=sigma1, lower.tail=FALSE))-
    log((1-epsilon)*pnorm(x, mean=mu0, sd=sigma0, lower.tail=FALSE)+epsilon*pnorm(x, mean=mu0, sd=sigma1, lower.tail=FALSE))
}

####
lr.neg<-function(x) -lr(x)
ex.pt1<-nlm(f = lr.neg, p = -3.5)
ex.pt1$extr.val<--ex.pt1$minimum

ex.pt2<-nlm(f = lr, p = -1.5)
ex.pt2$extr.val<-ex.pt2$minimum

ex.pt3<-nlm(f = lr.neg, p = 3)
ex.pt3$extr.val<--ex.pt3$minimum

ex.pt4<-nlm(f = lr, p = mu1+3)
ex.pt4$extr.val<-ex.pt4$minimum

ex.pt<-rbind(ex.pt1, ex.pt2, ex.pt3, ex.pt4)
ex.pt[, c(2,6)]

extr.pt<-matrix(as.numeric(ex.pt[,c(2,6)]), nrow = 4, ncol = 2)



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
    lr.y<-function(x) lr(x)-y
    if(y>extr.pt[2,2] && y<extr.pt[1,2]) 
    {
      temp1<-as.numeric(uniroot(f = lr.y, lower = -100, upper = extr.pt[1,1])$root)
      temp2<-as.numeric(uniroot(f = lr.y, lower = extr.pt[1,1], upper = extr.pt[2,1])$root)
      temp3<-as.numeric(uniroot(f = lr.y, lower = extr.pt[2,1], upper = extr.pt[3,1])$root)
      r.f[i]<-dist.cdf(temp1)+dist.cdf(temp3)-dist.cdf(temp2)
    }
    else if(y>extr.pt[4,2] && y<extr.pt[3,2])
    {
      temp1<-as.numeric(uniroot(f = lr.y, lower = extr.pt[2,1], upper = extr.pt[3,1])$root)
      temp2<-as.numeric(uniroot(f = lr.y, lower = extr.pt[3,1], upper = extr.pt[4,1])$root)
      temp3<-as.numeric(uniroot(f = lr.y, lower = extr.pt[4,1], upper = 100)$root)
      r.f[i]<-dist.cdf(temp1)+dist.cdf(temp3)-dist.cdf(temp2)
    }
    
    else
    {
      temp<-as.numeric(uniroot(f = lr.y, lower = -100, upper = 100)$root)
      r.f[i]<-dist.cdf(temp)
    }
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
hvector<-c(hvector,h)
arlvector<-c(arlvector,arls)

testshift<-c(0.5,0.7,0.9,1.2,1.5,1.8,2.1,2.5,3,4,5,6,7,8,9,10)
arlvector.cusum<-matrix(data=NA,nrow=length(testshift),ncol=2)
arlvector.cusum[,1]<-testshift

for(i in 1:length(testshift)){
  arlvector.cusum[i,2]<-arl.cusum(h=h,mu=testshift[i])
}
arlvector.cusum

#########Rlcusum#######

cont.cdf<-function(x,mu) (1-epsilon)*pnorm(x, mean = mu, sd = sigma0)+epsilon*pnorm(x, mean = mu, sd = sigma1)
y1<-seq(from = round(extr.pt[2,2],2)+0.01, to = round(extr.pt[1,2],2)-0.01, by = 0.0001)
j<-length(y1)
a<-rep(0,j); b<-rep(0,j); cutoff.l<-rep(0,j)

for(i in seq(1,j))
{
  lr.y1<-function(x) lr(x)-y1[i]
  a[i]<-as.numeric(uniroot(f = lr.y1, lower = -50, upper = extr.pt[1,1])$root)
  b[i]<-as.numeric(uniroot(f = lr.y1, lower = extr.pt[2,1], upper = extr.pt[3,1])$root)
  cutoff.l[i]<-log((cont.cdf(b[i],mu1)-cont.cdf(a[i],mu1))/(cont.cdf(b[i],mu0)-cont.cdf(a[i],mu0)))
}
id.l<-sort(abs(cutoff.l-y1), index.return = T)$ix[1]

cont.cdf<-function(x,mu) (1-epsilon)*pnorm(x, mean = mu, sd = sigma0)+epsilon*pnorm(x, mean = mu, sd = sigma1)
y2<-seq(from = round(extr.pt[4,2]+0.01,2), to = round(extr.pt[3,2],2)-0.01, by = 0.0001)
j<-length(y2)
c<-rep(0,j); d<-rep(0,j); cutoff.r<-rep(0,j)

for(i in seq(1,j))
{
  lr.y2<-function(x) lr(x)-y2[i]
  c[i]<-as.numeric(uniroot(f = lr.y2, lower = extr.pt[2,1], upper = extr.pt[3,1])$root)
  d[i]<-as.numeric(uniroot(f = lr.y2, lower = extr.pt[4,1], upper = 50)$root)
  cutoff.r[i]<-log((cont.cdf(d[i],mu1)-cont.cdf(c[i],mu1))/(cont.cdf(d[i],mu0)-cont.cdf(c[i],mu0)))
}
id.r<-sort(abs(cutoff.r-y2), index.return = T)$ix[1]

arl.mrlcusum<-function(h,mu)
{
  t<-1000; 
  half.w<-h/(2*t+1)
  r.mat<-matrix(0,nrow = t, ncol = t)
  r.f<-rep(0,2*t-1)
  
  dist.cdf<-function(x) (1-epsilon)*pnorm(x, mean = mu, sd = sigma0)+epsilon*pnorm(x, mean = mu, sd = sigma1)
  
  for(i in seq(from = 1, to = 2*t-1))
  {
    y<-(i-t)*2*half.w+half.w
    lr.y<-function(x) lr(x)-y
    if(y<y1[id.l])
      x<-as.numeric(uniroot(f = lr.y, lower = -100, upper = a[id.l])$root)
    else if(y<y2[id.r] && y1[id.l]<=y)
      x<-as.numeric(uniroot(f = lr.y, lower = b[id.l]-0.1, upper = c[id.r]+0.1)$root)
    else
      x<-as.numeric(uniroot(f = lr.y, lower = d[id.r], upper = 100)$root)
    r.f[i]<-dist.cdf(x)
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
count=0
arls<-arl.mrlcusum(h=h,mu=0)
while((abs(arls-1000))>0.5){
  if(arls<1000){minh=h
  h=(minh+maxh)/2
  arls<-arl.mrlcusum(h=h,mu=0)}
  else{maxh=h
  h=(minh+maxh)/2
  arls<-arl.mrlcusum(h=h,mu=0)
  }
  count=count+1
  if(count>50) break
  if((count>47)&(abs(arls-1000)<20)) break
  if((count>40)&(abs(arls-1000)<5)) break
  if((count>30)&(abs(arls-1000)<3)) break
  if((count>20)&(abs(arls-1000)<1)) break
}
print(c(mu1,h,arls))
hvector<-c(hvector,h)
arlvector<-c(arlvector,arls)


testshift<-c(0.5,0.7,0.9,1.2,1.5,1.8,2.1,2.5,3,4,5,6,7,8,9,10)
arlvector.rlcusum<-matrix(data=NA,nrow=length(testshift),ncol=2)
arlvector.rlcusum[,1]<-testshift

for(i in 1:length(testshift)){
  arlvector.rlcusum[i,2]<-arl.mrlcusum(h=h,mu=testshift[i])
}
arlvector.rlcusum


#arl.mrlcusum(h=h,mu=0)
#########Rlcusum_l0norm####

ally_p<-seq(extr.pt[4,2]+0.01,extr.pt[3,2]-0.01,0.0001)
ally_n<-seq(extr.pt[2,2]+0.01,extr.pt[1,2]-0.01,0.0001)

findy_p<-numeric(length(ally_p))
findy_n<-numeric(length(ally_n))

findroot<-numeric(3)
#findmax_y<-function(x) lr(x)-extr.pt[3,2]
#max_y<-uniroot(findmax_y,lower=extr.pt[4,1],upper)
#min_y<-

for(i in 1:length(ally_p)){
  newlr<-function(x) lr(x)-ally_p[i]
  findroot[1]<-uniroot(newlr,lower=as.numeric(extr.pt[2,1]),upper=as.numeric(extr.pt[3,1]))$root
  findroot[2]<-uniroot(newlr,lower=as.numeric(extr.pt[3,1]),upper=as.numeric(extr.pt[4,1]))$root
  findroot[3]<-uniroot(newlr,lower=as.numeric(extr.pt[4,1]),upper=100)$root
  findy_p[i]<-abs(findroot[3]-findroot[1])
}

for(i in 1:length(ally_n)){
  newlr<-function(x) lr(x)-ally_n[i]
  findroot[1]<-uniroot(newlr,lower=-100,upper=as.numeric(extr.pt[1,1]))$root
  findroot[2]<-uniroot(newlr,lower=as.numeric(extr.pt[1,1]),upper=as.numeric(extr.pt[2,1]))$root
  findroot[3]<-uniroot(newlr,lower=as.numeric(extr.pt[2,1]),upper=as.numeric(extr.pt[3,1]))$root
  findy_n[i]<-abs(findroot[3]-findroot[1])
}

ycut_p<-ally_p[which.min(findy_p)]
ycut_n<-ally_n[which.min(findy_n)]
root_p<-numeric(3)
root_n<-numeric(3)
newlr_p<-function(x) lr(x)-ycut_p
newlr_n<-function(x) lr(x)-ycut_n
root_p[1]<-uniroot(newlr_p,lower=as.numeric(extr.pt[2,1]),upper=as.numeric(extr.pt[3,1]))$root
root_p[2]<-uniroot(newlr_p,lower=as.numeric(extr.pt[3,1]),upper=as.numeric(extr.pt[4,1]))$root
root_p[3]<-uniroot(newlr_p,lower=as.numeric(extr.pt[4,1]),upper=100)$root
root_n[1]<-uniroot(newlr_n,lower=-100,upper=as.numeric(extr.pt[1,1]))$root
root_n[2]<-uniroot(newlr_n,lower=as.numeric(extr.pt[1,1]),upper=as.numeric(extr.pt[2,1]))$root
root_n[3]<-uniroot(newlr_n,lower=as.numeric(extr.pt[2,1]),upper=as.numeric(extr.pt[3,1]))$root


arl.mrlcusum.plus.l0<-function(h,mu)
{
  t<-1000; 
  half.w<-h/(2*t+1)
  r.mat<-matrix(0,nrow = t, ncol = t)
  r.f<-rep(0,2*t-1)
  
  dist.cdf<-function(x) (1-epsilon)*pnorm(x, mean = mu, sd = sigma0)+epsilon*pnorm(x, mean = mu, sd = sigma1)
  
  for(i in seq(from = 1, to = 2*t-1))
  {
    y<-(i-t)*2*half.w+half.w
    lr.y<-function(x) lr(x)-y
    if(y<lr(root_n[1]))
      x<-as.numeric(uniroot(f = lr.y, lower = -100, upper = extr.pt[1,1])$root)
    else if(y<lr(root_p[1]) && lr(root_n[1])<=y)
      x<-as.numeric(uniroot(f = lr.y, lower = extr.pt[2,1], upper = extr.pt[3,1])$root)
    else
      x<-as.numeric(uniroot(f = lr.y, lower = extr.pt[4,1], upper = 100)$root)
    r.f[i]<-dist.cdf(x)
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
count=0
arls<-arl.mrlcusum.plus.l0(h=h,mu=0)
while((abs(arls-1000))>0.5){
  if(arls<1000){minh=h
  h=(minh+maxh)/2
  arls<-arl.mrlcusum.plus.l0(h=h,mu=0)}
  else{maxh=h
  h=(minh+maxh)/2
  arls<-arl.mrlcusum.plus.l0(h=h,mu=0)
  }
  count=count+1
  if(count>50) break
  if((count>47)&(abs(arls-1000)<20)) break
  if((count>40)&(abs(arls-1000)<5)) break
  if((count>30)&(abs(arls-1000)<3)) break
  if((count>20)&(abs(arls-1000)<1)) break
}
print(c(mu1,h,arls))
hvector<-c(hvector,h)
arlvector<-c(arlvector,arls)


testshift<-c(0.5,0.7,0.9,1.2,1.5,1.8,2.1,2.5,3,4,5,6,7,8,9,10)
arlvector.rlcusum.l0<-matrix(data=NA,nrow=length(testshift),ncol=2)
arlvector.rlcusum.l0[,1]<-testshift

for(i in 1:length(testshift)){
  arlvector.rlcusum.l0[i,2]<-arl.mrlcusum.plus.l0(h=h,mu=testshift[i])
}
arlvector.rlcusum.l0

#arl.mrlcusum.plus.l0(h=4.867,mu=0)

#########Rlcusum_l1norm####

ally_p<-seq(extr.pt[4,2]+0.01,extr.pt[3,2]-0.01,0.0001)
ally_n<-seq(extr.pt[2,2]+0.01,extr.pt[1,2]-0.01,0.0001)

findy_p<-numeric(length(ally_p))
findy_n<-numeric(length(ally_n))

findroot<-numeric(3)
#findmax_y<-function(x) lr(x)-extr.pt[3,2]
#max_y<-uniroot(findmax_y,lower=extr.pt[4,1],upper)
#min_y<-

for(i in 1:length(ally_p)){
  newlr<-function(x) lr(x)-ally_p[i]
  findroot[1]<-uniroot(newlr,lower=as.numeric(extr.pt[2,1]),upper=as.numeric(extr.pt[3,1]))$root
  findroot[2]<-uniroot(newlr,lower=as.numeric(extr.pt[3,1]),upper=as.numeric(extr.pt[4,1]))$root
  findroot[3]<-uniroot(newlr,lower=as.numeric(extr.pt[4,1]),upper=100)$root
  findy_p[i]<-abs(findroot[3]-2*findroot[2]+findroot[1])
}

for(i in 1:length(ally_n)){
  newlr<-function(x) lr(x)-ally_n[i]
  findroot[1]<-uniroot(newlr,lower=-100,upper=as.numeric(extr.pt[1,1]))$root
  findroot[2]<-uniroot(newlr,lower=as.numeric(extr.pt[1,1]),upper=as.numeric(extr.pt[2,1]))$root
  findroot[3]<-uniroot(newlr,lower=as.numeric(extr.pt[2,1]),upper=as.numeric(extr.pt[3,1]))$root
  findy_n[i]<-abs(findroot[3]-2*findroot[2]+findroot[1])
}

ycut_p<-ally_p[which.min(findy_p)]
ycut_n<-ally_n[which.min(findy_n)]
root_p<-numeric(3)
root_n<-numeric(3)
newlr_p<-function(x) lr(x)-ycut_p
newlr_n<-function(x) lr(x)-ycut_n
root_p[1]<-uniroot(newlr_p,lower=as.numeric(extr.pt[2,1]),upper=as.numeric(extr.pt[3,1]))$root
root_p[2]<-uniroot(newlr_p,lower=as.numeric(extr.pt[3,1]),upper=as.numeric(extr.pt[4,1]))$root
root_p[3]<-uniroot(newlr_p,lower=as.numeric(extr.pt[4,1]),upper=100)$root
root_n[1]<-uniroot(newlr_n,lower=-100,upper=as.numeric(extr.pt[1,1]))$root
root_n[2]<-uniroot(newlr_n,lower=as.numeric(extr.pt[1,1]),upper=as.numeric(extr.pt[2,1]))$root
root_n[3]<-uniroot(newlr_n,lower=as.numeric(extr.pt[2,1]),upper=as.numeric(extr.pt[3,1]))$root


arl.mrlcusum.plus.l1<-function(h,mu)
{
  t<-1000; 
  half.w<-h/(2*t+1)
  r.mat<-matrix(0,nrow = t, ncol = t)
  r.f<-rep(0,2*t-1)
  
  dist.cdf<-function(x) (1-epsilon)*pnorm(x, mean = mu, sd = sigma0)+epsilon*pnorm(x, mean = mu, sd = sigma1)
  
  for(i in seq(from = 1, to = 2*t-1))
  {
    y<-(i-t)*2*half.w+half.w
    lr.y<-function(x) lr(x)-y
    if(y<lr(root_n[1]))
      x<-as.numeric(uniroot(f = lr.y, lower = -100, upper =extr.pt[1,1])$root)
    else if(y<lr(root_p[1]) && lr(root_n[1])<=y)
      x<-as.numeric(uniroot(f = lr.y, lower = extr.pt[2,1], upper = extr.pt[3,1])$root)
    else if(y>=lr(root_p[1]))
      x<-as.numeric(uniroot(f = lr.y, lower = extr.pt[4,1], upper = 100)$root)
    r.f[i]<-dist.cdf(x)
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
count=0
arls<-arl.mrlcusum.plus.l1(h=h,mu=0)
while((abs(arls-1000))>0.5){
  if(arls<1000){minh=h
  h=(minh+maxh)/2
  arls<-arl.mrlcusum.plus.l1(h=h,mu=0)}
  else{maxh=h
  h=(minh+maxh)/2
  arls<-arl.mrlcusum.plus.l1(h=h,mu=0)
  }
  count=count+1
  if(count>50) break
  if((count>47)&(abs(arls-1000)<20)) break
  if((count>40)&(abs(arls-1000)<5)) break
  if((count>30)&(abs(arls-1000)<3)) break
  if((count>20)&(abs(arls-1000)<1)) break
}
print(c(mu1,h,arls))
hvector<-c(hvector,h)
arlvector<-c(arlvector,arls)


testshift<-c(0.5,0.7,0.9,1.2,1.5,1.8,2.1,2.5,3,4,5,6,7,8,9,10)
arlvector.rlcusum.l1<-matrix(data=NA,nrow=length(testshift),ncol=2)
arlvector.rlcusum.l1[,1]<-testshift

for(i in 1:length(testshift)){
  arlvector.rlcusum.l1[i,2]<-arl.mrlcusum.plus.l1(h=h,mu=testshift[i])
}
arlvector.rlcusum.l1
#########Rlcusum_distancebasedond#####


deltavalue<-extr.pt[3,2]-extr.pt[4,2]

lrp<-function(x){lr(x)+deltavalue
}

lrn<-function(x){lr(x)-deltavalue
}

arl.mrlcusum.fd<-function(h,mu)
{
  t<-1000; 
  half.w<-h/(2*t+1)
  r.mat<-matrix(0,nrow = t, ncol = t)
  r.f<-rep(0,2*t-1)
  
  dist.cdf<-function(x) (1-epsilon)*pnorm(x, mean = mu, sd = sigma0)+epsilon*pnorm(x, mean = mu, sd = sigma1)
  
  for(i in seq(from = 1, to = 2*t-1))
  {
    y<-(i-t)*2*half.w+half.w
    
    if(y<extr.pt[2,2]){lr.y<-function(x) lrn(x)-y
    x<-as.numeric(uniroot(f = lr.y, lower = -100, upper = extr.pt[1,1])$root)}
    else if(y<extr.pt[3,2] && extr.pt[2,2]<=y){lr.y<-function(x) lr(x)-y
    x<-as.numeric(uniroot(f = lr.y, lower = extr.pt[2,1], upper = extr.pt[3,1])$root)}
    else{lr.y<-function(x) lrp(x)-y
    x<-as.numeric(uniroot(f = lr.y, lower = extr.pt[4,1], upper = 100)$root)}
    r.f[i]<-dist.cdf(x)
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
count=0
arls<-arl.mrlcusum.fd(h=h,mu=0)
while((abs(arls-1000))>0.5){
  if(arls<1000){minh=h
  h=(minh+maxh)/2
  arls<-arl.mrlcusum.fd(h=h,mu=0)}
  else{maxh=h
  h=(minh+maxh)/2
  arls<-arl.mrlcusum.fd(h=h,mu=0)
  }
  count=count+1
  if(count>50) break
  if((count>47)&(abs(arls-1000)<20)) break
  if((count>40)&(abs(arls-1000)<5)) break
  if((count>30)&(abs(arls-1000)<3)) break
  if((count>20)&(abs(arls-1000)<1)) break
}
print(c(mu1,h,arls))
hvector<-c(hvector,h)
arlvector<-c(arlvector,arls)


testshift<-c(0.5,0.7,0.9,1.2,1.5,1.8,2.1,2.5,3,4,5,6,7,8,9,10)
testshift<-seq(0.5,10,0.5)
arlvector.rlcusum.fd<-matrix(data=NA,nrow=length(testshift),ncol=2)
arlvector.rlcusum.fd[,1]<-testshift

for(i in 1:length(testshift)){
  arlvector.rlcusum.fd[i,2]<-arl.mrlcusum.fd(h=h,mu=testshift[i])
}
arlvector.rlcusum.fd


###########sdrl############

fdlr<-function(x){
  if(x<extr.pt[1,1]){return(lrn(x))}
  else if(x<extr.pt[2,1]){return(extr.pt[2,2])}
  else if(x<extr.pt[3,1]){return(lr(x))}
  else if(x<extr.pt[4,1]){return(extr.pt[3,2])}
  else {return(lrp(x))}
}
xxxx<-seq(-10,10,0.01)
yyyy<-numeric(length(xxx))
for(i in 1:length(xxxx)){yyyy[i]<-fdlr(xxx[i])}

plot(xxxx,yyyy,type='l',ylim=c(-4,4))

N=100000
truesigma=sqrt(1.8)
shiftvector<-seq(0.5,7,0.5)
for (shift in shiftvector){
  rl<-numeric(N)
  for (i in 1:N){
    S=0
    num=0
    while(S<=h){
      num=num+1
      znow<-rnorm(n=1, mean =shift,sd=sample(c(sigma0, sigma1), size = 1, prob = c(1-epsilon, epsilon), replace = T))
      S=max(0,S+fdlr(znow))
    }
    rl[i]=num
  }
  print(c(shift,mean(rl),sd(rl)))
}
##########

#plot
library(ggplot2)

CUSUMTYPE<-c('CUSUM'=1,'RLCUSUM'=2,'RLCUSUM0'=3,'RLCUSUM1'=4,'RLCUSUMfd'=5)
SHAPES<-c('CUSUM'=20,'RLCUSUM'=16,'RLCUSUM0'=17,'RLCUSUM1'=18,'RLCUSUMfd'=15)
COLORS<-c('CUSUM'='black','RLCUSUM'='red','RLCUSUM0'='blue','RLCUSUM1'='green','RLCUSUMfd'='orange')
F1<-ggplot(data=NULL,aes(testshift,arlvector.cusum[,2]))+
  geom_line(aes(linetype='CUSUM',color='CUSUM'))+
  geom_line(aes(x=testshift,y=arlvector.rlcusum[,2],linetype='RLCUSUM',color='RLCUSUM'))+
  geom_line(aes(x=testshift,y=arlvector.rlcusum.l0[,2],linetype='RLCUSUM0',color='RLCUSUM0'))+
  geom_line(aes(x=testshift,y=arlvector.rlcusum.l1[,2],linetype='RLCUSUM1',color='RLCUSUM1'))+
  geom_line(aes(x=testshift,y=arlvector.rlcusum.fd[,2],linetype='RLCUSUMfd',color='RLCUSUMfd'))+
  geom_point(aes(x=Stestshift,y=Sarlvector.cusum,shape='CUSUM',color='CUSUM'))+
  geom_point(aes(x=Stestshift,y=Sarlvector.rlcusum,shape='RLCUSUM',color='RLCUSUM'))+
  geom_point(aes(x=Stestshift,y=Sarlvector.rlcusum.l0,shape='RLCUSUM0',color='RLCUSUM0'))+
  geom_point(aes(x=Stestshift,y=Sarlvector.rlcusum.l1,shape='RLCUSUM1',color='RLCUSUM1'))+
  geom_point(aes(x=Stestshift,y=Sarlvector.rlcusum.fd,shape='RLCUSUMfd',color='RLCUSUMfd'))+
  #geom_point(aes(x=c(textr.pt[,1]),y=c(textr.pt[,2])))+
  #geom_text(aes(x=c(-2.5,3.5),y=c(-1.7,1.7),label=c('A','B')))+
  #scale_colour_manual(name="class",values=cols)+
  scale_linetype_manual(name=NULL,values = CUSUMTYPE,labels=c('CUSUM','WCUSUM',expression(DCUSUM[L[0]]),expression(DCUSUM[L[1]]),expression(DCUSUM[L^minute])))+
  scale_shape_manual(name=NULL,values = SHAPES,labels=c('CUSUM','WCUSUM',expression(DCUSUM[L[0]]),expression(DCUSUM[L[1]]),expression(DCUSUM[L^minute])))+
  scale_colour_manual(name=NULL,values = COLORS,labels=c('CUSUM','WCUSUM',expression(DCUSUM[L[0]]),expression(DCUSUM[L[1]]),expression(DCUSUM[L^minute])))+
  xlab(expression(delta))+
  ylab("ARL")+
  theme_bw()
F1<-F1+coord_cartesian(xlim=c(1,10),ylim=c(0,20))
F1+theme(legend.position=c(0.7,0.7),legend.key.size = unit(30, "pt"))


CUSUMTYPE<-c('CUSUM'=1,'RLCUSUM'=2,'RLCUSUM0'=3,'RLCUSUM1'=4,'RLCUSUMfd'=5)
SHAPES<-c('CUSUM'=20,'RLCUSUM'=16,'RLCUSUM0'=17,'RLCUSUM1'=18,'RLCUSUMfd'=15)
COLORS<-c('CUSUM'='black','RLCUSUM'='red','RLCUSUM0'='blue','RLCUSUM1'='green','RLCUSUMfd'='orange')
F1<-ggplot(data=NULL,aes(testshift,arlvector.cusum[,2]))+
  geom_line(aes(linetype='CUSUM',color='CUSUM'))+
  geom_line(aes(x=testshift,y=arlvector.rlcusum[,2],linetype='RLCUSUM',color='RLCUSUM'))+
  geom_line(aes(x=testshift,y=arlvector.rlcusum.l0[,2],linetype='RLCUSUM0',color='RLCUSUM0'))+
  geom_line(aes(x=testshift,y=arlvector.rlcusum.l1[,2],linetype='RLCUSUM1',color='RLCUSUM1'))+
  geom_line(aes(x=testshift,y=arlvector.rlcusum.fd[,2],linetype='RLCUSUMfd',color='RLCUSUMfd'))+
  geom_point(aes(x=Stestshift,y=Sarlvector.cusum,shape='CUSUM',color='CUSUM'))+
  geom_point(aes(x=Stestshift,y=Sarlvector.rlcusum,shape='RLCUSUM',color='RLCUSUM'))+
  geom_point(aes(x=Stestshift,y=Sarlvector.rlcusum.l0,shape='RLCUSUM0',color='RLCUSUM0'))+
  geom_point(aes(x=Stestshift,y=Sarlvector.rlcusum.l1,shape='RLCUSUM1',color='RLCUSUM1'))+
  geom_point(aes(x=Stestshift,y=Sarlvector.rlcusum.fd,shape='RLCUSUMfd',color='RLCUSUMfd'))+
  #geom_point(aes(x=c(textr.pt[,1]),y=c(textr.pt[,2])))+
  #geom_text(aes(x=c(-2.5,3.5),y=c(-1.7,1.7),label=c('A','B')))+
  #scale_colour_manual(name="class",values=cols)+
  scale_linetype_manual(name=NULL,values = CUSUMTYPE,labels=c('CUSUM','WCUSUM',expression(DCUSUM[L[0]]),expression(DCUSUM[L[1]]),expression(DCUSUM[L^minute])))+
  scale_shape_manual(name=NULL,values = SHAPES,labels=c('CUSUM','WCUSUM',expression(DCUSUM[L[0]]),expression(DCUSUM[L[1]]),expression(DCUSUM[L^minute])))+
  scale_colour_manual(name=NULL,values = COLORS,labels=c('CUSUM','WCUSUM',expression(DCUSUM[L[0]]),expression(DCUSUM[L[1]]),expression(DCUSUM[L^minute])))+
  xlab(expression(delta))+
  ylab("ARL")+
  theme_bw()
F1<-F1+coord_cartesian(xlim=c(1,10),ylim=c(0,20))
#F2+theme(legend.position=c(0.7,0.7),legend.key.size = unit(30, "pt"))

F1+theme(legend.position="none")

##############New complex distribution##############
mu0<-0; mu1<-2;sigma0<-1;sigma1<-3;epsilon1<-0.3;epsilon2<-0.05
v<-3;tsigma<-sqrt((v-2)/v)

#tlr<-function(x){(v+1)/2*(log(tsigma^2*v+(x-tmu0)^2)-log(tsigma^2*v+(x-tmu1)^2))}
clr<-function(x)
{
  log((1-epsilon1-epsilon2)*dnorm(x, mean=mu1, sd=sigma0)+epsilon1*dnorm(x, mean=mu1, sd=sigma1)+epsilon2*dt(x=(x-mu1)/tsigma,df=v))-
    log((1-epsilon1-epsilon2)*dnorm(x, mean=mu0, sd=sigma0)+epsilon1*dnorm(x, mean=mu0, sd=sigma1)+epsilon2*dt(x=(x-mu0)/tsigma,df=v))
}

plot(clr,xlim=c(-20,20))
abline(h=0,v=0,lty=2)
#############
mu0<-0; mu1<-1;sigma0<-1;sigma1<-2;epsilon1<-0.4;epsilon2<-0.2
v<-3;tsigma<-sqrt((v-2)/v)


lr<-function(x)
{
  log((1-epsilon1-epsilon2)*dnorm(x, mean=mu1, sd=sigma0)+epsilon1*dnorm(x, mean=mu1, sd=sigma1)+epsilon2*dt(x=(x-mu1)/tsigma,df=v))-
    log((1-epsilon1-epsilon2)*dnorm(x, mean=mu0, sd=sigma0)+epsilon1*dnorm(x, mean=mu0, sd=sigma1)+epsilon2*dt(x=(x-mu0)/tsigma,df=v))
}

fn.surv.lr <- function(x)
{
  log((1-epsilon1-epsilon2)*pnorm(x, mean=mu1, sd=sigma0, lower.tail=FALSE)+epsilon1*pnorm(x, mean=mu1, sd=sigma1, lower.tail=FALSE)+epsilon2*pt(q=(x-mu1)/tsigma,df=v,lower.tail=FALSE))-
    log((1-epsilon1-epsilon2)*pnorm(x, mean=mu0, sd=sigma0, lower.tail=FALSE)+epsilon1*pnorm(x, mean=mu0, sd=sigma1, lower.tail=FALSE)+epsilon2*pt(q=(x-mu0)/tsigma,df=v,lower.tail=FALSE))
}
plot(lr,xlim=c(-20,20),type='l')
xxx<-seq(-20,20,0.01)
yyy<-fn.surv.lr(xxx)
lines(xxx,yyy,lty=2)
abline(h=0,v=0,lty=3)

N=100000
h=4.645

rl<-numeric(N)
count=0
for(i in 1:N){
  S=0
  count=0
  while(S<=h){
    count=count+1
    flag<-sample(c(1,2,3),size=1,pro=c(1-epsilon1-epsilon2,epsilon1,epsilon2),replace=TRUE)
    if(flag==1){randomconx<-rnorm(n=1,mean=mu0,sd=sigma0)}
    else if(flag==2){randomconx<-rnorm(n=1,mean=mu0,sd=sigma1)}
    else{randomconx<-tsigma*rt(n=1,df=v)+mu0}
    z<-lr(randomconx)
    if(S+z>0){S=S+z}
    else{S=0}
  }
  rl[i]<-count
}
mean(rl)

N=500000
testshift<-seq(0.5,10,by=0.1)
rlvector<-numeric(length((testshift)))
for(j in 1:length(testshift)){
  rl1<-numeric(N)
  count=0
  for(i in 1:N){
    S=0
    count=0
    deltamu<-testshift[j]
    while(S<=h){
      count=count+1
      flag<-sample(c(1,2,3),size=1,pro=c(1-epsilon1-epsilon2,epsilon1,epsilon2),replace=TRUE)
      if(flag==1){randomconx<-rnorm(n=1,mean=deltamu,sd=sigma0)}
      else if(flag==2){randomconx<-rnorm(n=1,mean=deltamu,sd=sigma1)}
      else{randomconx<-tsigma*rt(n=1,df=v)+deltamu}
      z<-lr(randomconx)
      if(S+z>0){S=S+z}
      else{S=0}
    }
    rl1[i]<-count
  }
  rlvector[j]<-mean(rl1)
}
###############

lr.neg<-function(x) -lr(x)

ex.pt2<-nlm(f = lr, p = -2)
ex.pt2$extr.val<-ex.pt2$minimum


ex.pt1<-nlm(f = lr.neg, p = -3)
ex.pt1$extr.val<--ex.pt1$minimum

ex.pt5<-nlm(f = lr, p = -10)
ex.pt5$extr.val<-ex.pt5$minimum

ex.pt2<-nlm(f = lr, p = -1)
ex.pt2$extr.val<-ex.pt2$minimum

ex.pt3<-nlm(f = lr.neg, p = 2)
ex.pt3$extr.val<--ex.pt3$minimum

ex.pt4<-nlm(f = lr, p =4)
ex.pt4$extr.val<-ex.pt4$minimum

ex.pt6<-nlm(f = lr.neg, p = 10)
ex.pt6$extr.val<--ex.pt6$minimum

ex.pt<-rbind(ex.pt1, ex.pt2, ex.pt3, ex.pt4,ex.pt5,ex.pt6)
ex.pt[, c(2,6)]

extr.pt<-matrix(as.numeric(ex.pt[,c(2,6)]), nrow = 6, ncol = 2)

deltavalue<-extr.pt[3,2]-extr.pt[4,2]

lrd<-function(x){
  if(x<extr.pt[5,1]){return(extr.pt[2,2]+extr.pt[5,2]-extr.pt[1,2])}
  else if(x<extr.pt[1,1]){return(lr(x)+extr.pt[2,2]-extr.pt[1,2])}
  else if(x<extr.pt[2,1]){return(extr.pt[2,2])}
  else if(x<extr.pt[3,1]){return(lr(x))}
  else if(x<extr.pt[4,1]){return(extr.pt[3,2])}
  else if(x<extr.pt[6,1]){return(lr(x)+extr.pt[3,2]-extr.pt[4,2])}
  else{extr.pt[6,2]+extr.pt[3,2]-extr.pt[4,2]}
}

xxx2<-seq(-20,20,by=0.01)
yyy2<-seq(-20,20,by=0.01)
for(i in 1:length(xxx2)){yyy2[i]<-lrd(xxx2[i])}
lines(xxx2,yyy2,lty=3)

N=100000
h=4.705

rl<-numeric(N)
count=0
for(i in 1:N){
  S=0
  count=0
  while(S<=h){
    count=count+1
    flag<-sample(c(1,2,3),size=1,pro=c(1-epsilon1-epsilon2,epsilon1,epsilon2),replace=TRUE)
    if(flag==1){randomconx<-rnorm(n=1,mean=mu0,sd=sigma0)}
    else if(flag==2){randomconx<-rnorm(n=1,mean=mu0,sd=sigma1)}
    else{randomconx<-tsigma*rt(n=1,df=v)+mu0}
    z<-lrd(randomconx)
    if(S+z>0){S=S+z}
    else{S=0}
  }
  rl[i]<-count
}
mean(rl)

N=500000
testshift<-seq(0.5,10,by=0.1)
rlvectord<-numeric(length((testshift)))
for(j in 1:length(testshift)){
  rl1<-numeric(N)
  count=0
  for(i in 1:N){
    S=0
    count=0
    deltamu<-testshift[j]
    while(S<=h){
      count=count+1
      flag<-sample(c(1,2,3),size=1,pro=c(1-epsilon1-epsilon2,epsilon1,epsilon2),replace=TRUE)
      if(flag==1){randomconx<-rnorm(n=1,mean=deltamu,sd=sigma0)}
      else if(flag==2){randomconx<-rnorm(n=1,mean=deltamu,sd=sigma1)}
      else{randomconx<-tsigma*rt(n=1,df=v)+deltamu}
      z<-lrd(randomconx)
      if(S+z>0){S=S+z}
      else{S=0}
    }
    rl1[i]<-count
  }
  rlvectord[j]<-mean(rl1)
}

Stestshift<-testshift[c(1,3,5,8,11,15,20,25,30,36,46,56,66,76,86,96)]
Srlvector<-rlvector[c(1,3,5,8,11,15,20,25,30,36,46,56,66,76,86,96)]
Srlvectord<-rlvectord[c(1,3,5,8,11,15,20,25,30,36,46,56,66,76,86,96)]


CUSUMTYPE<-c('CUSUM'=1,'RLCUSUM'=2)
SHAPES<-c('CUSUM'=15,'RLCUSUM'=16)
COLORS<-c('CUSUM'='black','RLCUSUM'='orange')
F1<-ggplot(data=NULL,aes(testshift,rlvector))+
  geom_line(aes(linetype='CUSUM',color='CUSUM'))+
  geom_line(aes(x=testshift,y=rlvectord,linetype='RLCUSUM',color='RLCUSUM'))+
  geom_point(aes(x=Stestshift,y=Srlvector,shape='CUSUM',color='CUSUM'))+
  geom_point(aes(x=Stestshift,y=Srlvectord,shape='RLCUSUM',color='RLCUSUM'))+
  #geom_point(aes(x=c(textr.pt[,1]),y=c(textr.pt[,2])))+
  #geom_text(aes(x=c(-2.5,3.5),y=c(-1.7,1.7),label=c('A','B')))+
  #scale_colour_manual(name="class",values=cols)+
  scale_linetype_manual(name=NULL,values = CUSUMTYPE,labels=c('CUSUM','DCUSUM'))+
  scale_shape_manual(name=NULL,values = SHAPES,labels=c('CUSUM','DCUSUM'))+
  scale_colour_manual(name=NULL,values = COLORS,labels=c('CUSUM','DCUSUM'))+
  xlab(expression(delta))+
  ylab("ARL")+
  theme_bw()
F1<-F1+coord_cartesian(xlim=c(1,10),ylim=c(0,20))
F1+theme(legend.position=c(0.7,0.8),legend.key.size = unit(50, "pt"))










