library(lme4)
library(lmerTest)

delta<-0.8
n<-30
n.1<-n/2
n.2<-n/2
alpha<-0.05

seps<-1

rho<-seq(0.01,0.99,0.01)

asim<-1000

pow.mean1<-vector("numeric",length(rho))
pow.mean2<-vector("numeric",length(rho))
pow.mean3<-vector("numeric",length(rho))
pow.mean4<-vector("numeric",length(rho))

for(i in 1:length(rho))
{progress_bar <- txtProgressBar(1, length(rho), style = 3)
pow.1<-vector("numeric",asim)
pow.2<-vector("numeric",asim)
pow.3<-vector("numeric",asim)
pow.4<-vector("numeric",asim)
for(k in 1:asim)
{
y<-rnorm(n,0,1)
y[1:n.1]<-y[1:n.1]+delta
x<-c(rep(1,n.1),rep(0,n.2))
pow.1[k]<-summary(lm(y~x))$coef[2,4]<alpha

y2<-rnorm(n*2,0,1)
y2[1:(n.1*2)]<-y2[1:(n.1*2)]+delta
x2<-c(rep(1,(n.1*2)),rep(0,(n.2*2)))
pow.2[k]<-summary(lm(y2~x2))$coef[2,4]<alpha

alpac<-sqrt(rho[i]^2/(1-rho[i]^2)*seps)
y3<-rnorm(n,0,1)
y3.2u<-alpac*y3+rnorm(n)
y3.2<-y3.2u/sqrt(var(y3.2u))
y3[1:n.1]<-y3[1:n.1]+delta
y3.2[1:n.1]<-y3.2[1:n.1]+delta
y3o<-c(y3,y3.2)
x3<-c(rep(1,n.1),rep(0,n.2),rep(1,n.1),rep(0,n.2))
subject<-rep(1:n,2)
mm<-lmer(y3o ~ x3 + (1 | subject))    
pow.3[k]<-summary(mm)$coef[2,5]<alpha

y3m<-(y3+y3.2)/2
pow.4[k]<-summary(lm(y3m~x))$coef[2,4]<alpha

}
setTxtProgressBar(progress_bar,i)
pow.mean1[i]<-mean(pow.1)
pow.mean2[i]<-mean(pow.2)
pow.mean3[i]<-mean(pow.3)
pow.mean4[i]<-mean(pow.4)
}
close(progress_bar)


