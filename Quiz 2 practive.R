trials=1000
n1=10
n2=10
mu1=1
mu2=1
sig=1
alpha=0.05
#assuming mu1=mu2.
threshold=qf(1-alpha,2,n1+n2-2,(n1*mu1^2-n2*mu1^2)/sig^2)

Q1=rep(0,trials)
Q2=rep(0,trials)
F_stat=rep(0,trials)

for (i in 1:trials) {
  
  samp1=rnorm(n1,mu1,sig)
  samp2=rnorm(n2,mu2,sig)
  
  Q1[i]=(n1*mean(samp1)^2-n2*mean(samp2)^2)/sig^2
  Q2[i]=((n1-1)*var(samp1)+(n2-1)*var(samp2))/sig^2
  
  F_stat[i]=(Q1[i]/2)/(Q2[i]/(n1+n2-2))
}

length(F_stat[F_stat>threshold])/trials



trials=1000
n1=10
n2=10
mu1=-2
mu2=0
sig=1
alpha=0.05
#assuming mu1=mu2.
threshold=qf(1-alpha,1,n1-1,n1*mu2)

Q1=rep(0,trials)
Q2=rep(0,trials)
F_stat=rep(0,trials)

for (i in 1:trials) {
  
  samp1=rnorm(n1,mu1,sig)
  
  Q1[i]=n1*mean(samp1)^2/sig^2
  Q2[i]=(n1-1)*var(samp1)/sig^2
  
  F_stat[i]=Q1[i]/(Q2[i]/(n1-1))
}

length(F_stat[F_stat>threshold])/trials