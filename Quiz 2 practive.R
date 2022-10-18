trials=10000
n1=10
n2=10
mu1=1
mu2=1
sig=1
alpha=0.05

threshold=qf(1-alpha,1,n1+n2-2)

Q1=rep(0,trials)
Q2=rep(0,trials)
F_stat=rep(0,trials)

for (i in 1:trials) {
  
  samp1=rnorm(n1,mu1,sig)
  samp2=rnorm(n2,mu2,sig)
  
  X_bar=(sum(samp1)+sum(samp2))/(n1+n2)
  X_bar1=mean(samp1)
  X_bar2=mean(samp2)
  
  Q1[i]=n1*(X_bar1-X_bar)^2+n2*(X_bar2-X_bar)^2
  Q2[i]=sum((samp1-X_bar1)^2)+sum((samp2-X_bar2)^2)
  
  F_stat[i]=Q1[i]/(Q2[i]/(n1+n2-2))
}

length(F_stat[F_stat>threshold])/trials

