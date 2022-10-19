**Question 2b**

For this question, we want to verify the answer we got for 2a. We expect
that for any
*α*
we put in, the proportion of trials that will have a high enough
*F*
statistic to reject the null hypothesis that
*μ*<sub>1</sub> = *μ*<sub>2</sub>
.

    #Set seed
    set.seed(743358478)

    #Set parameters
    N=10000
    n1=10
    n2=10
    mu1=1
    mu2=1
    sig=1
    alpha=0.05

    #Matrix such that the top left is a n1xn1 submatrix of 1/n1 and the bottom left is n2xn2 submatrix of 1/n2.
    M=c(rep(1,n1),rep(0,n2))%*%t(c(rep(1,n1),rep(0,n2)))/n1+c(rep(0,n1),rep(1,n2))%*%t(c(rep(0,n1),rep(1,n2)))/n2

    #The A_i matrices in order to satisfy the Fisher-Cochran Theorem.
    A_1=matrix(1/(n1+n2),ncol = n1+n2,nrow = n1+n2)
    A_2=M-A_1
    A_3=diag(n1+n2)-M

    #Assume mu1 and mu2 are equal.
    MU=c(rep(mu1,n1+n2))

    #The threshold to reject the null hypothesis.
    threshold=qf(1-alpha,1,n1+n2-2,t(MU)%*%(A_2)%*%MU/sig^2)

    #Store each statistic into its own vector.
    Q1=rep(0,N)
    Q2=rep(0,N)
    Q3=rep(0,N)
    F_stat=rep(0,N)

    #Iterate process N times.
    for (i in 1:N) {
      
      #Create samples and stack the samples into one vector.
      samp1=rnorm(n1,mu1,sig)
      samp2=rnorm(n2,mu2,sig)
      Y=c(samp1,samp2)
      
      #Compute the values for the Q_i values.
      Q1[i]=t(Y) %*% A_1 %*% Y/sig^2
      Q2[i]=t(Y) %*% A_2 %*% Y/sig^2
      Q3[i]=t(Y) %*% A_3 %*% Y/sig^2

      #Store the F statistic
      F_stat[i]=Q2[i]/(Q3[i]/(n1+n2-2))
    }

    #This tells us which proportion of the amount of trials would have Type I error. We expect this to be around alpha.
    length(F_stat[F_stat>threshold])/N

    ## [1] 0.0501

This comes out to be **0.044**. This is very close to
*α* = 0.05
which is exactly what we expect.

**Question 2c**

Just to show my computation hold, I have made a script here instead of
writing it by hand.

    #Setting parameters
    n1=10
    n2=10
    mu1=1
    mu2 = seq(1.1, 2, by = 0.1)
    sig=1
    alpha=0.05

    #Store the calculated power for each mu2.
    power_prediction=rep(0,length(mu2))

    #Matrix such that the top left is a n1xn1 submatrix of 1/n1 and the bottom left is n2xn2 submatrix of 1/n2.
    M = c(rep(1, n1), rep(0, n2)) %*% t(c(rep(1, n1), rep(0, n2))) / n1 + c(rep(0, n1), rep(1, n2)) %*%
      t(c(rep(0, n1), rep(1, n2))) / n2

    #The A_i matrices in order to satisfy the Fisher-Cochran Theorem.
    A_1 = matrix(1 / (n1+n2), ncol = n1 + n2, nrow = n1 + n2)
    A_2 = M - A_1
    A_3 = diag(n1 + n2) - M

    #Assume mu1 and mu2 are equal.
    MU_0 = c(rep(mu1, n1 + n2))

    #The threshold to reject the null hypothesis if mu1=mu2.
    threshold = qf(1 - alpha, 1, n1 + n2 - 2, t(MU_0) %*% (A_2) %*% MU_0 / sig ^
                     2)

    #Testing for various mu2 values
    for (i in 1:length(mu2)) {
      
      #Recalculate mu but with mu1 and mu2 different.
      MU=c(rep(mu1, n1),rep(mu2[i],n2))
      
      #The computation for 2c.
      power_prediction[i]=1-pf(threshold,1,n1+n2-2, t(MU) %*% (A_2) %*% MU / sig ^
                     2)
    }

    data.frame(mu=mu2,power=power_prediction)

    ##     mu      power
    ## 1  1.1 0.05516129
    ## 2  1.2 0.07082135
    ## 3  1.3 0.09742459
    ## 4  1.4 0.13545156
    ## 5  1.5 0.18509566
    ## 6  1.6 0.24593206
    ## 7  1.7 0.31667051
    ## 8  1.8 0.39506921
    ## 9  1.9 0.47805041
    ## 10 2.0 0.56200665

Now we can compare to part 2d to see how close this is.

mu power 1.1 0.06628236  
1.2 0.08763217  
1.3 0.11500345  
1.4 0.14923924  
1.5 0.19093427  
1.6 0.24029097  
1.7 0.29699530  
1.8 0.36014233  
1.9 0.42823673  
2.0 0.49928175

**Question 2d**

For this question we want to find the power of this test for various
values of
*μ*<sub>2</sub>
. We should expect that as
*μ*<sub>2</sub>
gets further away from
*μ*<sub>1</sub>
, the power increases.

    #Set seed
    set.seed(1269075846)

    #Set parameters
    N = 10000
    n1 = 10
    n2 = 10
    mu1 = 1
    #Testing for multiple mu2 values.
    mu2 = seq(1.1, 2, by = 0.1)
    sig = 1
    alpha = 0.05

    #Matrix such that the top left is a n1xn1 submatrix of 1/n1 and the bottom left is n2xn2 submatrix of 1/n2.
    M = c(rep(1, n1), rep(0, n2)) %*% t(c(rep(1, n1), rep(0, n2))) / n1 + c(rep(0, n1), rep(1, n2)) %*%
      t(c(rep(0, n1), rep(1, n2))) / n2

    #The A_i matrices in order to satisfy the Fisher-Cochran Theorem.
    A_1 = matrix(1 / (n1+n2), ncol = n1 + n2, nrow = n1 + n2)
    A_2 = M - A_1
    A_3 = diag(n1 + n2) - M

    #Assume mu1 and mu2 are equal.
    MU = c(rep(mu1, n1 + n2))

    #The threshold to reject the null hypothesis.
    threshold = qf(1 - alpha, 1, n1 + n2 - 2, t(MU) %*% (A_2) %*% MU / sig ^
                     2)

    #Store the calculated power for each mu2.
    power=rep(0,length(mu2))

    #Store each statistic into its own vector.
    Q1 = rep(0, N)
    Q2 = rep(0, N)
    Q3 = rep(0, N)
    F_stat = rep(0, N)

    #Testing each mu2 at once.
    for (j in 1:length(mu2)) {
      for (i in 1:N) {
      
      #Create samples and stack the samples into one vector.
        samp1 = rnorm(n1, mu1, sig)
        samp2 = rnorm(n2, mu2[j], sig)
        Y = c(samp1, samp2)
        
        #Compute the values for the Q_i values.
        Q1[i] = t(Y) %*% A_1 %*% Y / sig ^ 2
        Q2[i] = t(Y) %*% A_2 %*% Y / sig ^ 2
        Q3[i] = t(Y) %*% A_3 %*% Y / sig ^ 2
        
      #Store the F statistic
        F_stat[i] = Q2[i] / (Q3[i] / (n1 + n2 - 2))
      }
      #Estimate the power after all N trials for each mu2.
    power[j]=length(F_stat[F_stat > threshold]) / N
    }

    data.frame(mu=mu2,power,estimation=power_prediction)

    ##     mu  power estimation
    ## 1  1.1 0.0528 0.05516129
    ## 2  1.2 0.0668 0.07082135
    ## 3  1.3 0.0977 0.09742459
    ## 4  1.4 0.1357 0.13545156
    ## 5  1.5 0.1800 0.18509566
    ## 6  1.6 0.2456 0.24593206
    ## 7  1.7 0.3214 0.31667051
    ## 8  1.8 0.3985 0.39506921
    ## 9  1.9 0.4758 0.47805041
    ## 10 2.0 0.5679 0.56200665

This is exactly in line with my prediction. It is very close to the
estimation from 2c.

mu power estimation 1.1 0.0599 0.06628236  
1.2 0.0830 0.08763217  
1.3 0.1108 0.11500345  
1.4 0.1499 0.14923924  
1.5 0.2016 0.19093427  
1.6 0.2443 0.24029097  
1.7 0.3152 0.29699530  
1.8 0.3799 0.36014233  
1.9 0.4450 0.42823673  
2.0 0.5266 0.49928175
