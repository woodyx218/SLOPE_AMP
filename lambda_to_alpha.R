library(SLOPE)
########## soft-thresholding operator
eta_soft=function(x,b){
  return((x-sign(x)*b)*(abs(x)>b))
}

########## alpha_min for Lasso as in 'The LASSO risk for gaussian matrices' (1.14)
amin= function(delta,length_precision=10000){
  A=seq(0,10,length.out = length_precision)
  return(A[which.min(abs((1+A^2)*pnorm(-A)-A*dnorm(A)-delta/2))])
}


################################## state evolution equation for SLOPE in (2.5)
F= function(tau,alpha,prior=x0,iteration=100,sigma=0){
  result=0
  for (i in 1:iteration){
  result=result+mean((prox_sorted_L1(prior+tau*rnorm(p),alpha*tau)-prior)^2)/iteration
  }
  return(sigma^2+result/delta)
}


################################## alpha to tau calibration for SLOPE in (2.5)
alpha_to_tau=function(alpha_seq,max_iter=100,prior=x0,sigma=0){
  tau=sqrt(sigma^2+mean(prior^2)/delta)
  record_tau=rep(0,max_iter) # initialize
  for (t in 1:max_iter){
    tau=sqrt(F(tau,alpha_seq,prior,sigma=sigma))
    record_tau[t]=tau #record each tau
  }
  #print(record_tau)
  return(mean(record_tau))
}

################################## alpha to lambda calibration for SLOPE in Lemma 2.2
alpha_to_lambda=function(alpha_seq,max_iter=100,prior=x0,second_iter=100,sigma=0){
  tau=sqrt(sigma^2+mean(prior^2)/delta)
  record_tau=rep(0,max_iter) # initialize
  
  for (t in 1:max_iter){
    tau=sqrt(F(tau,alpha_seq,prior,sigma=sigma))
    record_tau[t]=tau #record each tau
  }
  tau=mean(record_tau)

  E=0
  for (t in 1:second_iter){
    prox_solution=prox_sorted_L1(prior+tau*rnorm(p),alpha_seq*tau)
    E=E+length(unique(abs(prox_solution)))/p/delta/second_iter
  }
  lambda_seq=(1-E)*alpha_seq*tau
  
  return(lambda_seq)
}
################################## lambda to alpha calibration (implicitly) in Lemma 2.2
lambda_to_alpha=function(lambda_seq,tol=1e-2,prior=x0,sigma=0){
  # standardize lambda sequence
  lambda_seq=sort(lambda_seq,T)
  l=lambda_seq/lambda_seq[1]

  
  ### find interval that has larger and smaller at endpoints compared to lambda_seq[1] for bisection
  # starting guess
  alpha1=l/2
  alpha2=l*2
  
  lambda1=alpha_to_lambda(alpha1,prior=prior,sigma=sigma)
  lambda2=alpha_to_lambda(alpha2,prior=prior,sigma=sigma)

  while ((lambda1[1]<lambda_seq[1])*(lambda2[1]>lambda_seq[1])==0){
    if (lambda1[1]<lambda_seq[1]){
      alpha1=alpha1*2
      alpha2=alpha2*2
    } else{
      alpha1=alpha1/2
      alpha2=alpha2/2
    }
    lambda1=alpha_to_lambda(alpha1,prior=prior,sigma=sigma)
    lambda2=alpha_to_lambda(alpha2,prior=prior,sigma=sigma)
#    print(c(lambda1[1],lambda2[1]))
#    print(c(alpha1[1],alpha2[1]))
  }
  
  
  ### bisection to find the alpha_seq which is parallel to lambda_seq
  while ((alpha2[1]-alpha1[1])>tol){
    middle_alpha_seq=(alpha1+alpha2)/2
    middle_lambda=alpha_to_lambda(middle_alpha_seq,prior=prior,sigma=sigma)
    if (middle_lambda[1]>lambda_seq[1]){
      alpha2=middle_alpha_seq
    }else if (middle_lambda[1]<lambda_seq[1]){
      alpha1=middle_alpha_seq
    }else{
      break
    }
  }
  return(middle_alpha_seq)
}

### SLOPE proximal operator is equivalent to some separable soft-thresholding operator
### in the sense of https://arxiv.org/abs/1903.11582, Proposition 1

# Adaptive threshold
Ahat=function(x,piplusZ,A){
  A=sort(abs(A),decreasing = TRUE)
  piplusZ=sort(abs(piplusZ),decreasing = TRUE)
  Ahat=piplusZ-prox_sorted_L1(piplusZ,A)
  output=x*0
  for(i in 1:length(x)){
    index=which.min(abs(abs(x[i])-piplusZ))
    output[i]=Ahat[index]
  }
  return(output)
}

# Limiting scalar function
h=function(x,piplusZ,A){
  A=sort(abs(A),decreasing = TRUE)
  piplusZ=sort(abs(piplusZ),decreasing = TRUE)
  prox=prox_sorted_L1(piplusZ,A)
  output=x*0
  for(i in 1:length(x)){
    index=which.min(abs(abs(x[i])-piplusZ))
    output[i]=prox[index]*sign(x[i])
  }
  return(output)
}
