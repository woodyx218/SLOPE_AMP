library(SLOPE)
##########
eta_soft=function(x,b){
  return((x-sign(x)*b)*(abs(x)>b))
}

########## alpha range
amin= function(delta,length_precision=10000){
  A=seq(0,10,length.out = length_precision)
  return(A[which.min(abs((1+A^2)*pnorm(-A)-A*dnorm(A)-delta/2))])
}


################################## state evolution equation
F= function(tau,alpha,prior=x0,delta=delta,iteration=100,sigma=0){
  result=0
  for (i in 1:iteration){
  result=result+mean((prox_sorted_L1(prior+tau*rnorm(length(alpha)),alpha*tau)-prior)^2)/iteration
  }
  return(sigma^2+result/delta)
}


################################## alpha to tau calibration
alpha_to_tau=function(alpha_seq,prior=x0,delta=delta,max_iter=100,sigma=0,verbose=0){
  tau=sqrt(sigma^2+mean(prior^2)/delta)
  record_tau=rep(0,max_iter) # initialize
  
  for (t in 1:max_iter){
    tau=sqrt(F(tau,alpha_seq,prior,delta,sigma=sigma))
    record_tau[t]=tau #record each tau
  }
  if(verbose){print(record_tau)}
  return(mean(tail(record_tau,50)))
}

################################## alpha to lambda calibration
alpha_to_lambda=function(alpha_seq,prior=x0,delta=delta,max_iter=100,second_iter=100,sigma=0){
  tau=sqrt(sigma^2+mean(prior^2)/delta)
  record_tau=rep(0,max_iter) # initialize
  p=length(prior)
  
  for (t in 1:max_iter){
    tau=sqrt(F(tau,alpha_seq,prior,delta=delta,sigma=sigma))
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
################################## lambda to alpha calibration (implicitly)
lambda_to_alpha=function(lambda_seq,prior=x0,delta=delta,tol=1e-3,sigma=0){
  # standardize lambda sequence
  lambda_seq=sort(lambda_seq,T)
  l=lambda_seq/lambda_seq[1]

  
  ### find interval that has larger and smaller at endpoints compared to lambda_seq[1] for bisection
  # starting guess
  alpha1=l/2
  alpha2=l*2
  
  lambda1=alpha_to_lambda(alpha1,prior=prior,delta=delta,sigma=sigma)
  lambda2=alpha_to_lambda(alpha2,prior=prior,delta=delta,sigma=sigma)

  while ((lambda1[1]<lambda_seq[1])*(lambda2[1]>lambda_seq[1])==0){
    if (lambda1[1]<lambda_seq[1]){
      alpha1=alpha1*2
      alpha2=alpha2*2
    } else{
      alpha1=alpha1/2
      alpha2=alpha2/2
    }
    lambda1=alpha_to_lambda(alpha1,prior=prior,delta=delta,sigma=sigma)
    lambda2=alpha_to_lambda(alpha2,prior=prior,delta=delta,sigma=sigma)
#    print(c(lambda1[1],lambda2[1]))
#    print(c(alpha1[1],alpha2[1]))
  }
  
  
  ### bisection to find the alpha_seq which is parallel to lambda_seq
  while ((alpha2[1]-alpha1[1])>tol){
    middle_alpha_seq=(alpha1+alpha2)/2
    middle_lambda=alpha_to_lambda(middle_alpha_seq,prior=prior,delta=delta,sigma=sigma)
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

averaging=function(C,iter=100000){
  # aim C decreasing
  l=length(C)
  iii=1
  while(prod((C[2:l]-C[1:(l-1)])<=0)==0){
    # find the first non-decreasing sub-sequence
    start=0
    whether_start=0
    final=0
    for(i in 2:l){
      # start sub-sequence
      if((C[i]>C[i-1])&(whether_start==0)){
        start=i-1
        whether_start=1
      }
      # end sub-sequence
      if((C[i]<C[i-1])&(whether_start==1)){
        final=i-1
        break
      }
      if(i==l){final=i}
    }
    # average the sub-sequence
    C[start:final]=mean(C[start:final])
    
    ### allow fixed iteration break
    iii=iii+1
    if(iii>iter){break}
  }
  return(C)
}

## essential penalty of infinite dimension
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


# essential penalty vector of finite dimension
A_f=function(pseudo,a,A_seq,z_seq){
  output=pseudo*0
  for(i in 1:length(pseudo)){
    if(abs(pseudo[i])<a){output[i]=a}
    else{
      index=which.min(abs(abs(pseudo[i])-z_seq))
      output[i]=A_seq[index]
    }
  }
  return(output)
}

## limiting scalar function
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