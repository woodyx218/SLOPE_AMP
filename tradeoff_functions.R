library(quadprog)

# This code simulates the FDP-TPP curve, based heavily on https://github.com/wjsu/fdrlasso
lsandwich=function(t, delta, epsi){
  Lnume = (1-epsi)*(2*(1+t^2)*pnorm(-t) - 2*t*dnorm(t)) + epsi*(1+t^2) - delta;
  Ldeno = epsi*((1+t^2)*(1-2*pnorm(-t)) + 2*t*dnorm(t));
  L = Lnume/Ldeno;
  return(L)
}

rsandwich=function(t, u){
  return((1-u)/(1-2*pnorm(-t)))
}

# from (C.5), need careful thinking
plus_f = function(x,delta){return((delta+2*x*dnorm(x) - 2*(1+x^2)*pnorm(-x))/(1+x^2-2*(1+x^2)*pnorm(-x)+2*x*dnorm(x)))}

epsilonDT=function(delta,dx=0.01){
  if(delta==1||delta==0){return(delta)}
  epsi = plus_f(seq(0,8,dx),delta)
  return(max(epsi))
}

r_from_u=function(u,delta,epsilon,dx=0.001){
  tstar=thresfromu_upper(u_star(delta,epsilon,dx=dx),delta,epsilon,dx=dx)
  estar=epsilonDT(delta)
  return(qnorm((epsilon*(2-u)-estar)/2/(epsilon-estar))/tstar)
}

w_from_u=function(u,delta,epsilon,dx=0.001){
  tstar=thresfromu_upper(u_star(delta,epsilon,dx=dx),delta,epsilon,dx=dx)
  estar=epsilonDT(delta)
  r=r_from_u(u,delta,epsilon)
  return(estar+2*(1-estar)/(1-r)*(pnorm(-tstar)-r*pnorm(-r*tstar)-(dnorm(-tstar)-dnorm(-r*tstar))/tstar))
}

# Supplement to FDR Page 17 last row
u_star=function(delta, epsilon,dx=0.001){
  if (delta >= 1){return(1)};
  epsilon_star = epsilonDT(delta,dx);
  if (epsilon <= epsilon_star){return(1)}
  return((epsilon - epsilon_star)*(delta - epsilon_star)/epsilon/(1 - epsilon_star) + epsilon_star/epsilon)
}

# main function: q^* and q_*
# Lasso part of q^*
qfromu=function(u,delta,epsi,dx=0.01,verbose=0){
  if (u>u_star(delta,epsi,dx)){return(1)}
  if (u==0){return(0)}
  stepsize=dx # smaller better but slower
  tmax=50
  tmin=tmax-stepsize
  if (lsandwich(tmin,delta,epsi)<rsandwich(tmin,u)){print('tmax not large enough')}
  
  while(tmin>0){
    if (lsandwich(tmin,delta,epsi)<rsandwich(tmin,u)){break}
    tmax=tmin
    tmin=tmax-stepsize
  }
  
  diff = stepsize;
  while (diff > 1e-7){
    tmid = 0.5*tmax + 0.5*tmin;
    if (lsandwich(tmid, delta, epsi) > rsandwich(tmid, u)){tmax = tmid}
    else{tmin = tmid}
    diff = tmax - tmin
  }
  
  t = (tmax + tmin)/2;
  if(verbose==1){print(c('zero threshold: ',t))}
  return(2*(1-epsi)*pnorm(-t)/(2*(1-epsi)*pnorm(-t) + epsi*u));
}

# entire q^*
qfromu_upper=function(u_seq,delta,epsi,dx=0.01){
  output=u_seq*0
  ustar=u_star(delta,epsi,dx)
  estar=epsilonDT(delta,dx)
  for (i in 1:length(u_seq)){
    if(u_seq[i]<ustar){
      output[i]=qfromu(u_seq[i],delta,epsi,dx)
    }
    else{
      output[i]=(epsi*u_seq[i]*(1-epsi)-estar*(1-epsi))/(epsi*u_seq[i]*(1-estar)-estar*(1-epsi))
    }
  }
  return(output)
}

# this is fast implementation which requires decreasing input, also coarse; for finer result use thresfromu_lower
qfromu_lower=function(u_seq,delta,epsi,dx=0.01,dz=0.1,dt=0.1,tmax=12,zmax=15,verbose=0,trape=FALSE){
  if(any(diff(u_seq)>0)){return('Input sequence not decreasing')}
  output=u_seq*0
  c_start=dx
  t_seq=seq(0,tmax,dt)
  
  for(ind in 1:length(u_seq)){ # 16 seconds per u
    u=u_seq[ind]
    sta=Sys.time()
    c_seq=seq(c_start,10,dx)
    for(c in c_seq){
      z_seq=seq(c,zmax,dz)
      height=matrix(0,nrow=length(t_seq),ncol=length(t_seq))
      stop_sign=0
      
      for(i in 1:length(t_seq)){
        for(j in i:length(t_seq)){
          t1=t_seq[i];t2=t_seq[j]
          p1=(u-pnorm(t2-c)-pnorm(-t2-c))/((pnorm(t1-c)+pnorm(-t1-c))-pnorm(t2-c)-pnorm(-t2-c))
          if(p1<0 || p1>1){next}
          
          A=NULL
          try({A=optimAhat(epsi,c,u,t1,t2,dz=dz,zmax=max(z_seq))})
          
          if(is.null(A)==FALSE){
            height[i,j]=2*(1-epsi)*sum((z_seq-A)^2*dnorm(z_seq))*dz+epsi*p1*(t1^2*(pnorm(c-t1)-pnorm(-c-t1))+sum((z_seq-t1-A)^2*dnorm(z_seq-t1))*dz+sum((-z_seq-t1+A)^2*dnorm(-z_seq-t1))*dz)+epsi*(1-p1)*(t2^2*(pnorm(c-t2)-pnorm(-c-t2))+sum((z_seq-t2-A)^2*dnorm(z_seq-t2))*dz+sum((-z_seq-t2+A)^2*dnorm(-z_seq-t2))*dz)
            if(trape==TRUE){height[i,j]=height[i,j]-0.5*epsi*(p1*t1^2*(dnorm(c-t1)+dnorm(-c-t1))+(1-p1)*t2^2*(dnorm(c-t2)+dnorm(-c-t2)))*dz}
            if(height[i,j]<delta & height[i,j]>0){stop_sign=1;break}
            }
        }
        if(stop_sign){break}
        }

    if(any(height)==TRUE && min(height[height>0])>delta){
      if(verbose==1){
        print(Sys.time()-sta)
        print(paste0('TPP: ',u,'; FDP: ',2*(1-epsi)*pnorm(-c)/(2*(1-epsi)*pnorm(-c)+epsi*u),'; ',qfromu_upper(u,delta,epsi),' ; thresh: ',c))
      }
      output[ind]=2*(1-epsi)*pnorm(-c)/(2*(1-epsi)*pnorm(-c)+epsi*u)
      c_start=c
      break
    }
    }
  }
  return(output)
}

optimAhat=function(epsilon,zero_thresh,u,t1,t2,dz=0.1,zmax=8){
  z_seq=seq(zero_thresh,zmax,dz)
  if(t1!=t2){
  p1=(u-pnorm(t2-zero_thresh)-pnorm(-t2-zero_thresh))/((pnorm(t1-zero_thresh)+pnorm(-t1-zero_thresh))-pnorm(t2-zero_thresh)-pnorm(-t2-zero_thresh))
  p2=1-p1}
  else{p1=1;p2=0}
  
  Q=2*(1-epsilon)*diag(dnorm(z_seq))*dz+epsilon*diag(p1*(dnorm(z_seq-t1)+dnorm(-z_seq-t1))+p2*(dnorm(z_seq-t2)+dnorm(-z_seq-t2)))*dz
  d=2*(1-epsilon)*z_seq*dnorm(z_seq)*dz+epsilon*dz*(p1*((z_seq-t1)*dnorm(z_seq-t1)+(z_seq+t1)*dnorm(-z_seq-t1))+p2*((z_seq-t2)*dnorm(z_seq-t2)+(z_seq+t2)*dnorm(-z_seq-t2)))
  # compute Q
  #temp=p1*(dnorm(z_seq-t1)+dnorm(-z_seq-t1))+(1-p1)*(dnorm(z_seq-t2)+dnorm(-z_seq-t2))
  #Q=2*(1-epsilon)*diag(dnorm(z_seq))*dz+epsilon*diag(temp+1e-8)*dz
  
  Constraint=diag(length(z_seq))
  for(row in 2:nrow(Constraint)){Constraint[row,row-1]=-1}
  bvec=c(zero_thresh,rep(0,length(z_seq)-1))
  QP=solve.QP(Q,d,t(Constraint),bvec)
  return(QP$solution)
}


thresfromu_upper=function(u,delta,epsi,dx=0.0001,tol=1e-7){
  u_st=u_star(delta,epsi,dx)
  e_st=epsilonDT(delta,dx)
  if (u>u_st){return(-qnorm((1-(1-u)*epsi*(1-e_st)/(epsi-e_st)-epsi*u)/2/(1-epsi)))}
  if (u==0){return(Inf)}
  stepsize=dx # smaller better but slower
  tmax=20
  tmin=tmax-stepsize
  if (lsandwich(tmin,delta,epsi)<rsandwich(tmin,u)){print('tmax not large enough')}
  
  while(tmin>0){
    if (lsandwich(tmin,delta,epsi)<rsandwich(tmin,u)){break}
    tmax=tmin
    tmin=tmax-stepsize
  }
  diff = stepsize;
  while (diff > tol){
    tmid = 0.5*tmax + 0.5*tmin;
    if (lsandwich(tmid, delta, epsi) > rsandwich(tmid, u)){tmax = tmid}
    else{tmin = tmid}
    diff = tmax - tmin
  }
  
  return((tmax + tmin)/2);
}

ddnorm=function(x){return(-x*dnorm(x))}
K_f=function(z,t,epsi){return(4*(1-epsi)*dnorm(z)+2*epsi*(dnorm(z-t)+dnorm(z+t)))}
Kp_f=function(z,t,epsi){return(4*(1-epsi)*ddnorm(z)+2*epsi*(ddnorm(z-t)+ddnorm(z+t)))}
Kint_f=function(z,t,epsi){return(4*(1-epsi)*pnorm(z)+2*epsi*(pnorm(z-t)+pnorm(z+t)))}
ints=function(z){return(-z*dnorm(z)+pnorm(z))}
A_S=function(z,a,t,epsi){return(pmax(a,-Kp_f(z,t,epsi)/K_f(z,t,epsi)))}


## single-point prior given tpp and zero-threshold
tfromu=function(u,a){
  t_seq=seq(1e-6,10,1e-6)
  u_seq=pnorm(t_seq-a)+pnorm(-t_seq-a)
  return(t_seq[which.min(abs(u_seq-u))])
}

thresfromu_lower=function(u_seq,delta,epsi,dx=0.01,dz=0.1,dt=0.1,tmax=12,zmax=15,verbose=0,trape=TRUE,tol=1e-7){
  output=u_seq*0
  t_seq=seq(0,tmax,dt)
  
  for(ind in 1:length(u_seq)){ # 16 seconds per u
    u=u_seq[ind]
    sta=Sys.time()
    c_min=thresfromu_upper(u,delta,epsi,tol=tol)
    c_max=2*c_min
    while(c_max-c_min>tol){
      c_mid=c_max*0.5+c_min*0.5
      z_seq=seq(c_mid,zmax,dz)
      height=matrix(0,nrow=length(t_seq),ncol=length(t_seq))
      stop_sign=0
      
      for(i in 1:length(t_seq)){
        for(j in i:length(t_seq)){
          t1=t_seq[i];t2=t_seq[j]
          p1=(u-pnorm(t2-c_mid)-pnorm(-t2-c_mid))/((pnorm(t1-c_mid)+pnorm(-t1-c_mid))-pnorm(t2-c_mid)-pnorm(-t2-c_mid))
          if(p1<0 || p1>1){next}
          
          A=NULL
          try({A=optimAhat(epsi,c_mid,u,t1,t2,dz=dz,zmax=max(z_seq))})
          
          if(is.null(A)==FALSE){
            height[i,j]=2*(1-epsi)*sum((z_seq-A)^2*dnorm(z_seq))*dz+epsi*p1*(t1^2*(pnorm(c_mid-t1)-pnorm(-c_mid-t1))+sum((z_seq-t1-A)^2*dnorm(z_seq-t1))*dz+sum((-z_seq-t1+A)^2*dnorm(-z_seq-t1))*dz)+epsi*(1-p1)*(t2^2*(pnorm(c_mid-t2)-pnorm(-c_mid-t2))+sum((z_seq-t2-A)^2*dnorm(z_seq-t2))*dz+sum((-z_seq-t2+A)^2*dnorm(-z_seq-t2))*dz)
            if(trape==TRUE){height[i,j]=height[i,j]-0.5*epsi*(p1*t1^2*(dnorm(c_mid-t1)+dnorm(-c_mid-t1))+(1-p1)*t2^2*(dnorm(c_mid-t2)+dnorm(-c_mid-t2)))*dz}
            if(height[i,j]<delta & height[i,j]>0){stop_sign=1;c_min=c_mid;break}
          }
        }
        if(stop_sign){break}
      }
      
      if(any(height)==TRUE && min(height[height>0])>delta){
        c_max=c_mid
      }
      if(verbose==1){print(c_mid)}
    }
    output[ind]=c_mid
  }
  return(output)
}
