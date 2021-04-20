library(SLOPE)
library(latex2exp)
source('lambda_to_alpha.R')
source('tradeoff_functions.R')

########################################################
########### upper bound and lower bound for delta=0.3, epsilon=0.5
########################################################
d=0.3;e=0.5
par(mar=c(5,5,1,1))
x2=seq(u_star(d,e),0.999,0.001)
plot(x2,qfromu_upper(x2,d,e,dx=0.001),col=4,lwd=3,type='l',ylab='FDP',xlab='TPP',
     xlim=c(0,1),ylim=c(0,0.55),xaxs = "i", yaxs = "i",cex.lab=2,cex.axis=1.5)

x1=c(seq(0.002,u_star(d,e),0.002),u_star(d,e))
curve=qfromu_upper(x1,d,e,dx=0.001)
lines(x1,curve,lwd=3,col='chartreuse3')

#x3=c(seq(0.999,0.5,-0.05),seq(0.5,0.1,-0.02),seq(0.09,0.01,-0.01));
#thres=thresfromu_lower(x3,d,e,tol=1e-5,verbose=1)
#fdp_lower_save=2*(1-e)*pnorm(-thres)/(2*(1-e)*pnorm(-thres)+e*x3)
x3=c(seq(1,0.5,-0.05),seq(0.5,0.1,-0.02),seq(0.09,0.01,-0.01));
fdp_lower_save=qfromu_lower(x3,d,e,verbose=1,dx=0.002) #5 sec per point
lines(c(x3,0),c(fdp_lower_save,0),col=2,lwd=3,lty=2)

legend("bottomright",inset=0.02,legend=c('Upper (Lasso) ','Upper (Mobius) ','Lower'),col=c('chartreuse3',4,2),lwd=2.5,cex=1.5)

## zoom-in
par(mar=c(5,5,1,1))
plot(x2[x2<0.51],qfromu_upper(x2,d,e,dx=0.001)[x2<0.51],col=4,lwd=3,xlab='TPP',ylab='FDP',
     cex.lab=2,cex.axis=1.5,type='l',xlim=c(0.2,0.51),ylim=c(0.27,0.45),xaxs='i')
lines(x1[x1>0.2],curve[x1>0.2],lwd=3,col='chartreuse3')
lines(c(x3,0),c(fdp_lower_save,0),col=2,lwd=3,lty=2)
legend("bottomright",inset=0.02,legend=c('Upper (Lasso) ','Upper (Mobius) ','Lower'),col=c('chartreuse3',4,2),lwd=2.5,cex=1.5)

########################################################
############# complete diagram
########################################################
par(mar=c(5,5,1,1))
sss=seq(0,u_star(d,e),0.005)
curve=qfromu_upper(sss,d,e,0.005)
plot(c(sss,rep(u_star(d,e),2)),c(curve,tail(curve,1),1),type='l',ylab='FDP',xlab='TPP',xlim=c(0,1),ylim=c(0,1),
     lwd=3, xaxs = "i", yaxs = "i",cex.lab=2,cex.axis=1.5)
x=seq(0,u_star(d,e),0.001)
lines(x,1-(1-max(curve))/u_star(d,e)*x,lty=2,lwd=2)
tpp_rec=seq(1,0.01,-0.01)
fdp_upper=qfromu_upper(tpp_rec,d,e)
fdp_lower=qfromu_lower(tpp_rec,d,e,dx=0.005)

polygon(c(1,tpp_rec,0,0),c(1,rep(1-e,length(tpp_rec)+1),1), border=NA,col=rgb(1, 0, 0,0.3))#upper unachieve
polygon(c(tpp_rec,0,1),c(fdp_lower,0,0), border=NA,col=rgb(1, 0, 0,0.3))#lower unachieve
polygon(c(1,tpp_rec[tpp_rec>u_star(d,e)],u_star(d,e),e*u_star(d,e)/(1-max(curve))),
        c(1-e,fdp_upper[tpp_rec>u_star(d,e)],qfromu_upper(u_star(d,e),d,e),1-e), border=NA,col=rgb(0,1,0,0.3),lwd=2)#slope achieve
polygon(c(u_star(d,e),tpp_rec[tpp_rec<u_star(d,e)],0,0,e*u_star(d,e)/(1-max(curve))),
        c(qfromu_upper(u_star(d,e),d,e),fdp_upper[tpp_rec<u_star(d,e)],0,1-e,1-e), border=NA,col=rgb(0,0.5,1,0.3),lwd=2)#lasso achieve
polygon(c(1,tpp_rec,rev(tpp_rec)),c(1-e,fdp_upper,rev(fdp_lower)), border=NA,col=8)#possibly achieve by slope but not lasso
lines(c(sss,rep(u_star(d,e),2)),c(curve,tail(curve,1),1),lwd=2)
text(0.7,0.5-e/2,'Unachievable',cex=2)
text(0.5,1-e/2,'Unachievable',cex=2)



########################################################
######## upper bound different than lower
########################################################
d=0.3;e=0.2
par(mar=c(5,5,1,2))
x1=c(seq(0.5,u_star(d,e),0.002),u_star(d,e))
x2=seq(u_star(d,e),0.58,0.001)
curve=qfromu_upper(x1,d,e,dx=0.001)
plot(c(x1,u_star(d,e)),c(curve,1),type='l',ylab='FDP',xlab='TPP',xlim=c(0.5,0.58),ylim=c(0.545,0.66),
     lwd=4, xaxs = "i", yaxs = "i",cex.lab=2,cex.axis=1.5)
lines(x2,qfromu_upper(x2,d,e,dx=0.001),col=4,lwd=4)

start=Sys.time()
tpp_seq=seq(0.58,0.5,-0.004)
fdp_seq=tpp_seq*0
for(index in 1:length(tpp_seq)){
  u=tpp_seq[index]
  a_min=0;a_max=2
  while(a_max-a_min>0.0001){
    a=a_max/2+a_min/2
    t=tfromu(u,a)
    dz=0.01;z=seq(a,15,dz)
    Aopt=pmax(a,-Kp_f(z,t,e)/K_f(z,t,e))
    print(paste0('Monotonicy constraint satisfied: ',all(diff(Aopt)>=0)))
    if(any(diff(Aopt)<0)){break}
    
    int_temp=function(z){return(0.5*K_f(z,t,e)*A_S(z,a,t,e)^2+Kp_f(z,t,e)*A_S(z,a,t,e))}
    SE=integrate(int_temp,a,35)$value+2*(1-e)*ints(-a)+e*ints(-a+t)+e*ints(-a-t)+e*t^2*(pnorm(a-t)-pnorm(-a-t))
    if(SE>d){a_max=a}else{a_min=a}
    print(c(u,2*(1-e)*pnorm(-a)/(2*(1-e)*pnorm(-a)+e*u),a,SE))
  }
  fdp_seq[index]=2*(1-e)*pnorm(-a)/(2*(1-e)*pnorm(-a)+e*u)
}
print(Sys.time()-start) # 14 min

lines(tpp_seq[fdp_seq!=0],fdp_seq[fdp_seq!=0],lwd=4,col=rgb(1,0.5,0))
lines(rep(0.5283,2),c(0.5283,qfromu_upper(0.5283,d,e)),col=rgb(1,0.5,0),lty=2,lwd=2)
abline(v=0.5676,lty=2,lwd=2)
legend("topleft",legend=c(expression(Lasso~tradeoff~q^'*'),expression(Mobius~part~q^'*'),expression(SLOPE~with~constant~pi^'* ')),col=c(1,4,rgb(1,0.5,0)),lwd=3,cex=1.5)
text(0.531,0.555,expression(u^'\u2020'),cex=1.5)
text(0.571,0.555,expression(u[DT]^'*'),cex=1.5)       


########################################################
#### plots of trade-off-dividing areas without simulation
########################################################
d=0.6;e=0.25    # change epsilon
sss=seq(u_star(d,e)-1e-6,0,-0.002)
curve=qfromu_upper(sss,d,e)
estar=epsilonDT(d);


#### purely LASSO with shades of unachievable region
par(mar=c(5,5,1,1))
plot(c(u_star(d,e),sss,0),c(1,curve,0),type='l',ylab='FDP',xlab='TPP',xlim=c(0,1),ylim=c(0,1),
     lwd=2, xaxs = "i", yaxs = "i",cex.lab=2,cex.axis=1.5)
if(u_star(d,e)<1){polygon(c(u_star(d,e),sss,0,1,1),c(1,curve,0,0,1), border=NA,col=rgb(1, 0, 0,0.3))
  }else{polygon(c(sss,1),c(curve,0), border=NA,col=rgb(1, 0, 0,0.3))}
text(0.72,0.26-e/2,'Lasso Unachievable',cex=2)#0.65,0.45-e/2 for d=0.3,e=0.2 or 0.5;0.72,0.26-e/2 for d=0.6,e=0.25

#### four-area plot with shades of unachievable region
par(mar=c(5,5,1,1))
sss=seq(u_star(d,e)-1e-6,0,-0.002)
curve=qfromu_upper(sss,d,e)
plot(c(u_star(d,e),sss,0),c(1,curve,0),type='l',ylab='FDP',xlab='TPP',xlim=c(0,1),ylim=c(0,1),
     lwd=2, xaxs = "i", yaxs = "i",cex.lab=2,cex.axis=1.5)
tpp_rec=seq(0.9999,0.01,-0.02)
fdp_upper=qfromu_upper(tpp_rec,d,e)
fdp_lower=qfromu_lower(tpp_rec,d,e)
#polygon(c(1,1,tpp_rec[tpp_rec>u_star(d,e)],u_star(d,e),u_star(d,e)),c(1,1-e,fdp_upper[tpp_rec>u_star(d,e)],min(fdp_upper[tpp_rec>u_star(d,e)]),1), border=NA,col=rgb(0,0.5,1,0.5))
polygon(c(1,1,tpp_rec,0),c(0,1-e,fdp_lower,0), border=NA,col=rgb(1, 0, 0,0.3))
polygon(c(tpp_rec,1,rev(tpp_rec)),c(fdp_upper,min(fdp_lower),rev(fdp_lower)), border=NA,col=8)
text(0.72,0.26-e/2,'SLOPE Unachievable',cex=2)#0.65,0.45-e/2 for d=0.3,e=0.2 or 0.5
lines(c(u_star(d,e),sss,0),c(1,curve,0),type='l',ylab='FDP',xlab='TPP',xlim=c(0,1),ylim=c(0,1),
      lwd=2, xaxs = "i", yaxs = "i",cex.lab=2,cex.axis=1.5)

## whole hyperbola of Mobius part of trade-off upper bound
par(mar=c(5,5,1,1))
estar=epsilonDT(d)
x=seq(0.01,0.99,0.01)
plot(c(sss,rep(u_star(d,e),2)),c(curve,tail(curve,1),1),type='l',ylab='FDP',xlab='TPP',xlim=c(0,1),ylim=c(0,2),
     lwd=4, xaxs = "i", yaxs = "i",cex.lab=2,cex.axis=1.5)
lines(x,1-x*e/(1-(1-x)*e*(1-estar)/(e-estar)),col=2,lty=2,lwd=4)

x=seq(u_star(d,e),0.99,0.01)
phi2=(x-estar/e)/(1-estar/e)
lines(x,(1-e)*phi2/((1-e)*phi2+e*x),col=3,lwd=4)

########################################################
#### plots of trade-off upper/lower bound curves varying (delta,epsilon)
########################################################
tpp_rec=seq(0.995,0.01,-0.005)

par(mar=c(5,5,1,1))
fdp_upper=qfromu_upper(tpp_rec,0.2,0.2)
fdp_lower=qfromu_lower(tpp_rec,0.2,0.2)
plot(tpp_rec,fdp_upper,type='l',ylab='FDP',xlab='TPP',xlim=c(0,1),ylim=c(0,1),
     lwd=2, xaxs = "i", yaxs = "i",cex.lab=2,cex.axis=1.5,col=2)
lines(tpp_rec,fdp_lower,lwd=2,col=2,lty=2)

fdp_upper=qfromu_upper(tpp_rec,0.3,0.2)
lines(tpp_rec,fdp_upper,lwd=2,col=3)
fdp_lower=qfromu_lower(tpp_rec,0.3,0.2)
lines(tpp_rec,fdp_lower,lwd=2,col=3,lty=2)

fdp_upper=qfromu_upper(tpp_rec,0.5,0.2)
lines(tpp_rec,fdp_upper,lwd=2,col=4)
fdp_lower=qfromu_lower(tpp_rec,0.5,0.2)
lines(tpp_rec,fdp_lower,lwd=2,col=4,lty=2)

legend("topleft",legend=c(expression(paste(delta,"=0.2 ")),expression(paste(delta,"=0.3 ")),
                          expression(paste(delta,"=0.5 "))),
       col=c(2,3,4),lwd=2,cex=1.5)

###
par(mar=c(5,5,1,1))
fdp_upper=qfromu_upper(tpp_rec,0.3,0.1)
fdp_lower=qfromu_lower(tpp_rec,0.3,0.1)
plot(tpp_rec,fdp_upper,type='l',ylab='FDP',xlab='TPP',xlim=c(0,1),ylim=c(0,1),
     lwd=2, xaxs = "i", yaxs = "i",cex.lab=2,cex.axis=1.5,col=2)
lines(tpp_rec,fdp_lower,lwd=2,col=2,lty=2)

fdp_upper=qfromu_upper(tpp_rec,0.4,0.1)
lines(tpp_rec,fdp_upper,lwd=2,col=3)
fdp_lower=qfromu_lower(tpp_rec,0.4,0.1)
lines(tpp_rec,fdp_lower,lwd=2,col=3,lty=2)

fdp_upper=qfromu_upper(tpp_rec,0.5,0.1)
lines(tpp_rec,fdp_upper,lwd=2,col=4)
fdp_lower=qfromu_lower(tpp_rec,0.5,0.1)
lines(tpp_rec,fdp_lower,lwd=2,col=4,lty=2)

legend("topleft",legend=c(expression(paste(delta,"=0.3 ")),expression(paste(delta,"=0.4 ")),
                          expression(paste(delta,"=0.5 "))),
       col=c(2,3,4),lwd=2,cex=1.5)


###
par(mar=c(5,5,1,1))
fdp_upper=qfromu_upper(tpp_rec,0.9,0.1)
fdp_lower=qfromu_lower(tpp_rec,0.9,0.1)
plot(tpp_rec,fdp_upper,type='l',ylab='FDP',xlab='TPP',xlim=c(0,1),ylim=c(0,1),
     lwd=2, xaxs = "i", yaxs = "i",cex.lab=2,cex.axis=1.5,col=2)
lines(tpp_rec,fdp_lower,lwd=2,col=2,lty=2)

fdp_upper=qfromu_upper(tpp_rec,0.9,0.2)
lines(tpp_rec,fdp_upper,lwd=2,col=3)
fdp_lower=qfromu_lower(tpp_rec,0.9,0.2)
lines(tpp_rec,fdp_lower,lwd=2,col=3,lty=2)

fdp_upper=qfromu_upper(tpp_rec,0.9,0.5)
lines(tpp_rec,fdp_upper,lwd=2,col=4)
fdp_lower=qfromu_lower(tpp_rec,0.9,0.5)
lines(tpp_rec,fdp_lower,lwd=2,col=4,lty=2)

legend("topleft",legend=c(expression(paste(epsilon,"=0.1 ")),expression(paste(epsilon,"=0.2 ")),
                          expression(paste(epsilon,"=0.5 "))),
       col=c(2,3,4),lwd=2,cex=1.5)

###
par(mar=c(5,5,1,1))
fdp_upper=qfromu_upper(tpp_rec,0.1,0.1)
fdp_lower=qfromu_lower(tpp_rec,0.1,0.1)
plot(tpp_rec,fdp_upper,type='l',ylab='FDP',xlab='TPP',xlim=c(0,1),ylim=c(0,1),
     lwd=2, xaxs = "i", yaxs = "i",cex.lab=2,cex.axis=1.5,col=2)
lines(tpp_rec,fdp_lower,lwd=2,col=2,lty=2)

fdp_upper=qfromu_upper(tpp_rec,0.1,0.2)
lines(tpp_rec,fdp_upper,lwd=2,col=3)
fdp_lower=qfromu_lower(tpp_rec,0.1,0.2)
lines(tpp_rec,fdp_lower,lwd=2,col=3,lty=2)

fdp_upper=qfromu_upper(tpp_rec,0.1,0.5)
lines(tpp_rec,fdp_upper,lwd=2,col=4)
fdp_lower=qfromu_lower(tpp_rec,0.1,0.5)
lines(tpp_rec,fdp_lower,lwd=2,col=4,lty=2)

legend("bottomright",legend=c(expression(paste(epsilon,"=0.1 ")),expression(paste(epsilon,"=0.2 ")),
                              expression(paste(epsilon,"=0.5 "))),
       col=c(2,3,4),lwd=2,cex=1.5)
