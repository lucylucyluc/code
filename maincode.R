set1_exp = function(dat,par0,k,intebound){
  ## dat: data = list(Y,delta,X,U)
  ## par0: initial value, including parametric part theta0=(alpha0,beta0,nu0) and nonparametric part xi0 
  ## k :  number of inner knots of spline
  ## intebound: including the Upper and lower bounds of integral
  library(splines)
  library(survival)
  library(numDeriv)
  library(orthogonalsplinebasis)
  library(MASS)
  Y = dat[[1]]   ## survival time or censoring time
  delta = dat[[2]]  ## censoring index
  X = dat[[3]]   ##covariates 
  U = dat[[4]]   ## variable with threshold effect
  n = length(Y)
  px = dim(X)[2]
  
  ####  M,I spline 
  od <-3     # order of splines
  p <- k+od  # number of I-spline basis
  
  knotc<-expand.knots(quantile(Y,probs=seq(0,1,1/(k+1))),order=od)
  basc=SplineBasis(knotc,order=od)
  valc=evaluate(basc,Y)  ## Evaluates the basis functions and the points provided in x. Returns a matrix with length(x) rows and dim(object)[2] columns.
  dbasc=deriv(basc)      ###
  valdc=evaluate(dbasc,Y)
  Isp<-matrix(0,p,n)
  for(i in 1:p){
    for (j in 1:n){
      ac<-valc[j,]
      Isp[i,j]<-sum(ac[i:p])
    }
  }
  Msp<-matrix(0,p,n)
  for (i in 1:p){
    for (j in 1:n){
      acd<-valdc[j,]
      Msp[i,j]<-sum(acd[i:p])
    }
  }
  
  
  ## sieve loglikelihood function
  m=300
  inteup = intebound[2]
  intelow = intebound[1]
  grid=seq(intelow,inteup,length.out =m)
  
  seieveloglike_exp = function(par){
    alpha = par[1:px]
    beta  = par[(px+1):(2*px)]
    lam   = par[2*px+1]
    xi    = (par[(2*px+2):(2*px+1+p)])^2
    tt = matrix(rep(grid,n),nrow=m)#m*n
    t1 = matrix(rep(t(xi%*%Msp),m),nrow=n)#n*m
    t2 = matrix(rep(X%*%alpha,m),nrow=n)
    t3 = matrix(rep(X%*%beta,m),nrow=n)
    t4 = matrix(rep(U,m),nrow=n)#n*m
    t5 = t4>t(tt)
    t6 = exp(t2+(t3*t5))
    t7 = matrix(rep(t(xi%*%Isp),m),nrow=n)#n*m
    term1 = t1*t6
    term2 = exp(-t7*t6)
    term3 = lam*exp(-lam*(t(tt)))
    Delta = matrix(rep(delta,m),nrow=n)
    gridd = grid[2]-grid[1]
    Gd    = matrix(rep(gridd,n*m),nrow=n)
    int1  = apply(Delta*term1*term2*term3*Gd,1,sum)
    int2  = apply((1-Delta)*term2*term3*Gd,1,sum)
    lint  = log(int1+int2)
    loglf = -sum(lint)
    return(loglf)
  }
  result =optim(par0,seieveloglike_exp,hessian = TRUE)
  BIC  = result$value + log(n)*(length(par0))
  return(list(result = result,BIC = BIC))
}


set1_norm = function(dat,par0,k,intebound){
  ## dat: data = list(Y,delta,X,U)
  ## par0: initial value, including parametric part theta0=(alpha0,beta0,nu0) and nonparametric part xi0 
  ## k :  number of inner knots of spline
  ## intebound: including the lower and upper  bounds of integral
  library(splines)
  library(survival)
  library(numDeriv)
  library(orthogonalsplinebasis)
  library(MASS)
  
  Y = dat[[1]]   ## survival time or censoring time
  delta = dat[[2]]  ## censoring index
  X = dat[[3]]   ##covariates 
  U = dat[[4]]   ## variable with threshold effect
  n = length(Y)
  px = dim(X)[2]
  ##############-------------- M,I spline ---------------##################
  od <-3     # order of splines
  p <- k+od  # number of I-spline basis
  
  knotc<-expand.knots(quantile(Y,probs=seq(0,1,1/(k+1))),order=od)
  basc=SplineBasis(knotc,order=od)
  valc=evaluate(basc,Y)  ## Evaluates the basis functions and the points provided in x. Returns a matrix with length(x) rows and dim(object)[2] columns.
  dbasc=deriv(basc)      ###
  valdc=evaluate(dbasc,Y)
  
  Isp<-matrix(0,p,n)
  for(i in 1:p){
    for (j in 1:n){
      ac<-valc[j,]
      Isp[i,j]<-sum(ac[i:p])
    }
  }
  Msp<-matrix(0,p,n)
  for (i in 1:p){
    for (j in 1:n){
      acd<-valdc[j,]
      Msp[i,j]<-sum(acd[i:p])
    }
  }
 
  ## sieve loglikelihood function
  m=300
  inteup = intebound[2]
  intelow = intebound[1]
  grid=seq(intelow,inteup,length.out =m)
  
  seieveloglike_norm =function(par){
    
    alpha  = par[1:px]
    beta   = par[(px+1):(2*px)]
    muu    = par[2*px+1]
    sigma  = par[2*px+2]
    xi     = (par[(2*px+3):(2*px+2+p)])^2
    tt = matrix(rep(grid,n),nrow=m)#m*n
    t1 = matrix(rep(t(xi%*%Msp),m),nrow=n)#n*m
    t2 = matrix(rep(X%*%alpha,m),nrow=n)
    t3 = matrix(rep(X%*%beta,m),nrow=n)
    t4 = matrix(rep(U,m),nrow=n)#n*m
    t5 = t4>t(tt)
    t6 = exp(t2+(t3*t5))
    t7 = matrix(rep(t(xi%*%Isp),m),nrow=n)#n*m
    t9=-(t(tt)-muu)^2/(2*sigma^2)
    term1 = t1*t6
    term2 = exp(-t7*t6)
    term3 = exp(t9)/(sqrt(2*pi)*sigma)#lam*exp(-lam*(t(tt)))
    Delta = matrix(rep(delta,m),nrow=n)
    gridd = grid[2]-grid[1]
    Gd = matrix(rep(gridd,n*m),nrow=n)
    int1  = apply(Delta*term1*term2*term3*Gd,1,sum)
    int2  = apply((1-Delta)*term2*term3*Gd,1,sum)
    lint  = log(int1+int2)
    loglf = -sum(log(int1+int2))
    return(loglf)
  }
  result =optim(par0,seieveloglike_norm,hessian = TRUE,control = list(maxit=5000))
  BIC  = result$value + log(n)*(length(par0))
  return(list(result = result,BIC = BIC))
}


set2 = function(dat,par0,k,intebound){
  ## dat: data = list(Y,delta,Z,X,U)
  ## par0: initial value, including parametric part theta0=(alpha0,beta0,nu0) and nonparametric part xi0 
  ## k :  number of inner knots of spline
  ## intebound: including the Upper and lower bounds of integral
  library(splines)
  library(survival)
  library(numDeriv)
  library(orthogonalsplinebasis)
  library(MASS)
  
  Y = dat[[1]]   ## survival time or censoring time
  delta = dat[[2]]  ## censoring index
  Z = dat[[3]]   ## covariates with change point effects
  X = dat[[4]]   ##covariates 
  U = dat[[5]]   ## variable with threshold effect
  n = length(Y)
  px = dim(X)[2]
  pz = dim(Z)[2]
  
  ## M,I spline 
  od <-3     # order of splines
  p <- k+od  # number of I-spline basis
  knotc<-expand.knots(quantile(Y,probs=seq(0,1,1/(k+1))),order=od)
  basc=SplineBasis(knotc,order=od)
  valc=evaluate(basc,Y)  ## Evaluates the basis functions and the points provided in x. Returns a matrix with length(x) rows and dim(object)[2] columns.
  dbasc=deriv(basc)      ###
  valdc=evaluate(dbasc,Y)
  Isp<-matrix(0,p,n)
  for(i in 1:p){
    for (j in 1:n){
      ac<-valc[j,]
      Isp[i,j]<-sum(ac[i:p])
    }
  }
  Msp<-matrix(0,p,n)
  for (i in 1:p){
    for (j in 1:n){
      acd<-valdc[j,]
      Msp[i,j]<-sum(acd[i:p])
    }
  }
  
  ## sieve loglikelihood function
  m=300
  inteup = intebound[2]
  intelow = intebound[1]
  grid=seq(intelow,inteup,length.out =m)
  
  seieveloglike2_norm =function(par){
    omega = par[1:pz]
    alpha  = par[(pz+1):(pz+px)]
    beta   = par[(pz+px+1):(pz+pz*2)]
    muu    = par[pz+pz*2+1]
    sigma  = (par[pz+pz*2+2])^2
    xi     = (par[(pz+pz*2+3):(pz+pz*2+2+p)])^2
    tt = matrix(rep(grid,n),nrow=m)#m*n
    t0 = matrix(rep(Z%*%omega,m),nrow=n)
    t1 = matrix(rep(t(xi%*%Msp),m),nrow=n)#n*m
    t2 = matrix(rep(X%*%alpha,m),nrow=n)
    t3 = matrix(rep(X%*%beta,m),nrow=n)
    t4 = matrix(rep(U,m),nrow=n)#n*m
    t5 = t4>t(tt)
    t6 = exp(t0+t2+(t3*t5))
    t7 = matrix(rep(t(xi%*%Isp),m),nrow=n)#n*m
    t9=-(t(tt)-muu)^2/(2*sigma^2)
    term1 = t1*t6
    term2 = exp(-t7*t6)
    term3 = exp(t9)/(sqrt(2*pi)*sigma)#lam*exp(-lam*(t(tt)))
    Delta = matrix(rep(delta,m),nrow=n)
    gridd = grid[2]-grid[1]
    Gd = matrix(rep(gridd,n*m),nrow=n)
    int1  = apply(Delta*term1*term2*term3*Gd,1,sum)
    int2  = apply((1-Delta)*term2*term3*Gd,1,sum)
    lint  = log(int1+int2)
    loglf = -sum(lint)
    return(loglf)
  }
  result =optim(par0,seieveloglike2_norm,hessian = TRUE)
  BIC  = result$value + log(n)*(length(par0))
  return(list(result = result,BIC = BIC))
}


set3 = function(dat,par0,k,intebound){
  ## dat: data = list(Y,delta,X,U)
  ## par0: initial value, including parametric part theta0=(alpha0,beta0,nu0) and nonparametric part xi0 
  ## k :  number of inner knots of spline
  ## intebound: including the Upper and lower bounds of integral
  library(splines)
  library(survival)
  library(numDeriv)
  library(orthogonalsplinebasis)
  library(MASS)
  
  Y = dat[[1]]   ## survival time or censoring time
  delta = dat[[2]]  ## censoring index
  X = dat[[3]]   ##covariates 
  U = dat[[4]]   ## variable with threshold effect
  n = length(Y)
  px = dim(X)[2]
  
  ##M,I spline
  od <-3     # order of splines
  p <- k+od  # number of I-spline basis
  
  knotc<-expand.knots(quantile(Y,probs=seq(0,1,1/(k+1))),order=od)
  basc=SplineBasis(knotc,order=od)
  valc=evaluate(basc,Y)  ## Evaluates the basis functions and the points provided in x. Returns a matrix with length(x) rows and dim(object)[2] columns.
  dbasc=deriv(basc)      ###
  valdc=evaluate(dbasc,Y)
  Isp<-matrix(0,p,n)
  for(i in 1:p){
    for (j in 1:n){
      ac<-valc[j,]
      Isp[i,j]<-sum(ac[i:p])
    }
  }
  Msp<-matrix(0,p,n)
  for (i in 1:p){
    for (j in 1:n){
      acd<-valdc[j,]
      Msp[i,j]<-sum(acd[i:p])
    }
  }
  
  ## sieve loglikelihood function
  m=300
  inteup = intebound[2]
  intelow = intebound[1]
  grid=seq(intelow,inteup,length.out =m)
  
  seieveloglike =function(par){
    alpha=par[1:px]
    beta=par[(px+1):(2*px)]
    gamma=par[(2*px+1):(3*px)]
    sigma=(par[3*px+1])
    xi=(par[(3*px+2):(3*px+1+p)])^2
    tt=matrix(rep(grid,n),nrow=m)#m*n
    t1=matrix(rep(t(xi%*%Msp),m),nrow=n)#n*m
    t2=matrix(rep(X%*%alpha,m),nrow=n)
    t3=matrix(rep(X%*%beta,m),nrow=n)
    t4=matrix(rep(U,m),nrow=n)#n*m
    t5=t4>t(tt)
    t6=exp(t2+t3*t5)
    
    t7=matrix(rep(t(xi%*%Isp),m),nrow=n)#n*m
    t8=matrix(rep(X%*%gamma,m),nrow=n)
    t9=-(t(tt)-t8)^2/(2*sigma^2)
    term1=t1*t6
    term2=exp(-t7*t6)
    term3=exp(t9)/(sqrt(2*pi)*sigma)
    Delta=matrix(rep(delta,m),nrow=n)
    gridd=grid[2]-grid[1]
    Gd=matrix(rep(gridd,n*m),nrow=n)
    int1=apply(Delta*term1*term2*term3*Gd,1,sum)
    int2=apply((1-Delta)*term2*term3*Gd,1,sum)
    lint = int1+int2
    loglf=-sum(log(lint))
    return(loglf)
  }
  
  result =optim(par0,seieveloglike,hessian = TRUE)
  BIC  = result$value + log(n)*(length(par0))
  return(list(result = result,BIC = BIC))
}


realdata2 = function(dat,par0,k,intebound){
  ## dat: data = list(Y,delta,Z,X,U)
  ## par0: initial value, including parametric part theta0=(alpha0,beta0,nu0) and nonparametric part xi0 
  ## k :  number of inner knots of spline
  ## intebound: including the Upper and lower bounds of integral
  library(splines)
  library(survival)
  library(numDeriv)
  library(orthogonalsplinebasis)
  library(MASS)
  
  Y = dat[[1]]   ## survival time or censoring time
  delta = dat[[2]]  ## censoring index
  Z = dat[[3]]   ## covariates with change point effects
  X = dat[[4]]   ##covariates 
  U = dat[[5]]   ## variable with threshold effect
  n = length(Y)
  px = dim(X)[2]
  pz = dim(Z)[2]
  
  ##M,I spline 
  od <-3     # order of splines
  p <- k+od  # number of I-spline basis
  
  knotc<-expand.knots(quantile(Y,probs=seq(0,1,1/(k+1))),order=od)
  basc=SplineBasis(knotc,order=od)
  valc=evaluate(basc,Y)  ## Evaluates the basis functions and the points provided in x. Returns a matrix with length(x) rows and dim(object)[2] columns.
  dbasc=deriv(basc)      ###
  valdc=evaluate(dbasc,Y)
  Isp<-matrix(0,p,n)
  for(i in 1:p){
    for (j in 1:n){
      ac<-valc[j,]
      Isp[i,j]<-sum(ac[i:p])
    }
  }
  Msp<-matrix(0,p,n)
  for (i in 1:p){
    for (j in 1:n){
      acd<-valdc[j,]
      Msp[i,j]<-sum(acd[i:p])
    }
  }
  
  ## sieve loglikelihood function
  m=300
  inteup = intebound[2]
  intelow = intebound[1]
  grid=seq(intelow,inteup,length.out =m)
  
  seieveloglikehood1 = function(par){
    omega  = par[1:pz]
    alpha  = par[(pz+1):(pz+px)]
    beta   = par[(pz+px+1):(pz+2*px)]
    gamma  = par[(pz+2*px+1):(pz+2*px+pw)]
    sigma  = par[pz+2*px+pw+1]
    xi     = (par[(pz+2*px+pw+2):(pz+2*px+pw+8)])^2  
    tt = matrix(rep(grid,n),nrow=m)#m*n
    t1 = matrix(rep(t(xi%*%Msp),m),nrow=n)#n*m
    t0 = matrix(rep(Z%*%omega,m),nrow=n)#n*m
    t2 = matrix(rep(X%*%alpha,m),nrow=n)
    t3 = matrix(rep(X%*%beta,m),nrow=n)
    t4 = matrix(rep(U,m),nrow=n)#n*m
    t5 = t4>t(tt)
    t6 = exp(t0+t2+t3*t5)
    
    t7 = matrix(rep(t(xi%*%Isp),m),nrow=n)#n*m
    t8 = matrix(rep(W%*%gamma,m),nrow=n)
    t9 = -(t(tt)-t8)^2/(2*sigma^2)#
    term1 = t1*t6
    term2 = exp(-t7*t6)
    term3 = exp(t9)/(sqrt(2*pi)*sigma)
    Delta = matrix(rep(delta,m),nrow=n)
    gridd = grid[2]-grid[1]
    Gd    = matrix(rep(gridd,n*m),nrow=n)
    int1  = apply(Delta*term1*term2*term3*Gd,1,sum)
    int2  = apply((1-Delta)*term2*term3*Gd,1,sum)
    lint  = int1+int2
    loglf = -sum(log(lint))
    return(loglf)
  }
  result = optim(par0,seieveloglikehood1,hessian = TRUE,control = list(maxit=1000,abstol=10^-6,reltol=10^-6))
  BIC  = result$value + log(n)*(length(par0))
  return(list(result = result,BIC = BIC))
}



