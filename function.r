library(gglasso)
library(glmnet)

data.generate.k.fluctuate2 = function(d0,change.point,p,n,sigma,kappa){
  seed = 10
  true_change_point = change.point
  ncp = length(change.point)
  X = matrix(rnorm(p*n,0,1),p,n)
  y = matrix(0,n,1)
  nonzero.element.loc = c(1:d0)
  cp = c(0,true_change_point,n)
  beta = matrix(0,p,ncp+1)
  betafullmat = matrix(0,p,n)
  for (i in 1:(ncp+1)){
    if (i%%2 == 1){
      beta[nonzero.element.loc,i] = kappa/(2*sqrt(d0))
    }
    else{
      beta[nonzero.element.loc,i] = -kappa/(2*sqrt(d0))
    }
    y[(1+cp[i]):cp[i+1],] = rnorm(cp[i+1] - cp[i],t(X[,(1+cp[i]):cp[i+1]])%*%beta[,i],sigma)
    for (j in (1+cp[i]):cp[i+1]){
      betafullmat[,j] = beta[,i] 
    }
  }
  
  List = list(true.change.point = true_change_point, X = X, y = y, betafullmat = betafullmat)
  return(List)
}

DPR = function(gamma,lambda,dataX, datay,gridsize,delta){
  ## N must be the multiple of gridsize
  N = ncol(dataX)
  p = nrow(dataX)
  Bestvalue = c()
  partition = c()
  Bestvalue = c(Bestvalue,-gamma*log(max(N,p)))
  #gridsize = 1
  for (r in 1:N){
    b = Inf
    if (r %% gridsize == 0){
      b = sapply(1:r,function(l)inner.func(Bestvalue,gamma,lambda,dataX,datay,l,r,gridsize,delta))
      Bestvalue = c(Bestvalue,min(b))
      partition = c(partition, (which.min(b)/gridsize - 1))           
    }
  }              
  R = N/gridsize
  L = c()
  while (R > 0) {
    L = c(partition[R], L)
    R = partition[R]
  }
  L = L[-1]              
  return(L*gridsize)
}     


inner.func = function(Bestvalue,gamma,lambda,dataX,datay,l,r,gridsize,delta){
  btemp = Inf
  N = ncol(dataX)
  p = nrow(dataX)
  if (l %% gridsize == 0){
    btemp = Bestvalue[l/gridsize] + gamma*log(max(N,p)) + distanceR(l,r,dataX,datay,lambda,delta)
  } 
  return(btemp)
}
# residual
distanceR = function(s, e, dataX,datay,lambda,delta){
  n = ncol(dataX)
  p = nrow(dataX)
  if (abs(s-e) >delta){
    fit = glmnet(x=t(dataX[, s:e]), y=datay[s:e,], family=c("gaussian"),
                 alpha = 1, lambda=lambda*sqrt(max(log(max(n,p)),e-s))*sqrt(log(max(n,p)))/(e-s),intercept=F)
    coef_est = t(t(as.vector(fit$beta)))
    yhat = t(dataX[,s:e])%*%coef_est
    d = norm(datay[s:e] - yhat, type = "2")
  }
  else{
    d = Inf
  }
  return(d^2)
}                 

## divide X into [s,eta] and [eta+1,e]
convert.design.matrix.one.change=function(X,eta,s_ceiling){
  ee=ncol(X)
  xx1=t(X)
  t = eta - s_ceiling +1
  xx1[ (t+1):ee,]=0
  xx2=t(X)
  xx2[1:t,]=0
  xx=cbind(xx1/sqrt(t-1),xx2/sqrt(ee-t+1))
  return(xx)
}

LocalRefine = function(init_cp,X,y,zeta){
  w = 0.9
  partition = c(0)
  n = ncol(X)
  init_cp_ = c(0,init_cp,n)
  K = length(init_cp_) - 2
  for (k in 1:K){
    s = w*partition[k] + (1-w)*init_cp_[k+1]
    e = (1-w)*init_cp_[k+1] + w*init_cp_[k+2]
    b = c()
    lower = ceiling(s)+1
    upper = floor(e)-1
    b = sapply(lower:upper,function(eta)inner.func.lr(s,e,eta,X,y,zeta))
    partition = c(partition, ceiling(s) + which.min(b) )
  }
  return(partition[-1])
}

inner.func.lr = function(s,e,eta,X,y,zeta){
  p = nrow(X)
  group = rep(1:nrow(X),2)
  convertX = convert.design.matrix.one.change(X[,(ceiling(s)):(floor(e))],eta,ceiling(s))
  y_ = y[(ceiling(s)):(floor(e)),]
  lambda.LR = zeta*sqrt(log(max(ncol(X),nrow(X))))
  auxfit = gglasso(x=convertX,y=y_,group=group, loss="ls",
                   lambda=lambda.LR/(floor(e)-ceiling(s)+1),intercept = FALSE,eps =
                     0.001)
  coef = as.vector(auxfit$beta)
  coef1 = coef[1:p]
  coef2 = coef[(p+1):(2*p)]
  btemp = norm(y_ - convertX %*% coef, type = "2")^2 + lambda.LR*sum(sqrt(coef1^2 + coef2^2))
  return(btemp)
}


Hausdorff = function(vec1,vec2){
  dist = matrix(0,nrow = length(vec1),ncol = length(vec2))
  for (i in 1:nrow(dist)){
    for (j in 1:ncol(dist)){
      dist[i,j] = abs(vec1[i] - vec2[j])
    }
  }
  
  dH = max(max(apply(dist, 2, function(x) min(x) )), max(apply(dist, 1, function(x) min(x) )))
  return(dH)
}                                                             


# CV DP                                                               
distance.cv.dp = function(s, e, X,y, dataX,datay,lambda){
  n = ncol(dataX)
  p = nrow(dataX)
  lower = s-1
  upper = e
  fit = glmnet(x=t(dataX[, lower:upper]), y=datay[lower:upper,], family=c("gaussian"),alpha = 1, 
               lambda=lambda*sqrt(max(log(max(n,p)),upper-lower))*sqrt(log(max(n,p)))/(upper-lower),intercept=F)
  
  coef_est = t(t(as.vector(fit$beta)))
  yhat = t(dataX[,s:e])%*%coef_est
  d = norm(datay[s:e,] - yhat, type = "2")
  result = list("mse" = d^2, "beta" = coef_est)
  return(result)
}

residual = function(lower, upper,dataX,datay,beta_est){
  res = norm(datay[lower:upper] - t(dataX[,lower:upper])%*%beta_est, type = "2")^2
  return(res)
} 


cv.grid.dp = function(lambda,gamma,X,y,grid,delta){
  n = ncol(X)
  even_indexes = seq(2,n,2)
  odd_indexes = seq(1,n,2)
  train.X = X[,odd_indexes]
  train.y = y[odd_indexes,]
  train.y = as.matrix(train.y)
  validation.X = X[,even_indexes]
  validation.y = y[even_indexes,]
  colnames(train.X) = c()
  row.names(train.y) = c()
  init_cp_train = DPR(gamma,lambda,train.X, train.y,grid,delta)
  init_cp_train.long = c(0,init_cp_train,ncol(train.X))
  diff.point = diff(init_cp_train.long)
  if (length(which(diff.point == 1)) > 0){
      print(paste("gamma =", gamma,",", "lambda =", lambda, ".","Warning: Consecutive points detected. Try a larger gamma."))
      init_cp = c()
      for (i in init_cp_train){
          init_cp = c(init_cp,odd_indexes[i])
      }
      len = length(init_cp)
      valiloss = Inf
      trloss = Inf
      result = list("cp" = init_cp, "K" = len, "res" = valiloss,"trloss" = trloss)
      return(result)
  }
  else{
      init_cp = c()
      for (i in init_cp_train){
          init_cp = c(init_cp,odd_indexes[i])
      }
      #track on the process
      #print("error")
      #print(Hausdorff(init_cp,GT))
      #print(init_cp)
      #print(c(lambda,gamma))
      len = length(init_cp)
      init_cp_long = c(init_cp_train,n/2)
      interval = matrix(0,nrow = len + 1,ncol = 2)
      interval[1,] = c(1,init_cp_long[1])
      if (len > 0){
          for (j in 2:(1+len)){
              interval[j,] = c(init_cp_long[j-1]+1,init_cp_long[j])
          }
      }
      p = nrow(train.X)
      trainmat = sapply(1:(len+1),function(index) distance.cv.dp(interval[index,1], interval[index,2], X, y, train.X,train.y,lambda))
      betamat = matrix(0,nrow = p,ncol = len+1)
      training_loss = matrix(0,nrow = 1,ncol = len+1)
      for (col in 1:(len+1)){
          betamat[,col] = as.numeric(trainmat[2,col]$beta)
          training_loss[,col] = as.numeric(trainmat[1,col]$mse)
      }
      validationmat = sapply(1:(len+1),function(index) residual(interval[index,1], interval[index,2], validation.X,validation.y,betamat[,index]))
      result = list("cp" = init_cp, "K" = len, "res" = sum(validationmat),"trloss" = sum(training_loss))
      return(result)
  }
}   


cv.grid.search.dp = function(lambda.set,gamma.set,X,y,grid,delta){
  output = sapply(1:length(lambda.set), function(i) sapply(1:length(gamma.set), 
                                                           function(j) cv.grid.dp(lambda.set[i],gamma.set[j],X,y,grid,delta)))
  print(output)
  cp = output[seq(1,4*length(gamma.set),4),]## estimated change points
  len = output[seq(2,4*length(gamma.set),4),]## number of estimated change points
  ply = output[seq(3,4*length(gamma.set),4),]## validation loss
  trloss = output[seq(4,4*length(gamma.set),4),]## training loss                                                      
  result = list("changepoint" = cp, "K" = len, "res" = ply, "trloss" = trloss)
  return(result)
}                                                             



# CV DPLR                                                         


distance.cv.lr = function(s, e, X,y, dataX,datay,zeta){
  #dataX is training X not X
  n = ncol(dataX)
  p = nrow(dataX)
  lower = s-1
  upper = e
  group = c(1:nrow(dataX))
  convertX = convert.matrix.1.group(dataX[,(ceiling(s)):(floor(e))])
  y_ = datay[(ceiling(s)):(floor(e)),]
  lambda.LR = zeta*sqrt(log(max(n,p)))
  fit = glmnet(x=t(dataX[,lower:upper]), y=datay[lower:upper,], family=c("gaussian"),
               alpha = 1,lambda = lambda.LR,intercept=F)
  coef_est = t(t(as.vector(fit$beta)))
  yhat = t(dataX[,s:e])%*%coef_est
  d = norm(datay[s:e,] - yhat, type = "2")
  result = list("mse" = d^2, "beta" = coef_est)
  return(result)
}                                                            

convert.matrix.1.group=function(X){
  ee=ncol(X)
  xx1=t(X)
  xx=xx1#/sqrt(ee)
  return(xx)
}    

cv.grid.dp.lr = function(lambda,gamma,zeta,X,y,grid,delta){
  n = ncol(X)
  even_indexes = seq(2,n,2)
  odd_indexes = seq(1,n,2)
  train.X = X[,odd_indexes]
  train.y = y[odd_indexes,]
  train.y = as.matrix(train.y)
  validation.X = X[,even_indexes]
  validation.y = y[even_indexes,]
  colnames(train.X) = c()
  row.names(train.y) = c()
  init_cp_train.dp = DPR(gamma,lambda,train.X, train.y,grid,delta)
  if (length(init_cp_train.dp) != 0){
    init_cp.dp = c()
    for (i in init_cp_train.dp){
      init_cp.dp = c(init_cp.dp,odd_indexes[i])
    }
    init_cp = LocalRefine(init_cp.dp, X, y,zeta)
  }
  else{
    init_cp.dp = c()
    init_cp = c()
  }
  
  len = length(init_cp)
  init_cp_train = (1+init_cp)/2
  init_cp_long = c(init_cp_train,n/2)
  interval = matrix(0,nrow = len + 1,ncol = 2)
  interval[1,] = c(1,init_cp_long[1])
  if (len > 0){
    for (j in 2:(1+len)){
      interval[j,] = c(init_cp_long[j-1]+1,init_cp_long[j])
    }
  }
  p = nrow(train.X)
  trainmat = sapply(1:(len+1),function(index) distance.cv.lr(interval[index,1], interval[index,2], X,y,train.X,train.y,zeta))
  betamat = matrix(0,nrow = p,ncol = len+1)
  training_loss = matrix(0,nrow = 1,ncol = len+1)                
  for (col in 1:(len+1)){
    betamat[,col] = as.numeric(trainmat[2,col]$beta)
    training_loss[,col] = as.numeric(trainmat[1,col]$mse)
  }      
  validationmat = sapply(1:(len+1),function(index) residual(interval[index,1], interval[index,2], validation.X,validation.y,betamat[,index]))                       
  result = list("cp" = init_cp, "K" = len, "res" = sum(validationmat),"trloss" = sum(training_loss))                       
  return(result)
}   


cv.grid.search.dp.lr.lg = function(lambda.set,gamma.set,zeta,X,y,grid,delta){
  output = sapply(1:length(lambda.set), function(i) sapply(1:length(gamma.set), 
                                                           function(j) cv.grid.dp.lr(lambda.set[i],gamma.set[j],zeta, X,y,grid,delta)))
  print(output)
  cp = output[seq(1,4*length(gamma.set),4),]## estimated change points
  len = output[seq(2,4*length(gamma.set),4),]## number of estimated change points
  ply = output[seq(3,4*length(gamma.set),4),]## validation loss
  trloss = output[seq(4,4*length(gamma.set),4),]## training loss                                                         
  result = list("changepoint" = cp, "K" = len, "res" = ply, "trloss" = trloss)
  return(result)
}                           

cv.grid.search.dp.lr = function(lambda.set,gamma.set,zeta.set,X,y,grid,delta){
  output.2 = sapply(1:length(zeta.set), function(q) cv.grid.search.dp.lr.lg(lambda.set,gamma.set,zeta.set[q], X,y,grid,delta))
  print("output with zeta")
  print(output.2) 
  cp = output.2[1,]## estimated change points
  len = output.2[2,]## number of estimated change points
  ply = output.2[3,]## validation loss
  trloss = output.2[4,]## training loss   
  return(output.2)                  
}                                              
