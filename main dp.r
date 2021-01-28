source("./function.r")
library(gglasso)
library(glmnet)

set.seed(0)
start.all = Sys.time()
print('starttime')
print(start.all)
d0 = 15# the number of nonzero elements = p-sparsity
RR = 1## replication
grid = 1##gridsize to save computational cost for DP but might increase the error, default = 1
sigma = 1## variance of data
kappa = 5## spectral norm
delta = 5 ## minimal spacing to help DP more robust
n = 600## sample size, need to be the multiple of 2*gridsize
p = 200
change.point = c(120,220,350,450)
gamma.dp.set = c(1.125,5.625,11.250)
lambda.dp.set = c(0.195,0.385,1.935)
len.est = matrix(0,RR,1)
init.change.point = matrix(0,20,RR)
len.est = matrix(0,RR,1)
error = matrix(0,RR,1)
nrowmat = length(gamma.dp.set)
ncolmat = length(lambda.dp.set)
for (time in 1:RR){
  data = data.generate.k.fluctuate2(d0,change.point,p,n,sigma,kappa)
  X = data$X
  y = data$y
  true = data$true.change.point
  result = cv.grid.search.dp(lambda.dp.set,gamma.dp.set,X,y,grid,delta)
  ##output a table which contains estimated change points, number of estimated change points, 
  ##validation loss and training error.
  K = matrix(as.numeric(result$K),nrowmat,ncolmat)## number of estimated change points
  resmat = as.numeric(result$res)## validation loss
  cp = result$changepoint## estimated change points
  haus.index = which(matrix(resmat,nrowmat,ncolmat) == min(matrix(resmat,nrowmat,ncolmat)), arr.ind=TRUE)
  K.from.pair = c()
  for(i in 1:nrow(haus.index)){
    K.from.pair = c(K.from.pair, K[haus.index[i,1],haus.index[i,2]])
  }
  I = which.min(K.from.pair)
  haus.pair = c(gamma.dp.set[haus.index[I,1]],lambda.dp.set[haus.index[I,2]])
  # in case there are multiple min values
  selected.change.point = as.numeric(unlist(cp[haus.index[I,1],haus.index[I,2]]))
  print(" selected change points =")
  print(selected.change.point)
  init.change.point[1:length(selected.change.point), time] = selected.change.point
  #write.csv(init.change.point, file = "kappa_5_dp_d_0_15.csv")
  error[time,1] = Hausdorff(selected.change.point,change.point)
  #write.csv(error, file = "kappa_5_dp_error_d_0_15.csv")
  len.est[time,1] = K[haus.index[I,1],haus.index[I,2]]
  #write.csv(len.est, file = "kappa_5_dp_ncp_d_0_15.csv")
}
init.change.point ## estimated change points
error
len.est## number of estimated change points
end.all = Sys.time()
print('duration')
print(difftime(end.all, start.all, units = "secs")[[1]])
