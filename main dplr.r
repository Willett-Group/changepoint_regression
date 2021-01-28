source("./function.r")
library(gglasso)
library(glmnet)

set.seed(0)
start.all = Sys.time()
print('starttime')
print(start.all)
d0 = 15## the number of nonzero elements = p-sparsity
RR = 1##replication
grid = 1##gridsize to save computational cost for DP but might increase the error, default = 1
sigma = 1
kappa = 5
delta = 5 ## minimal spacing to help DP more robust
gamma.dp.set = c(1.125,5.625)
lambda.dp.set = c(0.195,0.385)
zeta.dp = c(0.01,0.05,0.1,1)#c(0.01,0.05,0.1)
len.est = matrix(0,RR,1)
lr.change.point = matrix(0,20,RR)
p = 200 
n = 600## sample size, need to be the multiple of 2*gridsize
change.point = c(120,220,350,450)
len.est = matrix(0,RR,1)
error = matrix(0,RR,1)
nrowmat = length(gamma.dp.set)*length(zeta.dp)
ncolmat = length(lambda.dp.set)
for (time in 1:RR){
  data = data.generate.k.fluctuate2(d0,change.point,p,n,sigma,kappa)
  X = data$X
  y = data$y
  true = data$true.change.point
  result = cv.grid.search.dp.lr(lambda.dp.set,gamma.dp.set,zeta.dp,X,y,grid,delta)
  ##output the table containing estimated change points, number of estimated change points, 
  ##validation loss and training error.
  resmat = c()
  cp = c()
  K = c()
  for (nzeta in 1:length(zeta.dp)){
    cp = rbind(cp,result[[(nzeta-1)*4+1]])## estimated change points
    resmat = rbind(resmat,result[[(nzeta-1)*4+3]])## validation loss
    K = rbind(K,result[[(nzeta-1)*4+2]])## number of estimated change points
  }
  haus.index = which(matrix(unlist(resmat),nrowmat,ncolmat) == min(matrix(unlist(resmat),nrowmat,ncolmat)), arr.ind=TRUE)[1,]
  # in case there are multiple min values
  selected.change.point = unlist(cp[haus.index[1],haus.index[2]])
  print(" selected change points =")
  print(selected.change.point)
  lr.change.point[1:length(selected.change.point), time] = selected.change.point
  #write.csv(lr.change.point, file = "1_dplr_d_0_15.csv")
  error[time,1] = Hausdorff(selected.change.point,change.point)
  #write.csv(error, file = "1_dplr_error_d_0_15.csv")
  len.est[time,1] = unlist(K[haus.index[1],haus.index[2]])
  #write.csv(len.est, file = "1_dplr_len_d_0_15.csv")
}
lr.change.point## estimated change points
error
len.est## number of estimated change points
end.all = Sys.time()
print('duration')
print(difftime(end.all, start.all, units = "secs")[[1]])
