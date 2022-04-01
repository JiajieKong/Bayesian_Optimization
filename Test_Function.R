##Function_any_dimension##
rm(list=ls())
library(GA)
library(tgp)
library(laGP)
library(mvtnorm)
library(ggplot2)
library(gridExtra)
library(CompModels)
library(scatterplot3d)
library(ggplot2)

bov <- function(y, end = length(y)){ # Calculates best overall value
  prog <- rep(min(y), end)
  prog[1:min(end, length(y))] <- y[1:min(end, length(y))]
  for(i in 2:end)
    if(is.na(prog[i]) || prog[i] > prog[i-1]) prog[i] <- prog[i-1]
  return(prog)
}

obj.UEI <- function(x, fmin, gpi, j, med, beta, pred=predGPsep){
  - UEI(gpi, x, fmin, j, med, beta)
} 

eps <- .Machine$double.eps 

UEI.search <- function(X, y, gpi, j, med, beta, A, pred=predGPsep, multi.start=5, tol=eps){
  
  l = dim(A)[1]
  m <- which.min(y)
  fmin <- y[m]
  start <- matrix(X[m,], nrow=1)
  if(multi.start > 1)
    start <- rbind(start, lhs(multi.start - 1, A))
  xnew <- matrix(NA, nrow = nrow(start), ncol = ncol(X) + 1)
  for(i in 1:nrow(start)) {
    if(UEI(gpi, start[i,], fmin, j, med, beta) <= tol) {
      #print(UEI(gpi, start[i,], fmin, j, med, beta))
      out <- list(value = -Inf); next }
    out <- optim(start[i,], obj.UEI, method = "L-BFGS-B",
                 lower = A[,1], 
                 upper = A[,2], 
                 gpi = gpi, pred = pred, fmin = fmin, j = j, med = med, beta = beta)
    xnew[i,] <- c(out$par, -out$value)
  }
  solns <- data.frame(cbind(start, xnew))
  names(solns) <- c(paste("s", 1:l, sep = ""), paste("x", 1:l, sep = ""), "val")
  #solns <- solns[solns$val > tol,]
  return(solns)
}

optim.UEI <- function(f, ninit, end, med, beta, A, seed){
  ## initialization
  l = dim(A)[1]
  set.seed(seed)
  X <- lhs(ninit, A)
  y <- f(X)
  gpi <- newGPsep(X, y, d = 0.1, g = 1e-6, dK = TRUE)
  da <- darg(list(mle = TRUE), lhs(1000, A))
  mleGPsep(gpi, param = "d", tmin = da$min, tmax = da$max, ab = da$ab)
  ## optimization loop of sequential acquisitions
  maxuei <- c()
  for(j in (ninit + 1):end) {
    solns <- UEI.search(X, y, gpi, j, med, beta, A)
    m <- which.max(solns$val)
    maxuei <- c(maxuei, solns$val[m])
    xnew <- as.matrix(solns[m, (l+1):(l+l)])
    ynew <- f(xnew)
    updateGPsep(gpi, xnew, ynew)
    mleGPsep(gpi, param = "d", tmin = da$min, tmax = da$max, ab = da$ab)
    X <- rbind(X, xnew)
    y <- c(y, ynew)
  }
  deleteGPsep(gpi)
  return(list(X = X, y = y, maxuei = maxuei))
}

UEI <- function(gpi, x, fmin, j, method, beta, pred = predGPsep){
  
  if(is.null(nrow(x))) x <- matrix(x, nrow=1)
  p <- pred(gpi, x, lite=TRUE)
  d <- fmin - p$mean
  sigma <- sqrt(p$s2)
  dn <- d / sigma
  ei <- d * pnorm(dn) + sigma * dnorm(dn)
  vi <- sigma^2 * ((dn^2 + 1) * pnorm(dn) + dn * dnorm(dn)) - ei^2
  vi <- ifelse(vi < 0, 0, vi) 
  if (vi == 0){
    sei <- 0
  }else{
    sei <- ei/sqrt(vi)
  }
  uei <- ei + beta * sqrt(vi)
  vei <- ei + beta * vi
  dei <- ei + 80/j * sqrt(vi)
  pei <- ei^2 + vi
  
  if (method == "uei"){
    return(uei)
  }else if (method == "sei"){
    return(sei)
  }else if (method == "dei"){
    return(dei)  
  }else if (method == "pei"){
    return(pei)
  }else{
    return(vei)
  }
}

#########Main#################
# test function
#1D
f.herbie <- function(X){
  ans <- sin(10 * pi * X) / (2 * X) + (X - 1)^4
  return(ans)
}
#2D
f.Rosenbrock <- function(X){
  ans <- (1-X[,1])^2+100*(X[,2]-X[,1]^2)^2
  return(ans)
}
f.HumpCamel <- function(X){
  ans <- -(cos((X[,1] - .1) * X[,2])^2) - X[,1] * sin (3 * X[,1] + X[,2])
  return(ans)
}
f.Ackley <- function(X){
  ans <- -20*exp(-0.2*sqrt(0.5*(X[,1]^2+X[,2]^2))) - exp(0.5*(cos(2*pi*X[,1])+cos(2*pi*X[,2]))) + exp(1) + 20
  return(ans)
}
f.rastrigin2 <- function(X){
  ans <- 10 * 2 + X[,1]^2 + X[,2]^2 - 10 * cos(2*pi*X[,1]) - 10 * cos(2*pi*X[,2])
  return(ans)
}
#6D
f.Hartmann <- function(X){
  f_tmp = function(X_matrix){(-sum(b* exp(-apply(A*(X_matrix - P)^2, MARGIN = 2, FUN = sum))))}
  b = c(1, 1.2, 3, 3.2)
  A = matrix(c(10, 3, 17, 3.5, 1.7, 8,
               0.05, 10, 17, 0.1, 8, 14,
               3, 3.5, 1.7, 10, 17, 8,
               17, 8, 0.05, 10, 0.1, 14), nrow = 6)
  P = 10^(-4)*matrix(c(1312, 1696, 5569, 124, 8283, 5886,
                        2329, 4135, 8307, 3736, 1004, 9991,
                        2348, 1451, 3522, 2883, 3047, 6650,
                        4047, 8828, 8732, 5743, 1091, 381), nrow = 6)
  
  X_array = array(rep(X, 4), dim = c(dim(X)[1],dim(X)[2],4))
  ans <- apply(X_array, MARGIN = 1, f_tmp)
  return(ans)
}
#2D
#f.GoldsteinPrice <- function(X){
  ans <- (1+(X[,1]+X[,2]+1)^2*(19-14*X[,1]+3*X[,1]^2-14*X[,2]+6*X[,1]*X[,2]+3*X[,2]^2))*
  (30+(2*X[,1]-3*X[,2])^2*(18-32*X[,1]+12*X[,1]^2+48*X[,2]-36*X[,1]*X[,2]+27*X[,2]^2))
  return(ans)
}
#4D
#f.shekel <- function(X){
  f_tmp = function(X_matrix){-sum(1/(apply((X_matrix - C)^2, MARGIN = 2, FUN = sum) + b))}
  b = c(1,2,2,4,4,6,3,7,5,5)/10
  C = matrix(c(4,1,8,6,3,2,5,8,6,7,
               4,1,8,6,7,9,3,1,2,3.6,
               4,1,8,6,3,2,5,8,6,7,
               4,1,8,6,7,9,3,1,2,3.6), byrow = TRUE, nrow = 4)
  X_array = array(rep(X, 10), dim = c(dim(X)[1],dim(X)[2],10))
  ans <- apply(X_array, MARGIN = 1, f_tmp)
  return(ans)
}
#10D
#f.rastrigin10 <- function(X){
  ans <- 10 * 10 + apply(X^2 - 10*cos(2*pi*X), MARGIN = 1, FUN = sum)
  return(ans) 
} ##

# Monte Carlo experiment setup
ninit <- 10
end <- 500
reps <- 100

main = function(f, method, beta, A, ninit = 10, end = 100, reps = 100){
  
  os <- list()
  prog<- matrix(NA, nrow=reps, ncol=end)
  
  for(r in 1:reps) {
    print(r)
    os[[r]] <- optim.UEI(f, ninit, end, method, beta, A, seed = r*35+140)
    prog[r,] <- bov(os[[r]]$y)
  }
  
  return(list(os = os, prog = prog))
}

res = list()

system.time({
  A = matrix(c(.5,2.5), ncol = 2)
  res[[1]] = list()
  res[[1]][[1]] = main(f.herbie ,"uei", 0, A, ninit, end, reps)
  res[[1]][[2]] = main(f.herbie, "uei", .5, A, ninit, end, reps)
  res[[1]][[3]] = main(f.herbie, "uei", 2, A, ninit, end, reps)
  res[[1]][[4]] = main(f.herbie, "uei", 100, A, ninit, end, reps)
  res[[1]][[5]] = main(f.herbie, "dei", 1, A, ninit, end, reps)
  res[[1]][[6]] = main(f.herbie, "vei", -.5, A, ninit, end, reps)
  res[[1]][[7]] = main(f.herbie, "sei", 1, A, ninit, end, reps)
  res[[1]][[8]] = main(f.herbie, "pei", 100, A, ninit, end, reps)
})
system.time({
  A = matrix(c(-2,-2, 2,2), ncol = 2)
  res[[2]] = list()
  res[[2]][[1]] = main(f.Rosenbrock ,"uei", 0, A, ninit, end, reps)
  res[[2]][[2]] = main(f.Rosenbrock, "uei", .5, A, ninit, end, reps)
  res[[2]][[3]] = main(f.Rosenbrock, "uei", 2, A, ninit, end, reps)
  res[[2]][[4]] = main(f.Rosenbrock, "uei", 100, A, ninit, end, reps)
  res[[2]][[5]] = main(f.Rosenbrock, "dei", 1, A, ninit, end, reps)
  res[[2]][[6]] = main(f.Rosenbrock, "vei", -.5, A, ninit, end, reps)
  res[[2]][[7]] = main(f.Rosenbrock, "sei", 1, A, ninit, end, reps)
  res[[2]][[8]] = main(f.Rosenbrock, "pei", 100, A, ninit, end, reps)
})
system.time({
  A = matrix(c(-2,-2, 2,2), ncol = 2)
  res[[3]] = list()
  res[[3]][[1]] = main(f.HumpCamel ,"uei", 0, A, ninit, end, reps)
  res[[3]][[2]] = main(f.HumpCamel, "uei", .5, A, ninit, end, reps)
  res[[3]][[3]] = main(f.HumpCamel, "uei", 2, A, ninit, end, reps)
  res[[3]][[4]] = main(f.HumpCamel, "uei", 100, A, ninit, end, reps)
  res[[3]][[5]] = main(f.HumpCamel, "dei", 1, A, ninit, end, reps)
  res[[3]][[6]] = main(f.HumpCamel, "vei", -.5, A, ninit, end, reps)
  res[[3]][[7]] = main(f.HumpCamel, "sei", 1, A, ninit, end, reps)
  res[[3]][[8]] = main(f.HumpCamel, "pei", 100, A, ninit, end, reps)
})
system.time({
  A = matrix(c(-2,-2, 2,2), ncol = 2)
  res[[4]] = list()
  res[[4]][[1]] = main(f.Ackley ,"uei", 0, A, ninit, end, reps)
  res[[4]][[2]] = main(f.Ackley, "uei", .5, A, ninit, end, reps)
  res[[4]][[3]] = main(f.Ackley, "uei", 2, A, ninit, end, reps)
  res[[4]][[4]] = main(f.Ackley, "uei", 100, A, ninit, end, reps)
  res[[4]][[5]] = main(f.Ackley, "dei", 1, A, ninit, end, reps)
  res[[4]][[6]] = main(f.Ackley, "vei", -.5, A, ninit, end, reps)
  res[[4]][[7]] = main(f.Ackley, "sei", 1, A, ninit, end, reps)
  res[[4]][[8]] = main(f.Ackley, "pei", 100, A, ninit, end, reps)
})
system.time({
  A = matrix(c(-2,-2, 2,2), ncol = 2)
  res[[5]] = list()
  res[[5]][[1]] = main(f.rastrigin2 ,"uei", 0, A, ninit, end, reps)
  res[[5]][[2]] = main(f.rastrigin2, "uei", .5, A, ninit, end, reps)
  res[[5]][[3]] = main(f.rastrigin2, "uei", 2, A, ninit, end, reps)
  res[[5]][[4]] = main(f.rastrigin2, "uei", 100, A, ninit, end, reps)
  res[[5]][[5]] = main(f.rastrigin2, "dei", 1, A, ninit, end, reps)
  res[[5]][[6]] = main(f.rastrigin2, "vei", -.5, A, ninit, end, reps)
  res[[5]][[7]] = main(f.rastrigin2, "sei", 1, A, ninit, end, reps)
  res[[5]][[8]] = main(f.rastrigin2, "pei", 100, A, ninit, end, reps)
}) 
system.time({
  A = matrix(c(rep(0, 6), rep(1, 6)), ncol = 2)
  res[[6]] = list()
  res[[6]][[1]] = main(f.Hartmann ,"uei", 0, A, ninit, end, reps)
  res[[6]][[2]] = main(f.Hartmann, "uei", .5, A, ninit, end, reps)
  res[[6]][[3]] = main(f.Hartmann, "uei", 2, A, ninit, end, reps)
  res[[6]][[4]] = main(f.Hartmann, "uei", 100, A, ninit, end, reps)
  res[[6]][[5]] = main(f.Hartmann, "dei", 1, A, ninit, end, reps)
  res[[6]][[6]] = main(f.Hartmann, "vei", -.5, A, ninit, end, reps)
  res[[6]][[7]] = main(f.Hartmann, "sei", 1, A, ninit, end, reps)
  res[[6]][[8]] = main(f.Hartmann, "pei", 100, A, ninit, end, reps)
}) 
system.time({
  A = matrix(c(0,0,0,0,10,10,10,10), ncol = 2)
  res[[7]] = list()
  res[[7]][[1]] = main(f.shekel ,"uei", 0, A, ninit, end, reps)
  res[[7]][[2]] = main(f.shekel, "uei", .5, A, ninit, end, reps)
  res[[7]][[3]] = main(f.shekel, "uei", 2, A, ninit, end, reps)
  res[[7]][[4]] = main(f.shekel, "uei", 100, A, ninit, end, reps)
  res[[7]][[5]] = main(f.shekel, "dei", 1, A, ninit, end, reps)
  res[[7]][[6]] = main(f.shekel, "vei", -.5, A, ninit, end, reps)
  res[[7]][[7]] = main(f.shekel, "sei", 1, A, ninit, end, reps)
  res[[7]][[8]] = main(f.shekel, "pei", 100, A, ninit, end, reps)
}) 
system.time({
  A = matrix(c(rep(-2, 10), rep(2,10)), ncol = 2)
  res[[8]] = list()
  res[[8]][[1]] = main(f.rastrigin10 ,"uei", 0, A, ninit, end, reps)
  res[[8]][[2]] = main(f.rastrigin10, "uei", .5, A, ninit, end, reps)
  res[[8]][[3]] = main(f.rastrigin10, "uei", 2, A, ninit, end, reps)
  res[[8]][[4]] = main(f.rastrigin10, "uei", 100, A, ninit, end, reps)
  res[[8]][[5]] = main(f.rastrigin10, "dei", 1, A, ninit, end, reps)
  res[[8]][[6]] = main(f.rastrigin10, "vei", -.5, A, ninit, end, reps)
  res[[8]][[7]] = main(f.rastrigin10, "sei", 1, A, ninit, end, reps)
  res[[8]][[8]] = main(f.rastrigin10, "pei", 100, A, ninit, end, reps)
}) 


#saveRDS(res2, file = "res_Apr1_Mac.rds")
#setwd("/Users/jiajiekong/Dropbox/Herbie")
#saveRDS(res, file = "res_Mar24.rds")
#res_win = readRDS("/Users/jiajiekong/Dropbox/Herbie/res_Mar31_Thinkpad.rds")


#function plot
ContourPlot2D = function(f){
  par(mfrow = c(1,1))
  f.2 <- function(x,y,f=f){X = cbind(x,y);return(f(X))}
  nn <- 200
  x <- seq(-2, 2, length.out = nn)
  y <- seq(-2, 2, length.out = nn)
  z <- outer(x , y , f.2, f=f)
  min_ind <- which(z == min(z), arr.ind = TRUE)
  filled.contour(x, y, z,
                 plot.axes={axis(1, seq(-2.5, 2.5, by = 1));
                   axis(2, seq(-2.5, 2.5, by = 1));
                   points(x[min_ind[1]],y[min_ind[2]],pch = 19, col = "blue")},
                 main = "2D function")
}
ContourPlot2D_Track = function(f, result, rep_num){
  os = result$os[[rep_num]]
  prog = result$prog
  f.2 <- function(x,y,f=f){X = cbind(x,y);return(f(X))}
  nn <- 200
  x <- seq(-2, 2, length.out = nn)
  y <- seq(-2, 2, length.out = nn)
  z <- outer(x , y , f.2, f=f)
  m <- which.min(os$y)
  filled.contour(x, y, z, 
                 main = paste("EI, f_min=",round(mean(prog[,end]),3),sep = ""),
                 plot.axes={points(os$X[,1],os$X[,2],pch = 19, cex = .5);
                   points(os$X[m,1],os$X[m,2],pch = 19, col = "red", cex = .7);
                   points(os$X[end,1], os$X[end,2], pch = 19, col = "blue", cex = .7);
                   axis(1, seq(-2.5, 2.5, by = 1));
                   axis(2, seq(-2.5, 2.5, by = 1))},
                 xlim = c(-2.5, 2.5),
                 ylim = c(-2.5, 2.5))
}
ContourPlot2D(f.Rosenbrock)
ContourPlot2D(f.HumpCamel)
ContourPlot2D(f.Ackley)
ContourPlot2D(f.rastrigin2)
ContourPlot2D_Track(f.Rosenbrock, res[[2]][[1]], 1)
ContourPlot2D_Track(f.HumpCamel, res[[3]][[1]], 1)
ContourPlot2D_Track(f.Ackley, res[[4]][[1]], 1)
ContourPlot2D_Track(f.rastrigin2, res[[5]][[7]], 3)

#Comparison with CI
fminPlot_CI = function(result,  CI = TRUE){
  n = dim(result[[1]]$prog)[2]
  y = c()
  for (j in 1:8){
    y = c(y, colMeans(result[[j]]$prog[,3:end]))
  }
  plot(c(0, n), c(max(y), min(y)), type="n", 
         xlab="black-box evaluations (n)", ylab="average best objective value")
  legend("topright", c("EI", "UEI.05", "UEI.2","UEI100", "UEI.Dym", "VEI","SEI","PEI"), lwd = 2, 
           col = 1:8)
  for (i in 1:8){
    lines(colMeans(result[[i]]$prog), col=i, lwd=2)
    if (CI == TRUE){
    lines(apply(result[[i]]$prog, MARGIN = 2, quantile, prob = 0.05), col = i, lty =2)
    lines(apply(result[[i]]$prog, MARGIN = 2, quantile, prob = 0.95), col = i, lty =2)
    }
  }
}
fminPlot_CI(res[[1]], FALSE)
fminPlot_CI(res[[2]]) #VEI need recalculate
fminPlot_CI(res[[3]], FALSE)
fminPlot_CI(res[[4]]) #miss
fminPlot_CI(res[[5]], FALSE) 
fminPlot_CI(res[[6]]) #Hartmann 4D
#fminPlot_CI(res[[7]]) #try 2D
#fminPlot_CI(res[[8]], FALSE) 

#Converge
obj.best = function(result){
  plot(c(0, end), c(max(result$prog[,3:end]), min(result$prog)), type="n", 
       xlab="black-box evaluations (n)", ylab="best objective value")
  for (j in 1:reps){
    lines(result$prog[j,], col=3, type="l")
  }
}
obj.best(res[[5]][[7]])

















