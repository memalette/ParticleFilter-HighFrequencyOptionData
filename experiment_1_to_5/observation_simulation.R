# Initializing vectors

X1          <- matrix(rnorm(n=S$n.traj.obs*S$n.step.obs*S$n.days, mean=0, sd=1), 
                      nrow=S$n.traj.obs, ncol=S$n.step.obs*S$n.days)
X2          <- matrix(rnorm(n=S$n.traj.obs*S$n.step.obs*S$n.days, mean=0, sd=1), 
                      nrow=S$n.traj.obs, ncol=S$n.step.obs*S$n.days)

V.start      <- rep(P$v0, S$n.traj.obs)
Y.start      <- rep(P$y0, S$n.traj.obs)

W.v          <- X1
Z.v          <- matrix(nrow=S$n.traj.obs, ncol=S$n.step.obs*S$n.days)

V_           <- matrix(nrow=S$n.traj.obs, ncol=S$n.step.obs*S$n.days+1)
V            <- matrix(nrow=S$n.traj.obs, ncol=S$n.step.obs*S$n.days+1)
V.close      <- matrix(nrow=S$n.traj.obs, ncol=S$n.days+1)
V_[,1]       <- V.start
V[,1]        <- V.start
V.close[,1]  <- V.start
sum.jump.v   <- matrix(ncol=S$n.days,nrow=S$n.traj.obs)


W.y          <- P$rho * X1 + sqrt(1-P$rho^2) * X2
Z.y          <- matrix(nrow=S$n.traj.obs, ncol=S$n.step.obs*S$n.days)

Y            <- matrix(nrow=S$n.traj.obs, ncol=S$n.step.obs*S$n.days+1)
Y_           <- matrix(nrow=S$n.traj.obs, ncol=S$n.step.obs*S$n.days+1)
Y.close      <- matrix(nrow=S$n.traj.obs, ncol=S$n.days+1)
Y.close[,1]  <- Y.start
Y[,1]        <- Y.start
Y_[,1]       <- Y.start
sum.jump.y   <- matrix(ncol=S$n.days,nrow=S$n.traj.obs)

dI           <- matrix(ncol=S$n.step.obs*S$n.days,nrow=S$n.traj.obs)
idx          <- matrix(ncol=S$n.days,nrow=S$n.traj.obs)
quad.var     <- matrix(ncol=S$n.days,nrow=S$n.traj.obs)
int.var      <- matrix(ncol=S$n.days,nrow=S$n.traj.obs)

day <- 1 


# SIMULATION

for (j in 1:(S$n.days*S$n.step.obs)){   
  
  
  # Variance jumps
  # ----------------------------------------------------------------------- #
  

  # is there a jump in the variance for this step?
  prob.jump.v        <- 1 - exp(-P$lambda.v0 * S$h.obs) 
  unif.jump.v        <- runif(S$n.traj.obs,0,1) 
  
  # for experiment #2 and experiment #3
  mat               <- cbind(unif.jump.v, prob.jump.v)
  idx.T             <- which(mat[,1] < mat[,2])
  ind.jump.v        <- rep(0, S$n.traj.obs)
  ind.jump.v[idx.T] <- 1
  
  # which trajectory jumped for this step? which did not?
  idx.T              <- which(ind.jump.v == 1)
  idx.F              <- which(ind.jump.v == 0)
  
  # Jump size: for trajectories that jumped, generate a exponential 
  Z.v[idx.T,j]       <- rexp(n=length(idx.T),rate=1/P$mu.v)
  Z.v[idx.F,j]       <- 0 
  
  
  # Variance
  # ----------------------------------------------------------------------- #
  
  # pre-jump variance
  V_[,j+1]        <- V[,j] + P$kappa * (P$theta - V[,j]) * S$h.obs +
                     P$sigma * sqrt(V[,j]) * W.v[,j] *sqrt(S$h.obs)
  
  # If pre-jump variance is negative, set to 0
  idx.T           <- which(V_[,j+1] < exp(-12))
  V_[idx.T,j+1]   <- exp(-12)
  
  # add jump
  V[,j+1]         <- V_[,j+1] + Z.v[,j]
  
  # Integrated variance
  # ----------------------------------------------------------------------- #
  
  dI[,j]  <- S$h.obs * (V[,j] + V_[,j+1])/2
  
  # Log price jumps
  # ----------------------------------------------------------------------- #
  
    
  if ( j == 1){
    lambda.y <- P$lambda.y0
  }else {
    if ( j == 2){
      lambda.y <- P$lambda.y1}
    else{
      lambda.y <- P$lambda.y0 + P$lambda.y1 * V[,j]}
  }
  
  # is there a jump in the variance for this step?
  prob.jump.y     <- 1 - exp(-lambda.y * S$h.obs)
  unif.jump.y     <- runif(S$n.traj.obs,0,1)
  
  # If unif < jump prob, jump indicator = 1
  mat               <- cbind(unif.jump.y, prob.jump.y)
  idx.T             <- which(mat[,1] < mat[,2])
  ind.jump.y        <- rep(0, S$n.traj.obs)
  ind.jump.y[idx.T] <- 1
  
  # which trajectory jumped for this step? which did not?
  idx.T           <- which(ind.jump.y == 1)
  idx.F           <- which(ind.jump.y == 0)
  
  # Jump size: for trajectories that jumped, generate a exponential 
  Z.y[idx.T,j]    <- rnorm(n=length(idx.T), mean=P$mu.y, sd=P$sigma.y)
  Z.y[idx.F,j]    <- 0
  
  
  # Log price
  # ----------------------------------------------------------------------- #
  
  mgfy     <- exp(P$mu.y * 1 + 0.5 * P$sigma.y^2 * 1^2)
  alpha    <- P$r - P$q + (RP$eta.y - 0.5)*V[,j]  + 
              ( RP$gamma.y.min - (mgfy - 1))*lambda.y
  Y_[,j+1] <- Y[,j] +  alpha *S$h.obs + sqrt(V[,j])*W.y[,j] *sqrt(S$h.obs)
  Y[,j+1]  <- Y_[,j+1] + Z.y[,j]
  
  if( (j %% S$n.step.obs ==  0) == TRUE){
    V.close[,day+1]  <- V[,j+1]
    Y.close[,day+1]  <- Y[,j+1]
    int.var[,day]    <- apply(dI[,(j-S$n.step.obs+1):j], 1, sum)
    sum.jump.y[,day] <- apply(Z.y[,(j-S$n.step.obs+1):j]^2, 1, sum)
    sum.jump.v[,day] <- apply(Z.v[,(j-S$n.step.obs+1):j]^2, 1, sum)
    quad.var[,day]   <- int.var[,day] + sum.jump.y[,day]
    day              <- day + 1
  }
  
  print(paste0("iteration: ",j))
  
}

# # saving in obs object
obs             <-list()
obs$Y           <- Y[1,]
obs$Y_          <- Y_[1,]
obs$V           <- V[1,]
obs$V_          <- V_[1,]
obs$Y.close     <- Y.close[1,]
obs$V.close     <- V.close[1,]
obs$int.var     <- int.var[1,]
obs$sum.jump.y  <- sum.jump.y[1,]
obs$sum.jump.v  <- sum.jump.v[1,]
obs$Z.y         <- Z.y[1,]
obs$Z.v         <- Z.v[1,]
obs$quad.var    <- quad.var[1,]
obs$dI          <- dI[1,]

par(mfrow=c(3,2))
plot(obs$Y, type="l", main = "Log-price", xlab = "step", ylab="")
plot(obs$V, type="l", main = "Variance", xlab = "step", ylab="")
plot(obs$quad.var, type="l", main = "Quadratic variation", xlab = "day", ylab="")
plot(obs$int.var, type="l", main = "Integrated variance", xlab = "day", ylab="")
plot(obs$Z.y, type="l", main = "Jumps Y", xlab = "step", ylab="")
plot(obs$Z.v, type="l", main = "Jumps V", xlab = "step", ylab="")

imgname   <- paste0("observations", ".png")
dev.copy(png, imgname)
dev.off()
  
