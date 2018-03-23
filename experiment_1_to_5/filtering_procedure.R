# --------------------------------------------------------------------------- #
# FILTETRING PROCEDURE
# --------------------------------------------------------------------------- #

# Initializing vectors
# Monte Carlo simulation starting values
V.start      <- rep(P$v0, S$n.traj)
Y.start      <- rep(P$y0, S$n.traj)

# Generate Brownian motions and V-jumps outside of loop (less costly)
W.v          <- matrix(rnorm(n=S$n.traj*S$n.step*S$n.days, mean=0, sd=1), 
                       nrow=S$n.traj, ncol=S$n.step*S$n.days)
Z.v          <- matrix(nrow=S$n.traj, ncol=S$n.step*S$n.days)

# Initialize Variance related vectors and matrices
V_           <- matrix(nrow=S$n.traj, ncol=S$n.step*S$n.days+1) # Var prior to jump (contains intraday values)
V            <- matrix(nrow=S$n.traj, ncol=S$n.step*S$n.days+1) # Var (containing intraday values)
V.close      <- matrix(nrow=S$n.traj, ncol=S$n.days+1) # end of day Var
V_[,1]       <- V.start # First column is set to starting var value
V[,1]        <- V.start # First column is set to starting var value
V.close[,1]  <- V.start # First column is set to starting var value
sum.jump.v   <- matrix(ncol=S$n.days,nrow=S$n.traj) # cumulative sum of V-jumps at the end of the day

# Generate Brownian motions and Y-jumps outside of loop (less costly)
W.y          <- matrix(rnorm(n=S$n.traj*S$n.step*S$n.days, mean=0, sd=1), 
                       nrow=S$n.traj, ncol=S$n.step*S$n.days)
Z.y          <- matrix(nrow=S$n.traj, ncol=S$n.step*S$n.days)

# Initialize Log-price related vectors and matrices
Y            <- matrix(nrow=S$n.traj, ncol=S$n.step*S$n.days+1)
Y_           <- matrix(nrow=S$n.traj, ncol=S$n.step*S$n.days+1)
Y.close      <- matrix(nrow=S$n.traj, ncol=S$n.days+1)
Y.close[,1]  <- Y.start
Y[,1]        <- Y.start
Y_[,1]       <- Y.start
sum.jump.y   <- matrix(ncol=S$n.days,nrow=S$n.traj)

# Initialize vectors and matrices related to integrated variance and quadratic variation
m            <- matrix(ncol=S$n.step*S$n.days,nrow=S$n.traj)
s            <- matrix(ncol=S$n.step*S$n.days,nrow=S$n.traj)
BK.approx    <- matrix(ncol=S$n.step*S$n.days,nrow=S$n.traj)
dI           <- matrix(ncol=S$n.step*S$n.days,nrow=S$n.traj)
idx          <- vector(length = S$n.traj)
quad.var     <- matrix(ncol=S$n.days,nrow=S$n.traj)
int.var      <- matrix(ncol=S$n.days,nrow=S$n.traj)
C1           <- (1 + exp(-P$kappa*S$h))/(1 - exp(-P$kappa*S$h)) # constants for int.var moments generation
C2           <- 4*exp(-P$kappa*S$h)/((1 - exp(-P$kappa*S$h))^2)

# Matrix that saves resampled var values at the end of each day 
V.SIR        <- matrix(nrow=S$n.traj, ncol=S$n.days)


# results
RMSE         <- data.frame(row.names = c("V", "QV", "I", "Z.y", "Z.v"))

# end of day iterator
day <-1 


# ========================================================================= #
#                       SIMULATION AND FILTER
# ========================================================================= #


for (j in 1:(S$n.days*S$n.step)){   
  
  
            # ================================================#
            #                 1 -  SIMULATION
            # =============================================== #
  
  
  # Variance jumps
  # ----------------------------------------------------------------------- #
  
  
  # is there a jump in the variance for this step?
  prob.jump.v        <- 1 - exp(-P$lambda.v0 * S$h) 
  unif.jump.v        <- runif(S$n.traj,0,1) 
  ind.jump.v         <- ifelse(unif.jump.v < prob.jump.v, 1, 0) 
  
  
  # which trajectory jumped for this step? which did not?
  idx.T              <- which(ind.jump.v == 1)
  idx.F              <- which(ind.jump.v == 0)
  
  
  # Jump size: for trajectories that jumped, generate a exponential 
  Z.v[idx.T,j]       <- rexp(n=length(idx.T),rate=1/P$mu.v)
  Z.v[idx.F,j]       <- 0 
  
  
  # Variance
  # ----------------------------------------------------------------------- #
  
  
  # pre-jump variance
  V_[,j+1]        <- V[,j] + P$kappa * (P$theta - V[,j]) * S$h +
                     P$sigma * sqrt(V[,j] * S$h) * W.v[,j] 
  
  # If pre-jump variance is negative, set to 0
  idx.T           <- which(V_[,j+1] < P$vmin)
  V_[idx.T,j+1]   <- P$vmin
  
  # add jump
  V[,j+1]         <- V_[,j+1] + Z.v[,j]
  
  V1 <- V[,j]
  V2 <- V_[,j+1]
  
  
  # Integrated variance
  # ----------------------------------------------------------------------- #
  
  averI <- interp1(P$cache.x,P$moments.grid[,1],sqrt(V1*V2)) + 
                    (V1+ V2)*(C1/P$kappa - 0.5*S$h*C2)
  
  varI  <- interp1(P$cache.x,P$moments.grid[,2],sqrt(V1*V2)) + 
                  (V1 + V2)*(P$sigma^2*C1/P$kappa^3 + 
                   P$sigma^2*S$h*C2/(2*P$kappa^2)-
                   P$sigma^2*S$h^2*C1*C2/(2*P$kappa))
  
  dI[,j]  <- rnorm(n=S$n.traj, mean = averI, sd = sqrt(varI))

  
  # If dI is negative, set to 0 - happens a lot
  idx.T           <- which(dI[,j] < 0 | is.na(dI[,j]))
  dI[idx.T,j]     <- 0
  
  BK.approx[,j]   <- rnorm(n=S$n.traj, mean = 0, sd = sqrt(dI[,j]))

  
  # Log price jumps
  # ----------------------------------------------------------------------- #
  
  
  if ( j == 1){ 
    lambda.y <- P$lambda.y0
  }else {   
    if ( j == 2){ 
      lambda.y <- P$lambda.y1}
    else{
      lambda.y <- P$lambda.y0 + P$lambda.y1 * V_[,j]}
  }
  
  
  # is there a jump in the variance for this step?
  prob.jump.y     <- 1 - exp(-lambda.y * S$h)
  unif.jump.y     <- runif(S$n.traj,0,1)
  ind.jump.y      <- ifelse(unif.jump.y < prob.jump.v, 1, 0)
  
  
  # which trajectory jumped for this step? which did not?
  idx.T           <- which(ind.jump.y == 1)
  idx.F           <- which(ind.jump.y == 0)
  
  
  # Jump size: for trajectories that jumped, generate a exponential 
  Z.y[idx.T,j]    <- rnorm(n=length(idx.T), mean=P$mu.y, sd=P$sigma.y)
  Z.y[idx.F,j]    <- 0
  
  # Log price
  # ----------------------------------------------------------------------- #
  
  Y_[,j+1] <- Y[,j] + P$c1 * S$h + 
              P$c2 * dI[,j]   + 
              (P$rho / P$sigma) * (V_[,j+1] - V[,j]) +
              sqrt(1 - P$rho^2) * BK.approx[,j] -
              (P$rho / P$sigma) * Z.v[,j] 
  
  Y[,j+1]  <- Y_[,j+1] + Z.y[,j]
  
  
            # ================================================#
            #                   2 - FILTER
            # =============================================== #
  # Theoretical measures are computed at the end of the day
  # This part is skipped by the loop if the step is intraday: j %% S$n.days != 0
  
  if( (j %% S$n.step ==  0) == TRUE){
    
    if(day != S$n.days){ 

      # Closing price and variance to calculate the return
      V.close[,day+1]  <- V[,j+1]
      Y.close[,day+1]  <- Y[,j+1]
      
      
      #Aggregated measures
      if(S$n.step == 1){
        # if only one step a day
        int.var[,day]    <- dI[,j]
        sum.jump.y[,day] <- Z.y[,j]^2
        sum.jump.v[,day] <- Z.v[,j]^2
      }else{ 
        # if more than one step a day
        int.var[,day]    <- apply(dI[,(j-S$n.step+1):j], 1, sum)
        sum.jump.y[,day] <- apply(Z.y[,(j-S$n.step+1):j]^2, 1, sum)
        sum.jump.v[,day] <- apply(Z.v[,(j-S$n.step+1):j]^2, 1, sum)
      }
      
      quad.var[,day]   <- int.var[,day] + sum.jump.y[,day]
      
          
                # ===================================== #
                #      2a - INDEX-RELATED OBSERVABLES
                # ===================================== #
      

    # Log price - weights
    # --------------------------------------------------------------------- #
      
      if(S$include$y == TRUE){
        return.obs   <- rep(obs$Y.close[day+1] - obs$Y.close[day], S$n.traj)
        return.sim   <- Y.close[,day+1] - Y.close[,day]
        var          <- (1-P$rho^2) * int.var[,day]
        z1           <- (return.obs - return.sim) / sqrt(var)
        weights.y    <- dnorm(x = z1)
      } else{
        weights.y    <- rep(1, S$n.traj)
      }
      
      
    # Realized variance  - weights
    # --------------------------------------------------------------------- #
      
      if(S$include$rv == TRUE){ 
        z2           <- (quad.var[,day] - obs$RV[day])/obs$RV[day]
        weights.rv   <- dnorm(x = z2,mean = 0,sd = S$sde$z2)
      } else{
        weights.rv   <- rep(1, S$n.traj)
      }
      
      
    # Bipower variation - weights
    # --------------------------------------------------------------------- # 
      
      if(S$include$bv == TRUE){ 
        z3           <- (int.var[,day] - obs$BV[day])/obs$BV[day]
        weights.bv   <- dnorm(x = z3,mean = 0,sd = S$sde$z3)
      } else{
        weights.bv   <- rep(1, S$n.traj)
      }
      
      
              # ======================================== #
              #      2b - OPTIONS-RELATED OBSERVABLES
              # ======================================== #
      
    # CASE #1: Filter set Y, RV, BV, IV
    # ==================================
    if(S$include$iv == TRUE & S$include$rov == FALSE & S$include$rivv == FALSE){
      
      # Get moneynesses of simulated paths with respect to observed strike prices K
      K           <- obs$K[day+1,]
      moneyness30 <- exp(rep(Y[,j+1], S$ld))/K[1:S$ld]
      moneyness90 <- exp(rep(Y[,j+1], S$ld))/K[(S$ld+1):(2*S$ld)]

      # Get the log variance of simulated paths
      lm          <- length(moneyness30)
      logV30      <- log(rep(V[,j+1], each= S$ld))
      logV90      <- log(rep(V[,j+1], each= S$ld))
      
      # define loc for interpolation of implied volatilities
      loc30       <- cbind( logV30, moneyness30)
      loc90       <- cbind( logV90, moneyness90)
    
      # interpolate IVs: S$obj.30.iv and  S$obj.90.iv contains an grid of implied volatilities 
      # log-variance and moneyness as the 2 defining factors of that implied volatility
      IV.sim30    <- matrix(interp.surface( S$obj.30.iv, loc30), 
                            ncol=S$ld, byrow = TRUE)
      IV.sim90    <- matrix(interp.surface( S$obj.90.iv, loc90), 
                            ncol=S$ld, byrow = TRUE)
      IV.sim      <- cbind(IV.sim30, IV.sim90)
      
    } 
    
      
    # CASE #2, #3, #4: Filter set Y, RV, BV, IV, ROV or 
    #                  Filter set Y, RV, BV, IV, RIVV or
    #                  Filter set Y, RV, BV, IV, ROV, RIVV
    # ========================================================================
    if(S$include$rov == TRUE | S$include$rivv == TRUE){
      
      # dimensions of data sets
      dim1        <- c(S$n.traj, S$n.step+1, S$ld)
      dim2        <- c(S$n.traj, S$n.step+1, 2*S$ld)

      # strike, underlying price, moneyness and log variance for interpolation
      K           <- array(rep(obs$K[day+1,],each=S$n.traj*(S$n.step+1)), 
                           dim=dim2)
      S0          <- exp(array(rep(Y[,(j-S$n.step+1):(j+1)], S$ld),dim=dim1))
      moneyness30 <- S0/K[,,1:S$ld]
      moneyness90 <- S0/K[,,(S$ld+1):(2*S$ld)]
      
      # to avoid NAs
      moneyness30[which(moneyness30 < 0.31)] <- 0.31
      moneyness30[which(moneyness30 > 1.99)] <- 1.99
      moneyness90[which(moneyness90 < 0.31)] <- 0.31
      moneyness90[which(moneyness90 > 1.99)] <- 1.99
      
      
      lm          <- length(moneyness30)
      V0          <- array(rep(V[,(j-S$n.step+1):(j+1)], S$ld),dim=dim1)
      logV        <- log(V0)
      
      # define loc for interpolation of implied volatilities
      loc30       <- cbind( logV, moneyness30)
      loc90       <- cbind( logV, moneyness90)
      
      # interpolate IVs, deltas and vegas for option pricing and option quadratic variaiton computation
      # log-variance and moneyness as the 2 defining factors of that implied volatility
      # I did this to make the option pricing faster (faster to interpolate then price options individually)
      IV.sim30    <- array(interp.surface( S$obj.30.iv, loc30), dim=dim1)
      IV.sim90    <- array(interp.surface( S$obj.90.iv, loc90), dim=dim1)
      
      delta30     <- array(interp.surface( S$obj.30.d, loc30), dim=dim1) * S0
      delta90     <- array(interp.surface( S$obj.90.d, loc90), dim=dim1) * S0
      
      vega30     <- array(interp.surface( S$obj.30.v, loc30), dim=dim1) * S0
      vega90     <- array(interp.surface( S$obj.90.v, loc90), dim=dim1) * S0
      
      # Combining results 
      IV.sim     <- abind(IV.sim30, IV.sim90, along=3)
      delta      <- abind(delta30, delta90, along=3)
      vega       <- abind(vega30, vega90, along=3)
      moneyness  <- abind(moneyness30, moneyness90, along=3)

      # convert IVs in price
      t        <- array(rep(S$tau, each=lm),dim=dim2)
      
      price    <- GetBlackScholesCall(S=abind(S0,S0,along=3), 
                                      K=K, 
                                      t=t, 
                                      r=P$r, 
                                      q=P$q, 
                                      sigma=IV.sim, 
                                      is.call=1)
      
      # in order to compute the OQV, we need the option prices prior 
      # to a jump in variance or log-price
        S0_          <- exp(array(rep(Y_[,(j-S$n.step+1):(j+1)], S$ld),
                            dim=dim1))
        moneyness30_ <- S0_/K[,,1:S$ld]
        moneyness90_ <- S0_/K[,,(S$ld+1):(2*S$ld)]
        
        # to avoid NAs
        moneyness30_[which(moneyness30_ < 0.31)] <- 0.31
        moneyness30_[which(moneyness30_ > 1.99)] <- 1.99
        moneyness90_[which(moneyness90_ < 0.31)] <- 0.31
        moneyness90_[which(moneyness90_ > 1.99)] <- 1.99
        
        lm           <- length(moneyness30_)
        V0_          <- array(rep(V_[,(j-S$n.step+1):(j+1)], S$ld),dim=dim1)
        logV_        <- log(V0_)
        
        loc30_ <- cbind( logV_, moneyness30_)
        loc90_ <- cbind( logV_, moneyness90_)
        
        # IV interpolation
        IV.sim30_ <- array(interp.surface( S$obj.30.iv, loc30_), dim=dim1)
        IV.sim90_ <- array(interp.surface( S$obj.90.iv, loc90_), dim=dim1)
        IV.sim_   <- abind(IV.sim30_, IV.sim90_, along=3)
        
        # convert IVs in option price
        price_    <- GetBlackScholesCall(S=abind(S0_,S0_,along=3), 
                                        K=K, 
                                        t=t, 
                                        r=P$r, 
                                        q=P$q, 
                                        sigma=IV.sim_, 
                                        is.call=1)
        


    }
      
    
    # implied volatility - weights
    # --------------------------------------------------------------------- #
      
    if(S$include$iv == TRUE){ 
      
        IV           <- matrix(rep(obs$IV.close[day+1,], S$n.traj), 
                           ncol=ncol(obs$IV.close), 
                           byrow=TRUE)
        
        if(S$include$rov == TRUE  | S$include$rivv == TRUE  ){
          IV.sim.eod <- IV.sim[,S$n.step+1,]
        }else{
          IV.sim.eod <- IV.sim
        }
        
        z4           <- (IV.sim.eod - IV)/IV
        weights.iv   <- dnorm(x = z4, mean = 0, sd = S$sde$z4)
        weights.iv   <- apply(X=weights.iv,MARGIN=1,FUN=prod)^(1/S$n.options)
        
      } else{
        weights.iv   <- rep(1, S$n.traj)
      }
      
      
    # Option realized variance  - weights
    # --------------------------------------------------------------------- #
      
    if(S$include$rov == TRUE ){
      
      if(S$n.step != 1){
        indexdI <- (j-S$n.step+1):j
        indexgr <- 1:S$n.step
        indexpr <- 2:(S$n.step+1)
        
        dI.a        <- array(data=rep(dI[,indexdI], S$n.options), 
                             dim=c(S$n.traj,S$n.step,S$n.options))
        
        mat1        <- (delta[,indexgr,]^2 +
                          (P$sigma * vega[,indexgr,])^2 +
                          2 * P$sigma * P$rho * delta[,indexgr,] * 
                          vega[,indexgr,]) * dI.a[,indexgr,] 
        
        term1       <- foreach(i=1:S$n.options,.combine = "cbind") %do%
                       apply(mat1[,,i],1,sum)
        
        mat2        <- (price[,indexpr,] -  price_[,indexpr,])^2
        
        term2       <- foreach(i=1:S$n.options, .combine = "cbind") %do%
                       apply(mat2[,,i], 1, sum)
        
      }else{ 
        indexdI <- j
        indexgr <- 1
        indexpr <- 2
        
        dI.a        <- array(data=rep(dI[,indexdI], S$n.options), 
                             dim=c(S$n.traj,S$n.step,S$n.options))
        
        term1        <- (delta[,indexgr,]^2 +
                          (P$sigma * vega[,indexgr,])^2 +
                          2 * P$sigma * P$rho * delta[,indexgr,] * 
                          vega[,indexgr,]) * dI.a[,indexgr,] 
        
        # term1       <- foreach(i=1:S$n.options,.combine = "cbind") %do%
        #   apply(mat1[,,i],1,sum)
        
        term2        <- (price[,indexpr,] -  price_[,indexpr,])^2
        
        # term2       <- foreach(i=1:S$n.options, .combine = "cbind") %do%
        #   apply(mat2[,,i], 1, sum)
        }
      
        # Compute OQV

        
        OQV         <- (term1 + term2)

        ROV         <- matrix(rep(obs$ROV[day,], S$n.traj), 
                        ncol=ncol(obs$ROV), byrow=TRUE)
        
        z5          <- (OQV - ROV)/ROV
        weights.rov <- dnorm(x = z5, mean = 0, sd = S$sde$z5)
        weights.rov <- apply(X=weights.rov,MARGIN=1,FUN=prod)^(1/S$n.options)

      } else {
        weights.rov <- rep(1, S$n.traj)
      }
     
    # Implied Volatility Realized Variance - weights
    # --------------------------------------------------------------------- #
      
    if(S$include$rivv == TRUE){
      
    
      # Computing IV-greeks using finite differences method -------------------
      
        S0.p.epsy   <- exp(log(S0)+S$epsy)
        S0.m.epsy   <- exp(log(S0)-S$epsy)
        logV.p      <- log(V0+S$epsv)
        V.m         <- V0-S$epsv
        V.m[which(V.m < P$vmin)] <- P$vmin
        logV.m      <- log(V.m)
      
        # dIV.y 30-days options 
        mon.p.epsy30 <- S0.p.epsy/K[,,1:S$ld]
        mon.m.epsy30 <- S0.m.epsy/K[,,1:S$ld]
        loc.p.epsy30 <- cbind( as.vector(logV), as.vector(mon.p.epsy30))
        loc.m.epsy30 <- cbind( as.vector(logV), as.vector(mon.m.epsy30))
        IV.p.epsy30  <- array(interp.surface( S$obj.30.iv, loc.p.epsy30), 
                              dim=dim1)
        IV.m.epsy30  <- array(interp.surface( S$obj.30.iv, loc.m.epsy30), 
                              dim=dim1)
        
        
        # dIV.y 90-days options 
        mon.p.epsy90 <- S0.p.epsy/K[,,(S$ld+1):(2*S$ld)]
        mon.m.epsy90 <- S0.m.epsy/K[,,(S$ld+1):(2*S$ld)]
        loc.p.epsy90 <- cbind( as.vector(logV), as.vector(mon.p.epsy90))
        loc.m.epsy90 <- cbind( as.vector(logV), as.vector(mon.m.epsy90))
        IV.p.epsy90  <- array(interp.surface( S$obj.90.iv, loc.p.epsy90), 
                                    dim=dim1)
        IV.m.epsy90  <- array(interp.surface( S$obj.90.iv, loc.m.epsy90), 
                                    dim=dim1)
        
        
        # dIV.y 
        dIV.y30      <- (IV.p.epsy30-IV.m.epsy30)/(2*S$epsy)
        dIV.y90      <- (IV.p.epsy90-IV.m.epsy90)/(2*S$epsy)
        
   
        
        # dIV.v 30-days options
        loc.p.epsv  <- cbind( as.vector(logV.p), as.vector(moneyness))
        loc.m.epsv  <- cbind( as.vector(logV.m), as.vector(moneyness))
        IV.p.epsv30 <- array(interp.surface( S$obj.30.iv, loc.p.epsv), 
                             dim=dim1)
        IV.m.epsv30 <- array(interp.surface( S$obj.30.iv, loc.m.epsv), 
                             dim=dim1)
        
        
        # dIV.v 90-days options
        IV.p.epsv90 <- array(interp.surface( S$obj.90.iv, loc.p.epsv), 
                             dim=dim1)
        IV.m.epsv90 <- array(interp.surface( S$obj.90.iv, loc.m.epsv), 
                             dim=dim1)
        dIV.v30     <- (IV.p.epsv30-IV.m.epsv30)/(2*S$epsv) 
        dIV.v90     <- (IV.p.epsv90-IV.m.epsv90)/(2*S$epsv) 
    
        # dIV.v
        dIV.y       <- abind(dIV.y30, dIV.y90, along=3)
        dIV.v       <- abind(dIV.v30, dIV.v90, along=3)
        

        # Computing IVQV 
        
        if(S$n.step != 1){
          indexdI <- (j-S$n.step+1):j
          indexgr <- 1:S$n.step
          indexpr <- 2:(S$n.step+1)
          
          dI.a       <- array(data=rep(dI[,indexdI], S$n.options), 
                              dim=c(S$n.traj,S$n.step,S$n.options))
          
          mat1       <- (dIV.y[,indexgr,]^2 + 
                           P$sigma^2 * dIV.v[,indexgr,]^2 +
                           2 * P$rho * P$sigma * dIV.v[,indexgr,] * 
                           dIV.y[,indexgr,]) * dI.a[,indexgr,]
          
          IV.jump    <- (IV.sim[,indexpr,] - IV.sim_[,indexpr,])^2
          
          
          term1      <- foreach(i=1:S$n.options, .combine = "cbind") %do%
                        apply(mat1[,,i], 1, sum)
          
          term2      <- foreach(i=1:S$n.options, .combine = "cbind") %do%
                        apply(IV.jump[,,i], 1, sum)
          
        }else{ 
          indexdI <- j
          indexgr <- 1
          indexpr <- 2
          
          dI.a       <- array(data=rep(dI[,indexdI], S$n.options), 
                              dim=c(S$n.traj,S$n.step,S$n.options))
          
         term1       <- (dIV.y[,indexgr,]^2 + 
                           P$sigma^2 * dIV.v[,indexgr,]^2 +
                           2 * P$rho * P$sigma * dIV.v[,indexgr,] * 
                           dIV.y[,indexgr,]) * dI.a[,indexgr,]
          
          term2    <- (IV.sim[,indexpr,] - IV.sim_[,indexpr,])^2
          
        
          }
        

        
        IVQV       <- term1 + term2
      
        RIVV       <- matrix(rep(as.numeric(obs$RIVV[day,]), S$n.traj),
                                                    ncol=ncol(obs$RIVV),
                                                    byrow=TRUE)
        
        z6           <- (IVQV - RIVV)/RIVV
        weights.rivv <- dnorm(x = z6, mean = 0, sd = S$sde$z6)
        weights.rivv <- apply(X=weights.rivv,MARGIN=1,FUN=prod)^(1/10)

      } else{
        weights.rivv <- rep(1, S$n.traj)
      }
      
      
                # ===================================== #
                #              2c - RESAMPLING
                # ===================================== #
      

    weights    <- weights.y * weights.rv * weights.bv * 
                  weights.iv * weights.rov * weights.rivv * 1/S$n.traj
    weights    <- weights/sum(weights)
    
    idx.na     <- which(is.na(weights)) 
    
    if(length(idx.na)>0){ 
    print(paste0("NAs in weights:", length(idx.na), "/", S$n.traj))}
    
    weights[idx.na] <- 0

    idx        <- resample(weights = weights,
                          num.samples = S$n.traj,
                          method = "multinomial",
                          engine = "C",
                          normalized = TRUE)$indices

      
    # Replace old paths with resampled ones 
    V           <- V[idx,]
    V_          <- V_[idx,]
    Y           <- Y[idx,]
    V.close     <- V.close[idx,]
    V.SIR[,day] <- V.close[,day+1]
    Y.close     <- Y.close[idx,]
    int.var     <- int.var[idx,]
    quad.var    <- quad.var[idx,]
    sum.jump.y  <- sum.jump.y[idx,]
    sum.jump.v  <- sum.jump.v[idx,]
    
                # ===================================== #
                #          3 -  RMSE Computation
                # ===================================== #
    
 
    RMSE["V",day]   <- sqrt(mean( (V.SIR[,day] - obs$V.close[day+1])^2 )   )
    RMSE["QV",day]  <- sqrt(mean( (quad.var[,day] - obs$quad.var[day])^2 ) )
    RMSE["I",day]   <- sqrt(mean( (int.var[,day] - obs$int.var[day])^2 ) )   
    RMSE["Z.y",day] <- sqrt(mean( (sum.jump.y[,day] - obs$sum.jump.y[day])^2 ) )   
    RMSE["Z.v",day] <- sqrt(mean( (sum.jump.v[,day] - obs$sum.jump.v[day])^2 ) )
    
    day <- day + 1
    
  } else {
    
    
                # =============================================== #
                #        4 - CLOSING SIMULATION (FINAL LOOP)
                # =============================================== #
    
    
    # No filtering, just save end of day measures
    V.close[,day+1]  <- V[,j+1]
    Y.close[,day+1]  <- Y[,j+1]
    V.SIR[,day]      <- V.close[,day+1]
    
    #Aggregated measures
    if(S$n.step == 1){
      int.var[,day]    <- dI[,j]
      sum.jump.y[,day] <- Z.y[,j]^2
      sum.jump.v[,day] <- Z.v[,j]^2
    }else{ 
      int.var[,day]    <- apply(dI[,(j-S$n.step+1):j], 1, sum)
      sum.jump.y[,day] <- apply(Z.y[,(j-S$n.step+1):j]^2, 1, sum)
      sum.jump.v[,day] <- apply(Z.v[,(j-S$n.step+1):j]^2, 1, sum)
    }
    quad.var[,day]   <- int.var[,day] + sum.jump.y[,day]
    

                  # ======================================= #
                  #   4a -  RMSE Computation (Final loop) 
                  # ======================================= #
    
    
    RMSE["V",day]   <- sqrt(mean( (V.SIR[,day] - obs$V.close[day+1])^2 )   )
    RMSE["QV",day]  <- sqrt(mean( (quad.var[,day] - obs$quad.var[day])^2 ) )
    RMSE["I",day]   <- sqrt(mean( (int.var[,day] - obs$int.var[day])^2 ) )   
    RMSE["Z.y",day] <- sqrt(mean( (sum.jump.y[,day] - obs$sum.jump.y[day])^2 ) )   
    RMSE["Z.v",day] <- sqrt(mean( (sum.jump.v[,day] - obs$sum.jump.v[day])^2 ) )
    
    } 
    

  } 

  print(paste0("iteration: ", j))
    
} 


