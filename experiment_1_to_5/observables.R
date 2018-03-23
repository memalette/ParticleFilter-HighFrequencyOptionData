# ========================================================================= #
#                              OBSERVABLES
# ========================================================================= #

            # ================================================#
            #                 INDEX-RELATED
            # =============================================== #

# Compute Realized variance and bipower variation
Y.mat  <- matrix(obs$Y[-1], nrow=S$n.step.obs)
obs$RV <- apply(Y.mat, 2, GetRVSimulation)
obs$BV <- apply(Y.mat, 2, GetBVSimulation)

            # ================================================#
            #                 OPTION-RELATED
            # =============================================== #

# Compute strike panel from Y, V time series
lv        <- length(obs$V.close)
moneyness <- matrix(nrow = lv, ncol = S$n.options)

for(m in 1:length(obs$V.close)){
  
  tau   <- rep(S$tau, each=5)
  delta <- rep(S$delta, 2)
  
  for(n in 1:S$n.options){
    
    fun             <- function(x) GetDelta(x, V=obs$V.close[m], t=tau[n],
                                            r=P$r, q=P$q) - delta[n]
    moneyness[m,n]  <- uniroot(f = fun, interval = c(0,2))$root
    
    
  }
  
  print(m)
  
}

obs$moneyness <- moneyness
obs$K         <- t(matrix(rep(exp(obs$Y.close), S$n.options),
                   nrow=S$n.options,byrow=TRUE)) / obs$moneyness

# plots saving
par(mfrow=c(1,1))
matplot(obs$K, type="l", main = "Strikes", xlab = "day", ylab="")
imgname   <- paste0("strikes", ".png")
dev.copy(png, imgname)
dev.off()


# Get all inputs for grid interpolation
S0          <- t(matrix(rep(exp(obs$Y), S$ld),
                        nrow=S$ld,byrow=TRUE))
SS          <- cbind(S0,S0)
K30         <- rbind(obs$K[1,1:S$ld], 
                     matrix(apply(obs$K[-1,1:S$ld],1,rep,S$n.step.obs),
                      ncol=S$ld,byrow = TRUE))
K90         <- rbind(obs$K[1,(S$ld+1):(2*S$ld)],
                     matrix(apply(obs$K[-1,(S$ld+1):(2*S$ld)],1,rep,S$n.step.obs),
                      ncol=S$ld, byrow = TRUE))
K           <- cbind(K30,K90)
moneyness30 <- S0/K30
moneyness90 <- S0/K90
moneyness   <- cbind(moneyness30,moneyness90)
V0          <- matrix(rep(obs$V, each = S$ld),ncol= S$ld,byrow = TRUE)
logV        <- log(V0)

# location of observations
loc30 <- cbind( as.vector(logV), as.vector(moneyness30))
loc90 <- cbind( as.vector(logV), as.vector(moneyness90))

# implied vol interpolation
IV30 <- matrix(interp.surface( S$obj.30.iv, loc30), ncol=S$ld)
IV90 <- matrix(interp.surface( S$obj.90.iv, loc90), ncol=S$ld)

# delta interpolation
delta30 <- matrix(interp.surface( S$obj.30.d, loc30), ncol=S$ld)*S0
delta90 <- matrix(interp.surface( S$obj.90.d, loc90), ncol=S$ld)*S0

# vega interpolation
vega30 <- matrix(interp.surface( S$obj.30.v, loc30), ncol=S$ld)*S0
vega90 <- matrix(interp.surface( S$obj.90.v, loc90), ncol=S$ld)*S0

# bind 2 maturities together
IV     <- cbind(IV30, IV90)
delta  <- cbind(delta30, delta90)
vega   <- cbind(vega30, vega90)

# Get all inputs for grid interpolation (for the price before jumps)

# inputs
S0_          <- t(matrix(rep(exp(obs$Y_), S$ld),
                        nrow=S$ld,byrow=TRUE))
SS_          <- cbind(S0_,S0_)
moneyness30_ <- S0_/K30
moneyness90_ <- S0_/K90
V0_          <- matrix(rep(obs$V_, each = S$ld),ncol= S$ld,byrow = TRUE)
logV_        <- log(V0_)

# location of observations
loc30_ <- cbind( as.vector(logV_), as.vector(moneyness30_))
loc90_ <- cbind( as.vector(logV_), as.vector(moneyness90_))

# implied vol interpolation
IV30_ <- matrix(interp.surface( S$obj.30.iv, loc30_), ncol=S$ld)
IV90_ <- matrix(interp.surface( S$obj.90.iv, loc90_), ncol=S$ld)
IV_   <- cbind(IV30_, IV90_)

# conversion of implied vol to price
lm       <- length(moneyness30)
t30      <- matrix(rep(S$tau[1], lm), ncol=S$ld)
t90      <- matrix(rep(S$tau[2], lm), ncol=S$ld)
t        <- cbind(t30,t90)
prices   <- GetBlackScholesCall(S=SS,
                                K=K ,
                                t=t,
                                r=P$r,
                                q=P$q,
                                sigma=IV,
                                is.call=1)

prices_    <- GetBlackScholesCall(S=SS_,
                                  K=K,
                                  t=t,
                                  r=P$r,
                                  q=P$q,
                                  sigma=IV_,
                                  is.call=1)


# saving results in obs list
obs$prices    <- prices
obs$prices_   <- prices_
obs$delta     <- delta
obs$vega      <- vega
obs$IV        <- IV
obs$IV_       <- IV_
eod.idx       <- seq(from=1, 
                     to=(S$n.days*S$n.step.obs+1), 
                     by=S$n.step.obs) # end of day index
obs$IV.close  <- obs$IV[eod.idx,]

        # ================================================#
        #            Realized option variance (ROV)
        # =============================================== #

obs$OQV <- matrix(nrow=S$n.days,ncol=S$n.options)

for(day in 1:S$n.days){
  
  i         <- day*S$n.step.obs
  
  dI.mat    <- matrix(rep(obs$dI[(i-S$n.step.obs+1):i],S$ld*S$lt),
                      nrow=S$n.step.obs, ncol=S$n.options)
  
  mat1      <- (obs$delta[(i-S$n.step.obs+1):(i),]^2 +
                  (P$sigma * obs$vega[(i-S$n.step.obs+1):(i),]^2 ) +
                  2 * P$sigma * P$rho * obs$delta[(i-S$n.step.obs+1):(i),] *
                  obs$vega[(i-S$n.step.obs+1):(i),] )  * dI.mat
  
  term1     <- apply(mat1 ,2,sum)
  
  mat2      <- (obs$prices[(i-S$n.step.obs+2):(i+1),] -
                  obs$prices_[(i-S$n.step.obs+2):(i+1),])^2
  
  term2     <- apply(mat2, 2, sum)
  
  obs$OQV[day,] <- (term1 + term2)
  
}

            # ================================================#
            #   Realized implied volatility  variance (RIVV)
            # =============================================== #

obs$IVQV         <- matrix(nrow=S$n.days,ncol=S$n.options)

# The derivatives of IV with respect to Y and V have to be computed by finite differences
S0.plus.epsy       <- exp(log(SS)+S$epsy)
S0.minus.epsy      <- exp(log(SS)-S$epsy)
mon.plus.epsy      <- S0.plus.epsy/K
mon.minus.epsy     <- S0.minus.epsy/K
loc.plus.epsy      <- cbind( as.vector(logV), as.vector(mon.plus.epsy))
loc.minus.epsy     <- cbind( as.vector(logV), as.vector(mon.minus.epsy))
IV.plus.epsy       <- matrix(interp.surface( S$obj.30.iv, loc.plus.epsy), 
                             ncol=S$n.options)
IV.minus.epsy      <- matrix(interp.surface( S$obj.30.iv, loc.minus.epsy), 
                             ncol=S$n.options)
dIV.y              <- (IV.plus.epsy-IV.minus.epsy)/(2*S$epsy)

# Get IV with V+eps/V-eps
logV.plus          <- log(V0+S$epsv)
V.minus            <- V0-S$epsv
V.minus[which(V.minus < P$vmin )] <- P$vmin
logV.minus         <- log(V.minus)
loc.plus.epsv      <- cbind( as.vector(logV.plus), as.vector(moneyness))
loc.minus.epsv     <- cbind( as.vector(logV.minus), as.vector(moneyness))
IV.plus.epsv       <- matrix(interp.surface( S$obj.30.iv, loc.plus.epsv), 
                             ncol=S$n.options)
IV.minus.epsv      <- matrix(interp.surface( S$obj.30.iv, loc.minus.epsv), 
                             ncol=S$n.options)
dIV.v              <- (IV.plus.epsv-IV.minus.epsv)/(2*S$epsv)


for(day in 1:S$n.days){
  
  i               <- day*S$n.step.obs
  dI.mat          <- matrix(rep(obs$dI[(i-S$n.step.obs+1):i],S$ld*S$lt),
                            nrow=S$n.step.obs, ncol=S$n.options)
  
  # aggregation
  mat1            <- (dIV.y[(i-S$n.step.obs+1):(i),]^2 + 
                        P$sigma^2 * dIV.v[(i-S$n.step.obs+1):(i),]^2 +
                        2 * P$rho * P$sigma * dIV.v[(i-S$n.step.obs+1):(i),] *
                        dIV.y[(i-S$n.step.obs+1):(i),]) * dI.mat
  IV.jump         <- (obs$IV[(i-S$n.step.obs+2):(i+1),] -
                        obs$IV_[(i-S$n.step.obs+2):(i+1),])^2
  term1           <- apply(mat1 ,2,sum)
  term2           <- apply(IV.jump, 2, sum)
  
  #IVQV
  obs$IVQV[day,]      <- (term1 + term2)
  
  
}

# Add error terms 
dimiv   <- dim(obs$IV.close)
dimoqv  <- dim(obs$OQV)
dimivqv <- dim(obs$IVQV)

error1 <- rnorm(n = length(obs$quad.var), 0, S$sde$z2)
error2 <- rnorm(n = length(obs$int.var), 0, S$sde$z3)
error3 <- matrix(rnorm(n = dimiv[1]*dimiv[2], 0, S$sde$z4), nrow=dimiv[1])
error4 <- matrix(rnorm(n = dimoqv[1]*dimoqv[2], 0, S$sde$z5), nrow=dimoqv[1])
error5 <- matrix(rnorm(n = dimivqv[1]*dimivqv[2], 0, S$sde$z6), nrow=dimivqv[1])

obs$RV   <- obs$quad.var * (1+error1)
obs$BV   <- obs$int.var * (1+error2)
obs$IV   <- obs$IV.close * (1+error3)
obs$ROV  <- obs$OQV * (1+error4)
obs$RIVV <- obs$IVQV * (1+error5)
