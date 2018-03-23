# =========================================================================== #
#                                EXPERIMENT 4
# =========================================================================== #

# Load all functions related to the experiement
source("functions.R")

#load simulation parameters
source("simulation_setup.R")

# Generated 500 1-day paths 
source("model_setup.R")
source("observation_simulation.R")

# SELECT WITH AND WITHOUT JUMPS + SAMPLE SIZE
S$jump    <- FALSE
l.sample  <- 500

# get indices of paths with either a jump in V or Y
ijyv      <- unique(which(Z.y != 0 | Z.v != 0, arr.ind = TRUE)[,1])
# get indices of paths without any jumps 
inoj      <- as.vector(1:S$n.traj.obs)[-ijyv]

# sample 500 paths with/without jumps
if(S$jump == FALSE){
    i        <- sample(x = inoj,size = l.sample,replace = FALSE)
    V        <- V[i,]
    Y        <- Y[i,]
    Z.y      <- Z.y[i,]
    Z.v      <- Z.v[i,]
    dI       <- dI[i,]
}else{ 
    i        <- sample(x = ijyv,size = l.sample,replace = FALSE)
    V        <- V[i,]
    Y        <- Y[i,]
    Z.y      <- Z.y[i,]
    Z.v      <- Z.v[i,]
    dI       <- dI[i,]
}


# option related parameters
BSdelta <- S$delta
is.call <- 1
tau     <- 30/252
M       <- c(1,2,3,5,10,1560)

lm <- length(M)
rmse.orv   <- matrix(ncol = S$ld, nrow = lm)
rrmse.orv  <- matrix(ncol = S$ld, nrow = lm)
r2.orv     <- matrix(ncol = S$ld, nrow = lm)
rmse.ivrv  <- matrix(ncol = S$ld, nrow = lm)
rrmse.ivrv <- matrix(ncol = S$ld, nrow = lm)
r2.ivrv    <- matrix(ncol = S$ld, nrow = lm)

for(o in 1:S$ld){ 
  
  tito   <- as.character(o)

  # Get end of day strikes for each path
    GetMoneyness <- function(V, t, r, q, is.call, BSdelta){ 
      fun                <- function(x) GetDelta(x, V=V, t=t,
                             r=r, q=q, is.call = is.call) - BSdelta
      moneyness          <- uniroot(f = fun, interval = c(0,2))$root
      return(moneyness)
    }
    
    moneyness            <- mapply(FUN=GetMoneyness, V=V[,ncol(V)], 
                                   MoreArgs = list(t=tau, 
                                                   r=P$r,
                                                   q=P$q, 
                                                   is.call=is.call,
                                                   BSdelta=BSdelta[o]))
    
    K                  <- exp(Y[,ncol(Y)])/moneyness
    K.mat              <- matrix(rep(K,ncol(Y)),nrow = nrow(Y))
    
    # info for grid interpolation 
    S0                 <- exp(Y)
    moneyness.intraday <- S0/K
    V0                 <- V
    logV               <- log(V0)
    loc                <- cbind( as.vector(logV), as.vector(moneyness.intraday))
    
    # Get IV and price through grid interpolation
    IV       <- matrix(interp.surface( S$obj.30.iv, loc), ncol=S$n.step.obs+1)
    prices   <- GetBlackScholesCall(S=S0,
                                    K=K.mat ,
                                    t=tau,
                                    r=P$r,
                                    q=P$q,
                                    sigma=IV,
                                    is.call=is.call)
    
    # Get Greeks through grid interpolation
    delta <- matrix(interp.surface( S$obj.30.d, loc), ncol=S$n.step.obs+1)*S0
    vega  <- matrix(interp.surface( S$obj.30.v, loc), ncol=S$n.step.obs+1)*S0
    
    # compute IVQV derivative 
    d                <- (log(S0/K) + (P$r - P$q + 0.5*IV^2) * tau)/
                        (IV * sqrt(tau))
    phi.d            <- dnorm(d, 0, 1)
    inverse.BS.vega  <- 1/(S0 * phi.d * sqrt(tau))
    dIV.y            <- delta * inverse.BS.vega
    dIV.v            <- vega * inverse.BS.vega
    
    # Get cumulative sum values 
    Z.y.cs  <- t(apply(Z.y, 1, cumsum))
    Z.v.cs  <- t(apply(Z.v, 1, cumsum))
    dI.cs   <- t(apply(dI, 1, cumsum))
    
    #For each M and each option delta, compute OQV/IVQV
    for(m in 1:lm){
      
      titm  <- as.character(M[m])
      title <- paste0("M = ", titm)
      
      # depending on M used, get beginning of period indices and end of period indices
      if(m == 1){
        ide <- 1560
        idb <- 1
      }else{
        ide <- seq(from=S$n.step.obs/M[m], to = S$n.step.obs, by = S$n.step.obs/M[m])
        idb <- seq(1, to = S$n.step.obs, by =  S$n.step.obs/M[m] )
      }
      
      Z.y.mat   <- cbind(rep(0, l.sample), Z.y.cs[,ide])
      Z.v.mat   <- cbind(rep(0, l.sample), Z.v.cs[,ide])
      dI.mat    <- cbind(rep(0, l.sample), dI.cs[,ide])
      dZ.y.mat  <- matrix(apply(Z.y.mat, 1, diff),nrow=l.sample, byrow=TRUE)
      dZ.v.mat  <- matrix(apply(Z.v.mat, 1, diff),nrow=l.sample, byrow=TRUE)
      I.mat     <- matrix(apply(dI.mat, 1, diff),nrow=l.sample, byrow=TRUE)
    
      # Get before jumps values for interpolation
      S0_                      <- exp(Y[,ide+1]-dZ.y.mat)
      moneyness.intraday_      <- S0_/K.mat[,ide+1]
      V0_                      <- V[,ide+1]-dZ.v.mat
      V0_[which(V0_ < P$vmin)] <- P$vmin
      logV_                    <- log(V0_)
      loc_                     <- cbind( as.vector(logV_), 
                                         as.vector(moneyness.intraday_))
      
      # Get IV and prices before jumps
      IV_       <- matrix(interp.surface( S$obj.30.iv, loc_), ncol=ncol(S0_))
      prices_   <- GetBlackScholesCall(S=S0_,
                                      K=K.mat[,ide+1] ,
                                      t=tau,
                                      r=P$r,
                                      q=P$q,
                                      sigma=IV_,
                                      is.call=is.call)

        #OQV
        mat1      <- (delta[,idb]^2 + (P$sigma * vega[,idb])^2 + 
                      2*P$sigma*P$rho*delta[,idb]*vega[,idb])* I.mat
        term1     <- apply(mat1 ,1,sum)
        mat2      <- matrix((prices[,ide+1] - prices_)^2, nrow=l.sample)
        term2     <- apply(mat2, 1, sum)
        OQV       <- (term1 + term2) 
        ORV       <- apply(prices, 1, GetRVSimulation)
        
        # Scatterplot of ORV/OQV
                max <- max(OQV)
                ran <- c(0,max)
                plot(x = ORV, y = OQV,xlab = "ORV", ylab = "OQV", main= title, 
                     xlim=ran, ylim=ran,pch=19)
                lines(x=ran,y=ran)
                imgname   <- paste0("scatterplot_ORV_",tito,"_", titm, ".png")
                dev.copy(png, imgname)
                dev.off()

       # Regression ORV/OQV
        rmse.orv[m,o]  <- sqrt (mean((OQV-ORV)^2))
        rrmse.orv[m,o] <- (sqrt(sum(OQV-ORV)^2))/ 
                          (length(ORV) * mean(ORV)) 
        oqv.fit        <- lm(OQV ~ ORV)
        r2.orv[m,o]    <- summary(oqv.fit)$r.squared
        
        # IVQV
        mat1      <- (dIV.y[,idb]^2 + (P$sigma * dIV.v[,idb])^2 + 
                      2*P$sigma*P$rho*dIV.y[,idb]*dIV.v[,idb])* I.mat
        term1     <- apply(mat1 ,1,sum)
        mat2      <- matrix((IV[,ide+1] - IV_)^2, nrow=l.sample)
        term2     <- apply(mat2, 1, sum)
      
        IVQV        <- (term1 + term2 )
        idx.ignore  <- which(is.nan(IVQV) | IVQV > 0.1) # when the option is too deep out of the money, inverse BS vega = Inf, and IVQV is NaN
        
        if(length(idx.ignore) > 0){
        IVQV            <- IVQV[-idx.ignore]
        IVRV            <- apply(IV, 1, GetRVSimulation)[-idx.ignore]
        }else{
          IVQV            <- IVQV
          IVRV            <- apply(IV, 1, GetRVSimulation)
        }
          
        # Scatterplot of IVRV/IVQV
                max <- max(IVQV)
                ran <- c(0,max)
                plot(x = IVRV, y = IVQV,xlab = "IVRV", ylab = "IVQV", main= title, 
                     xlim=ran, ylim=ran,pch=19)
                lines(x=ran,y=ran)
                titm <- as.character(M[m])
                imgname   <- paste0("scatterplot_IVRV_",tito,"_", titm, ".png")
                dev.copy(png, imgname)
                dev.off()
                
        # Regression IVRV/VQV
        rmse.ivrv[m,o]  <- sqrt (mean((IVQV-IVRV)^2))
        rrmse.ivrv[m,o] <- (sqrt(sum(IVQV-IVRV)^2))/ 
                           (length(IVRV) * mean(IVRV)) 
        ivqv.fit   <- lm(IVQV ~ IVRV)
        r2.ivrv[m,o]    <- summary(ivqv.fit)$r.squared
        
    }

    

}

colname <- as.character(S$delta)
rowname <- as.character(M)

colnames(rmse.ivrv) <- colname
colnames(rmse.orv)  <- colname
colnames(rrmse.ivrv)<- colname
colnames(rrmse.orv) <- colname
colnames(r2.ivrv)   <- colname
colnames(r2.orv)    <- colname

row.names(rmse.ivrv) <- rowname
row.names(rmse.orv)  <- rowname
row.names(rrmse.ivrv)<- rowname
row.names(rrmse.orv) <- rowname
row.names(r2.ivrv)   <- rowname
row.names(r2.orv)    <- rowname

results.ivrv <- rbind(rmse.ivrv,rrmse.ivrv,r2.ivrv)
results.orv <- rbind(rmse.orv,rrmse.orv,r2.orv)

namefile <- ifelse(S$jump == TRUE, "withjump.csv", "withoutjump.csv")

write.table(x= results.ivrv, file = paste0("Results_IVRV_", namefile ), sep = ",", col.names = colname)
write.table(x= results.orv, file = paste0("Results_ORV_", namefile ), sep = ",", col.names = colname)




