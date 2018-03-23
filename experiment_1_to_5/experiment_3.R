# =========================================================================== #
#                                EXPERIMENT 3
# =========================================================================== #

# Load all functions related to the experiement
source("functions.R")

#load simulation parameters
source("simulation_setup.R")

#load model parameters
source("model_setup.R")

# Generated 1000 1-day paths 
source("observation_simulation.R")

# SELECT WITH AND WITHOUT JUMPS + SAMPLE SIZE
S$jump    <- TRUE # true for paths with jumps
l.sample  <- 1000

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

n.type <- 5
RV     <- apply(Y, 1, GetRVSimulation) 
dN.y   <- apply(ifelse(Z.y !=0,1,0), 1, sum)
dN.v   <- apply(ifelse(Z.v !=0,1,0), 1, sum)

# option related parameters
BSdelta <- c(S$delta)
is.call <- rep(1, each=S$ld)
tau     <- 30/252
ORV     <- matrix(nrow=S$n.traj.obs, ncol=n.type)
IVRV    <- matrix(nrow=S$n.traj.obs, ncol=n.type)
coeff.orv  <- matrix(ncol=8,nrow=n.type)
coeff.ivrv <- matrix(ncol=8,nrow=n.type)
R2.orv     <- matrix(ncol=2,nrow=n.type)
R2.ivrv    <- matrix(ncol=2,nrow=n.type)
coeff.orv2  <- matrix(ncol=6,nrow=n.type)
coeff.ivrv2 <- matrix(ncol=6,nrow=n.type)
coeff.ivrv3 <- matrix(ncol=6,nrow=n.type)
R2.orv2     <- matrix(ncol=2,nrow=n.type)
R2.ivrv2    <- matrix(ncol=2,nrow=n.type)
R2.ivrv3    <- matrix(ncol=2,nrow=n.type)

for(o in 1:n.type){
  
  # option parameters
  fun                <- function(x) GetDelta(x, V=P$v0, t=tau,
                                 r=P$r, q=P$q, is.call = is.call[o]) - BSdelta[o]
  moneyness          <- uniroot(f = fun, interval = c(0,2))$root
  K                  <- exp(P$y0)/moneyness
  S0                 <- exp(Y)
  moneyness.intraday <- S0/K
  V0                 <- V
  logV               <- log(V0)
  
  # location of observations
  loc <- cbind( as.vector(logV), as.vector(moneyness.intraday))
  
  #  interpolation
  IV <- matrix(interp.surface( S$obj.30.iv, loc), ncol=S$n.step.obs+1)
  prices   <- GetBlackScholesCall(S=S0,
                                  K=K ,
                                  t=tau,
                                  r=P$r,
                                  q=P$q,
                                  sigma=IV,
                                  is.call=is.call[o])
  
  ORV[,o] <- apply(prices, 1, GetRVSimulation)
  IVRV[,o]<- apply(IV, 1, GetRVSimulation)
    
    # ORV Regressions with RV
    x1 <- log(RV)
    x2 <- dN.y
    x3 <- dN.v
    
    ORV.fit      <- lm(log(ORV[,o]) ~ x1 + x2 + x3)
    summ         <- summary(ORV.fit) # show results
    coeff.orv[o,]<- c(summ$coefficients[,1], summ$coefficients[,2])
    R2.orv[o,]   <- c(summ$r.squared, summ$adj.r.squared)
    
    IVRV.fit      <- lm(log(IVRV[,o]) ~ x1 + x2 + x3)
    summ2         <- summary(IVRV.fit) # show results
    coeff.ivrv[o,]<- c(summ2$coefficients[,1], summ2$coefficients[,2])
    R2.ivrv[o,]   <- c(summ2$r.squared, summ$adj.r.squared)
    
    ORV.fit       <- lm(log(ORV[,o]) ~ x1 + x2 )
    summ          <- summary(ORV.fit) # show results
    coeff.orv2[o,]<- c(summ$coefficients[,1], summ$coefficients[,2])
    R2.orv2[o,]   <- c(summ$r.squared, summ$adj.r.squared)
    
    IVRV.fit       <- lm(log(IVRV[,o]) ~ x1 + x2 )
    summ2          <- summary(IVRV.fit) # show results
    coeff.ivrv2[o,]<- c(summ2$coefficients[,1], summ2$coefficients[,2])
    R2.ivrv2[o,]   <- c(summ2$r.squared, summ$adj.r.squared)
    
    IVRV.fit       <- lm(log(IVRV[,o]) ~ x1 + x3)
    summ3          <- summary(IVRV.fit) # show results
    coeff.ivrv3[o,]<- c(summ3$coefficients[,1], summ3$coefficients[,2])
    R2.ivrv3[o,]   <- c(summ3$r.squared, summ$adj.r.squared)
  
    print(o)

}

# Coeff ORV / IVRV
colnames <- c("b0", "b1", "b2", "b3", "SE0", "SE1", "SE2", "SE3")
rownames <- c("call 0.2", "call 0.35", "call 0.50", "call 0.65", "call 0.80")

colnames(coeff.orv) <- colnames
rownames(coeff.orv) <- rownames

colnames(coeff.ivrv) <- colnames
rownames(coeff.ivrv) <- rownames

# Coeff ORV / IVRV 2
colnames2 <- c("b0", "b1", "b2", "SE0", "SE1", "SE2")
rownames2 <- c("call 0.2", "call 0.35", "call 0.50", "call 0.65", "call 0.80")

colnames(coeff.orv2) <- colnames2
rownames(coeff.orv2) <- rownames2

colnames(coeff.ivrv2) <- colnames2
rownames(coeff.ivrv2) <- rownames2

# r2 ORV / IVRV
colnames3 <- c("R2", "R2 adjusted")
rownames3 <- c("call 0.2", "call 0.35", "call 0.50", "call 0.65", "call 0.80")

colnames(R2.orv) <- colnames3
rownames(R2.orv) <- rownames3

colnames(R2.ivrv) <- colnames3
rownames(R2.ivrv) <- rownames3

colnames(R2.orv2) <- colnames3
rownames(R2.orv2) <- rownames3

colnames(R2.ivrv2) <- colnames3
rownames(R2.ivrv2) <- rownames3


write.table(x = coeff.orv, file = "coeff_orv_reg1.csv", sep = ",", col.names = colnames)
write.table(x = coeff.ivrv, file = "coeff_ivrv_reg1.csv", sep = ",", col.names = colnames)
write.table(x = coeff.orv2, file = "coeff_orv_reg2.csv", sep = ",", col.names = colnames2)
write.table(x = coeff.ivrv2, file = "coeff_ivrv_reg2.csv", sep = ",", col.names = colnames2)
write.table(x = coeff.ivrv3, file = "coeff_ivrv_reg3.csv", sep = ",", col.names = colnames2)

write.table(x = R2.orv, file = "R2_orv_reg1.csv", sep = ",", col.names = colnames3)
write.table(x = R2.ivrv, file = "R2_ivrv_reg1.csv", sep = ",", col.names = colnames3)
write.table(x = R2.orv2, file = "R2_orv_reg2.csv", sep = ",", col.names = colnames3)
write.table(x = R2.ivrv2, file = "R2_ivrv_reg2.csv", sep = ",", col.names = colnames3)
write.table(x = R2.ivrv3, file = "R2_ivrv_reg3.csv", sep = ",", col.names = colnames3)

