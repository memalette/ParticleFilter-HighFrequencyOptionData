# =========================================================================== #
#                                EXPERIMENT 2
# =========================================================================== #

# Load all functions related to the experiement
source("functions.R")

#load simulation parameters
source("simulation_setup.R")

#load model parameters
source("model_setup.R")

# Generated 1000 1-day paths 
source("observation_simulation.R")

# number of paths used in this experiment
l.sample  <- 1000

# get indices of paths with either a jump in V or Y
ijyv      <- unique(which(Z.y != 0 | Z.v != 0, arr.ind = TRUE)[,1])
# get indices of paths without any jumps 
inoj      <- as.vector(1:S$n.traj.obs)[-ijyv]
i        <- sample(x = inoj,size = l.sample,replace = FALSE)
V        <- V[i,]
V.start  <- V.start[i]
Y        <- Y[i,]
Y.start  <- Y.start[i]
Z.y      <- Z.y[i,]
Z.v      <- Z.v[i,]
dI       <- dI[i,]

# randomly generate Black-Scholes delta
BSdelta <- runif(n=l.sample,0.1,0.9)
tau     <- 30/252
moneyness <- vector()

# Get implied moneyness from deltas
for(o in 1:l.sample){
  
  fun             <- function(x) GetDelta(x, V=V.start[o], t=tau,
                                          r=P$r, q=P$q, is.call = 1) - BSdelta[o]

  if(is.na(fun(0.1)) | is.na(fun(2))){print(o)}
  moneyness[o]    <- uniroot(f = fun, interval = c(0.1,2))$root

}

S0                 <- exp(Y)
K                  <- exp(Y.start)/moneyness
K.mat              <- matrix(rep(K, S$n.step.obs+1),ncol=S$n.step.obs+1)
moneyness.intraday <- S0/K.mat
V0                 <- V
logV               <- log(V0)

# location of observations
loc <- cbind( as.vector(logV), as.vector(moneyness.intraday))

#  interpolation
IV <- matrix(interp.surface( S$obj.30.iv, loc), ncol=S$n.step.obs+1)
prices   <- GetBlackScholesCall(S=S0,
                                K=K.mat ,
                                t=tau,
                                r=P$r,
                                q=P$q,
                                sigma=IV,
                                is.call=1)

loc  <- cbind(logV[,1], moneyness.intraday[,1])
delta <- interp.surface( S$obj.30.d, loc)*S0[,1]
vega  <- interp.surface( S$obj.30.v, loc)*S0[,1]

ROV  <- apply(prices, 1, GetRVSimulation) 
RIVV <- apply(IV, 1, GetRVSimulation) 
RV   <- apply(Y, 1, GetRVSimulation) 

epsy=0.001
epsv=0.00001
# Get IV with Y+eps/Y-eps
S0.plus.epsy       <- exp(Y[,1]+epsy)
S0.minus.epsy      <- exp(Y[,1]-epsy)
mon.plus.epsy      <- S0.plus.epsy/K
mon.minus.epsy     <- S0.minus.epsy/K
loc.plus.epsy      <- cbind( as.vector(logV[,1]), as.vector(mon.plus.epsy))
loc.minus.epsy     <- cbind( as.vector(logV[,1]), as.vector(mon.minus.epsy))
IV.plus.epsy       <- interp.surface( S$obj.30.iv, loc.plus.epsy)
IV.minus.epsy      <- interp.surface( S$obj.30.iv, loc.minus.epsy)
dIV.y30              <- (IV.plus.epsy-IV.minus.epsy)/(2*epsy)

# Get IV with V+eps/V-eps
logV.plus          <- log(V[,1]+epsv)
logV.minus         <- log(V[,1]-epsv)
loc.plus.epsv      <- cbind( as.vector(logV.plus), as.vector(moneyness.intraday[,1]))
loc.minus.epsv     <- cbind( as.vector(logV.minus), as.vector(moneyness.intraday[,1]))
IV.plus.epsv       <- interp.surface( S$obj.30.iv, loc.plus.epsv)
IV.minus.epsv      <- interp.surface( S$obj.30.iv, loc.minus.epsv)
dIV.v30              <- (IV.plus.epsv-IV.minus.epsv)/(2*epsv)


# ROV Regressions with RV
x1 <- delta^2 * RV
x2 <- 2 * delta * vega * RV
x3 <- vega^2 * RV 

ORV.fit <- lm(ORV ~ 0 + x1 + x2 + x3)
summary(ORV.fit) # show results

# RIVV Regressions with RV
x1 <- dIV.y30^2 * RV
x2 <- 2 * dIV.y30 * dIV.v30 * RV
x3 <- dIV.v30^2 * RV

IVRV.fit <- lm(IVRV ~ 0 + x1 + x2 + x3)
summ <- summary(IVRV.fit) # show results

output <- rbind(summ$coefficients, c(summ$r.squared, 0, 0, 0))

write.table(x = output, file = "results_experiment_2.csv" )






