# =========================================================================== #
#                                EXPERIMENT 5
# =========================================================================== #

# Load all functions related to the experiement
source("functions.R")

#load simulation parameters
source("simulation_setup.R")

#load model parameters
source("model_setup.R")

# Generated 1000 1-day paths 
source("observation_simulation.R")

l.sample  <- 1000

# get indices of paths with either a jump in V or Y
ijyv      <- unique(which(Z.y != 0 | Z.v != 0, arr.ind = TRUE)[,1])
# get indices of paths without any jumps
inoj     <- as.vector(1:S$n.traj.obs)[-ijyv]
i        <- sample(x = inoj,size = l.sample,replace = FALSE)
# ijyv   <- unique(which(Z.y != 0 | Z.v != 0, arr.ind = TRUE)[,1])
# i      <- sample(x = ijyv,size = l.sample,replace = FALSE)
V        <- V[i,]
V.start  <- V.start[i]
Y        <- Y[i,]
Y.start  <- Y.start[i]
Z.y      <- Z.y[i,]
Z.v      <- Z.v[i,]
dI       <- dI[i,]


BSdelta <- runif(n=l.sample,0.6,0.9)
tau     <- 90/252
moneyness <- vector()

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

ORV  <- apply(prices, 1, GetRVSimulation) 
IVRV <- apply(IV, 1, GetRVSimulation) 
RV   <- apply(Y, 1, GetRVSimulation) 


fit <- lm(IVRV ~ ORV )
summary(fit) # show results





