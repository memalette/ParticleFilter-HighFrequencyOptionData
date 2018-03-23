# =========================================================================== #
#                               GRID COMPUTATION
# =========================================================================== #
# this code should be used to compute the IV, delta, vega grids, if the model parameters are changed

# Load all functions related to the experiement
source("functions.R")

#load simulation parameters
source("simulation_setup.R")

#load model parameters
source("model_setup.R")

# param for both grids
lmon          <- length(S$mon)
llogv         <- length(S$logV)

# ============================================================================= #
#                        Grid 1: 30-day maturity options
# ============================================================================= #

price.grid.30 <- matrix(ncol=lmon, nrow=llogv)
delta.grid.30 <- matrix(ncol=lmon, nrow=llogv)
vega.grid.30  <- matrix(ncol=lmon, nrow=llogv)
IV.grid.30    <- matrix(ncol=lmon, nrow=llogv)


for(i in 1:llogv){ #llogv
  
  V                 <- rep(exp(S$logV[i]), lmon)

  price.grid.30[i,] <- GetPriceVectorized(A1i,A1r ,A2i ,A2r ,B1i ,B1r ,B2i ,B2r ,
                                          V=V, Y=log(1), K=1/S$mon, 1, S)
  greeks            <- GetGreeksVectorized(A1i,A1r ,A2i ,A2r ,B1i ,B1r ,B2i ,B2r , 
                                           V=V, Y=log(1), K=1/S$mon, 1, S)
  delta.grid.30[i,] <- greeks$delta
  vega.grid.30[i,]  <- greeks$vega

  IV.grid.30[i,]    <- mapply(FUN = GBSVolatility, 
                              price = price.grid.30[i,],
                              S = 1,
                              X = 1/S$mon,
                              MoreArgs = list(TypeFlag = "c" ,
                                              r = P$r,
                                              b = P$q,
                                              Time = S$tau[1]))
  
  print(paste0("iteration: ", i))
}

save(IV.grid.30, file = "IV_grid_30.RData")
save(price.grid.30, file = "price_grid_30.RData")
save(delta.grid.30, file = "delta_grid_30.RData")
save(vega.grid.30, file = "vega_grid_30.RData")

# ============================================================================= #
#                     Grid 2: 90-day maturity options
# ============================================================================= #

price.grid.90 <- matrix(ncol=lmon, nrow=llogv)
delta.grid.90 <- matrix(ncol=lmon, nrow=llogv)
vega.grid.90  <- matrix(ncol=lmon, nrow=llogv)
IV.grid.90    <- matrix(ncol=lmon, nrow=llogv)


for(i in 1:llogv){ #llogv
  
  V                 <- rep(exp(S$logV[i]), lmon)
  
  price.grid.90[i,] <- GetPriceVectorized(A1i,A1r ,A2i ,A2r ,B1i ,B1r ,B2i ,B2r ,
                                          V=V, Y=log(1),K=1/S$mon,2, S)
  greeks            <- GetGreeksVectorized(A1i,A1r ,A2i ,A2r ,B1i ,B1r ,B2i ,B2r ,
                                           V=V, Y=log(1), K=1/S$mon,2, S)
  delta.grid.90[i,] <- greeks$delta
  vega.grid.90[i,]  <- greeks$vega
  
  IV.grid.90[i,]    <- mapply(FUN = GBSVolatility, 
                              price = price.grid.90[i,],
                              S = 1,
                              X = 1/S$mon,
                              MoreArgs = list(TypeFlag = "c" ,
                                              r = P$r,
                                              b = P$q,
                                              Time = S$tau[2]))
  print(paste0("iteration: ", i))
  
}

save(IV.grid.90, file = "IV_grid_90.RData")
save(price.grid.90, file = "price_grid_90.RData")
save(delta.grid.90, file = "delta_grid_90.RData")
save(vega.grid.90, file = "vega_grid_90.RData")
