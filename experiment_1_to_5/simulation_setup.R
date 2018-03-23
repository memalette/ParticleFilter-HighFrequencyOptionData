S <- list() # simulation parameters

# User-defined parameters:
S$include       <- list(y = T, rv=F, bv=F, iv=F, rov=F, rivv=F) # TRUE = included as a source of information
S$n.days        <- 252 # number of days for the simulation
S$n.step        <- 1 # number of steps in a day

# Fixed simulation parameters 
S$n.days.year   <- 252 # number of business days in a year
S$n.traj.obs    <- 2 # number of observation path to simulate
S$n.traj        <- 50000 # number of trajectories in simulation
S$n.step.obs    <- 1560 # number of steps for observation generation
S$dt            <- 1/S$n.days.year # time discretization: 1 day in years
S$h.obs         <- S$dt/S$n.step.obs # size of steps in years 
S$h             <- S$dt/S$n.step # size of steps in years 
S$sde           <- list(z2=0.05, z3=0.05, z4=0.05, z5=0.05, z6=0.05) # measurement errors of observables
S$delta         <- c(0.2, 0.35, 0.5, 0.65, 0.8) # option delta
S$tau           <- c(30/252, 90/252) # option maturities
S$x             <- exp(seq(from=-6, to=6, by=0.005)) # Used for option pricing
S$n.options     <- length(S$tau)*length(S$delta) # number of options in panel
S$ld            <- length(S$delta) # number of different moneynesses
S$lt            <- length(S$tau) # number of different maturities

# grids for option pricing interpolation
S$logV  <- seq(from=-12, to=1, by=0.05)
S$mon   <- seq(from=0.3, to=2, by=0.005)
S$epsy  <- 0.001
S$epsv  <- 0.001

load("IV_grid_30.RData")
load("IV_grid_90.RData")
load("delta_grid_30.RData")
load("delta_grid_90.RData")
load("vega_grid_30.RData")
load("vega_grid_90.RData")

S$IV.grid.30    <- IV.grid.30
S$IV.grid.90    <- IV.grid.90
S$delta.grid.30 <- delta.grid.30
S$delta.grid.90 <- delta.grid.90
S$vega.grid.30  <- vega.grid.30
S$vega.grid.90  <- vega.grid.90

S$obj.30.iv     <- list(x=S$logV, y=S$mon, z=S$IV.grid.30)
S$obj.90.iv     <- list(x=S$logV, y=S$mon, z=S$IV.grid.90)
S$obj.30.d      <- list(x=S$logV, y=S$mon, z=S$delta.grid.30)
S$obj.90.d      <- list(x=S$logV, y=S$mon, z=S$delta.grid.90)
S$obj.30.v      <- list(x=S$logV, y=S$mon, z=S$vega.grid.30)
S$obj.90.v      <- list(x=S$logV, y=S$mon, z=S$vega.grid.90)

pkgs <- c("Rcpp", "timeDate","timeSeries", "fBasics","pracma", "statmod",
          "ggplot2", "gridExtra","reshape2", "labeling","fOptions", "NMOF", "profvis", 
          "jsonlite", "yaml","htmltools", "htmlwidgets", "foreach", "bizdays", "rootSolve", 
          "snow", "doSNOW", "doParallel","rbenchmark","iterators", "rlist", "dplyr", "fields", 
          "abind", "statmod", "smcUtils")

ipak(pkgs)