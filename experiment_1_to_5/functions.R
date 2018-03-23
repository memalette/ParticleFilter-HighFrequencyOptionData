
# =========================================================================== #
#                 Packages necessary for the experiments
# =========================================================================== #

# do you want to install or just load the packages, 
# specify the path where the libraries are saved 

# ipak function: install and load multiple R packages.
# check to see if packages are installed. Install them if they are not, then load them into the R session.

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# =========================================================================== #
#                  Black-Scholes option pricing (vectorized) 
# =========================================================================== #

# S, K, t , r, q, price have to be the same length
GetBlackScholesCall <- function(S, K, t, r, q, sigma, is.call){
  d1    <- (log(S/K) + (q + 0.5*sigma^2)* sqrt(t))/(sigma*sqrt(t))
  d2    <- d1 - sigma*sqrt(t)
  call  <- S*exp((q-r)*t)*pnorm(d1) - K*exp(-r*t)*pnorm(d2)
  price <- is.call * call + (1-is.call)*(call + exp(-r*t)*K - S)
  price[which(price < 0)] <- 0.0001
  return(price)
}


# =========================================================================== #
#          ABG/HESTON option pricing of Affine jump-diffusion-model 
# =========================================================================== #

# A and B matrices are provided by "model_parameters.R", 
# V, Y, K, t are vectors of same length
# S is the list of simulation parameters from "simulation_parameters.R"

GetPriceVectorized <- function(A1i,A1r ,A2i ,A2r ,B1i ,B1r ,B2i ,B2r , V, Y,  K, t, S){  
  
  k         <- log( exp(Y) / K )
  ly        <- length(Y)
  integral1 <- integral2 <- rep(0, ly)
  
  A1i.v <- A1i[1, t]
  A1r.v <- A1r[1, t]
  A2i.v <- A2i[1, t]
  A2r.v <- A2r[1, t]
  B1i.v <- B1i[1, t]
  B1r.v <- B1r[1, t]
  B2i.v <- B2i[1, t]
  B2r.v <- B2r[1, t]
  
  u   <- rep(S$x[1], ly)
  
  # P1 / P2
  s1 <- sin( A1i.v + B1i.v * V + k * u )
  x1 <- exp( A1r.v + B1r.v * V ) * s1 / u
  
  s3 = sin( A2i.v + B2i.v * V  +  k * u )
  x3 = exp( A2r.v + B2r.v * V) * s3 / u
  
  
  for(j in 2:length(S$x)){
    
    A1i.v <- A1i[j, t]
    A1r.v <- A1r[j, t]
    A2i.v <- A2i[j, t]
    A2r.v <- A2r[j, t]
    B1i.v <- B1i[j, t]
    B1r.v <- B1r[j, t]
    B2i.v <- B2i[j, t]
    B2r.v <- B2r[j, t]
    
    u.1   <- rep(S$x[j], ly)
    
    # P1 / P2
    s2 <- sin( A1i.v + B1i.v * V + k * u.1 )
    x2 <- exp( A1r.v + B1r.v * V ) * s2 / u.1
    
    s4 <- sin( A2i.v + B2i.v * V + k * u.1 )
    x4 <- exp( A2r.v + B2r.v * V ) * s4 / u.1
    
    integral1 <- integral1 + (x2 + x1)*(u.1 - u)/2
    integral2 <- integral2 + (x3 + x4)*(u.1 - u)/2
    
    x1 <- x2
    x3 <- x4
    
    u <- u.1
    
  }
  
  # P1 / P2
  P2 <- 0.5 + integral2 / pi
  P1 <- 0.5 + integral1 / pi
  price <- exp(Y-P$q*S$tau[t]) * P1 - K * exp(-P$r*S$tau[t]) * P2
  price[which(price < 0)] <- 0
  
  return(price=price)
  
}

# =========================================================================== #
#                 Greeks from appendix section of ABG2017
# =========================================================================== #

# A and B matrices are provided by "model_parameters.R", 
# V, Y, K, t are vectors of same length
# S is the list of simulation parameters from "simulation_parameters.R"

GetGreeksVectorized <- function(A1i,A1r ,A2i ,A2r ,B1i ,B1r ,B2i ,B2r , V, Y,  K, t, S){  
  
  k         <- log( exp(Y) / K )
  ly        <- length(Y)
  
  integral0 <- rep(0, ly)
  integral1 <- integral0
  integral2 <- integral0
  integral3 <- integral0
  integral4 <- integral0
  integral5 <- integral0
  
  u   <- S$x
  
  c1 <- cos( A1i[1, t] + B1i[1, t] * V + k * u[1] )
  s1 <- sin( A1i[1, t] + B1i[1, t] * V + k * u[1] )
  x1 <- exp( A1r[1, t] + B1r[1, t] * V ) 
  
  c3 = cos( A2i[1, t] + B2i[1, t] * V  +  k * u[1] )
  s3 = sin( A2i[1, t] + B2i[1, t] * V  +  k * u[1] )
  x3 = exp( A2r[1, t] + B2r[1, t] * V) 
  
  
  for(j in 2:length(u)){
    
    
    c2 <- cos( A1i[j, t] + B1i[j, t] * V + k * u[j] )
    s2 <- sin( A1i[j, t] + B1i[j, t] * V + k * u[j] )
    x2 <- exp( A1r[j, t] + B1r[j, t] * V ) 
    
    c4 <- cos( A2i[j, t] + B2i[j, t] * V  +  k * u[j] )
    s4 <- sin( A2i[j, t] + B2i[j, t] * V  +  k * u[j] )
    x4 <- exp( A2r[j, t] + B2r[j, t] * V) 
    
    integral1 <- integral1 + (-x2*c2 + -x1*c1)*(u[j] - u[j-1])/2
    integral2 <- integral2 + (x2*s2/u[j] + x1*s2/u[j-1])*(u[j] - u[j-1])/2
    integral3 <- integral3 + (-x4*c4 + -x3*c3)*(u[j] - u[j-1])/2
    integral4 <- integral4 + ((x2*s2*B1r[j, t]/u[j] + x2*c2*B1i[j, t]/u[j]) + 
                                (x1*s1*B1r[j-1, t]/u[j-1] + x1*c1*B1i[j-1, t]/u[j-1]))*(u[j] - u[j-1])/2
    integral5 <- integral5 + ((x4*s4*B2r[j, t]/u[j] + x4*c4*B2i[j, t]/u[j]) + 
                                (x3*s3*B2r[j-1, t]/u[j-1] + x3*c3*B2i[j-1, t]/u[j-1]))*(u[j] - u[j-1])/2
    
    
    c1 <- c2
    s1 <- s2
    x1 <- x2
    c3 <- c4
    s3 <- s4
    x3 <- x4
    
    
  }
  
  
  delta <- exp(Y) * integral1/pi + exp(Y) * (0.5 + integral2/pi) - exp(-P$r*S$tau[t]) * K * integral3/pi 
  vega  <- exp(Y) * integral4/pi - exp(-P$r*S$tau[t]) * K * integral5/pi
  
  return(list(delta=delta, vega=vega))
  
}

# =========================================================================== #
#     Computes A, B matrices for ABG/Heston option and greeks computation
# =========================================================================== #

# Q is the list of model parameters under the risk neutral measure from "model_parameters.R"
# u is the integration variable
# maturity is a vector of maturities in DAYS for which the matrices need to be computed

GetAB2 <- function( Q, u, maturity){
  
  A      <- matrix(nrow = length(u), ncol = length(maturity))
  B      <- matrix(nrow = length(u), ncol = length(maturity))
  
  mgf.1  <- exp(Q$mu.y * 1 + 0.5 * Q$sigma.y^2 * 1^2)
  mgf.u  <- exp(Q$mu.y * u + 0.5 * Q$sigma.y^2 * u^2)
  
  C0 <-  u/2 + u*Q$lambda.y1*(mgf.1-1) - Q$lambda.y1*(mgf.u-1) - u^2/2
  C1 <-  Q$kappa - Q$sigma*Q$rho*u
  C2 <- -Q$sigma^2/2
  lambda <- sqrt(C1^2 - 4*C0*C2)
  
  
  D0 <-  Q$lambda.y0 *(u*(mgf.1-1) - (mgf.u-1))
  D1 <- -Q$kappa*Q$theta
  D2 <- -Q$lambda.v0
  
  for(d in 1:length(maturity)){
    
    dt <- 1/252
    t <- dt * maturity[d]
    
    a <- -(t * (lambda + C1) + 2 * log(2 * lambda) - 
             2 * log( (exp(t * lambda) + 1) * lambda + C1 * (exp(t * lambda) - 1) ) )/
      (2 * C2)
    
    b <- -2 * C0 * Q$mu.v * ( lambda * t + 2 * log( 2 * lambda) - 
                                2 * log( (exp(t * lambda) + 1) * lambda + C1 * (exp(t * lambda) - 1) + 
                                           2 * (exp(t * lambda) - 1) * C0 * Q$mu.v ) + 
                                t * C1 + 2 * t * C0 * Q$mu.v ) / 
      ( (-lambda + C1 + 2 * C0 * Q$mu.v) * ( lambda + C1 + 2 * C0 * Q$mu.v) )
    
    A[, d] <- - D0*t - D1*a - D2*b
    B[, d] <-  2*C0*(exp(-lambda*t)-1)/( lambda*(exp(-lambda*t)+1) - C1*(exp(-lambda*t)-1) )
    
  }
  
  colnames(A) <- as.character(maturity)
  colnames(B) <- as.character(maturity)
  
  Ar <- Re(A)
  Ai <- Im(A)
  Br <- Re(B)
  Bi <- Im(B)
  
  coeff <- list(Ar=Ar, Ai=Ai, Br=Br, Bi=Bi)
  
  return(coeff)
  
}

# =========================================================================== #
#             Computes Cache for the Integrated Variance Moments
# =========================================================================== #

# P is the list of model parameters under the physical measure: "model_parameters.R"
# cache.x is the vector of values for which moments are calculated (sqrt(V1*V2))
# delta is the time interval (1/252)/M, refered to as S$h in the rest of the code

GetCacheIntVar <- function(P, cache.x , delta){
  
  kappa <- P$kappa
  theta <- P$theta
  sigma <- P$sigma
  
  d    <-  4*theta*kappa/(sigma^2)
  nu   <- d/2-1
  z    <- 2*kappa*(sigma^2*sinh(kappa*delta/2))^(-1)*cache.x
  f    <- kappa*delta/2
  # C1   <- coth(f)
  C1   <- (1 + exp(-2*f))/(1 - exp(-2*f))
  C2   <- 4*exp(-2*f)/((1 - exp(-2*f))^2)
  # C2   <- csch(kappa*delta/2)^2
  D    <- 4*kappa*theta/sigma^2
  G <- (-8 + 2 * P$kappa * delta * C1 + P$kappa^2 * delta^2 * C2)
  # If dt is very small, G can be small negative and causes problems
  if(G < 0){  G <- 0 }
  
  En   <-  z*besselI(z, nu+1,expon.scaled = T)/(2*besselI(z, nu,expon.scaled = T))
  En2  <- z^2*besselI(z, nu+2,expon.scaled = T)/(4*besselI(z, nu,expon.scaled = T)) + En
  
  EX2  <- D*sigma^2*(-2+kappa*delta*C1)/(4*kappa^2)
  sX2  <- D*sigma^4*G/(8*kappa^4)
  EZ   <- 4*EX2/D
  sZ   <- 4*sX2/D
  
  averI  <- En*EZ + EX2
  varI   <- En*sZ + (En2-En^2)*EZ^2 + sX2
  
 
  return(cbind(averI, varI))
  
  
}


# =========================================================================== #
#                RV and BV, simple simulation (no subsampling)
# =========================================================================== #

# input a time series and it computes the RV/BV. 

GetRVSimulation <- function(my.ts){
  
  RV <- sum(diff(my.ts)^2)
  
  return(RV)
  
}

GetBVSimulation <- function(my.ts){
  
  diff.abs <- abs(diff(my.ts))
  BV       <- 0
  
  for(i in 2:length(diff.abs)){
    
    BV <- BV + diff.abs[i]*diff.abs[i-1]
    
  }
  
  return(BV)
  
}

# =========================================================================== #
#             Function that computes Delta (BS) of an option 
# =========================================================================== #

# all inputs without defaults have to be vectors of same length

GetDelta <- function(moneyness, V, t, r, q){
  
  d1    <- ( log(moneyness) + (r - q + 0.5*V) * t ) / ( sqrt(t*V) )
  delta <- pnorm(d1, 0, 1)
  
  return(delta)
  
}

# =========================================================================== #
#       Moneyness from fixed call-equivalent deltas (vectorized Bisection) 
# =========================================================================== #

# all inputs without defaults have to be vectors of same length

GetMoneynessFromDelta <- function(V, t, r, q, delta, low = -5, high=5, n.max = 1000, tol = 1e-4){ 
  
  # SETUP 
  
  i           <- 1
  lp          <- length(delta)
  
  # create data set with pricing parameters
  data        <- data.frame(id=1:lp,
                            V=V,
                            t=t,
                            r=r,
                            q=q,
                            delta = delta,
                            a=rep(low,lp),
                            b=rep(high,lp),
                            moneyness=rep(FALSE,lp))
  
  GetDelta <- function(moneyness, V, t, r, q){
    
    d1    <- ( log(moneyness) + (r - q + 0.5*V) * t ) / ( sqrt(t*V) )
    delta <- pnorm(d1, 0, 1)
    
    return(delta)
    
  }
  
  # create function for which we want to find roots, where the implied volatility is the only variable
  fun <- function(moneyness) GetDelta(moneyness, 
                                      V=data$V,
                                      t=data$t,
                                      r=data$r,
                                      q=data$q) 
  
  
  
  # Initial conditions: lower boundary must be lower than upper boundary
  # And, sign of function evaluated at the boundaries must be opposite
  fa         <- fun(data$a) - data$delta
  fb         <- fun(data$b) - data$delta
  
  # Return NA automatically if these conditions are satisfied
  idx.na       <- which(is.na(fa) |is.na(fb) | data$a > data$b)
  #idx.sign     <- which(sign(diff.at.low) == sign(diff.at.high))
  
  if(length(idx.na) > 0){
    new.data <- data[-idx.na,]
    data$moneyness[idx.na] <- NA
  }else{
    new.data <- data
  }
  
  while(i < n.max){
    #redefine fun for fixed new.data
    fun <- function(moneyness) GetDelta(moneyness, 
                                        V=new.data$V,
                                        t=new.data$t,
                                        r=new.data$r,
                                        q=new.data$q) 
    
    # Get midpoint and difference at midpoint
    new.data$c   <- (new.data$a + new.data$b) / 2
    new.data$fc  <- fun(new.data$c) - new.data$delta
    
    # This is the condition that will stop when we reach the implied volatility from the call price
    idx.success   <- which(new.data$fc == 0 | (new.data$b - new.data$a)/2 < tol )
    id.success    <- new.data$id[idx.success]
    if(length(id.success) > 0){
      data$moneyness[id.success] <- new.data$c[idx.success]
      if(sum(as.numeric(is.na(data$moneyness)|data$moneyness!=FALSE | data$moneyness==0))==lp) {
        return(data$moneyness)
      }else{
        new.data <- new.data[-idx.success,]
      }
    }
    
    if(isempty(new.data)){return(data$moneyness)}
    
    i <- i + 1
    
    # This is how we determine how to narrow down the range
    idx.a             <- which( fun(new.data$c)  < new.data$delta)
    idx.b             <- which( fun(new.data$c)  > new.data$delta)
    new.data$a[idx.a] <- new.data$c[idx.a]
    new.data$b[idx.b] <- new.data$c[idx.b]
    new.data
    data
    
  }
  
  print('Too many iterations')
  return(data$vol)
}


# =========================================================================== #
#                              Memory Management 
# =========================================================================== #

# improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing=FALSE, head=FALSE, n=5) {
  napply <- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.size <- napply(names, object.size)
  obj.dim <- t(napply(names, function(x)
    as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.dim)
  names(out) <- c("Type", "Size", "Rows", "Columns")
  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}
# shorthand
lsos <- function(..., n=10) {
  .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}
