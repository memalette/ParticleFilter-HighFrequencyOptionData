# Parameters under P
P  <- list()
# Log process
P$r           <- 0.01
P$q           <- 0.01
P$rho         <- -0.4336
P$lambda.y0   <- 2.1530
P$lambda.y1   <- 29.0744
P$mu.y        <- -0.005
P$sigma.y     <- 0.0163
P$y0          <- log(1000)

# variance process
P$kappa       <- 5.7808
P$theta       <- 0.0085
P$sigma       <- 0.8935
P$lambda.v0   <- 7.6498
P$mu.v        <- 0.0214
P$v0          <- 0.0118

# Risk premiums
RP <- list()
RP$eta.y       <- 0.7635
RP$gamma.y.min <- 0.0058
RP$eta.v       <- -0.2199
RP$gamma.v.maj <- 0.5811

# Moment generating functions
MGF   <- function(x) exp(P$mu.y * x + 0.5 * P$sigma.y^2 * x^2)
fun   <- function(y) (MGF(1) - 1 ) + MGF(y) - MGF(1+y) - RP$gamma.y.min
RP$gamma.y.maj  <- uniroot(f = fun, interval = c(-50,50))$root

mgf.y  <- MGF(RP$gamma.y.maj)
mgf.y1 <- MGF(1+RP$gamma.y.maj)
mgf.v  <- 1/(1 - RP$gamma.v.maj * P$mu.v)

# Simulation constants
P$c1     <- P$r - P$q + P$lambda.y0 * (mgf.y - mgf.y1) -
  (P$rho * P$kappa * P$theta / P$sigma)
P$c2     <- RP$eta.y - (1/2) + P$lambda.y1 * (mgf.y - mgf.y1) +
  (P$rho * P$kappa ) / P$sigma

# For integrated variance interpolation
P$C1    <- (1 + exp(-P$kappa*S$h))/(1 - exp(-P$kappa*S$h))
P$C2    <- 4*exp(-P$kappa*S$h)/((1 - exp(-P$kappa*S$h))^2)

# Parameters under Q
Q <- list()
# Log process
Q$r         <- P$r
Q$q         <- P$q
Q$rho       <- P$rho
Q$lambda.y0 <- mgf.y * P$lambda.y0 #
Q$lambda.y1 <- mgf.y * P$lambda.y1 #
Q$mu.y      <- P$mu.y + RP$gamma.y.maj * P$sigma.y^2 #
Q$sigma.y   <- P$sigma.y
Q$y0        <- P$y0

# variance process
Q$kappa     <- P$kappa + P$sigma * RP$eta.v #
Q$theta     <- (P$kappa * P$theta)/(Q$kappa + P$sigma * RP$eta.v) #
Q$sigma     <- P$sigma
Q$lambda.v0 <- mgf.v * P$lambda.v0 #
Q$mu.v      <- mgf.v * P$mu.v #
Q$v0        <- P$v0

# Coefficients for Gil Pelaez option pricing
coeff.i     <- GetAB2(Q=Q, u=S$x*1i, maturity=c(30,90))
coeff.i1    <- GetAB2(Q=Q, u=S$x*1i+1, maturity=c(30,90))

A1i <- coeff.i1$Ai
A1r <- coeff.i1$Ar
A2i <- coeff.i$Ai
A2r <- coeff.i$Ar

B1i <- coeff.i1$Bi
B1r <- coeff.i1$Br
B2i <- coeff.i$Bi
B2r <- coeff.i$Br

# For the integrated variance
P$vmin         <- exp(-12)
P$vmax         <- 8 * P$sigma
P$cache.x      <- seq(from = P$vmin, to=P$vmax, length.out = 10000)
P$moments.grid <- GetCacheIntVar(P, P$cache.x, S$h) 