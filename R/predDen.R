predDen <- 
  function(x, data, cens, par, 
           h1.fn = function(x, p) 1 / p * exp( - x / p),
           h2.fn = function(x, p) 1 / p * exp( - x / p),
           mu1.fn = function(x, p){
             exp(dweibull(x, shape = p[1], scale = p[2], log = TRUE) -
                   pweibull(x, shape = p[1], scale = p[2], lower.tail = FALSE, 
                            log.p = TRUE))
           },
           mu2.fn = function(x, p){
             exp(dweibull(x, shape = p[1], scale = p[2], log = TRUE) -
                   pweibull(x, shape = p[1], scale = p[2], lower.tail = FALSE, 
                            log.p = TRUE))
           },
           H1.fn = function(x, p) pexp(x, rate = 1 / p),
           H2.fn = function(x, p) pexp(x, rate = 1 / p),
           Mu1.fn = function(x, p){
             - pweibull(x, shape = p[1], scale = p[2], lower.tail = FALSE, 
                        log.p = TRUE)
           },
           Mu2.fn = function(x, p){
             - pweibull(x, shape = p[1], scale = p[2], lower.tail = FALSE, 
                        log.p = TRUE)
           }) {
  tms <- data[, 1]
  z <- data[, 2]
  n <- length(tms)
  tms1 <- tms[which(data[, 2] == 1)]
  tms2 <- tms[which(data[, 2] == 2)]
  
  p.mu1 <- par[1:2]
  p.mu2 <- par[3:4]
  p.h1 <- par[5]
  p.h2 <- par[6]
  eta11 <- par[7]
  eta12 <- par[8]
  eta21 <- par[9]
  eta22 <- par[10]
  
  ##auxiliary functions
  mu1 <- function(t) mu1.fn(t, p.mu1)
  mu2 <- function(t) mu2.fn(t, p.mu2)
  Mu1 <- function(t) Mu1.fn(t, p.mu1)
  Mu2 <- function(t) Mu2.fn(t, p.mu2)
  h1 <- function(t) h1.fn(t, p.h1)
  h2 <- function(t) h2.fn(t, p.h2)
  H1 <- function(t) H1.fn(t, p.h1)
  H2 <- function(t) H2.fn(t, p.h2)
  phi1 <- function(s) {
    eta11 * sum(h1(s - tms1[tms1 < s])) + 
      eta12 * sum(h1(s - tms2[tms2 < s]))
  }
  phi2 <- function(s) { 
    eta21 * sum(h2(s - tms1[tms1 < s])) + 
      eta22 * sum(h2(s - tms2[tms2 < s]))
  }  
  Phi <- function(s) {
    eta11 * sum(H1(s - tms1[tms1 < s])) + 
      eta12 * sum(H1(s - tms2[tms2 < s])) + 
      eta21 * sum(H2(s - tms1[tms1 < s])) + 
      eta22 * sum(H2(s - tms2[tms2 < s]))
  }
  
  Snp1 <- outer(exp( - Mu1(cens - c(0, tms)) + Mu1(tms[n] - c(0, tms))), 
                exp( - Mu2(cens - c(0, tms)) + Mu2(tms[n] - c(0, tms)))) * 
    exp( - Phi(cens) + Phi(tms[n])) 
  
  tmp <- mllMRH1(data, cens, par)
  pnp1 <- tmp$p
  den <- sum(Snp1 * pnp1)
  
  sapply(x, function(s){
    dnp1 <- 
      exp( - outer(Mu1(cens + s - c(0, tms)) - 
                     Mu1(tms[n] - c(0, tms)), 
                   Mu2(cens + s - c(0, tms)) - 
                     Mu2(tms[n] - c(0, tms)),
                   "+")
      )*
      exp( - Phi(cens + s) + Phi(tms[n])) *
      outer(mu1(cens + s - c(0, tms)) + phi1(cens + s), 
            mu2(cens + s - c(0, tms)) + phi2(cens + s), "+")
    
    sum( pnp1 * dnp1 / den) 
  })
}