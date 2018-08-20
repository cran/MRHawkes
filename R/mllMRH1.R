mllMRH1 <- 
  function(data, cens, par,
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
           }){
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
  
  if(n == 0) return(exp( - Mu1.fn(cens, p.mu1) - Mu2.fn(cens, p.mu2)))
  if(n == 1) return( if(z[1] == 1) { 
    exp( - Mu1.fn(tms[1], p.mu1) - Mu2.fn(tms[1], p.mu2)) * 
      mu1.fn(tms[1], p.mu1) * 
      exp( - Mu1.fn(cens - tms[1], p.mu1) - Mu2.fn(cens, p.mu2) + 
             Mu2.fn(tms[1], p.mu2))
  } else {
    exp( - Mu1.fn(tms[1], p.mu1) - Mu2.fn(tms[1], p.mu2)) * 
      mu2.fn(tms[1], p.mu2) * 
      exp( - Mu1.fn(cens, p.mu1) + Mu1.fn(tms[1], p.mu1) - 
             Mu2.fn(cens - tms[1], p.mu2))
  })
  
  ## auxiliary functions
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
  dPhi <- function(s, t){
    eta11 * sum(H1(s - tms1[tms1 < s])) + 
      eta12 * sum(H1(s - tms2[tms2 < s])) + 
      eta21 * sum(H2(s - tms1[tms1 < s])) +
      eta22 * sum(H2(s - tms2[tms2 < s])) -
      eta11 * sum(H1(t - tms1[tms1 < t])) - 
      eta12 * sum(H1(t - tms2[tms2 < t])) - 
      eta21 * sum(H2(t - tms1[tms1 < t])) -
      eta22 * sum(H2(t - tms2[tms2 < t]))
  }
  
  ## vectors to hold p(tau_i, z_i | tau_{1:i-1}, z_{1:i-1}),
  ## matrix to hold d_i(j,k), d_{i-1}(j,k), p_i(j,k), p_{i-1}(j,k)
  den <- numeric(n + 1);
  con.den <- con.denmo <- pi <- pimo <- matrix(0, nrow = n + 1, ncol = n + 1);
  
  den[1] <- exp( - Mu1(tms[1]) - Mu2(tms[1])) * 
    if(z[1] == 1) mu1(tms[1]) else mu2(tms[1]) ; 
  
  i <- 2;
  while(i <= n + 1){
    if(i == 2){
      pi[2, 1] <- if(z[1] == 1) 1 else 0 #p_2(1,0)
      pi[1, 2] <- if(z[1] == 2) 1 else 0 #p_2(0,1)
      if(i <= n){
        dPh <- dPhi(tms[2], tms[1])
        con.den[2, 1] <- exp( - Mu1(tms[2] - tms[1]) - Mu2(tms[2]) + 
                                Mu2(tms[1]) - dPh) * 
          if(z[2] == 1){
            mu1(tms[2] - tms[1]) + phi1(tms[2])
          }else{
            mu2(tms[2]) + phi2(tms[2])
          }
        con.den[1, 2] <- exp( - Mu1(tms[2]) - Mu2(tms[2] - tms[1]) + 
                                Mu1(tms[1]) - dPh) * 
          if(z[2] == 1){
            mu1(tms[2]) + phi1(tms[2])
          }else{
            mu2(tms[2] - tms[1]) + phi2(tms[2])
          }
      } else{ 
        dPh <- dPhi(cens, tms[1]) 
        con.den[2, 1] <- exp( - Mu1(cens - tms[1]) - 
                                Mu2(cens) + Mu2(tms[1]) - dPh) 
        con.den[1, 2] <- exp( - Mu1(cens) + Mu1(tms[1]) - 
                                Mu2(cens - tms[1]) - dPh) }
    } else{
      if(z[i - 1] == 1){
        ph1 <- phi1(tms[i - 1])
        m1 <- mu1(tms[i - 1] - c(0, tms[1:(i - 2)]))
        tmp <- ph1 / (m1 + ph1)
        mat <- con.denmo[1:(i - 1),1:(i - 1)] * pimo[1:(i - 1), 1:(i - 1)] / 
          den[i - 1]
        pi[1:(i - 1), 1:(i - 1)] <- mat1 <- tmp * mat 
        pi[i, 1:(i - 1)] <- colSums(mat-mat1)
      }
      
      if(z[i - 1] == 2){
        ph2 <- phi2(tms[i - 1])
        m2 <- mu2(tms[i - 1] - c(0, tms[1:(i - 2)]))
        tmp <- ph2 / (m2 + ph2) 
        mat <- t(con.denmo[1:(i - 1), 1:(i - 1)] * pimo[1:(i - 1), 1:(i - 1)]) / 
          den[i-1]
        pi[1:(i - 1), 1:(i - 1)] <- mat1 <- t(tmp * mat)
        pi[1:(i - 1), i] <- rowSums(t(mat - t(mat1))) 
      }
      
      if(i <= n){ 
        instant <- matrix(0, nrow = i - 1, ncol = i - 1) 
        if(z[i] == 1){ 
          instant[1:(i - 1), ] <- mu1(tms[i] - tms[1:(i - 1)]) + 
            phi1(tms[i])
        }else{ 
          instant[, 1:(i - 1)] <- mu2(tms[i] - tms[1:(i - 1)]) + 
            phi2(tms[i]); 
          instant <- t(instant)
        }
        dPh <- dPhi(tms[i], tms[i - 1])
        con.den[2:i, 2:i] <- 
          outer(exp( - Mu1(tms[i] - tms[1:(i-1)]) + 
                       Mu1(tms[i - 1] - tms[1:(i - 1)])), 
                exp( - Mu2(tms[i] - tms[1:(i-1)]) + 
                       Mu2(tms[i - 1] - tms[1:(i - 1)]))) * 
          exp( - dPh) * instant
        
        instant <- numeric(i - 1)
        if(z[i] == 1){
          instant[1:(i - 1)] <- mu1(tms[i]) + phi1(tms[i])
        }else{
          instant[1:(i - 1)] <- mu2(tms[i] - tms[1:(i - 1)]) + 
            phi2(tms[i])
        }  
        con.den[1, 2:i] <- exp( - Mu1(tms[i]) - 
                                  Mu2(tms[i] - tms[1:(i - 1)]) + 
                                  Mu1(tms[i - 1]) + 
                                  Mu2(tms[i - 1] - tms[1:(i - 1)]) - 
                                  dPh) * instant 
        if(z[i] == 1) { 
          instant[1:(i - 1)] <- mu1(tms[i] - tms[1:(i - 1)]) + 
            phi1(tms[i]) 
        } else{ 
          instant[1:(i - 1)] <- mu2(tms[i]) + phi2(tms[i]) 
        } 
        con.den[2:i, 1] <- exp( - Mu1(tms[i] - tms[1:(i - 1)]) - 
                                  Mu2(tms[i]) + 
                                  Mu1(tms[i - 1] - tms[1:(i - 1)]) + 
                                  Mu2(tms[i - 1]) - dPh) * 
          instant 
        diag(con.den) <- 0
      }else{
        dPh <- dPhi(cens, tms[n]) 
        con.den[2:(n + 1), 2:(n + 1)] <- 
          outer(exp( - Mu1(cens - tms[1:n]) + Mu1(tms[n] - tms[1:n])), 
                exp(- Mu2(cens - tms[1:n]) + Mu2(tms[n] - tms[1:n]))) * 
          exp( - dPh) 
        con.den[1, 2:i] <- 
          exp( - Mu1(cens) - Mu2(cens - tms[1:(i - 1)]) + 
                 Mu1(tms[n]) + Mu2(tms[n] - tms[1:(i - 1)]) - dPh) 
        con.den[2:i, 1] <- 
          exp( - Mu1(cens - tms[1:(i - 1)]) - Mu2(cens) + 
                 Mu1(tms[n] - tms[1:(i - 1)]) + Mu2(tms[n]) - dPh) 
        diag(con.den) <- 0
      }
    }
    
    tmp <- con.den[1:i, 1:i] * pi[1:i, 1:i]
    den[i] <- sum(tmp)
    pimo[1:i, 1:i] <- pi[1:i, 1:i]
    con.denmo[1:i, 1:i] <- con.den[1:i, 1:i] 
    i <- i + 1
  }
  mll <- - sum(log(den[1:(n + 1)]))
  pi <- pi # store the last immigrant probabilities at the censoring time
  list(mll = mll, p = pi, n = n)
}