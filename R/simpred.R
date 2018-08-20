simpred <- 
  function(data, par, cens, cens.tilde, 
           re.dist1 = rweibull, 
           par.redist1 = list(shape = par[1], scale = par[2]),
           re.dist2 = rweibull, 
           par.redist2 = list(shape = par[3], scale = par[4]),
           h1.fn = function(x, p.h1) 1 / p.h1 * exp( - x / p.h1),
           h2.fn = function(x, p.h2) 1 / p.h2 * exp( - x / p.h2),
           p.h1 = par[5], p.h2 = par[6],
           eta11 = par[7], eta12 = par[8], eta21 = par[9], eta22 = par[10],
           B = 10, B0 = 50, pnp1 = NULL,
           max.h1 = max(optimize(h1.fn, c(0, cens.tilde - cens), maximum = TRUE,
                                 p = p.h1)$obj, h1.fn(0, p.h1), 
                        h1.fn(cens.tilde - cens, p.h1)) * 1.1,
           max.h2 = max(optimize(h2.fn, c(0, cens.tilde - cens), maximum = TRUE,
                                 p = p.h2)$obj, h2.fn(0, p.h2), 
                        h2.fn(cens.tilde - cens, p.h2)) * 1.1
  ){ 
  tms <- data[, 1]
  z <- data[, 2]
  n <- length(tms)
  tms1 <- tms[which(data[, 2] == 1)]
  tms2 <- tms[which(data[, 2] == 2)]
  
  ## auxillary functions
  h1 <- function(t) h1.fn(t, p.h1)
  h2 <- function(t) h2.fn(t, p.h2)
  phi1<-function(t)sapply(t,function(s)eta11*sum(h1(s-tms1[tms1<s]))+
                            eta12*sum(h1(s-tms2[tms2<s]))) 
  phi2<-function(t)sapply(t,function(s)eta21*sum(h2(s-tms1[tms1<s]))+
                            eta22*sum(h2(s-tms2[tms2<s]))) 
  
  ## Step 1: Simulate offspring of individuals already in the 
  ## population at the censoring time
  tmp <- simNSMHP(TT = cens.tilde - cens,
                  nu1 = function(t) phi1(cens + t),
                  nu2 = function(t) phi2(cens + t),
                  g11 = function(t) eta11*h1(t),
                  g21 = function(t) eta21*h2(t),
                  g12 = function(t) eta12*h1(t),
                  g22 = function(t) eta22*h2(t))
  off1 <- cens + unlist(tmp$offspr1); off2 <- cens + unlist(tmp$offspr2)
  
  ## Step 2: Simulate the last immigrants from I(T):
  if(length(pnp1) == 0){
    fit.mod <- mllMRH1(data, cens, par)
    pnp1 <- fit.mod$p
  }
  last.im1 <- sample(x = 0:n, 1, prob = rowSums(pnp1))
  last.im2 <- sample(x = 0:n, 1, prob = colSums(pnp1))
  last.im.tm1 <- if(last.im1 == 0){ 0 } else { tms[last.im1] } ; 
  last.im.tm2 <- if(last.im2 == 0){ 0 } else { tms[last.im2] } ; 
  
  ## Step 3: Simulate the next immigrant given the waiting time
  ## is greater than cens - last.im.tm for m = 1,2
  
  ## For m = 1
  ## wt.fni1 <- rweibull(1, shape = par[1], scale = par[2])
  wt.fni1 <- do.call(re.dist1, args = c(n = 1, par.redist1))
  while(wt.fni1 <= cens - last.im.tm1){
    wt.fni1 <- do.call(re.dist1, args = c(n = 1, par.redist1))
  }
  re.tm1 <- last.im.tm1 + wt.fni1
  
  ## For m = 2
  wt.fni2 <- do.call(re.dist2, args = c(n = 1, par.redist2))
  while(wt.fni2 <= cens - last.im.tm2){
    wt.fni2 <- do.call(re.dist2, args = c(n = 1, par.redist2))
  }
  re.tm2 <- last.im.tm2 + wt.fni2
  
  ## Step 4: Simulate the remaining immigrants
  wtms1 <- do.call(re.dist1, args = c(n = B0, par.redist1))
  wtms2 <- do.call(re.dist2, args = c(n = B0, par.redist2))
  n1 <- ceiling(((sqrt(4 * var(wtms1) + 4*mean(wtms1) * (cens.tilde-cens)) - 
                    2 * sd(wtms1)) / (2 * mean(wtms1))) ^ 2)
  n2 <- ceiling(((sqrt(4 * var(wtms2) + 4*mean(wtms2) * (cens.tilde-cens)) - 
                    2 * sd(wtms2)) / (2 * mean(wtms2))) ^ 2)
  wtms1 <- c(wtms1, do.call(re.dist1, args = c(n = n1, par.redist1)))
  wtms2 <- c(wtms2, do.call(re.dist2, args = c(n = n2, par.redist2)))
  
  last.i.tm1 <- re.tm1 + sum(wtms1)
  while(last.i.tm1 <= cens.tilde){
    tmp1 <- do.call(re.dist1, args = list(n = B, unlist(par.redist1)))
    wtms1 <- c(wtms1, tmp1)
    last.i.tm1 <- last.i.tm1 + sum(tmp1)
  }
  tmp1 <- re.tm1 + cumsum(wtms1)
  i.tms1 <- c(re.tm1, tmp1)
  i.tms1 <- i.tms1[i.tms1 < cens.tilde]
  
  last.i.tm2 <- re.tm2 + sum(wtms2)
  while(last.i.tm2 <= cens.tilde){
    tmp2 <- do.call(re.dist2, args = list(n = B, unlist(par.redist2)))
    wtms2 <- c(wtms2, tmp2)
    last.i.tm2 <- last.i.tm2 + sum(tmp2)
  }
  tmp2 <- re.tm2 + cumsum(wtms2)
  i.tms2 <- c(re.tm2, tmp2)
  i.tms2 <- i.tms2[i.tms2 < cens.tilde]
  
  ## Step 5: Simulate the offspring for each immigrant simmulated above
  ## function to simulate the offspring from an immigrant given its 
  ## event type
  ## Now simulate the offspring for every immigrant
  n1 <- length(i.tms1); n2 <- length(i.tms2)
  off11 <- off21 <- vector('list', n1)
  off12 <- off22 <- vector('list', n2)
  if(n1 > 0){
    for(i in 1:n1){
      sim <- simNSMHP(TT = cens.tilde - i.tms1[i],
                      nu1 = function(t) eta11*h1.fn(t, p.h1),
                      nu2 = function(t) eta21*h2.fn(t, p.h2),
                      g11 = function(t) eta11*h1.fn(t, p.h1),
                      g21 = function(t) eta21*h2.fn(t, p.h2),
                      g12 = function(t) eta12*h1.fn(t, p.h1),
                      g22 = function(t) eta22*h2.fn(t, p.h2))
      off11[[i]] <- i.tms1[i] + unlist(sim$offspr1)
      off21[[i]] <- i.tms1[i] + unlist(sim$offspr2)
    }
  }
  if(n2 > 0){
    for(i in 1:n2){
      sim <- simNSMHP(TT = cens.tilde - i.tms2[i],
                      nu1 = function(t) eta12*h1.fn(t, p.h1),
                      nu2 = function(t) eta22*h2.fn(t, p.h2),
                      g11 = function(t) eta11*h1.fn(t, p.h1),
                      g21 = function(t) eta21*h2.fn(t, p.h2),
                      g12 = function(t) eta12*h1.fn(t, p.h1),
                      g22 = function(t) eta22*h2.fn(t, p.h2))
      off12[[i]] <- i.tms2[i] + unlist(sim$offspr1)
      off22[[i]] <- i.tms2[i] + unlist(sim$offspr2)
    }
  }
  
  tms <- c(off1, i.tms1, unlist(off11), unlist(off12),
           off2, i.tms2, unlist(off21), unlist(off22))
  ones <- n1 + length(off1) + length(unlist(off11)) +
    length(unlist(off12))
  twos <- n2 + length(off2) + length(unlist(off21)) + 
    length(unlist(off22))
  zs <- c(rep(1, times = ones), rep(2, time = twos)) 
  o <- order(tms)
  cbind(tms[o], zs[o])
}
