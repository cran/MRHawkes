simMRHawkes <- 
  function(re.dist1 = rweibull, par.redist1 = list(shape = 3, scale = 1.2), 
           re.dist2 = rweibull, par.redist2 = list(shape = 1/3, scale = 0.2), 
           h1.fn = function(x, p.h1) 1/p.h1 * exp(-x/p.h1), 
           h2.fn = function(x, p.h2) 1/p.h2 * exp(-x/p.h2), 
           p.h1 = 1, p.h2 = 1, 
           eta11 = 0.3, eta12 = 0.1, eta21 = 0.1, eta22 = 0.3, cens = 100, 
           B = 10, B0 = 50, 
           max.h1 = max(optimize(h1.fn, c(0, cens), maximum = TRUE, 
                                 p = p.h1)$obj, 
                        h1.fn(0, p.h1), h1.fn(cens, p.h1)) * 1.1, 
           max.h2 = max(optimize(h2.fn, c(0, cens), maximum = TRUE, 
                                 p = p.h2)$obj, 
                        h2.fn(0, p.h2), h2.fn(cens, p.h2)) * 1.1
  ){
  wtms1 <- do.call(re.dist1, args = c(n = B0, par.redist1))
  wtms2 <- do.call(re.dist2, args = c(n = B0, par.redist2))
  n1 <- ceiling(((sqrt(4 * var(wtms1) + 4 * mean(wtms1) * cens) - 
                    2 * sd(wtms1))/(2 * mean(wtms1)))^2)
  n2 <- ceiling(((sqrt(4 * var(wtms2) + 4 * mean(wtms2) * cens) - 
                    2 * sd(wtms2))/(2 * mean(wtms2)))^2)
  wtms1 <- c(wtms1, do.call(re.dist1, args = c(n = n1, par.redist1)))
  wtms2 <- c(wtms2, do.call(re.dist2, args = c(n = n2, par.redist2)))
  last.i.tm1 <- sum(wtms1)
  while (last.i.tm1 <= cens) {
    tmp1 <- do.call(re.dist1, args = list(n = B, unlist(par.redist1)))
    wtms1 <- c(wtms1, tmp1)
    last.i.tm1 <- last.i.tm1 + sum(tmp1)
  }
  tmp1 <- cumsum(wtms1)
  i.tms1 <- tmp1[tmp1 <= cens]
  last.i.tm2 <- sum(wtms2)
  while (last.i.tm2 <= cens) {
    tmp2 <- do.call(re.dist2, args = list(n = B, unlist(par.redist2)))
    wtms2 <- c(wtms2, tmp2)
    last.i.tm2 <- last.i.tm2 + sum(tmp2)
  }
  tmp2 <- cumsum(wtms2)
  i.tms2 <- tmp2[tmp2 <= cens]
  n1 <- length(i.tms1)
  n2 <- length(i.tms2)
  
  ## Now simulate the offspring for every immigrant
  off11 <- off21 <- vector('list', n1)
  off12 <- off22 <- vector('list', n2)
  if(n1 > 0){
    for(i in 1:n1){
      sim <- simNSMHP(TT = cens - i.tms1[i],
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
      sim <- simNSMHP(TT = cens - i.tms2[i],
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
  
  ## combine the times and labels
  tms <- c(i.tms1, unlist(off11), unlist(off12), i.tms2, 
           unlist(off21), unlist(off22))
  zs <- c(rep(1, times = n1 + length(c(unlist(off11), unlist(off12)))), 
          rep(2, times = n2 + length(c(unlist(off21), unlist(off22)))))
  o <- order(tms)
  cbind(tms[o], zs[o])
}