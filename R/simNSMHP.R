simNSMHP <- 
  function(TT = 100,
           nu1 = function(t) 0.6*exp(-t),
           nu2 = function(t) 0.2*exp(-t),
           g11 = function(t) 0.6*exp(-t),
           g12 = function(t) 0.2*exp(-t),
           g21 = function(t) 0.1*exp(-t),
           g22 = function(t) 0.5*exp(-t)
  ){
  gen <- 1
  if(nu1(0) > 0) o1 <- simPois(int = nu1, cens = TT) else o1 <- numeric(0)
  if(nu2(0) > 0) o2 <- simPois(int = nu2, cens = TT) else o2 <- numeric(0)
  
  offtms1 <- offtms2 <- list()
  offtms1[[gen]] <- unlist(o1); offtms2[[gen]] <- unlist(o2)
  l1 <- length(offtms1[[gen]]); l2 <- length(offtms2[[gen]])
  
  while(l1 > 0 | l2 > 0){
    if(g11(0) > 0 & l1 > 0){
      o11 <- sapply(offtms1[[gen]],
                    function(t) t + simPois(int = g11, cens = TT - t))
    } else o11 <- numeric(0)
    if(g21(0) > 0 & l1 > 0){
      o21 <- sapply(offtms1[[gen]],
                    function(t) t + simPois(int = g21, cens = TT - t))
    } else o21 <- numeric(0)
    if(g12(0) > 0 & l2 > 0){
      o12 <- sapply(offtms2[[gen]],
                    function(t) t + simPois(int = g12, cens = TT - t))
    } else o12 <- numeric(0)
    if(g22(0) > 0 & l2 > 0){
      o22 <- sapply(offtms2[[gen]],
                    function(t) t + simPois(int = g22, cens = TT - t))
    } else o22 <- numeric(0)
    
    gen <- gen + 1
    offtms1[[gen]] <- c(unlist(o11), unlist(o12))
    offtms2[[gen]] <- c(unlist(o21), unlist(o22))
    l1 <- length(offtms1[[gen]])
    l2 <- length(offtms2[[gen]])
  }
  list(offspr1 = offtms1, offspr2 = offtms2)
} 

