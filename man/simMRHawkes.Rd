\name{simMRHawkes}
\alias{simMRHawkes}
\title{
  Simulate an (bivariate) renewal Hawkes (MRHawkes) process
}
\description{
  Simulate an (bivairate) renewal Hawkes (MRHawkes) process with given renewal 
  immigration distribution functions \eqn{\mu}, offspring density functions 
  \eqn{h} and branching ratios \eqn{\eta} using the cascading structure of the 
  process.
}
\usage{
simMRHawkes(re.dist1 = rweibull, par.redist1 = list(shape = 3, scale = 1.2), 
            re.dist2 = rweibull, par.redist2 = list(shape = 1/3, scale = 0.2), 
            h1.fn = function(x, p.h1) 1/p.h1 * exp(-x/p.h1), 
            h2.fn = function(x, p.h2) 1/p.h2 * exp(-x/p.h2), 
            p.h1 = 1, p.h2 = 1, 
            eta11 = 0.3, eta12 = 0.1, eta21 = 0.1, eta22 = 0.3, cens = 100, 
            B = 10, B0 = 50, 
            max.h1 = max(optimize(h1.fn, c(0, cens), maximum = TRUE, p = p.h1)$obj, 
                         h1.fn(0, p.h1), h1.fn(cens, p.h1)) * 1.1, 
            max.h2 = max(optimize(h2.fn, c(0, cens), maximum = TRUE, p = p.h2)$obj, 
                         h2.fn(0, p.h2), h2.fn(cens, p.h2)) * 1.1)
}
\arguments{
 \item{re.dist1}{
    The renewal distribution for type one events.
  }
  \item{re.dist2}{
    The renewal distribution for type two events.
  }  
  \item{par.redist1}{
    A numeric list. The parameters of the renewal distribution for type one 
    events.
  } 
  \item{par.redist2}{
    A numeric list. The parameters of the renewal distribution for type two 
    events.
  }   
  \item{h1.fn}{
    A (vectorized) function. The offspring density function for type one 
    events.    
  }
  \item{h2.fn}{
    A (vectorized) function. The offspring density function for type two 
    events.
  }   
  \item{p.h1}{
    A numeric vector. The paramters of the offspring density for type one 
    events.
  }  
  \item{p.h2}{
    A numeric vector. The paramters of the offspring density for type two 
    events.
  }    
 \item{eta11}{
  A numeric scalar. The self-exciting branching ratio for type one events.
  }
  \item{eta12}{
  A numeric scalar. The cross-exciting branching ratio for type one events due 
  to the effects of a type two event.
  }
  \item{eta21}{
  A numeric scalar. The cross-exciting branching ratio for type two events due 
  to the effects of a type one event.
  }
  \item{eta22}{
  A numeric scalar. The self-exciting branching ratio for type two events.
  }    
  \item{cens}{
    A scalar. The censoring time. 
  }
  \item{B}{
  A numeric scalar. Tuning parameter
  }
  \item{B0}{
  A numeric scalar. Tuning parameter
  } 
  \item{max.h1}{
    A numeric scalar. The maximum value of the offspring density for 
    type one events.
  }  
  \item{max.h2}{
    A numeric scalar. The maximum value of the offspring density for 
    type two events.
  }   
 
}
\details{
  The function works by simulating the arrival times of immigrants
  accoridng to the respective renewal immigration distribution for each event 
  type. The birth times ofoffspring from each immigrant are then simulated 
  according to an non-stationary multivariate Hawkes Process (NSMHP) with 
  appropriate baseline and excitation functions.
}
\value{
  A numeric matrix. The row coloumn contains the event times in ascending 
  order while the second coloumn contains the corresponding event type.
}

\author{
  Tom Stindl <t.stindl@unsw.edu.au>
  Feng Chen <feng.chen@unsw.edu.au>
}

\examples{
  B <- 10; i <- 0;
  data <- replicate(B, 
                    {cat(i<<-i+1,'\n'); 
                    simMRHawkes(re.dist1 = rweibull, 
                                par.redist1 = list(shape = 3, scale = 1.2),
                                re.dist2 = rweibull, 
                                par.redist2 = list(shape = 1 / 3, scale = 0.2),
                                p.h1 = 1, p.h2 = 1,
                                eta11 = 0.3, eta12 = 0.1, 
                                eta21 = 0.1, eta22 = 0.3,
                                cens = 100)
                    })
}
\keyword{renewal process}
\keyword{self-exciting}
