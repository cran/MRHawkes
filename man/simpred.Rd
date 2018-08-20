\name{simpred}
\alias{simpred}
\title{
  Simulate a fitted (bivariate) MRHawkes process model
}
\description{
  Simulate a fitted bivariate MRHawkes process model after the censoring time 
  \code{cens} to a future time point \code{cens.tilde} using the cascading 
  structure of the process.
}
\usage{
  simpred(data, par, cens, cens.tilde, 
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
                      h2.fn(cens.tilde - cens, p.h2)) * 1.1)
}
\arguments{
  \item{data}{
    A two column matrix. The first column contains the event times sorted in 
    ascending order. The second column contains the corresponding event type
    with the label one or two.
  }
  \item{par}{
    A numeric vector. Contains the ten parameters of the model, in order of 
    the immigration parameters \eqn{\mu(.)} for the two renewal distributions, 
    the two offspring parameters \eqn{h(.)} and lastly the four branching 
    ratios \eqn{\eta}.
  }
  \item{cens}{
    A scalar. The censoring time.
  }
  \item{cens.tilde}{
    A scalar. The time that the simulation run uptil. 
  }
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
  \item{B}{
  A numeric scalar. Tuning parameter
  }
  \item{B0}{
  A numeric scalar. Tuning parameter
  }
  \item{pnp1}{
  A numeric square matrix. The joint last immigrant probabilities. 
  }
  \item{max.h1}{
    A numeric scalar. The maximum value of the offspring density for 
    type one events.
  }  
  \item{max.h2}{
    A numeric scalar. The maximum value of the offspring density for 
    for type two events.
  }  
} 

\value{
  A numeric matrix that contains the simulated event times from censoring time 
  \code{cens} up until \code{cens.tilde} and the corresponding event types.
}

\author{
  Tom Stindl <t.stindl@unsw.edu.au>
  Feng Chen <feng.chen@unsw.edu.au>
}


\examples{
  \donttest{
    ## Magnitude 5.5 or greater earthquakes over the 25 year period from 
    ## 01/01/1991 to 31/12/2015.
    data(fivaqks); 
    near.fiji <- grep("Fiji", fivaqks$place)
    near.vanuatu <- grep("Vanuatu", fivaqks$place)
    t.beg <- strptime("1991-01-01 00:00:00", "\%Y-\%m-\%d \%H:\%M:\%S", tz = "UTC")
    t.end <- strptime("2015-12-31 23:59:59", "\%Y-\%m-\%d \%H:\%M:\%S", tz = "UTC")
    t0 <- 0
    t1 <- as.numeric(t.end - t.beg)
    tms <- strptime(fivaqks$time, "\%Y-\%m-\%dT\%H:\%M:\%OSZ", tz = "UTC")
    ts <- as.numeric(tms[-1] - t.beg)
    ts <- c(as.numeric(tms[1] - t.beg)/24, ts)
    ts.fi <- ts[near.fiji]; ts.fi <- ts.fi[ts.fi >= 0 & ts.fi <= t1]
    ts.va <- ts[near.vanuatu]; ts.va <- ts.va[ts.va >=0 & ts.va <= t1]
    ts.c <- c(ts.fi, ts.va)
    z.c <- c(rep(1, times = length(ts.fi)), rep(2, times = length(ts.va)))
    o <- order(ts.c)
    data <- cbind(ts.c[o], z.c[o])
    # simulate future event time based on MLE fitted Rhawkes model
    N <- 100; i <- 0;
    data.pred <- replicate(N, 
                          {cat(i<<-i+1,'\n'); 
                          simpred(data = data,
                                    par = c(0.488, 20.10, 0.347, 9.53, 
                                            461, 720,
                                            0.472, 0.293, 0.399, -0.0774), 
                                            cens = t1, cens.tilde = t1 + 1000)
                                            })
  }
}
\keyword{ predict }
\keyword{ simulations }
\keyword{ forecast }