\name{TmllMRH}
\alias{TmllMRH}
\title{
  Minus loglikelihood of an (bivariate) MRHawkes model with truncated most 
  recent immigrant probabilities
}
\description{
  Calculates the minus loglikelihood of an (bivariate) MRHawkes model with 
  given immigration hazard functions \eqn{\mu}, offspring density functions 
  \eqn{h} and bracnhing ratios \eqn{\eta} for event times and event types 
  \code{data} on interval \eqn{[0,cens]} with truncated most recent immigrant 
  probabilities looking back at only the previous \code{B} event times.
}
\usage{
TmllMRH(data, cens, par, B = 25,
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
        })
}

\arguments{
  \item{data}{
    A two column matrix. The first column contains the event times sorted in 
    ascending order. The second column contains the corresponding event type
    with the label one or two.
  }
  \item{cens}{
    A scalar. The censoring time.
  }
  \item{par}{
    A numeric vector. Contains the ten parameters of the model, in order of 
    the immigration parameters \eqn{\mu(.)} for the two renewal distributions, 
    the two offspring parameters \eqn{h(.)} and lastly the four branching 
    ratios \eqn{\eta}.
  }
  \item{B}{
    A scalar. The number of previous events that are considered to be possible
    last immigrant arrivals.
  }
  \item{h1.fn}{
    A (vectorized) function. The offspring density function for type one 
    events. 
  }
    \item{h2.fn}{
    A (vectorized) function. The offspring density function for type two 
    events. 
  }
  \item{mu1.fn}{
    A (vectorized) function. The immigration hazard function for events of type
    one. 
  }
  \item{mu2.fn}{
    A (vectorized) function. The immigration hazard function for events of type 
    two. 
  }
  \item{H1.fn}{
    A (vectorized) function. Its value at \code{t} gives the integral of
    the offspring density function from 0 to \code{t} for type one events. 
  }
  \item{H2.fn}{
    A (vectorized) function. Its value at \code{t} gives the integral of
    the offspring density function from 0 to \code{t} for type two events. 
  }
  \item{Mu1.fn}{
    A (vectorized) function. Its value at \code{t} gives the integral of
    the immigrant hazard function from 0 to \code{t} for type one events. 
  }
  \item{Mu2.fn}{
    A (vectorized) function. Its value at \code{t} gives the integral of
    the immigrant hazard function from 0 to \code{t} for type two events. 
  }
}

\value{
  Value of the negative log-liklihood. 
}

\author{
  Tom Stindl <t.stindl@unsw.edu.au>
  Feng Chen <feng.chen@unsw.edu.au>
}

\examples{
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
  
  ## calculate the minus loglikelihood of an (bivariate) MRHawkes with some 
  ## parameters the default hazard functions and density functions are Weibull 
  ## and exponential respectivley
  TmllMRH(data, cens = t1, par = c(0.488, 20.10, 0.347, 9.53, 461, 720,
                                  0.472, 0.293, 0.399, -0.0774), B = 200)
                                  
  ## calculate the MLE for the parameter assuming known parametric forms
  ## of the immigrant hazard function and offspring density functions.  
  \donttest{
    system.time(est <- optim(c(0.488, 20.10, 0.347, 9.53, 461, 720,
                                  0.472, 0.293, 0.399, -0.0774), 
                             TmllMRH, data = data, cens = t1, B = 200,
                             control = list(maxit = 5000, trace = TRUE),
                             hessian = TRUE)
    )
    ## point estimate by MLE
    est$par
    ## standard error estimates:
    diag(solve(est$hessian))^0.5
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ MRHawkes }
\keyword{ likelihood }
\keyword{ MLE }