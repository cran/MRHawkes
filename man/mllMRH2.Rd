\name{mllMRH2}
\alias{mllMRH2}
\title{
  Minus loglikelihood of an (bivariate) MRHawkes model with Rosenblatt 
  residuals
  }
  
\description{
  Calculates the minus loglikelihood of an (bivariate) RHawkes model with 
  given immigration hazard functions \eqn{\mu}, common offspring density 
  functions \eqn{h} and bracnhing ratios \eqn{\eta} for event times and 
  event types \code{data} on interval \eqn{[0,cens]}. The same as 
  \code{mllMRH} although this version also returns the Rosenblatt residuals 
  for goodness-of-fit assessment of the event times.
}

\usage{
mllMRH2(data, cens, par,
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

\details{
  Calculate the MRHawkes point process Rosenblatt residuals
}

\value{
  \item{mll}{minus log-likelihood}
  \item{W}{Rosenblatt residuals of observed event times}
}

\author{
  Tom Stindl <t.stindl@unsw.edu.au>
  Feng Chen <feng.chen@unsw.edu.au>
}

\seealso{
  \code{mllMRH}}
  
\examples{
\donttest{
  n <- 1000
  data <- cbind(sort(runif(n,0,1000)), 
                sample(1:2, size = n, replace = TRUE))
  tmp <- mllMRH2(data = data, cens = 1001, 
                 par = c(1,1,1,1,1,1,0.5,0.2,0.2,0.3))              
  pp <- ppoints(n)
  par(mfrow=c(1,2))
  plot(quantile(tmp$W,prob=pp),pp,type="l",
       main="Uniform QQ plot",
       xlab="Sample quantiles",ylab="Theoretical quantiles")
  abline(a = 0, b = 1, col = 2)
  a <- acf(tmp$W, main = "ACF Plot")
  ks.test(tmp$W,"punif")
  Box.test(tmp$W,lag=tail(a$lag,1))
  }
}

\keyword{ residual }
\keyword{ point process }
