\name{simNSMHP}
\alias{simNSMHP}
\title{
  Simulate a (bivariate) non-stationary multivariate Hawkes process (NSMHP)
}
\description{
  Simulate a bivariate non-stationary multivariate Hawkes process (NSMHP) with 
  given given baseline intensity functions and self-excitation functions using 
  the cascading structure of the process.
}
\usage{
simNSMHP(TT = 100,
        nu1 = function(t) 0.6*exp(-t),
        nu2 = function(t) 0.2*exp(-t),
        g11 = function(t) 0.6*exp(-t),
        g12 = function(t) 0.2*exp(-t),
        g21 = function(t) 0.1*exp(-t),
        g22 = function(t) 0.5*exp(-t))
}
\arguments{
  \item{TT}{
    A scalar. The censoring time. 
  }
 \item{nu1}{
    Basline intensity function for type one events.
  }
  \item{nu2}{
    Basline intensity function for type two events.
  }
  \item{g11}{
  Self-exciting function for type one events given the parent is a type two 
  event.
  }
  \item{g12}{
  Cross-exciting function for type one events given the parent is a type two
  event.
  }
  \item{g21}{
  Cross-exciting function for type two events given the parent is a type one
  event.
  }
  \item{g22}{
  Self-exciting function for type two events given the parent is a type two
  event.
  }
}

\details{
  The function works by simulating generation 0 events according to independent
  Poisson processes with the baseline intensity functions; then keep simulating 
  future generation events as long as the number of the previous generation 
  events of any type is non-zero. For each event type, we simulate these events
  according to M independent Poisson processes with the appropriate excitation 
  intensity. When this recursive process stops, return events of all generations 
  with their respective type labels as the events of the NSMHP on the interval 
  (0,T].
}

\value{
  \item{offspr1}{All offspring events of type one}
  \item{offspr2}{All offspring events of type two}
}

\author{
  Tom Stindl <t.stindl@unsw.edu.au>
  Feng Chen <feng.chen@unsw.edu.au>
}


\examples{
  B <- 10; i <- 0;
  data <- replicate(B, 
                    {cat(i<<-i+1,'\n'); 
                    simNSMHP(TT = 100)
                    })
}
\keyword{self-exciting}
