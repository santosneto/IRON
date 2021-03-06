% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RCodes.R
\name{IRON}
\alias{IRON}
\alias{diron}
\alias{qiron}
\alias{piron}
\alias{riron}
\title{The IRON distribution.}
\usage{
diron(x, family = "norm", tau, shape, scale, kappa = NULL, log = FALSE, ...)

qiron(p, shape, scale, tau, family = "norm", ...)

piron(q, shape, scale, tau, family = "norm", log.p = FALSE, ...)

riron(n, shape, scale, tau, family = "normal", df = 4, kappa = NULL, ...)
}
\arguments{
\item{x, q}{vetor of quantiles.}

\item{family}{kernel used.}

\item{tau}{quantile value.}

\item{shape}{shape parameter.}

\item{scale}{scale parameter.}

\item{kappa}{shape parameter.}

\item{log, log.p}{logical; if TRUE, the log-density or log(p) is used.}

\item{...}{additional arguments to be passed.}

\item{p}{vector of probabilities.}

\item{n}{number of observations.}

\item{df}{degrees of freedom.}
}
\value{
diron gives the density, piron gives the distribution function, qiron gives the quantile function and riron generate pseudo-random numbers.
}
\description{
Density function, distribution function, quantile function and random generation for the IRON distribution.
}
\details{
The density function of the IRON distribution is

\deqn{g(t|\lambda,\beta,\alpha) = \alpha f(a_t)[F(a_t)]^{\alpha-1} \vdot \frac{t^{-3/2}(t+\beta)}{2\lambda\sqrt{\beta}}\qc t,\lambda,\beta,\alpha,}
where \eqn{f(.)} and \eqn{F(.)} are, respectively, the density function and cumulative distribution of a symmetric distribution.
}
\examples{
diron(x = 10, tau = 0.5, shape = 1, scale = 1) #Birnbaum-Saunders distribution
piron(q = 10, tau = 0.5, shape = 1, scale = 1) #Birnbaum-Saunders distribution
qiron(p = 0.5, tau = 0.5, shape = 1, scale = 1) #Birnbaum-Saunders distribution
riron(n = 10, tau = 0.5, shape = 1, scale = 1) #Birnbaum-Saunders distribution
}
\author{
Manoel Santos-Neto \email{manoel.ferreira at professor.ufcg.edu.br} and Diego I. Gallardo \email{diego.gallardo.mateluna at gmail.com}
}
