% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optionstrat.R
\name{vertical}
\alias{vertical}
\title{Vertical Spread Analytics}
\usage{
vertical(options = c("call", "put"), s, x1, x2, t, r, sigma,
  sigma2 = sigma, vol = sigma, d = 0)
}
\arguments{
\item{options}{Character string. Either "call", or "put"}

\item{s}{Spot price of the underlying asset}

\item{x1}{Strike price of the short option}

\item{x2}{Strike price of the long option}

\item{t}{Time to expiration in years}

\item{r}{Annual continuously compounded risk-free rate}

\item{sigma}{Implied volatility of the short option (annualized)}

\item{sigma2}{Implied volatility of the long option (annualized)}

\item{vol}{Manual over-ride for the volatility of the underlying asset (annualized)}

\item{d}{Annual continuously compounded dividend yield}
}
\value{
Returns a data.frame
}
\description{
Calculates the key analytics of a vertical spread
}
\examples{
vertical("call", s = 100, x1 = 90, x2 = 110, t = (45/365), r =  0.025, sigma = 0.20, vol = 0.25)
}
