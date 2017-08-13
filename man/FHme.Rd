\name{FHme}
\alias{FHme}
\title{
    Small area estimation when the covariates subject to sampling error
}
\description{
    - Calculate the estimated value of Y with auxiliary variable X where the auxiliary variable X also measured with error\cr
    - Calculate the MSE of estimated Y with jackknife  method
}
\usage{
    FHme(y, x, vardir, C, iter)
}
\arguments{
    \item{y}{vector of direct estimation for each area}
    \item{x}{matrix of auxiliary variable from another survey. It contains two or more columns. first column is vector of 1 (intercept). It should be like \code{cbind(1, aux_var)}. See \strong{example}.}
    \item{vardir}{vector of design variance of variable y}
    \item{C}{matrix of variance of variable x. It has same column number like \var{x}. The first column is vector of 0 (intercept variance). It should be like \code{cbind(0, variance_x)}. See \strong{example}}
    \item{iter}{Number of iteration for finding w. The dafault is 2. See \strong{reference} for more details.}
}

\value{
    \item{y_me}{small area estimation for each area}
    \item{gamma}{estimation of gamma for each area}
    \item{sigma}{estimation of sigma for each area}
    \item{beta}{estimation of beta coefficient}
    \item{psi}{MSE of direct estimator}
    \item{mse}{MSE of small area estimation}
    \item{call}{function call}
}
\references{
    Ybarra, L. M., & Lohr, S. L. (2008). Small area estimation when auxiliary information is measured with error. Biometrika, 919-931.
}
\author{
    Amin Fathullah (\email{amin_fathullah@yahoo.com})
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

    \examples{
data(bmi)

#calculate SAE with measurement error
sae_me <- FHme(y = bmi$y, x = cbind(1, bmi$x), vardir = bmi$mse_y, C = cbind(0, bmi$mse_x))

sae_me
summary(sae_me)
plot(sae_me)
}
\keyword{sae}
\keyword{covariate}
\keyword{measuremnt}
\keyword{error}
