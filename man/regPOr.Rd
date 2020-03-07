\name{regPOr}
\alias{regPOr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Regression analysis of right censored data under the proportional odds model
}
\description{
To fit the semi-parametric maximum likelihood estimator for right censored survival data using the EM algorithm.
}
\usage{
regPOr <- function(time,event,X,order = 3,equal_space = T,nknot, myknots,conv_cri = 1e-9,c_initial = F, initial)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{This is the follow up time.}
  \item{event}{The status indicator, normally 0=alive, 1=dead. }
  \item{X}{Covariate matrix, an n x p matrix.}
  \item{order}{The order of I-spline basis function.}
  \item{equal_space}{Logical, if TRUE knots are placed evenly across the observed intervals based on the input data set; if FALSE knots should be specified by the user. see myknots. The default is TRUE.}
  \item{nknot}{The total number of knots.}
  \item{myknots}{Knots specified by the user.}
  \item{conv_cri}{The convergence criterion of the EM algorithm.}
  \item{c_initial}{Logical, if FALSE 1 will be used as initial values for all the parameters ; if TRUE costumized initial values should be specified by the user. The default value for c_initial is FALSE.}
  \item{initial}{Customized initial values for the parameters to be optimized over}
}
\details{
The EM algorithm converges when the sum of square of the difference in the parameter estimates is less than conv_cri.
}
\value{
\item{beta}{Estimates of regression coefficients.}
\item{beta_se}{The standard errors of regression coefficients estimate.}
\item{CI}{Confidence intervals for the regression coefficients.}
\item{spline_coef}{Estimates of spline coefficients.}
\item{knots}{The knots for I-spline.}
\item{AIC}{The Akaike information criterion. }
\item{BIC}{The Bayesian information criterion.}
\item{Baseline_Surv}{Plot of the estimated baseline survival curve.}
\item{Baseline_hazard}{Plot of the estimated baseline hazard rate curve.}
\item{Baseline_odds}{Plot of the estimated baseline odds curve.}
}
\references{
Lu Wang & Lianming Wang. An EM algorithm for analyzing right-censored survival data under the semiparametric proportional odds model
}
\author{
Wang,L. and Wang,L.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(survival)
data(ovarian)
ovarian[,'resid.ds'] = ovarian[,'resid.ds'] - 1
fitpo <- regPOr(time = ovarian$futime, event = ovarian$fustat,X = ovarian$resid.ds,nknot = 10)
fitpo
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
