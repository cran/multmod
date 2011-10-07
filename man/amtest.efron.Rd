\name{amtest.efron}

\Rdversion{1.1}

\alias{amtest.efron}

\title{
Adjusted p-values for multiple tests
}

\description{
  Provides adjusted p-values for multiple z-tests utilizing correlation between pairs of z-tests
}

\usage{
  amtest.efron(modelList, varName, sig.level = 0.05, display = TRUE)
}

\arguments{
  \item{modelList}{
  list of model fits. The fits may be obtained using \code{\link{lm}}, \code{\link{glm}}, or 
  \code{\link[survival]{coxph}}. The fits should all be based on exactly the same dataset (data frame) and missing
  values should be handled using \code{\link{na.omit}}. All models should contain the explanatory variable of
  interest.
}
  \item{varName}{string specifying the name of the explanatory variable for which the effect is to be assessed. 
}
  \item{sig.level}{
  numeric specifying the nominal significance level or type I error rate (default is 0.05).
}
  \item{display}{logical indicating whether or not the results should be shown at the command line.
}
}

\details{
  The function calculates adjusted p-values based on an approximation provided in Efron (1997) to bound the type I error. An algorithm described in Jensen et al (2011) is utilized to pair the test statistics so as to ensure an optimal bound based on this approximation. Correlation between test statistics is obtained by means of the i.i.d. decomposition as detailed in Pipper et al (2011). 
If test statistics are correlated the implemented method is less conservative than the Bonferroni correction. The function is designed to be applicable to large scale multiple testing.      
  
}

\value{
  Efron adjusted p-values, Sleipan adjusted p-values, original p-values. 
  
  These values are also returned invisibly in a list.
}

\references{

Efron, B. (1997). 
The Length Heuristic for Simultaneous Hypothesis Tests.
\emph{Biometrika, 84, 143-157.}

Jensen, S.M., Ritz, C., Pipper, C.B. (2011).
An efficient single-step adjustment in large-scale multiple testing problems.
\emph{Submitted}.

Pipper, C. B., Ritz, C., Bisgaard, H. (2011).
A versatile method for confirmatory evaluation of the effects of a covariate in multiple models.
\emph{To appear in JRSSC.}

}

\author{
  Signe Marie Jensen, Christian Ritz, and Christian Pipper
}

%\note{}

%\seealso{}

\examples{


##Example from: McCullagh, P. and Nelder, J.A. (1989,p. 239). Generalized Linear Models. 
##Second Edition. Chapman & Hall/CRC.


model <- glm(y2~group, family=binomial,na.action=na.omit,data=mn6.9)


##Multiple testing in case of perfectly correlated test statistics


amtest.efron(list(model,model),"group2")

}


\keyword{models}
\keyword{htest}

\concept{type I error}
\concept{Family wise error rate}
\concept{Adjusted p-value}

