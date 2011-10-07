\name{amtest}

\Rdversion{1.1}

\alias{amtest}

\title{
Testing multiple outcomes
}

\description{
  Testing multiple endpoints or outcomes measures by providing a corrected significance level to assess the 
  uncorrected, marginal p-values against. 
}

\usage{
  amtest(modelList, varName, vcov. = c("sandwich", "model-based"), 
         sig.level = 0.05, display = TRUE, adjp=FALSE)
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
  \item{vcov.}{
  string specifying the type of estimated variance-covariance matrix to use for calculating the
  test statistics. Two options are available model-based estimated standard errors or robust estimated sandwich
  standard errors, which is the default.
} 
  \item{sig.level}{
  numeric specifying the nominal significance level or type I error rate (default is 0.05).
}
  \item{display}{logical indicating whether or not the results should be shown at the command line.
}
  \item{adjp}{logical indicating whether or not the correction should be reported in terms of adjusted p-values.
}
}

\details{
  The function calculates a corrected overall significance level that asymptotically reaches the desired type I
  error using i.i.d. decompositions of parameter estimates by means of score components. The details are
  provided by Pipper and Ritz (2011).
  
  The implemented method is less conservative than the practically identical Bonferroni/Slepian corrections and it flexibly adapts to the 
  correlation structure between outcomes. Different types of outcome measures such as continuous, binomial, 
  and event times, as well as missing values  in outcome measures and explanatory variables can be handled 
  by the procedure.
}

\value{
  The corrected significance level as well as the Slepian correction  
 \eqn{1-(1-Nominal significance level)^(1/Number of tests)} and the marginal
  p-values are shown by default (can be switched off using the argument 'display').
  
  These values are also returned invisibly in a list.
}

\references{

Pipper, C. B., Ritz, C., Bisgaard, H. (2011).
A versatile method for confirmatory evaluation of the effects of a covariate in multiple models.
\emph{To appear in JRSSC.}

}

\author{
  Christian B. Pipper and Christian Ritz
}

%\note{}

%\seealso{}

\examples{


##Example from: McCullagh, P. and Nelder, J.A. (1989,p. 239). Generalized Linear Models. 
##Second Edition. Chapman & Hall/CRC.


model1 <- glm(y1~group, family=binomial,na.action=na.omit,data=mn6.9)

model2 <- glm(y2~group, family=binomial,na.action=na.omit,data=mn6.9)

model3 <- glm(y3~group, family=binomial,na.action=na.omit,data=mn6.9)

model4 <- glm(y4~group, family=binomial,na.action=na.omit,data=mn6.9)

amtest(list(model1,model2,model3,model4),"group2",vcov.="model-based")

amtest(list(model1,model2,model3,model4),"group2")




##Multiple testing in case of perfectly correlated outcomes


amtest(list(model2,model2),"group2",vcov.="model-based")

amtest(list(model2,model2),"group2")



## Multiple testing in a subgroup analysis and overall analysis


subgroup<-mn6.9

set.seed(2982)

subgroup[sample(1:2982,1500),]<-NA

model1<- glm(y2~group, family=binomial,na.action=na.omit,data=mn6.9)

model2<- glm(y2~group, family=binomial,na.action=na.omit,data=subgroup)

amtest(list(model1,model2),"group2",vcov.="model-based")

amtest(list(model1,model2),"group2")

}

\keyword{models}
\keyword{htest}

\concept{type I error}
\concept{Family wise error rate}

