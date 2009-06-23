\name{mn6.9}
\alias{mn6.9}
\docType{data}
\title{I.Q. and attitude towards science}
\description{
  Responses given by 2982 New Jersey high-school seniors on 4 questions concerning attitude towards science. Also recorded was whether students had a high or low I.Q. 
}
\usage{data(mn6.9)}
\format{
  A data frame with 2982 observations on the following 5 variables.
  \describe{
    \item{\code{y1}}{Agree=1/disagree=0 to "The development of new ideas is the scientist's greatest source of satisfaction"}
    \item{\code{y2}}{Agree=1/disagree=0 to "Scientists and engineers should be eliminated form the military draft"}
    \item{\code{y3}}{Agree=1/disagree=0 to "The scientist will make his maximum contribution to society when he has freedom to work on problems that interest him"}
    \item{\code{y4}}{Agree=1/disagree=0 to "The monetary compensation of a Nobel Prize-winner in physics should be at least equal to that given to popular entertainers"}
    \item{\code{group}}{I.Q.levels: 1=low, 2=high}
  }
}

\source{
  McCullagh, P. and Nelder, J.A. (1989, p. 239). \emph{Generalized Linear Models}. Second Edition. Chapman \& Hall/CRC.  
}

\examples{

##Marginal assessment of effect of I.Q. group on response probability

model1 <- glm(y1~group, family=binomial,data=mn6.9)

model2 <- glm(y2~group, family=binomial,data=mn6.9)

model3 <- glm(y3~group, family=binomial,data=mn6.9)

model4 <- glm(y4~group, family=binomial,data=mn6.9)



}
\keyword{datasets}
