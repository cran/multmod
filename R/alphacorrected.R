"alphacorrected" <- function(K, alpha, sigma, error = .Machine$double.eps^0.25)  # argument 'error' not used
{
    K <- dim(sigma)[1]

#   appFct <- function(upper)
#    {
#        pmvnorm(upper = upper, mean = rep(0, K), sigma = sigma)
#    }
#    pointMat0 <- t(permutations(2, K, c(-1, 1), repeats.allowed = TRUE))
#    signVec <- apply(pointMat0, 2, function(x){prod(sign(x))})

    ## Defining helper function
#    "probFct" <- function(alphatilde)
#    {
#        atHalf <- alphatilde / 2
##        pointMat <- t(permutations(2, K, c(qnorm(atHalf), qnorm(1 - atHalf)), repeats.allowed = TRUE))
##        pointMat <- qnorm(1 - atHalf) * pointMat0
        
#        retVec <- 1 - sum(signVec * apply(qnorm(1 - atHalf) * pointMat0, 2, appFct))
##        print(retVec)
#        retVec
#    } 
    
    probFct<-function(alphatilde)
    {
    
        atHalf<-alphatilde/2
        1-pmvnorm(lower=rep(qnorm(atHalf),K),upper=rep(qnorm(1-atHalf),K),mean=rep(0,K),sigma=sigma)
    }
    
    
    ## Finding corrected overall significance level
    uniroot(function(alphatilde){probFct(alphatilde) - alpha}, c(0, alpha), tol = error)$root
}