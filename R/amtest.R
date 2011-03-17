"amtest" <- function(
modelList, varName, vcov. = c("sandwich", "model-based"), sig.level = 0.05, display = TRUE, adjp=FALSE)
{
    vcov. <- match.arg(vcov.)  # only "sandwich" and "model-based" acceptable values
    
    ## Loading packages that are used
    require(gtools, quietly = TRUE)
    require(mvtnorm, quietly = TRUE)
    require(sandwich, quietly = TRUE)
    
    ## Obtaining iid decompositions
    makeIIDdecomp <- function(modelObject)
    {
        ## Retrieving number of observations used in the fit
        numObsUsed <- ifelse(inherits(modelObject, "coxph"), modelObject$"n", nrow(modelObject$"model"))
     
        ## Using functions from the 'sandwich' package
        iidVec0 <- bread(modelObject)[varName, , drop = FALSE] %*% t(estfun(modelObject))
        moNAac <- modelObject$"na.action"
        numObs <- numObsUsed + length(moNAac) 
        iidVec <- rep(0, numObs)        
        if (!is.null(moNAac))
        {
            iidVec[-moNAac] <- sqrt(numObs / numObsUsed) * iidVec0
        } else {
            iidVec <- iidVec0
        }
        list(iidVec = iidVec, numObsUsed = numObsUsed, numObs = numObs)
    }    
    iidList <- lapply(modelList, makeIIDdecomp)
    numModels <- length(modelList)   
    iidresp <- matrix(as.vector(unlist(lapply(iidList, function(listElt){listElt[[1]]}))), 
    nrow = numModels, byrow = TRUE)

    ## Defining function for retrieving values (e.g. estimates or standard errors) from the individual fits    
    pickFct <- function(modelObject, matchStrings)
    {
        as.vector(na.omit((coef(summary(modelObject))[varName, ])[matchStrings]))
    }
    
    ## Retrieving variance estimates from the individual fits    
    numObsUsed <- as.vector(unlist(lapply(iidList, function(listElt){listElt[[2]]})))
    varest <- diag(sqrt(1 / (numObsUsed * (as.vector(unlist(lapply(modelList, pickFct, 
    matchStrings = c("Std. Error", "se(coef)")))))^2)))
    # "Std. Error" and "se(coef)" used in glm(), lm() and coxph() summary output, respectively 

    ##  Retrieving parameter estimates from the individual fits
    estVec <- as.vector(unlist(lapply(modelList, pickFct, matchStrings = c("Estimate", "coef"))))    
    # "Estimate", "coef" used in glm(), lm() and coxph() summary output, respectively 
    
    ## Defining estimated variance-covariance matrix
    covar <- (iidresp %*% t(iidresp)) / iidList[[1]]$"numObs"

    ## Calculating test statistics
    zVec <- sqrt(numObsUsed) * estVec

    ## Calculating estimated variance-covariance matrix for test statistics
    if (identical(vcov., "model-based"))
    {  
        ## Using model-based estimated standard errors
        varestalt <- diag(1 / sqrt(diag(covar)))
        vcMat <- varestalt %*% covar %*% varestalt
        zVec <- zVec * diag(varest)
    } else {
        ## Using iid decomposition based standard errors  
        varestalt <- diag(1 / sqrt(diag(covar)))
        zVec <- zVec *diag(varestalt)
        vcMat <- varestalt %*% covar %*% varestalt
    }
    pvals <- 1 - pchisq(zVec * zVec, 1)

 probFct<-function(alphatilde)
    {
    
        atHalf<-alphatilde/2
        1-pmvnorm(lower=rep(qnorm(atHalf),numModels),upper=rep(qnorm(1-atHalf),numModels),mean=rep(0,numModels),sigma=vcMat)
    }

    ## Calculating asymptotic correction
    asympCorr <- alphacorrected(numModels, sig.level, vcMat)
    sleCorr <- 1-(1-sig.level)^(1 / numModels)     #Slepian correction
  
    ##Calculating adjusted p-values
    if (adjp){
     slepvals<-1-(1-pvals)^(numModels)
     alphacorrVec <- Vectorize(probFct, "alphatilde")
     asymppvals<-alphacorrVec(pvals)

    }
       

    ## Showing the results 
    if (display)
    {  
        numDigits <- 4  # hardcoded at the moment
        cat(paste("Nominal level:", format(sig.level, digits = numDigits), "\n\n")) 
        if (!adjp)
        {
        cat(paste("Slepian-corrected level:", format(sleCorr, digits = numDigits), "\n\n"))    
        cat(paste("Asymptotically corrected level:", format(asympCorr, digits = numDigits), "\n\n"))
        cat("P-values:\n")
        cat(format(pvals, digits = numDigits))
        cat("\n")
        }
        if(adjp){
        cat("Slepian-corrected P-values:\n")
        cat(format(slepvals, digits = numDigits))
        cat("\n\n")
        cat("Asymptotically corrected P-values:\n")
        cat(format(asymppvals, digits = numDigits))
        cat("\n\n")
        warning("Not to be interpreted as ordinary p-values",call.=FALSE)
        
        }
        
    }

    if (adjp) 

    {
         return(invisible(list(nominal = sig.level, slepian = sleCorr, asymptotic = asympCorr, p.values = pvals,
                   asympadj.pvalues = asymppvals, slepianadj.pvalues = slepvals)))
 
    }


    if (!adjp) 

    {
         return(invisible(list(nominal = sig.level, slepian = sleCorr, asymptotic = asympCorr, p.values = pvals)))
 
    }

}
