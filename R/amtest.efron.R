

"amtest.efron"<-
function (modelList, varName, vcov. = c("sandwich", "model-based"), 
    sig.level = 0.05, display = TRUE) 
{
    vcov. <- match.arg(vcov.)
    require(gtools, quietly = TRUE)
    require(mvtnorm, quietly = TRUE)
    require(sandwich, quietly = TRUE)
    makeIIDdecomp <- function(modelObject) {
        numObsUsed <- ifelse(inherits(modelObject, "coxph"), 
            modelObject$n, nrow(modelObject$model))
        iidVec0 <- bread(modelObject)[varName, , drop = FALSE] %*% 
            t(estfun(modelObject))
        moNAac <- modelObject$na.action
        numObs <- numObsUsed + length(moNAac)
        iidVec <- rep(0, numObs)
        if (!is.null(moNAac)) {
            iidVec[-moNAac] <- sqrt(numObs/numObsUsed) * iidVec0
        }
        else {
            iidVec <- iidVec0
        }
        list(iidVec = iidVec, numObsUsed = numObsUsed, numObs = numObs)
    }
    iidList <- lapply(modelList, makeIIDdecomp)
    numModels <- length(modelList)
    iidresp <- matrix(as.vector(unlist(lapply(iidList, function(listElt) {
        listElt[[1]]
    }))), nrow = numModels, byrow = TRUE)
    pickFct <- function(modelObject, matchStrings) {
        as.vector(na.omit((coef(summary(modelObject))[varName, 
            ])[matchStrings]))
    }
    numObsUsed <- as.vector(unlist(lapply(iidList, function(listElt) {
        listElt[[2]]
    })))
    estVec <- as.vector(unlist(lapply(modelList, pickFct, matchStrings = c("Estimate", 
        "coef"))))
    zVec1 <- sqrt(numObsUsed) * estVec
    zVec<-zvecfun(iidresp,iidList[[1]]$numObs,zVec1)
    pvals <- 1 - pchisq(zVec * zVec, 1)
    efronsum<-efron(iidresp,iidList[[1]]$numObs,nrow(iidresp),abs(zVec)) 
    efronpvals<-rep(0,length(zVec))
    for(i in 1:length(zVec)){
  efronpvals[i]<-min(1,pvals[i] + efronsum[i])
    	}
	
    sleCorr <- 1 - (1 - sig.level)^(1/numModels)

	slepvals <- 1 - (1 - pvals)^(numModels)
    
    if (display) {
        numDigits <- 4
        cat(paste("Nominal level:", format(sig.level, digits = numDigits), 
            "\n\n"))
            cat("Slepian-corrected P-values:\n")
            cat(format(slepvals, digits = numDigits))
            cat("\n\n")
            cat("Efron corrected P-values:\n")
            cat(format(efronpvals, digits = numDigits))
            cat("\n\n")
            cat("P-values:\n")
            cat(format(pvals, digits = numDigits))
                
    }
	return(invisible(list(nominal = sig.level, slepian = sleCorr, 
            p.efron = efronpvals, p.values = pvals, p.slepian = slepvals)))
}
