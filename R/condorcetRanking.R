condorcetRanking <- function (scoreMatrix, mutationMatrix, pen = 0.85, parallel = NULL, par = T) {
    pretrunc <- apply(mutationMatrix, 1, sum)
    trunc <- names(pretrunc[pretrunc > 0])
    trunc <- intersect(trunc, rownames(scoreMatrix))
    pwc <- function(p1, p2, scoreMatrix, mutationMatrix) {
        voters <- colnames(mutationMatrix)[(mutationMatrix[p1,
                                                           ] != 0) | (mutationMatrix[p2, ] != 0)]
        fightmat <- as.matrix(scoreMatrix[c(p1, p2), voters])
        p1wins <- sum((fightmat[1, ] - fightmat[2, ]) > 0)
        p2wins <- ncol(fightmat) - p1wins
        return(c(p1wins, p2wins))
    }
    scoreMatrix[mutationMatrix == 0] <- pen * scoreMatrix[mutationMatrix ==
                                                          0]
    scoreMatrix <- scoreMatrix[trunc, ]
    cmat <- matrix(0, nrow = length(trunc), ncol = length(trunc),
                   dimnames = list(trunc, trunc))
    
    if (par) {
        
        parFun <- function(i, trunc, scoreMatrix, mutationMatrix) {
            colTmp <- rowTmp <- numeric(length(trunc))
            for (j in 1:length(trunc)) {
                if (j > i) {
                    vres <- pwc(trunc[i], trunc[j], scoreMatrix,
                                mutationMatrix)
                    rowTmp[j] <- vres[1]
                    colTmp[j] <- vres[2]
                }
            }
            return(list(r = rowTmp, c = colTmp))
        }

        if (!is.null(parallel)) {
            
            library(snowfall)

            sfInit(parallel = T, cpus = parallel)

            sfExport("trunc", "scoreMatrix", "mutationMatrix", "parFun")
            
            tmp <- sfLapply(1:nrow(cmat), parFun, trunc, scoreMatrix, mutationMatrix)
            
            sfStop()
            
        } else {
            tmp <- lapply(1:nrow(cmat), parFun, trunc, scoreMatrix, mutationMatrix)
        }
        
        tmp2 <- do.call("cbind", unlist(tmp, recursive = F))
        cmatW <- t(tmp2[, (2*(1:ncol(cmat)) - 1)])
        cmatL <- tmp2[, (2*(1:ncol(cmat)))]
        colnames(cmatW) <- colnames(cmatL) <- rownames(cmatW) <- rownames(cmatL) <- colnames(cmat)
        cmat <- cmatW + cmatL
        
    } else {
        for (i in 1:nrow(cmat)) {
            for (j in 1:ncol(cmat)) {
                if (j > i) {
                    vres <- pwc(trunc[i], trunc[j], scoreMatrix,
                                mutationMatrix)
                    cmat[i, j] <- vres[1]
                    cmat[j, i] <- vres[2]
                }
            }
        }
    }
    wins <- apply(cmat, 1, sum)
    losses <- apply(cmat, 2, sum)
    copelandcriterion <- sort(wins/(wins + losses), decreasing = TRUE)
    return(list(cmat, copelandcriterion))
}
