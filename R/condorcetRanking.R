#' Aggregate gene ranking over all patients
#' 
#' /code{condorcetRanking} Returns aggregate gene list over all patients.
#' @param scoreMatrixa matrix containing all the given DawnRank scores per patient. Rows are genes, columns are  patients
#' @param mutationMatrix the mutation matrix. Rows are genes, columns are patients
#' @param pen the penality parameter in the condorcet algorithm for missing data. Default 0.85
#' @param parallel number of cores for parallel calculations. Default NULL.
#' @return the ranks. A list of 2 including a [[1]] a matrix of all pairwise comparisons, [[2]] the final rankings based on the Condorcet score
#' @examples
#' ###using a small subset of the TCGA dataset and a small KEGG 
#' ###gene interaction network,
#' ###We will show how to get DawnRank Results
#' 
#' library(DawnRank)
#' 
#' #load the mutation data
#' data(brcaExampleMutation)
#' 
#' #load the tumor expression data
#' data(brcaExampleTumorExpression)
#' 
#' #load the normal expression data
#' data(brcaExampleNormalExpression)
#' 
#' #load the pathway data
#' data(brcaExamplePathway)
#' 
#' #load the gold standard
#' data(goldStandard)
#' 
#' #normalize the tumor and normal data to get the differential expression
#' normalizedDawn<-DawnNormalize(tumorMat=brcaExampleTumorExpression,
#' normalMat=brcaExampleNormalExpression)
#' 
#' #get the DawnRank Score
#' dawnRankScore<-DawnRank(adjMatrix=brcaExamplePathway,
#' mutationMatrix=brcaExampleMutation,expressionMatrix=normalizedDawn, 
#' mu=3,goldStandard=goldStandard)
#' 
#' #look at the DawnRank scores for a few patients
#' dawnRankFrame<-dawnRankScore[[3]]
#' head(dawnRankFrame)
#' 
#' #get the aggregate DawnRank scores
#' aggregateDawnRankScore<-condorcetRanking(scoreMatrix=dawnRankScore[[2]],
#' mutationMatrix=brcaExampleMutation)
#' 
#' #look at top 10 ranked genes
#' top10<-aggregateDawnRankScore[[2]][1:10]
#' top10
#' 
#' #get the individual cutoff for patient TCGA-A2-A04P
#' dawnRankFrame$isCGC<-dawnRankFrame$isGoldStandard
#' library(maxstat)
#' significance<-patspeccutoff(patient="TCGA-A2-A04P",ms=dawnRankFrame,
#' default=95)
#' 
#' @export
condorcetRanking <- function (scoreMatrix, mutationMatrix, pen = 0.85, parallel = NULL, par = T, dovec = T) {
    pretrunc <- apply(mutationMatrix, 1, sum)
    trunc <- names(pretrunc[pretrunc > 0])
    trunc <- intersect(trunc, rownames(scoreMatrix))
    pwc <- function(p1, p2, scoreMatrix, mutationMatrix) {
        voters <- colnames(mutationMatrix)[(mutationMatrix[p1,
                                                           ] != 0) | (mutationMatrix[p2, ] != 0)]
        fightmat <- as.matrix(scoreMatrix[c(p1, p2), voters])
        p1wins <- sum((fightmat[1, ] - fightmat[2, ]) > 0)
        p2wins <- ncol(fightmat) - p1wins # does that make sense? so even if one is mutated and the other is not, the non mutated can win? and if both are not mutated it basically automatically loses (the entry remains 0)?
        return(c(p1wins, p2wins))
    }
    scoreMatrix[mutationMatrix == 0] <- pen * scoreMatrix[mutationMatrix ==
                                                          0]
    scoreMatrix <- scoreMatrix[trunc, ]
    mutationMatrix <- mutationMatrix[trunc, ]
    cmat <- matrix(0, nrow = length(trunc), ncol = length(trunc),
                   dimnames = list(trunc, trunc))
    
    if (par) {
        
        parFun <- function(i, trunc, scoreMatrix, mutationMatrix) {
            rowSums2 <- function(x, ...) {
                if (!is.null(dim(x))) {
                    return(rowSums(x, ...))
                } else {
                    return(t(x))
                }
            }
            colTmp <- rowTmp <- numeric(length(trunc))
            if (dovec) {
                ## try vector instead of for:
                if (i != nrow(scoreMatrix)) {
                    votemat <- mutationMatrix # does this take more memory ?
                    fightmat <- scoreMatrix
                    fightvec <- fightmat[i, ]
                    fightres <- t(t(fightmat) - fightvec)*(-1) # win should be positive
                    fightres[which(fightres > 0)] <- 1
                    fightres[which(fightres < 0)] <- 0
                    testmut <- which(votemat == 0 & arrayInd(1:length(votemat), dim(votemat))[, 2] %in% which(votemat[i, ] == 0))
                    if (length(testmut) > 0) { # maybe do this next with try()
                        fightres[testmut] <- 0
                        rowTmp[-i] <- rowSums2(fightres[-i, ])
                        fightres <- 1 - fightres
                        fightres[testmut] <- 0
                        colTmp[-i] <- rowSums2(fightres[-i, ])
                    } else {
                        rowTmp[-i] <- rowSums2(fightres[-i, ])
                        fightres <- 1 - fightres
                        colTmp[-i] <- rowSums2(fightres[-i, ])
                    }
                    rowTmp[1:i] <- colTmp[1:i] <- 0
                }
            } else {
                for (j in 1:length(trunc)) {
                    if (j > i) {
                        vres <- pwc(trunc[i], trunc[j], scoreMatrix,
                                    mutationMatrix)
                        rowTmp[j] <- vres[1]
                        colTmp[j] <- vres[2]
                    }
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
