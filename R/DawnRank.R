DawnRank <- function(adjMatrix, expressionMatrix, mutationMatrix, mu = 20, 
    maxit = 100, epsilon = 1e-04, goldStandard = NULL, parallel = NULL) {
    damping <- dawnDamping(adjMatrix, mu)
    dawnMatrix<-DawnMatrix(adjMatrix)
    if (!is.null(parallel)) {
        
        library(snowfall)
        
        sfInit(parallel = T, cpus = parallel)

        parFun <- function(i) {
            Dawni <- Dawn(dawnMatrix, expressionMatrix[, i], mutationMatrix[, i], damping = damping, maxit = maxit, epsilon = epsilon, patientTag = colnames(expressionMatrix)[i], goldStandard = goldStandard)
            allPercentiles <- Dawni[[1]][, 2]
            allRanks <- Dawni[[1]][, 1]
            allMutated <- Dawni[[2]]
            iterations <- Dawni[[3]]
            return(list(p = allPercentiles, r = allRanks, m = allMutated, i = iterations))
        }

        Dawn <- get("Dawn", en = asNamespace("DawnRank"))

        sfExport("dawnMatrix", "expressionMatrix", "mutationMatrix", "damping", "maxit", "epsilon", "goldStandard", "parFun", "Dawn")

        tmp <- sfLapply(1:ncol(expressionMatrix), parFun)

        sfStop()
        
        tmp <- unlist(tmp, recursive = F)
        
        allPercentiles <- do.call("cbind", tmp[((1:ncol(expressionMatrix))*4 - 3)])

        allRanks <- do.call("cbind", tmp[((1:ncol(expressionMatrix))*4 - 2)])

        allMutated <- do.call("rbind", tmp[((1:ncol(expressionMatrix))*4 - 1)])

        iterations <- do.call("c", tmp[((1:ncol(expressionMatrix))*4 - 0)])
        
    } else {
        allRanks <- logical(0)
        allPercentiles <- logical(0)
        allMutated <- logical(0)
        iterations <- logical(0)
        for (i in 1:ncol(expressionMatrix)) {
            Dawni <- Dawn(dawnMatrix, expressionMatrix[, i], mutationMatrix[, 
                                                                           i], damping = damping, maxit = maxit, epsilon = epsilon, 
                          patientTag = colnames(expressionMatrix)[i], goldStandard = goldStandard)
            allPercentiles <- cbind(allPercentiles, Dawni[[1]][, 
                                                               2])
            allRanks <- cbind(allRanks, Dawni[[1]][, 1])
            allMutated <- rbind(allMutated, Dawni[[2]])
            iterations <- append(iterations, Dawni[[3]])
        }
    }
    rownames(allPercentiles) <- rownames(expressionMatrix)
    colnames(allPercentiles) <- colnames(expressionMatrix)
    rownames(allRanks) <- rownames(expressionMatrix)
    colnames(allRanks) <- colnames(expressionMatrix)
    allMutated$Frequency <- table(allMutated[, 1])[allMutated[, 
        1]]
    averageRank <- apply(allRanks, 1, mean)
    sdRank <- apply(allRanks, 1, sd)
    allMutated$Deviation <- (allMutated$Rank - averageRank[paste(allMutated$Gene)])/sdRank[paste(allMutated$Gene)]
    return(list(AllRanks = allRanks, AllPercentiles = allPercentiles, 
        AllMutated = allMutated, iterations = iterations))
}
