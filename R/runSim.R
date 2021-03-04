# Function run conan based from a correlation matric and a set of parameters specifying differences in modular connectivity and background connectivity within groups

## Arguments
### netConn = Output of simAdjMat()
### diff = Difference in modular connectivity
### diffbg = Difference in background connectivity
### nIter = Number of iterations to run ConAn
### nCores = Number of cores to run ConAn
### seed = Seed for simulating data and running ConAn

## Value
### Numeric value of community or background connectivity

# Function to run simulations
runSim <- function(netConn, diff, diffbg, n1, n0, nIter, nCores, seed, N_genes, iter_bg) {

    # Get correlation matrices
    corMat <- netConn$adjMat

    # Shrink modular matrices
    cor_modRed <- corMat
    cor1 <- cor0 <- corMat
    if(diff != 0) {
        target <- max(netConn$conns) - abs(diff)
        whCom <- names(netConn$conns)[which.max(netConn$conns)]

        ## Get shrinking factors for each module
        corModx <- corMat[netConn$coms == whCom, netConn$coms == whCom]
        if(target > 0) {
            shrFac <- findShrinkFactor(target, corModx)
        } else {
            shrFac <- 1
        }

        for(i in levels(netConn$coms)) {

            cor_modRed[netConn$coms == i, netConn$coms == i] <- shrinkCor(corMat[netConn$coms == i, netConn$coms == i], shrFac)

        }

        # Fix positive definiteness
        cor_modRed <- nearPD(cor_modRed, corr = TRUE)$mat

        # Fix background connectivity
        bgRed <- getComConn("bg", cor_modRed, netConn$coms)
        bgOrig <- getComConn("bg", corMat, netConn$coms)
        if(bgRed < bgOrig) {
            mdShrink <- findShrinkFactor(bgRed, corMat)
            corMat <- shrinkCor(corMat, mdShrink)
        } else {
            mdShrink <- findShrinkFactor(bgOrig, cor_modRed)
            cor_modRed <- shrinkCor(cor_modRed, mdShrink)
        }

        # Assign correlation matrix based on differences
        if(diff > 0) {
            cor1 <- corMat
            cor0 <- cor_modRed
        } else {
            cor0 <- corMat
            cor1 <- cor_modRed
        }

    }
    conn1 <- getComConn(levels(netConn$coms), cor1, netConn$coms)
    conn0 <- getComConn(levels(netConn$coms), cor0, netConn$coms)

    Pres <- NULL
    if((diff != 0 & !identical(conn1, conn0)) | diff == 0) {
        # Shrink background matrices
        if(diffbg != 0) {
            bgtarg <- netConn$bg - abs(diffbg)

            corComRem <- corMat
            for(i in levels(netConn$coms)) corComRem[netConn$coms == i, netConn$coms == i] <- NA

            bgshrinkfactor <- findShrinkFactor(bgtarg, corComRem)

            if(diffbg > 0) {
                cor0 <- shrinkCor(cor0, bgshrinkfactor)
            }

            if(diffbg < 0) {
                cor1 <- shrinkCor(cor1, bgshrinkfactor)
            }

        }

        ## Simulate expression set
        eSet <- eSetSim(n0, n1, cor0, cor1)

        ## Create module list
        modList <- lapply(levels(netConn$coms), function(x) rownames(eSet)[netConn$coms == x])
        names(modList) <- levels(netConn$coms)

        ## Run conan
        Pres <- conanWrapperPval(eSet, modList, nIter, nCores, seed, N_genes, iter_bg)

        ## Add results
        Pres$conn1 <- conn1
        Pres$conn0 <- conn0

        ## Add BG cor
        Pres$bg1 <- getComConn("bg", cor1, netConn$coms)
        Pres$bg0 <- getComConn("bg", cor0, netConn$coms)

        Pres$diff <- Pres$conn1 - Pres$conn0
        Pres$diffbg <- Pres$bg1 - Pres$bg0

        Pres$diffInput <- diff
        Pres$diffbgInput <- diffbg

    }

    return(Pres)

}
