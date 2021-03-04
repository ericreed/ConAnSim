## R packages to run simulation functions
require(magrittr)
library(reticulate)
library(Matrix)
library(BDgraph)
require(corpcor)
require(MASS)
require(Biobase)
require(ConAn)
require(igraph)

## Python modules
nx <- import("networkx")

## Directories
baseDir <- ".." # If running from the "./scripts" directory
rDir <- file.path(baseDir, "R")
pyDir <- file.path(baseDir, "python")
resDir <- file.path(baseDir, "results")

## Create output directory
time <- gsub("-| |:", "_", Sys.time())
outDir <- file.path(resDir, paste("sim_results", time, sample(seq(1000, 9999), 1), sep = "_"))
dir.create(outDir)

## Files to write out
simResFile <- file.path(outDir, "simResults.rds")

## Read in functions

## Read in R functions
rDir %>%
    list.files(., full.names = TRUE) %>%
    lapply(., source) %>%
    invisible # This keeps lapply from printing NULL

## Read in python functions
source_python(file.path(pyDir, "LFR_Communities.py"))

## Read in parameters
source("parameters.R")

simOut <- NULL

# Initialize timer
startTime <- Sys.time()
runTime <- 0

while(difftime(Sys.time(), startTime, units = "hours")[[1]] < 280) {

    # Set flag
    check <- TRUE

    while(check) {

        # Set seed ID
        seed <- sample(seq(100000000, 999999999), 1)

        # Simulate network
        netComs <- simNetComs(seed=seed,
                              n=nFeats,
                              tau1=3,
                              tau2=2,
                              mu=mu,
                              average_degree=aveDeg,
                              min_community=minComSize,
                              max_community=maxComSize)

        ## Use D matrix in correlation estimates
        set.seed(seed)
        netConn <- simAdjMat(netComs, D = Dval)

        ## Check if background connectivies are large enough
        check <- netConn$bg < minBG | max(netConn$conns) < minConn

        ## Print background connectivity
        cat("Background Connectivity =", netConn$bg, "\n")
        cat("Modular Connectivities =", netConn$conns, "\n")

    }

    # Get modular connectivities
    modConns <- netConn$conns

    # Floor modular connectivities
    modFloor <- floor(modConns/0.05)*0.05
    bgFloor <- floor(netConn$bg/.025)*.025

    # Get grid of connectivities to test
    modMax <- max(modFloor)
    connTest <- seq(-modMax, modMax, by = 0.05)
    connTest <- connTest[abs(connTest) < 0.51]
    bgTest <- round(seq(-bgFloor, bgFloor, by = 0.025), 3)
    bgTest <- bgTest[bgTest %in% c(-0.2, -0.1, -0.05, -0.025, 0, 0.025, 0.05, 0.1, 0.2)]
    testGrid <- expand.grid(connTest, bgTest)

    # Run simulations
    simRes <- do.call(rbind, lapply(seq(nrow(testGrid)), function(x) {
        sRes <- try(runSim(netConn, testGrid[x, 1], testGrid[x, 2], n0, n1, nIter, nCores, seed, N_genes, iter_bg), silent = TRUE)
        if(is(sRes, "try-error")) sRes <- NULL
        return(sRes)
    }))

    # Add module sizes
    sizeTab <- as.data.frame(table(netConn$coms))
    colnames(sizeTab) <- c("mod", "modsize")
    sizeTab$mod <- factor(sizeTab$mod, levels = levels(simRes$mod))
    simRes <- merge(simRes, sizeTab)

    # Add seed as identifier
    simRes$seedID <- seed

    if(is.null(simOut)) {
        simOut <- simRes
    } else {
        simOut <- rbind(simOut, simRes)
    }

    # Save results as rds file
    cat("Writing results for seed:", seed, "\n")
    saveRDS(simOut, simResFile)

}
