# Function to create grid of specifications

paramGrid <- function(connbins, diffs, diffbgs) {

    resCounts <- expand.grid(connbins, diffs, diffbgs)
    for(i in seq(ncol(resCounts))) resCounts[, i] <- round(resCounts[, i], 2)
    colnames(resCounts) <- c("conn0bin", "diff", "diffbg")
    resCounts <- resCounts[!(resCounts$conn0bin + resCounts$diff > 1 | resCounts$conn0bin + resCounts$diff < 0) ,]
    for(i in seq(ncol(resCounts))) resCounts[, i] <- as.character(resCounts[, i])
    resCounts$count <- 0
    rownames(resCounts) <- paste(resCounts$conn0bin, resCounts$diff, resCounts$diffbg, sep = "_")
    return(resCounts)

}
