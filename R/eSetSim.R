# Function to generate expression set from a set of paramters using the mvrnorm function

## Arguments
### n0 = Number of samples in group 0
### n1 = Number of samples in group 1
### cor0 = Correlation matrix from group 0
### cor1 = Correlation matrix from group 1

## Value
### An ExpressionSet object

eSetSim <- function(n0, n1, cor0, cor1) {
    
    ### Data matrix
    dat1 <- mvrnorm(n0, rep(0, ncol(cor0)), cor0)
    dat2 <- mvrnorm(n1, rep(0, ncol(cor1)), cor1)
    
    ### Convert to expression set
    eMat <- cbind(t(dat1), t(dat2))
    colnames(eMat) <- paste0("s", 1:ncol(eMat))
    rownames(eMat) <- paste0("f", 1:nrow(eMat))
    eSet <- ExpressionSet(eMat)
    pData(eSet) <- data.frame(group = paste0("g", c(rep(0, n0), rep(1, n1))), row.names = colnames(eMat))
    
    return(eSet)
}
