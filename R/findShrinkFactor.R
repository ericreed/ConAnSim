# Function to find shrinking factor using a correlation matrix and specified target

## Arguments
### target = The target connectivity
### adjMat = A correlation matrix

## Value
### An n by p matrix of multivariate normal data

findShrinkFactor <- function(target, adjMat) {
    
    if(target >= 1 | target <= 0 ) stop("Target must be between 0 and 1." )
    optimize(optFuncShrink, interval = c(0, 1), targ = target, cm = adjMat)$minimum
    
}

optFuncShrink <- function(sf, targ, cm) {
    
    zTrans <- atanh(cm) * sf
    Est <- zTrans[upper.tri(zTrans)] %>%
        tanh %>%
        `^`(2) %>%
        mean(., na.rm = TRUE)
    abs(targ - Est)
    
}
