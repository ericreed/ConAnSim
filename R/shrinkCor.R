# Function to shrink correlation by a shrinking factor between 0 and 1

## Arguments
### adjMat = correlation matrix
### shrFac = shrinking factor (number between 0 and 1)

## Value
### A correlation matrix

shrinkCor <- function(adjMat, shrFac) {
    
    adjMat %>%
        atanh %>%
        `*`(shrFac) %>%
        tanh

}