# Wrapper function for simulating data from the output of simAdjMat(). 

## Arguments
### netConn = The output of simAdjMat()
### n = Number of samples to simulate

## Value
### An n by p matrix of multivariate normal data

simAdj2data <- function(netConn, n) {
    
    # Generate data
    mvrnorm(n, rep(0, ncol(netConn$adjMat)), netConn$adjMat)
    
}
