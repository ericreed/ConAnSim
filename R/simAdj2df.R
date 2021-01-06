# Function to create matrix of simulated adjacencies

## Arguments
### netConn = The output of simAdjMat()

## Value
### Data frame with 3 columns
#### cors = Upper triangle of correlation matrix
#### edges = Upper triangle of adjacency matrix of original graph
#### coms = Upper triangle of matrix defining community connections ("0" if not in the same community)

simAdj2df <- function(netConn) {
    
    # Get correlations
    cors <- netConn$adjMat[upper.tri(netConn$adjMat)]
    
    # Get network edges
    edges <- as.matrix(netConn$adjNet)[upper.tri(as.matrix(netConn$adjNet))]
    
    # Get community connections
    comConns <- outer(netConn$coms, netConn$coms, function(x, y) ifelse(x == y, x, "0"))
    comConns <- comConns[upper.tri(comConns)]
    
    return(data.frame(cors = cors, edges = edges, coms = comConns))
    
}
