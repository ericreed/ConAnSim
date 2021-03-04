# Function to simulate an adjacency matrix from simulate.lfr() output

## Arguments
### netComs = The output of simNetComs()
### D = A p dimension symetric matrix. Controls shape of distribution of simulated connectivities (See BDgraph::bdgraph.sim)
### b = "The degree of freedom for G-Wishart distribution, W_G(b, D)." (See BDgraph::bdgraph.sim)

## Values
### Named list with 5 elements
#### adjMat = Simulated correlation matrix
#### conns = Named numeric vector of community connectivies
#### bg = Background connectivity
#### coms = Factor vector containing the community assignment for each node (Same as netComs$coms)
#### adjNet = Sparse matrix of converted from original graph input, i.e. 0's and 1's.

simAdjMat <- function(netComs, D = NULL, b = 3) {

    # Create adjacency matrix
    adjNet <- netComs$ig %>%
        as_adjacency_matrix %>%
        as.matrix

    # Simulate covariance matrix (Adapted from BDgraph::bdgraph.sim)
    diag(adjNet) <- 0
    p <- nrow(adjNet)
    threshold <- 1e-8 # This is also hardcoded in BDgraph::bdgraph.sim

    Dmat <- diag(p)
    if(!is.null(D)) {
        Dmat[outer(netComs$coms, netComs$coms, "==")] <- D
        diag(Dmat) <- 1
        Dmat <- nearPD(Dmat, corr = TRUE)$mat
    }

    result <- .C("rgwish_c",
                 as.integer(adjNet),
                 as.double(chol(solve(Dmat))),
                 K = as.double(matrix(0,p,p)),
                 as.integer(b),
                 as.integer(p),
                 as.double(threshold),
                 PACKAGE = "BDgraph")

    K <- matrix(result$K, p, p)
    covMat = solve(K)

    # Convert covariance matrix to correlation matrix
    adjMat <- cov2cor(covMat)

    # Get community connectivities
    coms <- netComs$coms
    conns <- getComConn(levels(coms), adjMat, coms)

    # Get background connectivity
    bg <- getComConn("bg", adjMat, coms)

    # Convert adjNet to sparse matrix
    adjNet <- Matrix(adjNet, sparse = TRUE)

    # Return list of adjacency matrix, connectivities, and communities
    return(list(
        adjMat = adjMat,
        conns = conns,
        bg = bg,
        coms = coms,
        adjNet = adjNet
    ))

}
