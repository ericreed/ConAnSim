# # A basic R-wrapper for calling the python LFR function from 'LFR_Comminities.py'.

## Arguments
### See (https://networkx.org/documentation/stable/reference/generated/networkx.generators.community.LFR_benchmark_graph.html)

## Values
### Named list with two elements
#### ig = 'igraph' object of scale-free network
#### coms = Factor vector containing the community assignment for each node

# A basic R-wrapper for calling the python LFR funciton
simNetComs <- function(n, tau1, tau2, mu, average_degree, min_community, max_community, ...) {
    
    lfrOut <- LFR(n=n, 
               tau1=tau1, 
               tau2=tau2, 
               mu=mu, 
               average_degree=average_degree, 
               min_community=min_community, 
               max_community=max_community,
               ...)
    
    # Create adjacency matrix
    adj <- nx$adjacency_matrix(lfrOut)
    ig <- igraph::graph_from_adjacency_matrix(adj, mode="undirected", diag=FALSE)
    ig <- igSimplify(ig)
    
    coms <- lapply(seq(0, ncol(adj)-1), function(x, lfrOut) {
                as.character(lfrOut$nodes[[x]]$community)
            }, lfrOut) %>%
                unlist %>%
                as.factor %>%
                as.numeric %>%
                as.factor
    
    return(list(ig = ig, coms = coms))
    
}

# Removes loops, multi-edges, and isolated vertices
igSimplify <- function(ig) {
    ig <- igraph::simplify(ig, remove.multiple=TRUE, remove.loops=TRUE)
    return(ig)
}