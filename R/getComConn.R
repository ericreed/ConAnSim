# Function to get connectivity for a given community from an adjacency matrix

## Arguments
### x = A single value or vector community identifiers or "bg" to calculate background connectivity
### adjMat = Adjacency matrix (p x p)
### coms = A vector of length p of community members

## Value
### Numeric value of community or background connectivity

getComConn <- function(x, adjMat, coms) {
    
    if(length(x) == 1 && x == "bg") {
        adjSub <- adjMat
        for(i in levels(coms)) adjSub[coms == i, coms == i] <- NA
        getConn(adjSub)
    } else {
        sapply(x, function(y, adjMat) {
            getConn(adjMat[coms == y, coms == y])
        }, adjMat)
    }
}

getConn <- function(adjSub) {
    
    adjSub %>%
        .[upper.tri(.)] %>%
        `^`(2) %>%
        mean(., na.rm = TRUE)
    
}