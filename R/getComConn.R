# Function to get connectivity for a given community from an adjacency matrix

## Arguments
### x = The community identifier or "bg" to calculate background connectivity
### adjMat = Adjacency matrix (p x p)
### coms = A vector of length p of community members

## Value
### Numeric value of community or background connectivity

getComConn <- function(x, adjMat, coms) {
    
    if(x == "bg") {
        adjSub <- adjMat
        for(i in levels(coms)) adjSub[coms == i, coms == i] <- NA
    } else {
        adjSub <- adjMat[coms == x, coms == x]
    }
    
    adjSub %>%
        .[upper.tri(.)] %>%
        `^`(2) %>%
        mean(., na.rm = TRUE)
}