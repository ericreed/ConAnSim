# Function to estimate moving discovery rates for differences

## Arguments
### conn = True connectivity of control group (vector)
### conn = True difference in connectivity (vector)
### sig = Discovery (logical)
### widthConn = Width of moving average window for connectivity
### widthDiff = Width of moving average window for difference
### gapConn = Gap between connectivity windows
### gapDiff = Gap between difference windows
### minN = Minimum number of samples to calculate DR


## Value
### An ExpressionSet object
movingDR2 <- function(conn, diff, sig, widthConn = 0.05, widthDiff = 0.05, gapConn = 0.01, gapDiff = 0.01, minN = 20) {
    
    # Create data frame of results
    df <- data.frame(conn, diff, sig)
    
    # Create vector of bins
    
    ## Connectivity
    binsConn <- data.frame(lower = seq(0, 1-widthConn, by = gapConn),
                           upper = seq(widthConn, 1, by = gapConn))
    
    ## Diff
    binsDiff <- data.frame(lower = round(seq(-1, 1-widthDiff, by = gapDiff), 4),
                           upper = round(seq(-1 + widthDiff, 1, by = gapDiff), 4))
    binsDiff <- binsDiff[sign(binsDiff$lower * binsDiff$upper) != -1,]
    
    # Get moving average of discovery rate
    mvDF <- do.call(rbind, apply(binsConn, 1, function(x, df, binsDiff) {
        dfSub <- df[df$conn > x[1] & df$conn < x[2],]
        
        do.call(rbind, apply(binsDiff, 1, function(y, x, dfSub) {
            dfSubDiff <- dfSub[dfSub$diff > y[1] & dfSub$diff < y[2],]
            
            # Get sample size of this bin
            n <- nrow(dfSubDiff)
            
            # Calculate discovery rate
            dr <- NA
            if(n >= minN) {
                dr <- mean(dfSubDiff$sig)
            }
            
            # Calculate standard error
            se <- sqrt((dr * (1-dr)) / n)
            
            # Calculate CIs
            uCI <- dr + qnorm(0.975) * se
            lCI <- dr - qnorm(0.975) * se
            
            return(data.frame(lowerConn = x[1],
                              upperConn = x[2],
                              lowerDiff = y[1],
                              upperDiff = y[2],
                              dr = dr,
                              se = sqrt((dr * (1-dr)) / n),
                              uCI = uCI,
                              lCI = lCI,
                              n = n))
        }, x, dfSub))
    }, df, binsDiff))
    
    mvDF <- mvDF[!is.na(mvDF$dr),]
    
    return(mvDF)
}
