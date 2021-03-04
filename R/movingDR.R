# Function to estimate moving discovery rates

## Arguments
### conn = True connectivity of each group (vector)
### sig = Discovery (logical)
### width = Width of moving average windeo
### gap = Distance between windown
### minN = Minimum number of samples to calculate DR

## Value
### An ExpressionSet object
movingDR <- function(conn, sig, width = 0.05, gap = 0.01, minN = 20) {
    
    # Create data frame of results
    df <- data.frame(conn, sig)
    
    # Create vector of bins
    bins <- data.frame(lower = seq(0, 1-width, by = gap), 
                       upper = seq(width, 1, by = gap))
    
    # Get moving average of discovery rate
    mvDF <- do.call(rbind, apply(bins, 1, function(x, df) {
        dfSub <- df[df$conn > x[1] & df$conn < x[2],]
        
        # Get sample size of this bin
        n <- nrow(dfSub)
        
        # Calculate discovery rate
        dr <- NA
        if(n >= minN) {
            dr <- mean(dfSub$sig)
        }
        
        # Calculate standard error
        se <-sqrt((dr * (1-dr)) / n)
        
        # Calculate CIs
        uCI <- dr + qnorm(0.975) * se
        lCI <- dr - qnorm(0.975) * se
        
        return(data.frame(lower = x[1],
                          upper = x[2],
                          dr = dr,
                          se = sqrt((dr * (1-dr)) / n),
                          uCI = uCI,
                          lCI = lCI,
                          n = n))
    }, df))
    
    return(mvDF)
}
