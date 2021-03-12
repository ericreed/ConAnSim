require(ggplot2)
require(parallel)

## Directories
### /restricted/projectnb/montilab-p/personal/eric/ConAnSim/scripts
baseDir <- ".." # If running from the "./scripts" directory
rDir <- file.path(baseDir, "R")
pyDir <- file.path(baseDir, "python")
resDir <- file.path(baseDir, "results")
plotDir <- file.path(resDir, "plots")

## Read in R functions
rDir %>%
    list.files(., full.names = TRUE) %>%
    lapply(., source) %>%
    invisible # This keeps lapply from printing NULL


# Get results files directories
simResDirs <- list.files(resDir, pattern = "sim_results_", full.names = TRUE)

## Files to read in
simResFiles <- file.path(simResDirs, "simResults.rds")

## Files to write out
TPRout <- file.path(plotDir, "Power_ConAn_BG.png")

## Read in results file
res <- do.call(rbind, lapply(simResFiles, function(file) readRDS(file)))

## Add discovery logical (pval < 0.05)
res$sig <- res$pval < 0.05

## Subset for diff !=0
res0 <- res[res$diffInput != 0,]

## Get unique parameters
pars <- unique(res0[, c("sim_type", "mean_correct", "mdc_type", "diffbgInput")])

## Create moving average of discovery rates for each
maPars <- do.call(rbind, mclapply(seq(nrow(pars)), function(i) {
    
    # Subset for par estimates
    res0pars <- res0[res0$sim_type == pars$sim_type[i] &
                     res0$mean_correct == pars$mean_correct[i] &
                     res0$mdc_type == pars$mdc_type[i] &
                     res0$diffbgInput == pars$diffbgInput[i],
                     ]
    
    # Get moving average discovery rate for set of parameters
    drMA <- movingDR2(res0pars$conn0, res0pars$diff, res0pars$sig)
    
    # Add parameter info
    drMA <- data.frame(drMA, pars[i,])
    
    # Return
    return(drMA)
}, mc.cores = 4))

# Add group breaks
gBreaks <- c(-Inf, 0.05, 0.1, 0.25, 0.50, 0.75, 0.9, Inf)
maPars$drBr <- cut(maPars$dr,
                      breaks = gBreaks,
                      labels = c("0.00", "0.05", "0.10", "0.25", "0.50", "0.75", "0.90"))

# Add color pallets
drCols <- colorRampPalette(c("grey70", "darkgreen", "red"))(7)
names(drCols) <- levels(maPars$drBr)

# Fix upper and lower 7 limits to be centered around 0
maPars$yDiff <- apply(maPars[, c("lowerDiff", "upperDiff")], 1, function(x) x[which.min(abs(x))])

## Subset for differences in background connectivity <= 10
maParsBG10 <- maPars[abs(maPars$diffbgInput) <= 0.1,]

## Get min and max difference connectivity
minDiff <- min(maParsBG10$yDiff)
maxDiff <- max(maParsBG10$yDiff)

invisible(lapply(c(TRUE, FALSE), function(bg) {
    
    # Subset for bgDiff = 0 and no mean correct
    maParsSub <- maParsBG10[maParsBG10$mean_correct == bg,]
    
    # Change filname
    if(bg) {
        TPRoutBG <- sub("_BG", "_BGcorrTRUE", TPRout)
    } else {
        TPRoutBG <- sub("_BG", "_BGcorrFALSE", TPRout)
    }
    
    # Plot
    png(TPRoutBG, height = 1000, width = 1100)
    print(ggplot(maParsSub, aes(x = lowerConn, y = yDiff, color = drBr)) +
        
        # Add points
        geom_point(shape = 15) +
        
        # Add horizontal line
        geom_hline(yintercept = 0, color = "white", size = 0.7) +
        
        # Create grid
        facet_grid(as.factor(diffbgInput) ~ sim_type + mdc_type) +
        
        # Assign colors
        scale_color_manual(name = "Power", values = drCols) +
        
        # Axis scales
        scale_y_continuous(name = "Difference") +
        scale_x_continuous(name = "Connectivity") +
        
        # Plot dimensions
        coord_cartesian(ylim = c(minDiff, maxDiff), xlim = c(0, 0.8), expand = TRUE) +
            
        guides(color = guide_legend(override.aes = list(size=4))) +
        
        # Theme
        theme_bw() +
        theme(legend.position = "left",
              strip.background = element_rect(fill = "white"),
              strip.text = element_text(size = 20),
              legend.text = element_text(size = 15),
              legend.title = element_text(size = 20),
              axis.text = element_text(size = 15),
              axis.title = element_text(size = 22)))
    dev.off()
    
}))

