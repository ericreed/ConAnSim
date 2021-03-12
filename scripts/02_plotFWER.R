require(tidyr)
require(ggplot2)

## Directories
### /restricted/projectnb/montilab-p/personal/eric/ConAnSim/scripts
baseDir <- ".."
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
FDRout <- file.path(plotDir, "FWER_ConAn.png")

## Read in results file
res <- do.call(rbind, lapply(simResFiles, function(file) readRDS(file)))

## Add discovery logical (pval < 0.05)
res$sig <- res$pval < 0.05

## Subset for diff = 0
res0 <- res[res$diffInput == 0,]

## Get unique parameters
pars <- unique(res0[, c("sim_type", "mean_correct", "mdc_type", "diffbgInput")])

## Create moving average of discovery rates for each
maPars <- do.call(rbind, lapply(seq(nrow(pars)), function(i) {
    
    # Subset for par estimates
    res0pars <- res0[res0$sim_type == pars$sim_type[i] &
                     res0$mean_correct == pars$mean_correct[i] &
                     res0$mdc_type == pars$mdc_type[i] &
                     res0$diffbgInput == pars$diffbgInput[i],
                     ]
    
    # Get moving average discovery rate for set of parameters
    drMA <- movingDR(res0pars$conn0, res0pars$sig)
    
    # Add parameter info
    drMA <- data.frame(drMA, pars[i,])
    
    # Return
    return(drMA)
}))

# Create tranformation function for y axis
minLCI <- -(min(maPars$lCI, na.rm = TRUE) - 1E-04)
trans_y <- scales::trans_new(name = "trans_y",
                        transform = function(x) abs(x)^(1/2)*sign(x),
                        inverse = function(x) x^2 * sign(x))


## Subset for differences in background connectivity <= 10
maParsBG10 <- maPars[abs(maPars$diffbgInput) <= 0.1,]

# Plot
png(FDRout, height = 1000, width = 1300)
ggplot(maParsBG10, aes(x = lower, y = dr, group = mean_correct, color = mean_correct)) +
    
    # Target DR
    geom_hline(yintercept = 0.05, alpha = 0.5, size = 2) +
    geom_hline(yintercept = 0) +
    
    # Discovery rate
    geom_line(size = 1.3) +
    
    # CIs
    geom_segment(aes(x = lower, xend = lower, y = lCI, yend = uCI), size = 1.7, alpha = 0.3) +
    geom_line(aes(x = lower, y = uCI), alpha = 0.2, size = 1) +
    geom_line(aes(x = lower, y = lCI), alpha = 0.2, size = 1) +
    
    # Create grid
    facet_grid(as.factor(diffbgInput) ~ sim_type + mdc_type) +
    
    # Assign colors
    scale_color_manual(name = "Background\nCorrection", values = c("TRUE" = "red4", "FALSE" = "blue")) +
    
    # Axis scales
    scale_y_continuous(name = "Family-wise Error Rate", trans = trans_y, breaks = c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
    scale_x_continuous(name = "Connectivity", breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)) +
    
    # Plot dimensions
    coord_cartesian(ylim = c(0, 1), xlim = c(0, 0.69), expand = TRUE) +
    
    # Theme
    theme_bw() +
    theme(legend.position = "left",
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 20),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 22))
dev.off()

