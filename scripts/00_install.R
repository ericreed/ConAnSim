# Run on the SCC with R v3.6.0

# K2Taxonomer package
devtools::install_github("sac4vf/ConAn")

# Functions for running simulations
install.packages(c("reticulate", "Matrix", "BDgraph", "corpcor", "MASS", "igraph"), repos='http://cran.us.r-project.org')
