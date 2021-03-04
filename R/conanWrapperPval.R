# Function to run conan for each way to run it

## Arguments
### eset = ExpressionSet
### mod_list = Named list of vectors of feature IDs for each module
### iter = Number of iterations to run
### cores = Number of cores to use
### seed = seed to use

## Value
### A data frame with p-values for each module and combination of conan parameters

conanWrapperPval <- function(eset, mod_list, iter, cores, seed, N_genes, iter_bg) {

    # Create grid of parameters
    parComs <- t(expand.grid(c("bootstrap","permutation"),
                           c(FALSE, TRUE),
                           c("fraction", "difference")))

    # Run conan on grid of parameters
    do.call(rbind, lapply(seq(ncol(parComs)), function(i, eset, mod_list, iter, cores, seed, N_genes, iter_bg) {

        # Get set of parameters
        x <- parComs[,i]

        RNGkind("L'Ecuyer-CMRG")
        set.seed(seed)
        res <- quiet(conan(eset,
              mod_list,
              "group", "g0", "g1",
              sim_type = x[1],
              iter = iter,
              mean_correct = as.logical(x[2]),
              cores = cores,
              N_genes = N_genes,
              iter_bg = iter_bg,
              mdc_type = x[3]))

        return(suppressWarnings(data.frame(sim_type = x[1],
                          mean_correct = as.logical(x[2]),
                          mdc_type = x[3],
                          mod = names(mod_list),
                          pval = res$significance$mdc_pval,
                          cv0 = res$stat$mods_mc_r,
                          cv1 = res$stat$mods_mc_t)))

    }, eset, mod_list, iter, cores, seed, N_genes, iter_bg))

}

quiet <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
}
