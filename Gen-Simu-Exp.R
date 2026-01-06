#' Simulate True Genetic and Environmental Effects for Multi-Environment Trials
#' 
#' @description
#' This function simulates true genetic effects ($g$) and residual environmental effects ($e$) 
#' based on a Multi-Environment Genomic Prediction model. It uses the Kronecker product 
#' of environmental covariance matrices and the genomic relationship matrix ($K$) to 
#' generate multivariate normal distributions across different GxE scenarios.
#' 
#' @param K Matrix. Genomic Relationship Matrix ($nc \times nc$).
#' @param env Integer. Number of environments/trials.
#' @param OmegaG_fixed Numeric. Baseline common genetic covariance between environments.
#' @param OmegaG_first Numeric. Baseline genetic variance for the first $env-1$ environments.
#' @param OmegaG_last Vector. A range of genetic variance values to test for the final environment.
#' @param h2 Numeric. Heritability ($h^2$), used to scale the residual variance relative to genetic variance.
#' @param n Integer. Number of simulation replicates (iterations).
#' @param DatasetName Character. Name of the dataset for file organization.
#' @param core Numeric. Ratio of CPU cores to utilize (defaulting to 0.6).
#' @param file_path Character. The file name/path for saving the generated effects.
#' 
#' @return List. A nested list containing:
#'   \itemize{
#'     \item \code{[[1]]}: Simulated true genetic effects (\code{g_true}) for each scenario.
#'     \item \code{[[2]]}: Simulated true environmental effects (\code{e_true}) for each scenario.
#'   }
#' 
#' @export
RandomEffects = function(K, env = 3, OmegaG_fixed = 10, OmegaG_first = 10, OmegaG_last = c(0, 10, 20), h2 = 0.5,
                         n = 2000, DatasetName, core = 0.6, file_path = "/RandomEffetcts.RData"){
  
  # Initialize lists to store simulated results
  g_true = e_true = list()
  
  # Load necessary libraries for simulation and parallelization
  suppressPackageStartupMessages({
    library(MASS)      # For mvrnorm (multivariate normal distribution)
    library(parallel)  # For parallel computing
    library(pbapply)   # For progress bars in apply functions
  })
  
  ## 1. Parallel Computing Setup
  myCoreNums = detectCores()
  cl <- makeCluster(14) # Configure cluster with a specific number of workers
  
  # Initialize worker nodes with required packages and functions
  invisible(clusterEvalQ(cl, suppressPackageStartupMessages({
    library(parallel); library(tidyverse); library(magrittr); 
    library(MASS); library(pbapply); library(compiler)
  })))
  
  ## 2. Simulation Loop for Different Variance Scenarios
  # Iterate through each variance value provided in OmegaG_last
  for (r in seq(OmegaG_last) ){
    
    # Console output for tracking current simulation status
    cat("### ", DatasetName, " ; OmegaG_last = ", (OmegaG_last[r]+OmegaG_fixed[1]), " ###\n", sep = "")
    
    # Construct the Genetic Covariance Matrix (Omega_G) for environments
    # Combines baseline fixed covariance and environment-specific variances
    Omega_G = matrix(OmegaG_fixed, env, env) + diag(c( rep(OmegaG_first, env-1) , OmegaG_last[r]))
    
    # Calculate Residual Covariance Matrix (omega_e) based on defined heritability
    # Formula derived from h^2 = Var(G) / (Var(G) + Var(E))
    omega_e = (diag(OmegaG_fixed, env) + diag(c( rep(OmegaG_first, env-1) , OmegaG_last[r]))) * (1-h2)/h2
    
    # Compute full Kronecker covariance matrices
    # G: (nc*env) x (nc*env) genetic covariance
    # E: (nc*env) x (nc*env) residual/environmental covariance
    G = Omega_G %x% K
    E = omega_e %x% diag(nrow(K))
    
    # Export current covariance matrices to the parallel cluster
    clusterExport(cl, c("G", "E"), envir = environment())
    
    ## 3. Simulating Genetic Effects (g_true)
    # Sampling from Multivariate Normal: g ~ MVN(0, Omega_G ⊗ K)
    cat("    Simulated Genetic Effects\n")
    g_true[[r]] = pbsapply(seq(n), cl = cl, function(i) {
      mvrnorm(n = 1, mu = rep(0, nrow(G)), Sigma = G)
    }) %>% as.data.frame()
    
    # Intermediate save to prevent data loss
    save(g_true, e_true, file = paste0("數據/", DatasetName, file_path))
    
    ## 4. Simulating Environmental Effects (e_true)
    # Sampling from Multivariate Normal: e ~ MVN(0, omega_e ⊗ I)
    cat("    Simulated Environmental Effects\n")
    e_true[[r]] = pbsapply(seq(n), cl = cl, function(i) {
      mvrnorm(n = 1, mu = rep(0, nrow(E)), Sigma = E)
    }) %>% as.data.frame()
    
    # Final save for the current scenario
    save(g_true, e_true, file = paste0("數據/", DatasetName, file_path))
  }
  
  # Clean up: stop the cluster and return results
  stopCluster(cl)
  return(list(g_true, e_true))
  
} %>% cmpfun()