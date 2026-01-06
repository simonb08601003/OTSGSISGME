#' Calculate CD(mean) for Multi-Environment Trials
#' 
#' @description
#' This function calculates the Coefficient of Determination (CD) mean for a Multi-Environment 
#' Genomic Prediction model. It accounts for both main genetic effects and Genotype-by-Environment 
#' (GxE) interactions. It supports two different calculation indices (v2 and MET) and can 
#' target either the entire population or a specific trial.
#' 
#' @param K Matrix. Genomic Relationship Matrix (nc x nc) for the whole candidate population.
#' @param train List. A list where each element contains the indices of individuals selected for 
#'   the training set in that specific environment (e.g., \code{train[[1]]} for Env 1).
#' @param index Character. The CD criterion to calculate: "CD_mean(v2)" or "CD_mean.MET".
#' @param var_G Numeric. Variance component for the genetic main effect ($\sigma_G^2$).
#' @param var_ax1 Vector. Specific genetic variances for each environment ($\sigma_{G \times E}^2$).
#' @param var_E Vector. Residual variances for each environment ($\sigma_E^2$). 
#'   If NULL, defaults to the diagonal of \code{Omega_G}.
#' @param trial Integer. The specific trial index to calculate CD for. Set to 0 for all trials.
#' 
#' @return Numeric. The calculated mean CD value.
#' 
#' @export
CD = function(K, train, index = "CD_mean(v2)",
              var_G = 10, 
              var_ax1 = c(rep(10, length(train)-1), 20), 
              var_E = NULL, 
              trial = 0){
  
  #' 1. Initialize parameters and Variance-Covariance structures
  env = length(train) # Number of environments/trials
  
  # Construct the Genetic Covariance Matrix between environments (Omega_G)
  # This includes the common genetic variance and environment-specific variances
  Omega_G = matrix(var_G, env, env) 
  diag(Omega_G) = diag(Omega_G) + var_ax1 
  
  # Set default residual variance if not provided
  if( is.null(var_E) ){var_E = diag(Omega_G)} 
  
  nc = nrow(K) # Total number of candidate individuals
  n = c()      # Vector to store training set size per environment
  for(e in seq(env)){n[e] = length(train[[e]])}
  
  #' 2. Initialize matrices Gt, Gct, and Mt
  # Gt: Genetic covariance within the training set
  # Mt: Projection/Centering matrix for the training set
  # Gct: Genetic covariance between the whole population and the training set
  Gt = Mt = matrix(0, nrow = sum(n), ncol = sum(n))
  Gct = matrix(0, nrow = nc*env, ncol = sum(n))
  
  # Fill in the elements of the matrices using block-wise logic
  for(i in seq(env)){
    # Calculate row/column index range for the current environment i in the training set
    idx_i = sum(n[seq(i)]) - n[i] + seq(n[i])
    
    for(j in seq(env)){
      # Calculate index range for environment j
      idx_j = sum(n[seq(j)]) - n[j] + seq(n[j])
      # Index range for the whole population in environment i
      idx_all = nc*(i-1) + seq(nc)
      
      # Fill Genetic matrices based on the GxE Kronecker structure
      Gt[idx_i, idx_j] = Omega_G[i,j] * K[train[[i]], train[[j]]] 
      Gct[idx_all, idx_j] = Omega_G[i,j] * K[, train[[j]]] 
    }
    
    # Fill Mt with the centering operator (I - 1/n) weighted by residual variance
    # This step effectively accounts for fixed environmental effects
    ni = n[i]
    Mt[idx_i, idx_i] = 1/var_E[i] * (diag(1, ni, ni) - matrix(1/ni, ni, ni)) 
  }
  
  #' 3. Calculate Variance of Predictions (A) and True Genetic Variance (B)
  # Matrix A: Var(g_hat) derived from the Mixed Model Equations logic
  A = Gct %*% solve(Mt %*% Gt + diag(1, sum(n), sum(n))) %*% Mt %*% t(Gct)
  # Matrix B: Var(g) = Omega_G %x% K
  B = Omega_G %x% K 
  
  #' 4. Return the CD value based on the requested index and trial scope
  if (trial == 0) { # Mean across all trials
    
    if(index == "CD_mean(v2)"){ 
      # Standard CD mean: average of the ratio of diagonal elements
      CD_value = mean(diag(A) / diag(B)) 
      
    } else if (index == "CD_mean.MET") { 
      # MET CD mean: aggregates variances of the same variety across environments
      A_all = B_all = 0
      for (i in seq(env)) {
        i_all = nc*(i-1) + seq(nc)
        for (j in seq(env)) {
          j_all = nc*(j-1) + seq(nc)
          
          # Sum diagonal elements of the sub-matrices for variety-wise aggregation
          A_sub = diag(A[i_all, j_all])
          B_sub = diag(B[i_all, j_all])
          
          A_all = A_all + A_sub
          B_all = B_all + B_sub
        }
      }
      CD_value = mean(A_all / B_all)
    }
    
  } else if (trial %in% seq(env)) { 
    # Calculate CD mean for a specific single trial
    idx_trial = seq(nc) + nc*(trial-1)
    CD_value = mean( c(diag(A) / diag(B))[idx_trial] )
  }
  
  return(CD_value)
} %>% cmpfun() # Pre-compile function for performance optimization


#' Generate Optimal Training Set via Genetic Algorithm for Multi-Environment Trials
#' 
#' @description
#' This function implements a Genetic Algorithm (GA) to optimize the selection of a training 
#' set for Genomic Prediction. It aims to maximize the Coefficient of Determination (CD) 
#' across multiple environments by iteratively evolving a population of potential 
#' training set solutions through evaluation, selection, crossover, and mutation.
#' 
#' @param K Matrix. Genomic Relationship Matrix (nc x nc).
#' @param DatasetName Character. Name of the dataset for logging and file naming.
#' @param n Integer Vector. Target training set sizes for each environment.
#' @param option Character. "Same" for identical training sets across trials; 
#'   "Diff" for environment-specific optimization.
#' @param n_iter Integer. Maximum number of iterations for the GA.
#' @param index Character. The CD criterion to maximize (e.g., "CD_mean(v2)").
#' @param new Logical. If TRUE, starts a new GA; if FALSE, resumes from a saved .RData file.
#' @param file_path Character. Directory path for saving progress and results.
#' 
#' @return List. Contains:
#'   \itemize{
#'     \item \code{Opt.set}: Data frame of the optimized training set indices.
#'     \item \code{CD}: Data frame recording the maximum CD scores per iteration.
#'   }
#' 
#' @export
GA = function(K, DatasetName, n = c(50, 100), option = "Diff", 
              n_iter = 13000, index, new = T, 
              file_path = paste0("數據/", DatasetName, "/", index, " ", "Number = ", paste(n, collapse = " & "))){
  
  # Log initial parameters to console
  cat("\n'", "Dataset :", DatasetName, ";",
      "Training Set Size =", paste(n, collapse = " & "), ";", 
      "Index =", index, "'\n")
  
  ## 1. Initialization and Parameter Setup
  N = nrow(K)          # Population size
  cand = seq(N)        # Candidate indices
  env = length(n)      # Number of environments
  nt = mean(n)         # Average training size
  iter = 0             # Iteration counter
  
  # Determine if solutions are trial-specific or global
  if (option == "Same" & length(unique(n)) == 1) { 
    list_len = 1 
  } else if (option == "Diff") {
    list_len = seq(env) 
  } else {
    stop("Invalid option. Use 'Same' or 'Diff'.")
  }
  
  # Define Solution Size (Population Size for GA)
  if (nt < 50){ ss = 50 } else if (nt >= 50 && nt <= 200){ ss = nt } else { ss = 200 }
  
  # Initialize or Resume GA Process
  if(new == T){ 
    # Generate initial random binary solutions (1 = selected, 0 = not selected)
    sol = lapply(seq(ss), function(i) {
      sapply(list_len, function(j) {
        sample(c(rep(1, n[j]), rep(0, N - n[j]))) 
      }) %>% as.data.frame()
    })
    max.score = data.frame(matrix(NA, 1, env + 1))
    colnames(max.score) = c("Overall", paste0("Trial", seq_len(env)))
    
  } else { 
    # Resume from existing progress file
    load(paste0(file_path, ".RData"))
  }
  
  # Placeholders for tracking changes across iterations
  prev_sol = lapply(seq(ss), function(i) {
    sapply(list_len, function(j) {rep(0, N)}) %>% as.data.frame()
  })
  prev_score = matrix(0, nrow = ss, ncol = env+1)
  
  ## 2. Parallel Computing Configuration
  suppressPackageStartupMessages({
    library(httr)
    library(progress)
    library(parallel)
  })
  
  # Initialize progress bar
  pb <- progress_bar$new(format = " [:bar] :current/:total :percent ; Time: :elapsedfull ; Rate: :tick_rate iter/sec", 
                         total = n_iter, clear = FALSE, width = 100)
  
  # Setup parallel cluster (utilizing 14 cores)
  cl <- makeCluster(14) 
  invisible(clusterEvalQ(cl, suppressPackageStartupMessages({
    library(parallel); library(tidyverse); library(magrittr); library(compiler)
  })))
  
  # Export static objects to worker nodes
  clusterExport(cl, c("CD", "ss", "env", "cand", "K", "N", "n", "list_len", "index"), envir = environment())
  
  ### Start Main Genetic Algorithm Loop
  
  
  stop_flag = F 
  while(stop_flag == F){
    iter = iter + 1 
    
    # Identify which solutions changed in the last iteration to avoid redundant calculations
    changed_sol = which(sapply(seq(ss), function(i) !all(sol[[i]] == prev_sol[[i]])))
    
    # Export dynamic objects for current iteration
    clusterExport(cl, c("sol", "ss", "changed_sol", "prev_score", "max.score"), envir = environment())
    
    ## Evaluation Phase: Calculate CD values
    score = parLapply(cl, seq(ss), function(i) {
      if (i %in% changed_sol) { 
        # Convert binary vector to candidate indices
        train = lapply(seq(env), function(e) {
          cand[which(sol[[i]][, (e - 1) %% ncol(sol[[i]]) + 1] == 1)]
        })
        # Calculate CD for overall and environment-specific trials
        sapply(0:env, function(trial_id) {
          CD(K = K, train = train, index = index, trial = trial_id)
        })
      } else {
        prev_score[i, ] 
      }
    }) %>% do.call(rbind, .) %>% as.data.frame()
    
    colnames(score) = c("Overall", paste0("Trial", 1:env))
    
    # Update historical tracking
    prev_sol = sol
    prev_score = score 
    max.score[iter, ] = score[which.max(score$Overall), ] 
    
    # Selection Phase: Identify top 10% (Elite) and bottom 60% (To be replaced)
    elite = which(rank(-score$Overall, ties.method = "max") <= ceiling(ss * 0.1))
    del = which(rank(-score$Overall, ties.method = "min") >= ceiling(ss * 0.6))
    
    ### Termination Criteria
    th = 0.0001 # Convergence threshold
    if(iter >= n_iter){ stop_flag = T }
    
    # Stop if improvement plateaus over 500 iterations
    if(iter > 12000){
      gains <- as.numeric(max.score[iter, ] - max.score[iter - 500, ])
      if (all(gains >= 0 & gains < th)) { stop_flag <- TRUE }
    }
    
    # Sanity check: Ensure max score does not decrease (Optimization integrity)
    if(iter > 2 && (max.score$Overall[iter] < max.score$Overall[iter-1])){
      stop("Max score drops! Check logic.")
    }
    
    ## 3. Crossover Phase (Recombination)
    # Replace the "delete" group by crossing over traits from the "elite" group
    for(i in del){ 
      repeat {
        parents = sample(seq(ss)[elite], 2) 
        split_point = sample(seq(N - 1), 1) 
        for(e in seq(list_len)){ 
          sol[[i]][,e] = c(sol[[ parents[1] ]][1:split_point, e], 
                           sol[[ parents[2] ]][(split_point+1):N, e])
        }
        # Ensure offspring maintains correct training set size
        n_check = colSums(sol[[i]])
        if (all(n_check == n)) { break } 
      }
    }
    
    ## 4. Mutation Phase
    # Introduce random variation to maintain genetic diversity
    sol.new = sol 
    mut_count = ceiling(0.04 * n) # 4% mutation rate
    
    for(i in seq(ss)){
      for(e in seq(list_len)){
        # Swap 'r' number of 0s with 1s and vice versa
        pos = c(sample(which(sol[[i]][,e] == 0), mut_count[e]), 
                sample(which(sol[[i]][,e] == 1), mut_count[e])) 
        sol.new[[i]][pos, e] = abs(sol.new[[i]][pos, e] - 1) 
      }
    }
    
    # Elitism logic: Only keep mutations in the elite group if they improve the score
    clusterExport(cl, c("ss", "sol.new", "sol", "elite", "score", "n"), envir = environment())   
    sol = parLapply(cl, seq(ss), function(i){
      if (!(i %in% elite)) {
        sol.new[[i]] 
      } else {
        old_val = score$Overall[i]
        train.new = lapply(seq(env), function(e) {
          cand[which(sol.new[[i]][, (e - 1) %% ncol(sol.new[[i]]) + 1] == 1)]
        })
        new_val = CD(K = K, train = train.new, index = index, trial = 0)
        if (new_val > old_val) { sol.new[[i]] } else { sol[[i]] } 
      }
    })
    
    # Progress feedback and checkpointing
    pb$tick()
    save(iter, sol, max.score, file = paste0(file_path, ".RData")) 
    
    if(stop_flag == T) { cat("\nGenetic Algorithm has ended!\n\n") }
  } 
  
  ## 5. Final Output Processing
  stopCluster(cl) 
  if (!pb$finished) { pb$terminate() }
  
  # Select the best solution found in the final iteration
  best_idx = which.max(score$Overall)
  opt.sol = sol[[best_idx[1]]]
  
  # Format the optimized training set for export
  opt.set = sapply(seq(env), function(e) {
    c(cand[which(opt.sol[, (e - 1) %% ncol(opt.sol) + 1] == 1)], rep(NA, max(n) - n[e]))
  }) %>% as.data.frame()
  
  return(list(Opt.set = opt.set, CD = max.score))
  
} %>% cmpfun()