#' Fit Multi-Environment Genomic Prediction Model
#' 
#' @description
#' This function fits a Mixed Model for Genomic Evaluation (MGE) using the \code{sommer} 
#' package. It accounts for fixed environmental effects and decomposes random genetic 
#' effects into a main effect (across environments) and a GxE interaction effect 
#' (environment-specific).
#' 
#' @param K Matrix. Genomic Relationship Matrix (nc x nc).
#' @param train_set Data Frame. The training data containing columns for 
#'   TBV (True Breeding Value), Obs (Phenotypes), Env (Environment), and Taxa (Individuals).
#' 
#' @return An object of class \code{mmes} containing fitted model parameters, 
#'   variance components, and predicted random effects (BLUPs).
#' 
#' @export
mmes_MGE = function(K, train_set){
  # Load the sommer package for mixed model solver
  library(sommer) %>% suppressPackageStartupMessages()
  
  # output: Fits the MME using the sommer solver
  output = mmes(
    # Fixed Effect: Accounting for the mean of each environment
    Obs ~ Env,  
    # Random Effects: 
    # 1. Main genetic effect (ism(Taxa) with relationship matrix K)
    # 2. GxE interaction effect (dsm(Env) allows environment-specific genetic variance)
    random = ~ vsm(ism(Taxa), Gu = K) + vsm(dsm(Env), ism(Taxa), Gu = K),  
    # Residual structure: Diagonal residual matrix per environment
    rcov = ~ vsm(dsm(Env), ism(units)),  
    data = train_set,  
    verbose = FALSE,
    nIters = 10000,
    dateWarning = F
  )
  return(output)
} %>% cmpfun()

#' Extract and Process Predictions from MGE Model
#' 
#' @description
#' This wrapper processes the output of \code{mmes_MGE}. It calculates the GEBVs 
#' (Genomic Estimated Breeding Values) by combining fixed effects (means) and 
#' random effects (BLUPs). It organizes results by environment and calculates 
#' rankings for selection accuracy evaluation.
#' 
#' @param K Matrix. Genomic Relationship Matrix.
#' @param Nc Integer. Total number of individuals in the population.
#' @param output Object. The fitted model object from \code{mmes_MGE}.
#' @param taxa Character Vector. Names of the individuals/taxa.
#' @param env Integer. Number of environments/trials.
#' @param Trait List. Contains the population metadata, including TBV and Env indicators.
#' @param g_true List. Simulated true genetic effects for validation.
#' @param r Integer. Index for the specific variance scenario (from RandomEffects).
#' @param iter Integer. Current simulation iteration.
#' 
#' @return List. A nested list where each element represents an environment (and one 
#'   for the "Overall" mean) containing TBV, GEBV, simulated g, and predicted g.
#' 
#' @export
mmes_output = function(K, Nc, output, taxa, env, Trait, g_true, r, iter){
  
  # Helper function: Calculates the average of a specific column across all environments
  mean_col = function(col) { Reduce(`+`, lapply(Env, function(df) df[[col]])) / length(Env) }
  
  #' 1. Extract Fixed Effects and calculate Environmental Means
  # b[1] is the intercept (Env1), others are deviations from Env1
  b = output[["b"]] 
  mean_test = b[1]
  for(e in 2:env) { mean_test[e] = b[1] + b[e] }
  
  #' 2. Extract and Reconstruct Random Effects (BLUPs)
  # gh extracts the random effects 'u' from the sommer output
  # Column 1: Main genetic effect; Columns 2 to (env+1): GxE effects
  gh = matrix(output[["u"]], nrow = Nc, ncol = env+1, byrow = F) %>% as.data.frame()
  
  # g_pred: Predicted genetic effect (Main + GxE interaction)
  g_pred = sapply(seq(env), function(e){ gh[,1] + gh[,e+1] }) %>% data.frame()
  
  # GEBV: Genomic Estimated Breeding Value (Predicted Genetic Effect + Environmental Mean)
  GEBV = sapply(seq(env), function(e){ g_pred[,e] + mean_test[e] }) %>% data.frame()
  
  # Assign row names for identification
  row.names(gh) = row.names(g_pred) = row.names(GEBV) = taxa
  
  #' 3. Organize results into a list for each environment
  Env = lapply(seq(env), function(e) { 
    df = data.frame(
      TBV = Trait$TBV[which(Trait$Env == e)], # True Breeding Value from data
      GEBV = GEBV[, e],                       # Final prediction (g + mean)
      g = g_true[[r]][Nc * (e - 1) + seq(Nc), iter], # Simulated true g
      gh = g_pred[, e],                       # Predicted random genetic effect
      Rank_T = rank(-Trait$TBV[which(Trait$Env == e)]), # Rank based on TBV
      Rank_G = rank(-GEBV[, e])                # Rank based on GEBV
    )
    row.names(df) = taxa
    return(df)
  })
  
  #' 4. Aggregate results for the "Overall" performance across environments
  Env$Overall = data.frame(
    TBV = mean_col("TBV"),
    GEBV = mean_col("GEBV"),
    g = mean_col("g"),
    gh = mean_col("gh"),
    Rank_T = rank(- mean_col("TBV") ),
    Rank_G = rank(- mean_col("GEBV") )
  )
  row.names(Env$Overall) = taxa
  
  return(Env)
} %>% cmpfun()

#' Calculate Normalized Discounted Cumulative Gain (NDCG)
#' 
#' @description
#' This function calculates the NDCG to evaluate the ranking accuracy of a selection process.
#' It compares the relevance of individuals ranked by GEBV (Genomic Estimated Breeding Values) 
#' against the ideal ranking based on TBV (True Breeding Values). A logarithmic 
#' discount is applied to penalize lower-ranked individuals in the top selection set.
#' 
#' @param input Data Frame. A data frame where the first column is the "True" value 
#'   (e.g., TBV or simulated genetic effect) and a column named \code{GEBV} 
#'   contains the predicted values.
#' @param top Integer. The number of top individuals to consider (e.g., top 5% of the population).
#' 
#' @return Numeric. The NDCG score ranging from 0 to 1, where 1 indicates an ideal ranking.
#' 
#' @export
ndcg = function(input, top){
  
  #' 1. Extract the actual selection based on predicted ranking (GEBV)
  # Sort individuals by GEBV and select the top 'n' individuals
  f = input[order(input$GEBV, decreasing = TRUE), ] %>% .[seq(top), ]
  
  #' 2. Extract the ideal selection based on true ranking (TBV)
  # Sort individuals by the true genetic values (first column) to get the ideal DCG
  f0 = input[order(input[, 1], decreasing = TRUE), ] %>% .[seq(top), ]
  
  #' 3. Calculate the Logarithmic Discount Factor
  # The discount factor decreases as the rank position increases
  # Formula: 1 / log2(rank + 1)
  d0 = 1/log(seq(top) + 1, base = 2) 
  
  #' 4. Compute Normalized DCG
  # DCG (Discounted Cumulative Gain) is the sum of true values weighted by the discount
  # NDCG = DCG / IDCG (Ideal Discounted Cumulative Gain)
  NDCG = sum(f[, 1] * d0) / sum(f0[, 1] * d0)
  
  return(NDCG)
} %>% cmpfun()

#' Evaluate Methods via Simulation
#' 
#' @description
#' This function performs a comprehensive simulation study to evaluate different 
#' training set selection methods. It simulates phenotypes based on pre-generated 
#' genetic and environmental effects, fits Multi-Environment Mixed Models (MGE), 
#' and records various performance indicators including predictive correlation, 
#' NDCG, and CD values.
#' 
#' @param K Matrix. Genomic Relationship Matrix (nc x nc).
#' @param datasets List. A pre-initialized nested list structure to store result metrics.
#' @param g_true List. Simulated true genetic effects from \code{RandomEffects()}.
#' @param e_true List. Simulated true environmental effects from \code{RandomEffects()}.
#' @param mean Numeric Vector. Population means for each environment.
#' @param iter Integer Vector. Sequence of simulation replicates to run.
#' @param sigma_G Numeric. Variance of the genetic main effect.
#' @param sigma_G1 Numeric. Genetic variance for environments 1 to (env-1).
#' @param sigma_G2 Numeric Vector. A range of genetic variances for the final environment.
#' @param num_train Matrix. Rows are environments, columns are different training size scenarios.
#' @param DatasetName Character. Name of the dataset for saving results.
#' @param method Character Vector. Names of the selection methods being compared.
#' @param file_path Character. File name for the output .RData file.
#' @param iter_n,iter_r,iter_met Integer Vectors. Indices for iterating through training 
#'   sizes, variance scenarios, and selection methods.
#' 
#' @return List. The updated \code{datasets} list containing all recorded performance metrics.
#' 
#' @export
simulation = function(K, datasets, g_true, e_true, mean = c(100, 150, 200), iter = seq(2000),
                      sigma_G = 10, sigma_G1 = 10, sigma_G2 = c(0, 10, 20), num_train, DatasetName,
                      method =  c("Random", "CD_mean_V2", "CD_mean_MET"), file_path = "Simulation_Results", 
                      iter_n = seq(ncol(num_train)), iter_r = seq(rho), iter_met = seq(method)){
  
  # 1. Initialize Simulation Parameters
  env = nrow(num_train)
  Nc = nrow(K)
  mu = rep(mean[seq(env)], each = Nc) # Expand environment means to population size
  taxa = row.names(K)
  
  # Setup progress monitoring
  suppressPackageStartupMessages({
    library(httr) 
    library(progress)
  })
  n_values = apply(num_train, 2, paste, collapse = " & ")
  
  # 2. Main Simulation Loops
  for (n in iter_n){       # Loop through Training Set Sizes
    for (r in iter_r){     # Loop through Variance Scenarios (GxE levels)
      
      cat("\n\n'", "Dataset : ", DatasetName, 
          "    Training Set Size = ", n_values[n], 
          "    sigma_G2 = ", (sigma_G2[r] + sigma_G),
          "'\n", sep = "")
      flush.console()
      
      pb = progress_bar$new(format = " [:bar] :current/:total :percent ; Time: :elapsedfull", 
                            clear = FALSE, width = 100, total = length(iter))
      
      for (g in iter){     # Loop through Simulation Replicates
        used_samples = c() 
        attempts = 0       
        
        repeat {
          if (attempts >= 10) { break }
          attempts = attempts + 1
          
          tryCatch({
            # Randomly pair a genetic effect with an environmental effect sample
            sample_idx = sample(setdiff(seq(2000), used_samples), 1)
            used_samples = c(used_samples, sample_idx)
            
            # Construct Phenotypic Data (Trait)
            # Obs = mu + g + e
            Trait = data.frame(TBV = mu + g_true[[r]][,g] ,
                               Obs = mu + g_true[[r]][,g] + e_true[[r]][, sample_idx] ,
                               Env = rep( seq(env), each = Nc ) %>% as.factor() ,
                               Taxa = rep(taxa, times = env) %>% as.factor())
            
            # Biological constraint: ensure no negative phenotypes (if applicable)
            if (any(Trait$Obs < 0)) { next }
            
            for (met in iter_met){ # Loop through Selection Methods
              
              # Step A: Select individuals for the Training Set based on method
              if(met == 1){ # Random Selection
                train_taxa = lapply(seq(env), function(e){ sample(taxa, size = num_train[e, n]) %>% sort() })  
              } else if(met == 2){ # Optimized via CD_mean_v2
                train_taxa = lapply(seq(env), function(e){ taxa[ Train_CD_mean_v2[[n]][[1]][,e] ] %>% sort() })
              } else if(met == 3){ # Optimized via CD_mean_MET
                train_taxa = lapply(seq(env), function(e){ taxa[ Train_CD_mean_MET[[n]][[1]][,e] ] %>% sort() })
              }
              
              # Subsetting training data from the population
              train_set = lapply(seq(env), function(e){ 
                Trait[ which(Trait$Taxa %in% train_taxa[[e]] & Trait$Env == e) , ]
              }) %>% do.call(rbind, .)
              train_set$Taxa = as.factor(train_set$Taxa)
              
              # Step B: Model Fitting and Prediction
              MGE = mmes_MGE(K, train_set) # Fit Mixed Model
              Env = mmes_output(K, Nc, MGE, taxa, env, Trait, g_true, r, g) # Extract BLUPs
              
              # Step C: Record Performance Indicators
              # 1. CD (Coefficient of Determination) Indicators
              for (e in seq(env)){ 
                datasets$`CD[mean(v2)]`[[e]][[n]][[r]][g, met] <- 
                  CD(K, train_taxa, index = "CD_mean(v2)", var_G = sigma_G, 
                     var_ax1 = c(rep(sigma_G1, (env-1)), sigma_G2[r]), trial = e) 
                datasets$`CD[mean.MET]`[[e]][[n]][[r]][g, met] <- 
                  CD(K, train_taxa, index = "CD_mean.MET", var_G = sigma_G, 
                     var_ax1 = c(rep(sigma_G1, (env-1)), sigma_G2[r]), trial = e) 
              }
              
              # 2. Accuracy and Ranking Metrics
              for (e in seq(env+1)){ # Trials + Overall
                top = Nc * 0.05 # Define selection threshold (Top 5%)
                
                # Predictive Squared Correlation (Reliability)
                datasets$`r^2` [[e]][[n]][[r]][g, met] <- cor(Env[[e]]$g, Env[[e]]$gh)^2
                # Normalized Discounted Cumulative Gain (Ranking Quality)
                datasets$NDCG  [[e]][[n]][[r]][g, met] <- ndcg(Env[[e]], top)
                # Predictive Correlation (Accuracy)
                datasets$` `    [[e]][[n]][[r]][g, met] <- cor(Env[[e]]$TBV, Env[[e]]$GEBV)
                
                # Top-Slice Metrics (Metrics focusing only on the selected top 5%)
                GEBV_top = Env[[e]][order(Env[[e]]$GEBV, decreasing = TRUE), ] %>% .[seq(top),]
                datasets$Cor   [[e]][[n]][[r]][g, met] <- cor(GEBV_top$TBV, GEBV_top$GEBV)
                datasets$SRC   [[e]][[n]][[r]][g, met] <- cor(GEBV_top$Rank_T, GEBV_top$Rank_G)
                datasets$`RS[ratio]`[[e]][[n]][[r]][g, met] <- sum(GEBV_top$Rank_G)/sum(GEBV_top$Rank_T)
              }
            }  
            # Break repeat loop if model converged successfully
            if ( !all(MGE$u == 0) ) {break}
          },
          error = function(e) { }) # Catch potential matrix singularity errors
        }
        
        pb$tick() 
        # Periodic save to prevent data loss during long simulations
        save(datasets, file = paste0("數據/", DatasetName, "/", file_path, ".RData"))
      }
      if (!pb$finished) { pb$terminate() }
    }
  }  
  
  return(datasets)
} %>% cmpfun()