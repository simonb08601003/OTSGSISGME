#' Calculate CD(mean) for Multi-Environment Trials
#' 
#' @description
#' This function calculates the Coefficient of Determinations (CD) 
#' for a specific set of individuals. It is designed for optimizing the training set in 
#' Genomic Prediction.
#' 
#' The calculations derives the Prediction Error Variance (PEV) based on Mixed Effect Model.
#' It returns the mean CD for the **test set** (individuals NOT included in \code{soln}).
#' 
#' @param soln Integer Vector. The indices of individuals selected for the **Training Set**.
#' @param Data A named List containing the necessary matrices:
#'  \itemize{
#'    \item \code {K}: Genomic Relationship Matrix (nc x nc).
#'    \item \code {G}: The covariance matrix for true genotypic values of whole candidate ((nc*env)x(nc*env)).
#'    \item \code {R}: The covariance matrix for residuals.
#'    \item \code {E}: Design matrix for environmental effect.
#'  }
#'  
#'  @return Numeric. The mean CD value of the test set (individuals not in \code{soln}).
#'  
#'  @export
CDMEANOPTwithEnvME<-function(soln, Data){
  
  #' 1. Extract matrices
  G<-Data[["G"]]
  R<-Data[["R"]]
  E<-Data[["E"]] # Design matrices for fixed effects
  
  #' 2. Compute Inverse of the Phenotypic Variance Matrix (V) for the training set
  Vinv<-solve(G[soln,soln]+R[soln,soln])
  
  #' 3. Calculate CD Matrix
  outmat<-(G[,soln]%*%(Vinv-Vinv%*%E%*%solve(t(E)%*%Vinv%*%E)%*%t(E)%*%Vinv)%*%G[soln,])/G
  
  #' 4. Return mean CD of the test set
  return(mean(diag(outmat[-soln,-soln])))
}

#' Wrapper for CDMEANOPTwithEnvME 
#' 
#' @description
#' A wrapper function that handles the subsetting of the design matrix (E)
#' to match the selected training set (soln) before calling the core calculation function.
#' 
#' @inheritParams CDMEANOPTwithEnvME
#' @export
myCDMEANOPTwithEnvME <- function(soln, Data){
  DataCopy <- Data
  # Critical Step: Subset E to match the dimensions of Vinv in the core function.
  DataCopy[["E"]] <- Data[["E"]][soln, , drop = FALSE]
  return(CDMEANOPTwithEnvME(soln = soln, Data = DataCopy))
}

#' Calculate CD(mean) version 2 for Multi-Environment Trials
#' 
#' @description
#' This function calculates the CD(mean) for a selected
#' training set in a multi-environment context. 
#' 
#' Unlike the standard mixed effect model approach (which inverts V = G + R), 
#' this version explicitly constructs a projection matrix (\code{Mt}) based on the 
#' environmental variances and the design of the training set. 
#' This approach is often computationally more efficient when \eqn{nt \ll Nc} or when
#' specific assumption about the fixed effects structure are made.
#' 
#' @param soln Integer Vector. The **row indices** of the selected individuals for the training set.
#'  \strong{Note:} In a stacked data structure, these indices must correspond to the 
#'   specific rows in the stacked matrices (Genotype + Environment combinations).
#' @param Data A named List containing the necessary matrices and parameters:
#'   \itemize{
#'     \item \code{OmegaG}: Genetic covariance matrix between environments (Env x Env).
#'     \item \code{VarE}: Vector of residual variances for each environment (Length = Env).
#'     \item \code{K}: Genomic Relationship Matrix (Nc x Nc).
#'     \item \code{G}: The full genetic covariance matrix ((Nc*Env) x (Nc*Env)).
#'     \item \code{R}: The full residual covariance matrix (unused in this specific logic but good to have).
#'     \item \code{E}: The design matrix for environmental effects ((Nc*Env) x Env). 
#'           Columns should represent environments (0/1 coding).
#'   }
#'   @return Numeric. The mean CD value across the entire population (or test set),
#'   calculated as the mean ratio of \eqn{diag(Var(\hat{g})) / diag(Var(g))}.
#'   
#'   @export
CDmeanv2 <- function(soln, Data){
  OmegaG <- Data[["OmegaG"]]
  VarE <- Data[["VarE"]]
  K <- Data[["K"]] # Nc*Nc
  G <- Data[["G"]] # (Nc*Env)*(Nc*Env)
  R <- Data[["R"]] # (Nc*Env)*(Nc*Env)
  E <- Data[["E"]] # (Nc*Env) * Env
  
  Nc <- nrow(K) # The number of whole candidate population
  n <- c() # The number of of training set varieties
  env <- ncol(E) # The number of environments 
  
  soln <- sort(soln)
  Esub <- E[soln, , drop=FALSE]
  n <- colSums(Esub)
  
  # The empty matrices of Mt
  Mt <- matrix(0, nrow = length(soln), ncol = length(soln))
  
  # Fill in the elements of the matrices Mt
  current_idx <- 0
  for(i in seq_len(env)){
    ni <- n[i]
    if(ni > 0){
      idx_range <- (current_idx+1):(current_idx+ni)
      SigmaSq <- VarE[i]
      
      Mt[idx_range, idx_range] <- (1/SigmaSq)*(diag(1, ni, ni)-matrix(1/ni, ni, ni))
      current_idx <- current_idx+ni
    }
  }
  
  # Fill in the Gt and Gct
  Gt <- G[soln, soln]; Gct <- G[, soln]
  
  # matrices A and B 
  A <- Gct %*% solve(Mt %*% Gt + diag(1, sum(n), sum(n))) %*% Mt %*% t(Gct) # A: variance matrix of g_hat
  B <- kronecker(OmegaG, K, make.dimnames = TRUE) # B: variance matrix of g
  
  # Calculate CDmean(v2)
  CDvalue <- mean(diag(A)/diag(B))
  return(CDvalue)
}