#' RSAVS: A package for robust subgroup analysis and variable selection.
#' 
#' The RSAVS can do subgroup analysis and variable selection simultaneously and 
#' it supports multiple choices of loss and penalty functions.
#' It implements the computation in a parallel manner.
#' 
#' @section RSAVS functions:
#' RSAVS_Compute_BIC
#' RSAVS_Determine_Mu
#' RSAVS_Further_Improve
#' RSAVS_Generate_D_Matrix
#' loss_function: RSAVS_Huber, RSAVS_L1, RSAVS_L2
#' RSAVS_LargeN
#' RSAVS_Mu_to_Mat
#' RSAVS_RI
#' RSAVS_Summary_Iteration
#' RSAVS_S_to_Groups
#' 
#' @section RSAVS datesets:
#' Student performance dataset: full_df, mat_df, por_df
#' 
#' @docType  package
#' @name RSAVS
NULL