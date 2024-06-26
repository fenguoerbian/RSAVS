# This file contains the pre-define functions for Robust Subgroup Analysis

#' Generate the pair-wise different matrix
#'
#' This function generate the pairwise difference matrix(D matrix in the paper).
#' By default it returns a sparse matrix(matrix.csr) from the package SparseM.
#'
#' In Original RSAVS via ADMM algorithm. This functions generates the pairwise
#'   difference matrix with \code{w_vec = 1}.
#'   
#' In RSI, the original objective function is converted to a sparse median regression
#'   via LLA. The pairwise difference matrix is coupled with weights extracted from
#'   the penalty functions. This is supplied into this function through \code{w_vec}.
#'   
#' @param n number of observations.
#' @param dense logical, whether the return type should be in dense matrix or not
#' @param w_vec a scalar value of a length \eqn{n * (n - 1) / 2} vector representing
#'   the weights for left part in each pair(e.g. \eqn{\mu_i} in \eqn{\mu_i - \mu_j}). 
#'   If \code{w_vec} is a scalar, it will be recycled for every pair. 
#'   Defaults to 1.
#' @return a difference matrix with size \eqn{(n * (n - 1) / 2) \times n}
#' @examples 
#' RSAVS:::RSAVS_Generate_D_Matrix(5)
#' RSAVS:::RSAVS_Generate_D_Matrix(4, dense = TRUE)
#' RSAVS:::RSAVS_Generate_D_Matrix(3, w_vec = 2)
#' RSAVS:::RSAVS_Generate_D_Matrix(4, dense = TRUE, w_vec = c(1, -2, 3, -4, 5, -6))
RSAVS_Generate_D_Matrix <- function(n, dense = FALSE, w_vec = 1){
  # This function generate the pairwise difference matrix(D matrix in the paper)
  # By default it returns a sparse matrix(matrix.csr) from the package SparseM
  # Args: n: number of observations. 
  #       dense: logical, whether the return type should be in dense matrix or not
  # Returns: res: a difference matrix with size (n * (n - 1) / 2) \times n
  
  # a sparse matrix with correct dimensions but non elements
  res <- as.matrix.csr(0, nrow = n * (n - 1) / 2, ncol = n)
  # compute the elements
  # Refer to the pdf manual of SparseM for detailed definitions of ra, ia and ja
  
  # ra: the non-zero elements, in row-major pattern
  # ra <- rep(c(1, -1), times = n * (n - 1) / 2)    # original design
  if(length(w_vec) == 1){
    w <- w_vec
    ra <- rep(c(w, -w), times = n * (n - 1) / 2)
  }else{
    if(length(w_vec) == (n * (n - 1) / 2)){
      tmp <- - w_vec
      ra <- rbind(w_vec, tmp)
      ra <- as.vector(ra)
    } else{
      stop("w_vec must be a scalar or a n * (n - 1) / 2 length vector!")
    }
  }
  
  # ia: starting point of each row in ra and ja
  ia <- (1 : (n * (n - 1) / 2 + 1)) * 2 - 1
  
  # ja: column index, according to ra
  positive_column_indicator <- rep(1 : (n - 1), times = (n - 1) : 1)
  negative_column_indicator <- unlist(lapply(2 : n, function(x) x : n))
  ja_matrix <- t(cbind(positive_column_indicator, negative_column_indicator))
  ja <- as.vector(ja_matrix)
  
  # construct the resulting matrix
  res@ra <- ra
  res@ja <- ja
  res@ia <- as.integer(ia)
  # check whether the dense format is needed
  if(dense){
    res <- as.matrix(res)
  }
  return(res)
}

#' Generate the intercept term matrix according to mu_vec.
#' 
#' This function generate the group id matrix according to the subgroup effect vector mu_vec.
#' 
#' \itemize{
#' \item n = length(mu_vec) is the number of observations.
#' \item p = length(unique(mu_vec)) is the number of subgroups.
#' }
#' 
#' @section Note: 
#' In the result, 1st observation of mu_vec will always be in the 1st subgroup.
#' 
#' @param mu_vec a length-n vector of the subgroup effect.
#' @return res: a (n * p) matrix. Each row is for one observation and res[i, j] = 1 if i \eqn{\in} group_j and res[i, j] = 0 o.w.
#' @examples 
#' set.seed(1024)
#' mu_vec <- sample(c(0.5, 0, -0.5), size = 10, replace = TRUE)
#' group_id_mat <- RSAVS:::RSAVS_Mu_to_Mat(mu_vec)
RSAVS_Mu_to_Mat <- function(mu_vec){
    # Generate the intercept term matrix according to mu_vec
    # n = length(mu_vec), p = length(unique(mu_vec))
    # 1st observation will always be in the 1st group
    # Args: mu_vec: a vector of the subgroup effect
    # Return: res: a (n * p) matrix.
    #              each row is for one observation 
    #              and res[i, j] = 1 if i \in group_j and res[i, j] = 0 if o.w.
    
    n <- length(mu_vec)
    mu_unique <- unique(mu_vec)
    p <- length(mu_unique)
    
    res <- matrix(0, nrow = n, ncol = p)
    for(i in 1 : n){
      j <- which(mu_unique == mu_vec[i])
      res[i, j] <- 1
    }
    return(res)
}


#' Convert pair-wise difference vector to subgroup structure
#' 
#' This function converts the augmented pair-wise difference vector(the s vector in the algorithm) 
#' to the subgrouping result.
#' 
#' Detailed definition of the s vector(length \code{n * (n - 1) / 2}) can be found in the paper.
#' 
#' @section Note:
#'   It's possible that there's logical contradiction in the s vector, 
#'   e.g. \code{s_{12} = s_{13} = 0}, but \code{s_{23} != 0}, especially when the algorithm is
#'   configured with a loose tolerance. In this function, \code{i} and \code{j} would be put in the
#'   same subgroup as long as there is a path that connects them. In previous example, 1st, 2nd and
#'   3rd observations will be classified as in the same subgroup even though \code{s_{23} != 0}. 
#'   
#'   This strategy would presumably provide a more concentrate results, hence less number of subgroups. 
#'   But of course it tends to merge subgroups into big ones and leaving some alone observations.
#'   
#'   For large scale data, especially when the number observation is big. It will be difficult to save
#'   all \code{s_vec}s during the algorithm considering the algorithm have to search over the 
#'   \code{lam1_length * lam2_length} grid. During the algorithm, the \code{s_vec} is utilized to
#'   improve the subgrouping results then discarded.
#'   
#' @param s_vec  the s vector(pair-wise difference vector), length \code{n * (n - 1) / 2} and 
#'   \eqn{s_{ij} = \mu_i - \mu_j}.
#' @param n number of observations.
#' @return a list containing the grouping result. 
#' @examples 
#' n <- 10    # number of observations
#' group_center <- c(-1, 0, 1)    # three group centers
#' # subgroup effect vector    
#' alpha_true <- sample(group_center, size = n, replace = TRUE)
#' d_mat <- RSAVS:::RSAVS_Generate_D_Matrix(n)    # pair-wise difference matrix
#' s_vec <- d_mat %*% alpha_true
#' RSAVS:::RSAVS_S_to_Groups(s_vec, n)
#' 
#' # if there's contradiction in s_vec
#' s_vec[1] <- 0    # s_{12} = 0
#' s_vec[2] <- 0    # s_{13} = 0
#' s_vec[n] <- 1    # s_{23} != 0
#' # 1, 2, 3 will still be in the same subgroup
#' RSAVS:::RSAVS_S_to_Groups(s_vec, n)
RSAVS_S_to_Groups <- function(s_vec, n){
  # This function converts the s vector from the algorithm to the subgrouping result
  # Definition of the s vector(length n * (n -1) / 2) is in the paper
  # Note: It's possible that there's logical contradiction in the s vector, eg s_{12} = s_{13} = 0, but s_{23} != 0
  #       In the algorithm, i and j would be put in the same subgroup as long there is a path that can connect them.
  #       ie: s_{12} = s_{13} = 0, but s_{23} != 0 would still results these 3 points in the same subgroup
  #       This strategy would presumably provide a more concentrate results, hence less number of subgroups
  # Args: s_vec: the s vector, length n * (n - 1) / 2 and s_{ij} = \mu_i - \mu_j
  #       n: number of observation
  # Return: res: a list containing the grouping result
  res <- list()
  res[[1]] <- c(1)    # 1st observation is always in the 1st subgroup
  
  # analyse the grouping result between 1st observation and others, in order to construct the basic of all observations
  for(i in 1 : (n - 1)){
    if(s_vec[i] == 0){    # (i+1)-th observation in the same group as the 1st observation
      res[[1]] <- c(res[[1]], i + 1)
    } else{    # (i+1)-th observation is not in the same group as the 1st observation
      res[[length(res) + 1]] <- i + 1
    }
  }
  
  positive_obs_vec <- rep(1 : (n - 1), times = (n - 1) : 1)
  negative_obs_vec <- unlist(lapply(2 : n, function(x) x : n))
  for(i in n : length(s_vec)){
    positive_obs <- positive_obs_vec[i]
    negative_obs <- negative_obs_vec[i]
    if(s_vec[i] == 0){
      # if s_vec[i] == 0, which means positive_obs and negative_obs are in the same subgroup, 
      # then we should put them in the same subgroup
      # otherwise nothing should be done
      
      # find the current subgroup IDs of positive and negative observations
      positive_groupid <- which(sapply(X = res, FUN = function(x, id){
        return(is.element(id, x))
      }, id = positive_obs))
      negative_groupid <- which(sapply(X = res, FUN = function(x, id){
        return(is.element(id, x))
      }, id = negative_obs))
      
      # if they are already in the same subgroup, then nothing to be done. Otherwise merge these 2 groups
      if(positive_groupid != negative_groupid){
        # currently merge them all into the positive group and remove the negative group
        res[[positive_groupid]] <- sort(unique(c(res[[positive_groupid]], res[[negative_groupid]])))
        res[[negative_groupid]] <- NULL
      }
    }
  }
  
  return(res)
}

#' Determine the subgroup effects
#' 
#' This function determines the final mu vector given the grouping results.
#' 
#' @param mu_vec numerical vector. \code{n = length(mu_vec)} is the number of observations.
#'   Presumably it comes from the ADMM algorithm and is not a very good grouping result.
#' @param group_res a list, containing the grouping results.
#'   Each element of group_res is a list containing the indecs from the same subgroup.
#'   See \code{\link{RSAVS_S_to_Groups}} for the output structure.
#' @param round_digits a positive integer digit. If \code{group_res} is missing and \code{round_digits}
#'   is provided, then the function will directly round the \code{mu_vec} to \code{round_digits}
#'   to provide the final grouping results.
#' @param klim a length-3 integer vector, storing `kmin`, `kmax`, `dudahart_kmax` 
#'   for specifying settings for clustering algorithm.
#' @param usepam boolen, whether to use `cluster::pam()` or instead `cluster::clara()` 
#'   to estimate center of clusters.
#' @return a length-n vector of subgroup effects. Observations belonging to the same
#'   subgroup should have the same effects.
#' @details 
#'   Currently, the resulting subgroup effects are the average value of \code{mu_vec}
#'   for those belonging to the same subgroup according to \code{group_res}.
#'   
#'   If \code{group_res} is missing and \code{round_digits} is \code{NULL}, 
#'   then a simple K-means is performed on \code{mu_vec}
#'   to get the estimated grouping result.
#' @examples 
#' # basic settings
#' set.seed(1024)
#' n <- 100    # number of observations
#' group_center <- c(-1, 0, 1)    # three group centers
#' # subgroup effect vector    
#' alpha_true <- sample(group_center, size = n, replace = TRUE)
#' alpha_est <- alpha_true + rnorm(n, sd = 0.25)    # add noise to group effect
#' table(alpha_true)
#' 
#' # Use `group_res` to determine estimated group effect from `alpha_est`
#' d_mat <- RSAVS:::RSAVS_Generate_D_Matrix(n)    # pair-wise difference matrix
#' # Normally, s_vec should come from the algorithm
#' #   Here we just use a simple estimate, which might lead to not so satisfing results
#' s_vec <- round(d_mat %*% alpha_est, 0)    
#' group_res <- RSAVS:::RSAVS_S_to_Groups(s_vec, n)
#' table(RSAVS:::RSAVS_Determine_Mu(alpha_est, group_res = group_res))
#' 
#' # Use `pamk` to determine estiamted group effect
#' table(RSAVS:::RSAVS_Determine_Mu(alpha_est))
#' 
#' # Use directly rounding to determine estimated group effect
#' table(RSAVS:::RSAVS_Determine_Mu(alpha_est, round_digits = 1))
RSAVS_Determine_Mu <- function(mu_vec, group_res, klim = c(2, 7, 4), usepam = (length(mu_vec) < 2000), round_digits = NULL){
  # This function determines the final mu vector given the grouping results
  # Args: mu_vec: The given mu vector, length n, probability comes from the ADMM algorithm and not a very good grouping result
  #       group_res: A list, containing the grouping results. 
  #                  Each element of group_res is a list containing the indecs from the same subgroup
  #                  You can refer to RSAVS_S_to_Groups for the output structure
  # Returns: res: a new mu vector. 
  #               Current strategy is taking average value of mu_vec for those belong to the same subgroup in group_res

  if(missing(group_res)){    # new version, use cluster method
    if(is.null(round_digits)){
      # no `round_digits` provided, use `pamk`
      n <- length(mu_vec)
      if(n < 2){
        stop("There are only 2 observations, not meaning for clustering.")
      }
      if(n < 10){
        warning("number of observations(n) to small, pamk might be unstable.")
        
      }
      kmax <- min(klim[2], 10)
      kmin <- klim[1]
      dudahart_kmax <- klim[3]
      if(kmin > 1){
        pamk_res <- try(fpc::pamk(mu_vec, krange = kmin : kmax, usepam = usepam), silent = T)
        if(!inherits(pamk_res, "try_error")){
          group_num <- pamk_res$nc
        }else{
          print("inner pamk fail!")
          group_num <- 1
        }
      }else{
        # try to determine k = 1 (only one cluster) is suitable for this data
        p_vec <- rep(0, dudahart_kmax - 1)
        for(k in 2 : dudahart_kmax){
          if(usepam){
            tmp2 <- cluster::pam(mu_vec, k = k)
          }else{
            tmp2 <- cluster::clara(mu_vec, k = k)
          }
          
          # print(paste("------ k = ", k, " ------", sep = ""))
          all_p_na <- TRUE
          for(i in 1 : (k - 1)){
            for(j in (i + 1) : k){
              idx <- union(which(tmp2$clustering == i), 
                           which(tmp2$clustering == j))
              current_p <- fpc::dudahart2(mu_vec[idx], tmp2$clustering[idx])$p.value
              all_p_na <- all_p_na & is.na(current_p)
              if(!is.na(current_p)){
                if(current_p > p_vec[k - 1]){
                  p_vec[k - 1] <- current_p
                }
              }
            }
          }
          if(all_p_na){
            p_vec[k - 1] <- 1
          }
        }
        # print(p_vec)
        double_check1 <- any(p_vec < 0.0001)
        if(!double_check1){
          group_num <- 1
        }else{
          pamk_res <- try(fpc::pamk(mu_vec, krange = 2 : kmax, usepam = usepam), silent = T)
          if(!inherits(pamk_res, "try_error")){
            group_num <- pamk_res$nc
          }else{
            print("inner pamk fail!")
            group_num <- 1
          }
        }
      }
      
      if(group_num == 1){
        res <- rep(mean(mu_vec), n)
      }else{
        res <- pamk_res$pamobject$medoids[pamk_res$pamobject$clustering]
      }
      return(res)
    }else{
      # simply round `mu_vec` based on `round_digits`
      res <- round(mu_vec, digits = round_digits)
      return(res)
    }

  }else{    # old version, use my result from RSAVS_S_to_Groups
    # The result is not very good because it tends to merge all observation into one big group
    n <- length(mu_vec)
    res <- mu_vec
    group_num <- length(group_res)
    for(i in 1 : group_num){
      group_ids <- group_res[[i]]
      res[group_ids] <- mean(mu_vec[group_ids])
    }
    return(res)
  }
}

#' Compute Modified BIC
#' 
#' This function computes the modified BIC(mBIC), given a specific solution.
#' 
#' @param y_vec numerical vector of response. \code{n = length(y_vec)} is the number of observations.
#' @param x_mat numerical matrix of covariates. \code{p = ncol(x_mat)} is the number of covariates.
#' @param beta_vec numerical vector of covariate effects
#' @param mu_vec numerical vector of subgroup effect(intercept term) for each observations
#' @param loss_type,loss_param type and parameters of the loss function
#'   \itemize{
#'     \item \code{loss_type = "L1"}: L-1 loss, no actual parameter is needed.
#'     \item \code{loss_type = "L2"}: L-2 loss, no actual parameter is needed.
#'     \item \code{loss_type = "Huber"}, Huber loss, and \code{c = loss_param[0]}, where \code{c}
#'       is the parameter needed by huber loss. A popular choice is \code{c = 1.345}.
#'   }
#' @param phi a positive constant. 
#' @param a a scalar controlling the regression part of mBIC. If not provided, then the number of observation(\code{n}) is used.
#' @param double_log_lik boolen, whether to use log logLik(double log likelihood) in the regression part of mBIC. Defaults to \code{TRUE}.
#' @return the mBIC value.
#' @details 
#' \deqn{BIC = log(1 / a * sum(loss(y - mu - x * beta)) + |S| * Phi ,}
#' where
#'   \itemize{
#'     \item mu is the intercept term of each observation. And the number of subgroups is
#'       \code{K = length(unique(mu_vec))}.
#'     \item beta is the covariate effect vector. And the number of active covariates is
#'       \code{Q = sum(beta_vec != 0)}.
#'     \item the loss function is determined by \code{loss_type} and \code{loss_param}.
#'     \item |S| is the complexity of the model and \eqn{|S| = K + Q}.
#'     \item Phi is a constant and \code{Phi = phi * log(log(n + p)) * log(n) / n}.
#'   }
#' Note that in most cases, \code{sum(loss_fun)} corresponds to log likelihood. 
RSAVS_Compute_BIC <- function(y_vec, x_mat, beta_vec, mu_vec, loss_type, loss_param, phi, a, double_log_lik = TRUE){
  # This function computes the BIC, given a specific solution.
  # BIC = log(1 / n * sum(loss(y - mu - x * beta)) + |S| * Phi
  # where 1. mu is the intercept term of each observation. 
  #          And the number of subgroups = length(unique(mu_vec))
  #       2. beta is the covariate effect vector
  #          And the number of active covariates = sum(beta_vec != 0)
  #       3. loss function is determined by loss_type and loss_param
  #          And |S| is the number of subgroups and active covariates
  #       4. Phi is a constant, Phi = phi * loglog(n + p) * log(n) / n
  # Args: y_vec, x_mat: the data
  #       beta_vec: vector of estimated covariate effects 
  #       mu_vec: vector of intercept effect of each observation
  #       loss_type, loss_param: type and parameters of the loss function
  #                              loss_type = "L1", L1 loss, no actual parameter is needed
  #                              loss_type = "L2", L2 loss, no actual parameter is needed
  #                              loss_type = "Huber", Huber loss, and c = loss_param[0]
  #       phi: a constant
  # Returns: bic.
  
  # prepare the data
  y_vec <- matrix(y_vec, ncol = 1)    # make sure y_vec is column vector
  n <- length(y_vec)
  x_mat <- matrix(x_mat, nrow = n)    # make sure x_mat is a matrix, even if it has only one column
  p <- ncol(x_mat)
  
  if(missing(a)){
    a <- n
  }
  
  # compute PHI
  # phi <- phi * log(log(n + p)) * log(n + p) / n
  phi <- phi * log(log(n + p)) * log(n) / n    # another version
  
  # Find loss function
  loss_fun <- RSAVS_L1
  if(loss_type == "Huber"){
    loss_fun <- RSAVS_Huber
  } else{
    if(loss_type == "L2"){
      loss_fun <- RSAVS_L2
    }
  }
  
  # Find the grouping results
  # Use my function to determine subgrouping results
  # group_res <- RSAVS_S_to_Groups(s_vec, n)
  # group_num <- length(group_res)
  # mu_vec <- RSAVS_Determine_Mu(mu_vec, group_res)
  
  # Use a cluster method to determine subgrouping results
  group_num <- length(unique(mu_vec))
  
  # active number of covariates
  active_beta_num <- sum(beta_vec != 0)
  
  # a thresholding for complexsity
  # if(active_beta_num + group_num > 15){
  #   return(10^4)
  # }
  
  # compute bic
  # bic_p1 <- log(1 / a * sum(loss_fun(y_vec - mu_vec - x_mat %*% beta_vec, loss_param)))
  bic_p1 <- 1 / a * sum(loss_fun(y_vec - mu_vec - x_mat %*% beta_vec, loss_param))
  if(double_log_lik){
    bic_p1 <- log(bic_p1)
  }
  bic_p2 <- (group_num + active_beta_num) * phi
  bic <- bic_p1 + bic_p2
  
  return(bic)
}

#' Summary the iterations during ADMM algorithm
#' 
#' This function is designed to summary and improve the resutls during the iteration of ADMM algorithm.
#' 
#' This function has two purposes:
#'   \itemize{
#'     \item Determine and improve beta_vec and mu_vec, if possible.
#'     \item Compute BIC.
#'   }
#' Since for large scale data set, especially with big number of observations, it's impossible to first 
#' save all the variables during the ADMM algorithm over the \code{lam1_length * lam2_length} grid points
#' of lambdas, then pick a best solution with mBIC. For \code{s_vec} alone, this means we have to save a
#' matrix with \code{n * (n - 1) / 2} rows and \code{lam1_length * lam2_length} columns, which is hard
#' for a single computer. Instead, we summarise each iteration during the algorithm. Then there's no need
#' for storing so many data.
#' @param y_vec numerical response vector. \code{n = length(y_vec)} is the number of observations.
#' @param x_mat numerical covariate matrix. \code{p = ncol(x_mat)} is the number of covariates.
#' @param beta_vec covariate effect vector during the ADMM algorithm.
#' @param mu_vec subgroup effect vector during the ADMM algorithm
#' @param s_vec augmented vector for pair-wise difference of \code{mu_vec} in ADMM algorithm.
#' @param w_vec augmented vector for \code{beta_vec} in ADMM algorithm.
#' @param loss_type character string indicating type of loss function.
#' @param loss_param numerical vector for necessary parameters in loss function.
#' @param phi a parameter needed in mBIC. It controls how strong mBIC penalizes the complexity of 
#'   the candidate model.
#' @return a list, containing:
#'   \itemize{
#'     \item \code{bic}: the bic value.
#'     \item \code{mu_vec}: the improved mu vector.
#'     \item \code{group_num}: number of subgroups in the improved \code{mu_vec}.
#'     \item \code{active_num}: number of active covariates in the \code{beta_vec}.
#'   }
#' @details In the ADMM algorithm, it will take many iterations to reach a sharp tolerance.
#'   But one can stop the algorithm early stage by setting a small \code{max_iter}. This is
#'   equivalent to setting a loose tolerance. Then the \code{mu_vec} and \code{beta_vec} will
#'   not be close to their augmented counterparts \code{s_vec} and \code{w_vec}. But these
#'   counterparts actually provides the sparsity information during the algorithm, therefore
#'   \itemize{
#'     \item \code{w_vec} will provide the estimate of covariate effect while \code{beta_vec}
#'       is just a intermediate variable.
#'     \item \code{mu_vec} is also the intermediate variable. Improvement is needed like forming 
#'       a reasonalbe subgroup structure. One possible solution is to utilize \code{s_vec} to 
#'       improve \code{mu_vec}. Another is to apply some cluster methods on \code{mu_vec}. See
#'       \code{\link{RSAVS_Determine_Mu}} for more details.
#'   }
#' @export
RSAVS_Summary_Iteration <- function(y_vec, x_mat, beta_vec, mu_vec, s_vec, w_vec, loss_type, loss_param, phi){
  # This function is designed to summary and improve the resutls during the iteration of ADMM algorithm
  # It does: 1. Determing and improving beta_vec and mu_vec, if possible
  #          2. Compute BIC.
  # Args:
  # Returns: res, a list containing:
  #          1. bic: the bic value
  #          2. mu_vec: the improved mu vector
  #          3. group_num: the number of subgroups in the improved mu vector
  #          4. active_num: the number of active covariates in the w vector
  # Note: In the ADMM algorithm:
  #       1. w_vec will provide the estimate of covariate effect while beta_vec is just a intermediate variable
  #       2. mu_vec is also the intermediate variable, improvement is needed
  #          One possible solution is to utilize s_Vec to improve mu_vec.
  #          Another is to apply some cluster methods on mu_vec
  #
  #
  
  # 1. Improve the mu vector
  mu_vec <- RSAVS_Determine_Mu(mu_vec)
  group_num <- length(unique(mu_vec))
  active_num <- sum(w_vec != 0)
  
  # 2. Compute BIC
  # We need to simplify RSAVS_Compute_BIC function.
  bic <- RSAVS_Compute_BIC(y_vec = y_vec, x_mat = x_mat, 
                           beta_vec = w_vec, mu_vec = mu_vec, 
                           loss_type = loss_type, loss_param = loss_param, 
                           phi = phi)
  
  # 3. Construct the res list
  res <- list(bic = bic, 
              mu_vec = mu_vec,
              group_num = group_num, 
              active_num = active_num)
  return(res)
}

#' Post-selection estimation
#' 
#' This function carries out the post-selection estimation given the model structure.
#' 
#' @param y_vec numeric vector for response. 
#'   \code{n = length(y_vec)} is the number of observations.
#' @param x_mat numeric matrix for covariates. Each row is for one observation.
#'   \code{p = ncol(x_mat)} is the number of covariates.
#' @param l_type character, type of loss function.
#'   \itemize{
#'     \item "L1": L-1 loss.
#'     \item "L2": L-2 loss.
#'     \item "Huber": Huber loss.
#'   }
#'   The default value is "1".
#' @param l_param vector of parameters needed by the loss function. 
#'   For Huber loss, \code{c = l_param[1]} and a popular choice is c = 1.345.
#'   The Default value is \code{NULL} since no additional parameters are
#'   needed for the default L-1 loss.
#' @param mu_vec a length-n vector for subgroup effect.
#'   The function uses this to determine the subgroup structure.
#' @param beta_vec a length-p vector for covariate effect.
#'   The function uses this to find the active covariates.
#' @return a list containing the improved mu(mu_vec) and beta(beta_vec) vector.
#' @details This function uses \code{mu_vec} to determine the subgroup structure matrix
#'   via \code{\link{RSAVS_Mu_to_Mat}}. The active covariates are those with non-zero
#'   \code{beta_vec} entries.
#'   
#'   Since this is the post-selection estimation. It's a regular estimate without
#'   any penalties. One should make sure the given model structure is identifiable.
#' @seealso \code{\link{RSAVS_Mu_to_Mat}} for getting the subgroup index matrix
#'   from the subgroup effect vector.
#' @export
RSAVS_Further_Improve <- function(y_vec, x_mat, l_type = "L1", l_param = NULL, mu_vec, beta_vec){
    # This function is designed for further improving the estimating of mu and beta
    # after the mu vector has been improved by the clustering method.
    # The idea is simple, the clustering method will provide a good estimate of subgroups
    # but generally it will increase the loss function.
    # Therefore we fit the loss function again, given the grouping and variable selection results.
    # Args: y_vec: response vector, length(y_vec) = n
    #       x_mat: covariate matrix, nrow(x_mat) = n, ncol(x_mat) = p
    #       l_type, l_param: type of loss function and the necessary parameters.
    #                        l_type = "1"(default) for L1 loss; "2" for L2 loss; "H" for Huber loss.
    #                        Currently only support L1 and L2 loss, and l_param is not used
    #                        Future plan is to support Huber loss, and c = l_param[1](c is the k in rlm's psi.huber)
    #       mu_vec: the intercept vector. This function uses this to find the grouping result.
    #       beta_vec: the covariate effect vector. This function uses this to find the active covariates.
    # Return: res: a list containing the improved mu(mu_vec) and beta(beta_vec) vector
    y_vec <- as.matrix(y_vec, ncol = 1)
    n <- length(y_vec)
    x_mat <- as.matrix(x_mat, nrow = n)
    p <- ncol(x_mat)
    
    # Generate the intercept term matrix
    mu_mat <- RSAVS_Mu_to_Mat(mu_vec)
    group_num <- length(unique(mu_vec))
    
    # Generate the active covariate matrix
    active_id <- which(beta_vec != 0)
    active_x <- x_mat[, active_id, drop = F]
    
    # put together the covariate matrix
    newx_mat <- cbind(mu_mat, active_x)
    
    # Check if the new design matrix is over-parameterized
    if(ncol(newx_mat) >= n){
      # For over-parametrized model, 
      # there is no meaning in further improve the estimation
      res <- list(mu_vec = mu_vec, beta_vec = beta_vec)
    } else{
      # re-fit the model
      if(l_type == "L1"){
          tmp <- rq(y_vec ~ newx_mat - 1, tau = 0.5)
          new_mu <- tmp$coefficients[1 : group_num]
          new_beta <- tmp$coefficients[(1 : length(active_id)) + group_num]
      } else{
        if(l_type == "L2"){
            tmp <- lm(y_vec ~ newx_mat - 1)
            new_mu <- tmp$coefficients[1 : group_num]
            new_beta <- tmp$coefficients[(1 : length(active_id)) + group_num]
        } else{
          if(l_type == "Huber"){
            tmp <- rlm(y_vec ~ newx_mat - 1, k = l_param[1])
            new_mu <- tmp$coefficients[1 : group_num]
            new_beta <- tmp$coefficients[(1 : length(active_id)) + group_num]
          }
        }
      }  
      mu_vec <- mu_mat %*% new_mu
      beta_vec[active_id] <- new_beta
      res <- list(mu_vec = mu_vec, beta_vec = beta_vec)      
    }
 
    return(res)   
}

#' Robust subgroup analysis and variable selection simultaneously for large scale dataset
#' 
#' This function carries out robust subgroup analysis and variable selection simultaneously. 
#' It's implemented in a parallel fashion. It supports different types of loss functions and penalties.
#' 
#' @param y_vec numerical vector of response. n = length(y_vec) is the number of observations.
#' @param x_mat numerical matrix of covariates. Each row for one observation and 
#'   \code{p = ncol(x_mat)} is the number of covariates.
#' @param l_type character string, type of loss function.
#'   \itemize{
#'     \item "L1": l-1 loss(absolute value loss)
#'     \item "L2": l-2 loss(squared error loss)
#'     \item "Huber": Huber loss. Its parameter is given in l_param.
#'   }
#'   Default value is "L1".
#' @param l_param numeric vector containing necessary parameters of the corresponding loss function. 
#'   The default value is NULL.
#' @param p1_type,p2_type a character indicating the penalty types for subgroup identification and variable selection.
#'   \itemize{
#'     \item "S": SCAD
#'     \item "M": MCP
#'     \item "L": Lasso
#'   }
#'   Default values for both parameters are "S".
#' @param p1_param,p2_param numerical vectors providing necessary parameters for the corresponding penalties.
#'   \itemize{
#'     \item For Lasso, lam = p_param[1]
#'     \item For SCAD and MCP, lam = p_param[1], gamma = p_param[2]
#'   }
#'   Default values for both parameters are \code{c(2, 3.7)}. 
#'   Note: This function searches the whole lam1_vec * lam2_vec grid for the best solution. 
#'   Hence the \code{lambda}s provided in these parameters will be ignored and overwritten.
#' @param lam1_vec,lam2_vec numerical vectors of customized lambda vectors. 
#'   For \code{lam1_vec}, it's preferred to be in the order from small to big.
#' @param min_lam_ratio the ratio between the minimal and maximal lambda, equals to (minimal lambda) / (maximal lambda).
#'   The default value is 0.03.
#' @param lam1_len,lam2_len integers, lengths of the auto-generated lambda vectors.
#' @param initial_vec list of vector, providing initial values for the algorithm. 
#'   \code{mu_initial = initial_vec$mu}
#'   and 
#'   \code{beta_initial = initial_vec$beta}.
#'   Note: Current the function uses auto-generated lambda vectors 
#'   and the corresponding initial values are all just shrinked to 0.
#'   Hence this parameter is ignored.
#' @param phi numerical variable. A parameter needed for mBIC.
#' @param tol numerical, convergence tolerance for the algorithm.
#' @param max_iter integer, max number of iteration during the algorithm.
#' @param subgroup_benchmark bool. Whether this call should be taken as a benchmark of subgroup identification. 
#'   If \code{TRUE}, then the penalty for variable selection will be surpressed to a minimal value.
#' @examples
#' # a toy example
#' # first we generate data
#' n <- 200    # number of observations
#' q <- 5    # number of active covariates
#' p <- 50    # number of total covariates
#' k <- 2    # number of subgroups
#' 
#' # k subgroup effect, centered at 0
#' group_center <- seq(from = 0, to = 2 * (k - 1), by = 2) - (k - 1)
#' # covariate effect vector
#' beta_true <- c(rep(1, q), rep(0, p - q))
#' # subgroup effect vector    
#' alpha_true <- sample(group_center, size = n, replace = TRUE)    
#' x_mat <- matrix(rnorm(n * p), nrow = n, ncol = p)    # covariate matrix
#' err_vec <- rnorm(n, sd = 0.1)    # error term
#' y_vec <- alpha_true + x_mat %*% beta_true + err_vec    # response vector
#' 
#' # a simple analysis using default loss and penalties
#' res <- RSAVS_LargeN(y_vec = y_vec, x_mat = x_mat, 
#'                     lam1_len = 50, lam2_len = 40, 
#'                     phi = 5)
#' 
#' # you can choose different loss and penalties
#' res_huber <- RSAVS_LargeN(y_vec = y_vec, x_mat = x_mat, 
#'                           l_type = "Huber", l_param = 1.345,
#'                           p1_type = "M", p2_type = "L",
#'                           lam1_len = 40, lam2_len = 30, 
#'                           phi = 5)
#' # you can do post-selection estimation by
#' ind <- res$best_ind    # pick an id of the solution
#' res2 <- RSAVS_Further_Improve(y_vec = y_vec, x_mat = x_mat, 
#'                               mu_vec = res$mu_improve_mat[ind, ], 
#'                               beta_vec = res$w_mat[ind, ])                           
#' @export
RSAVS_LargeN <- function(y_vec, x_mat, l_type = "L1", l_param = NULL, 
                         p1_type = "S", p1_param = c(2, 3.7), p2_type = "S", p2_param = c(2, 3.7), 
                         lam1_vec, lam2_vec, 
                         min_lam1_ratio = 0.03, min_lam2_ratio = 0.03, 
                         lam1_len, lam2_len, 
                         const_r123, const_abc = rep(1, 3), 
                         initial_values, phi = 1.0, tol = 0.001, max_iter = 10, 
                         cd_max_iter = 1, cd_tol = 0.001, 
                         subgroup_benchmark = FALSE){
  ### preparation ###
  # preparation for x and y #
  y_vec <- matrix(y_vec, ncol = 1)    # make sure y_vec is column vector
  n <- length(y_vec)
  x_mat <- matrix(x_mat, nrow = n)    # make sure x_mat is a matrix, even if it has only one column
  p <- ncol(x_mat)
  
  # preparation for loss function #
  if(l_type == "L2"){
    # l_type = "2"
    l_param = 0
  } else{
    if(l_type == "L1"){
      # l_type = "1"
      l_param = 0
    } else {
      if(l_type == "Huber"){
        # l_type = "H"
        huber_c = l_param[1]
      }else{
        stop("Error of Loss type!")
      }
    }
  }
  
  
  
  # preparation for penalties #
  if(p1_type == "L"){
    p1_lam = p1_param[1]
  } else{
    if((p1_type == "S") || (p1_type == "M")){
      p1_lam = p1_param[1]
      p1_gamma = p1_param[2]
    } else {
      stop("Error of Penalty 1 type!")
    }
  }
  if(p2_type == "L"){
    p2_lam = p2_param[1]
  } else {
    if((p2_type == "S") || (p2_type == "M")){
      p2_lam = p2_param[1]
      p2_gamma = p2_param[2]
    } else {
      stop("Error of Penalty 2 type!")
    }
  }
  
  
  # preparation for r1, r2 and r3
  if(missing(const_r123)){
    message("`const_r123` is missing. Use default settings!")
    r1 <- 2
    
    if(p1_type == "L"){
      r2 <- 2
    }
    if(p1_type == "S"){
      r2 <- max(1.3 * 1 / (p1_gamma - 1), 1 / (p1_gamma - 1) + 1)
    }
    if(p1_type == "M"){
      r2 <- max(1.3 * 1 / p1_gamma, 1 / p1_gamma + 1)
    }
    
    if(p2_type == "L"){
      r3 <- 2
    }
    if(p2_type == "S"){
      r3 <- max(1.3 * 1 / (p2_gamma - 1), 1 / (p2_gamma - 1) + 1)
    }
    if(p2_type == "M"){
      r3 <- max(1.3 * 1 / p2_gamma, 1 / p2_gamma + 1)
    }
  }else{
    const_r123 <- abs(const_r123)
    r1 <- const_r123[1]
    r2 <- const_r123[2]
    r3 <- const_r123[3]
    
    # check for r2
    if(p1_type == "S"){
      r2min <- max(1.3 * 1 / (p1_param[2] - 1), 1 / (p1_param[2] - 1) + 1)
    }
    if(p1_type == "M"){ 
      r2min <- max(1.3 * 1 / p1_param[2], 1 / p1_param[2] + 1)
    }
    
    if(r2 < r2min){
      message("r2 = ", r2, "< r2min = ", r2min, ". Modify value of r2!")
      r2 <- r2min
      const_r123[2] <- r2
    }
    # check for r3
    if(p2_type == "S"){
      r3min <- max(1.3 * 1 / (p2_param[2] - 1), 1 / (p2_param[2] - 1) + 1)
    }
    if(p2_type == "M"){
      r3min <- max(1.3 * 1 / p2_param[2], 1 / p2_param[2] + 1)
    }
    if(r3 < r3min){
      message("r3 = ", r3, "< r3min = ", r3min, ". Modify value of r3!")
      r3 <- r3min
      const_r123[3] <- r3
    }
  }
  
  # check const_abc
  const_abc <- abs(const_abc)
  # const_a <- const_abc[1]
  # const_b <- const_abc[2]
  # const_c <- const_abc[3]
  
  # preparation for lam_vec
  if(missing(lam1_vec)){
    # Old version, currently only designed for L-1 loss
    # lam1_max <- 1 / (n - 1) 
    # lam1_max <- 1 / (2 - 1)    #use 1 / 2 \sum\loss, not 1 / n \sum\loss
    
    # new version
    # d_mat <- RSAVS_Generate_D_Matrix(n = n)
    # if(l_type == "1"){    # L1 loss
    #   y_median <- median(y_vec)
    #   d_vec <- (-1) * sign(y_vec - y_median)
    # }
    # if(l_type == "2"){    # L2 loss
    #   y_mean <- mean(y_vec)
    #   d_vec <- (-2) * (y_vec - y_mean)
    # }
    # if(l_type == "H"){    # Huber loss
    #   tmp <- rlm(y_vec ~ 1, k = l_param[1])
    #   y_huber <- tmp$coefficients[1]
    #   d_vec <- (-1) * RSAVS_Huber(y_vec - y_huber, param = l_param[1], derivative = T)
    #   rm(tmp)
    # }
    # d_vec <- d_vec / n    # use 1 / n sum loss
    # # d_vec <- d_vec / 2    # use 1 / 2 sum loss
    # 
    # if(n <= 500){
    #   d_inv <- ginv(as.matrix(t(d_mat) %*% d_mat))
    #   lam1_max <- max(abs(as.matrix(d_mat %*% (d_inv %*% d_vec))))
    # }else{
    #   lam1_max <- 2 / n * max(abs(d_vec))
    # }
    # 
    # lam1_min <- lam1_max * min_lam_ratio
    # lam1_vec <- exp(seq(from = log(lam1_max), to  = log(lam1_min), length.out = lam1_length))
    
    # newer version
    message("lam1_vec is missing, use default values...")
    lam1_max <- RSAVS_Get_Lam_Max(y_vec = y_vec, l_type = l_type, l_param = l_param, 
                                  lam_ind = 1, 
                                  const_abc = const_abc, eps = 10^(-6))
    lam1_min <- lam1_max * min_lam1_ratio
    if(p1_type != "L"){
      lam1_max <- lam1_max * 100    # safe guard for non-convex penalties
    }
    
    # lam1_vec <- exp(seq(from = log(lam1_max), to = log(lam1_min), length.out = lam1_len))
    lam1_vec <- exp(seq(from = log(lam1_min), to = log(lam1_max), length.out = lam1_len))
  }else{
    # lam1_vec <- sort(abs(lam1_vec), decreasing = T)
    # lam1_length = length(lam1_vec)
    lam1_len <- length(lam1_vec)
    lam1_vec <- sort(abs(lam1_vec), decreasing = FALSE)
  }
  
  if(missing(lam2_vec)){
    # if(l_type == "1"){    # L-1 loss
    #   d_vec <- (-1) * sign(y_vec - median(y_vec))
    # }
    # if(l_type == "2"){    # L-2 loss
    #   d_vec <- (-2) * (y_vec - mean(y_vec))
    # }
    # if(l_type == "H"){    # Huber loss
    #   tmp <- rlm(y_vec ~ 1, k = l_param[1])
    #   y_huber <- tmp$coefficients[1]
    #   d_vec <- (-1) * RSAVS_Huber(y_vec - y_huber, param = l_param[1], derivative = T)
    #   rm(tmp)
    # }    
    # d_vec <- d_vec / n    # use 1 / n \sum\loss
    # # d_vec <- d_vec / 2    # use 1 / 2 \sum\loss
    # lam2_max <- max(abs(t(x_mat) %*% d_vec))
    # lam2_min <- lam2_max * min_lam_ratio
    # lam2_vec <- exp(seq(from = log(lam2_max), to = log(lam2_min), length.out = lam2_length)) 
    
    # newer version
    message("lam2_vec is missing, use default values...")
    lam2_max <- RSAVS_Get_Lam_Max(y_vec = y_vec, x_mat = x_mat, 
                                  l_type = l_type, l_param = l_param, 
                                  lam_ind = 2, 
                                  const_abc = const_abc, eps = 10^(-6))
    lam2_min <- lam2_max * min_lam2_ratio
    lam2_vec <- exp(seq(from = log(lam2_max), to = log(lam2_min), length.out = lam2_len))
  }else{
    lam2_vec <- sort(abs(lam2_vec), decreasing = TRUE)
    lam2_len <- length(lam2_vec)
  }
  
  if(subgroup_benchmark){
    # supress lambda for variable selection if this is subgroup benchmark
    lam2_vec <- lam2_vec * 0.05
  }
  
  # --- prepare other values ---
  message("prepare intermediate variables needed by the algorithm")
  # intermediate variables needed for the algorithm
  d_mat <- RSAVS_Generate_D_Matrix(n)    # pairwise difference matrix
  
  beta_lhs <- NA    # left part for updating beta
  if(n >= p){
    beta_lhs <- solve(r1 * t(x_mat) %*% x_mat + r3 * diag(nrow = p))
  }else{
    beta_lhs <- 1.0 / r3 * (diag(nrow = p) - r1 * t(x_mat) %*% solve(r1 * x_mat %*% t(x_mat) + r3 * diag(nrow = n)) %*% x_mat)
  }
  
  # left part for updating mu
  mu_lhs <- solve(SparseM::as.matrix(r1 * diag(nrow = n) + r2 * SparseM::t(d_mat) %*% d_mat))
  
  mu_beta_lhs <- NA    # left part for updating beta and mu together
  if(cd_max_iter == 0){
    mu_beta_lhs <- matrix(0, nrow = n + p, ncol = n + p)
    mu_beta_lhs[1 : n, 1 : n] <- SparseM::as.matrix(r1 * diag(nrow = n) + r2 * SparseM::t(d_mat) %*% d_mat) 
    mu_beta_lhs[n + (1 : p), 1 : n] <- r1 * t(x_mat)
    mu_beta_lhs[1 : n, n + (1 : p)] <- r1 * x_mat
    mu_beta_lhs[n + (1 : p), n + (1 : p)] <- r1 * t(x_mat) %*% x_mat + r3 * diag(nrow = p)
    mu_beta_lhs <- solve(mu_beta_lhs)
  }
  
  additional <- list(d_mat = d_mat,  
                     beta_lhs = beta_lhs, 
                     mu_lhs = mu_lhs, 
                     mu_beta_lhs = mu_beta_lhs)
  
  if(l_type == "L2"){
    # beta_left_inv <- solve(t(x_mat) %*% x_mat / n + r3 / 2 * diag(nrow = p))
    # d_mat <- RSAVS_Generate_D_Matrix(n)
    # mu_left_inv <- solve(diag(nrow = n) / n + r2 / 2 * as.matrix(t(d_mat) %*% d_mat))
    res <- RSAVS_LargeN_L2_Rcpp(x_mat, y_vec, n, p, p1_type, p1_param, p2_type, p2_param, lam1_vec, lam2_vec, r2, r3, phi, tol, max_iter)
  } else{
    res <- RSAVS_LargeN_Rcpp(x_mat, y_vec, n, p, l_type, l_param, p1_type, p1_param, p2_type, p2_param, lam1_vec, lam2_vec, r1, r2, r3, phi, tol, max_iter)
  }
  
  # Idealy, the result of c(lam1[1], lam2[1]) should be mu being median and beta being 0
  # So the result is directly set to this, without actually computing using ADMM
  # But in actuality, our derivation of lam1[1] and lam2[1] is not that accurate, 
  # Hence this direct setting is not that accurate
  # Given this situation, when best_ind being 1
  # We re-choose the best result
  if(res$best_ind == 1){
    min_bic_id <- which.min(res$bic_mat[-1])
    res$best_ind <- min_bic_id + 1
    res$best_i <- (min_bic_id + 1 - 1) %% length(lam1_vec) + 1
    res$best_j <- (min_bic_id + 1 - 1) %/% length(lam1_vec) + 1
    res$repick <- TRUE
  }
  return(res)
}

#' Rand Index
#' 
#' This function computes the Rand Index(RI) of an estimated mu vector(mu_est) 
#' against the target mu vector(mu_target).
#' 
#' @param mu_est numerical vector of the estimated subgroup effect.
#'   \code{n = length(mu_est)}.
#' @param mu_target numerical vector of the target(real) subgroup effect.
#'   \code{n = length(mu_target)}.
#' @param detail whether or not should the function computes the details(TP, TN, FP, FN) 
#'   of these 2 identifications(\code{mu_est} and \code{mu_target}).
#'   Default value is \code{FALSE}.
#' @return numerical value of Rand Index if \code{detail = FALSE}. 
#'   Otherwise a list containing \code{RI, TP, TN, FP, FN}.
#' @section Note:
#'   \code{mu_est} and \code{mu_target} must have the same length.
#' @section Details:
#'   For a subgroup identification of \code{n} observations, 
#'   there are \eqn{n * (n - 1) / 2} pairs of different observations.
#'   The Rand Index describes how these two subgroup identifications are
#'   based on these observation pairs and it is given by
#'   \deqn{RI = (#TP + #TN) / (#TP + #TN + #FP + #FN),}
#'   where for any two observation \code{i} and \code{j},
#'   \itemize{
#'     \item TP(True Positive): 
#'       these two observations belong to the same subgroup in both these two identification.
#'     \item TN(True Negative): 
#'       these two observations belong to different subgroups in both these two identification.
#'     \item FP(False Positive): 
#'       these two observations belong to the same subgropu in \code{mu_est}.
#'       But they belong to different subgroups in \code{mu_target}.
#'     \item FN(False Negative): 
#'       these two observations belong to different subgroups in \code{mu_est}.
#'       But they belong to the same subgroup in \code{mu_target}.   
#'   }
#'   Therefore RI ranges from 0 to 1 and RI = 1 if and only if these two identifications are identical.
#'   
#'   Note: if \code{detail = TRUE}, the function also returns the number of TP, TN, FP and FN.
#' @examples 
#' n <- 50    # number of observations
#' k <- 2    # k subgroup effect, centered at 0
#' group_center <- seq(from = 0, to = 2 * (k - 1), by = 2) - (k - 1)
#' # subgroup effect vector    
#' alpha1 <- sample(group_center, size = n, replace = TRUE) 
#' # another subgroup effect vector
#' alpha2 <- sample(group_center, size = n, replace = TRUE)
#' RSAVS_RI(alpha1, alpha2)
#' RSAVS_RI(alpha1, alpha2, detail = TRUE)
#' @export
RSAVS_RI <- function(mu_est, mu_target, detail = FALSE){
  # Compute the Rand Index(RI) of an estimated mu vector(mu_est) against the target mu vector(mu_target)
  # RI = (TP + TN) / (TP + TN + FP + FN)
  # mu_est and mu_target must have the same length(n), then TP + TN + FP + FN = n * (n - 1) / 2
  # Args: mu_est: a length n vector of the estimated group effect
  #       mu_target: a length n vector of the target group effect
  #       detail: indicator, whether or not should the function compute the detail of these 2 groups(TP, TN, FP, FN)
  # Return: res: the Rand Index
  if(length(mu_est) != length(mu_target)){
    stop("mu_est and mu_target must have same length!")
  }else{
    n <- length(mu_est)
  }
  
  # convert these 2 group effect vector to grouping results
  d_mat <- RSAVS_Generate_D_Matrix(n)
  mu_est <- d_mat %*% mu_est
  mu_target <- d_mat %*% mu_target
  # simplify the results:
  # 0: i, j in the same group. 1: i, j in the different group.
  mu_est <- as.numeric(mu_est != 0)
  mu_target <- as.numeric(mu_target != 0)
  
  # Find TP, TN and compute rand index
  TP <- sum((mu_est == 0) & (mu_target == 0))
  TN <- sum((mu_est != 0) & (mu_target != 0))
  # FP <- sum((mu_est == 0) & (mu_target != 0))
  # FN <- sum((mu_est != 0) & (mu_target == 0))
  res <- (TP + TN) / (n * (n - 1) / 2)
  
  # another way to do this
  # FP_FN_sum <- sum(xor(mu_est, mu_target))
  # res <- 1 - FP_FN_sum / (n * (n - 1) / 2)
  if(detail == TRUE){
    RI <- res
    FP <- sum((mu_est == 0) & (mu_target != 0))
    FN <- sum((mu_est != 0) & (mu_target == 0))
    res <- list(RI = RI, 
                TP = TP, 
                TN = TN, 
                FP = FP, 
                FN = FN)
  }
  
  return(res)
}
