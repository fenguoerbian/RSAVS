# This file contains the pre-define functions for Robust Subgroup Analysis


#' Built-in loss functions
#'
#' These are built-in loss functions.
#'
#' @aliases RSAVS_L2 RSAVS_L1 RSAVS_Huber
#' @usage RSAVS_L1(x, param)
#' RSAVS_L2(x, param)
#' RSAVS_Huber(x, param, derivative)
#' @param x input numeric vector
#' @param param parameters needed for the function, takes the form of 
#' numeric vector. Unused for L1 and L2.
#' @param derivative logical, whether the return is the loss value or the
#' derivative.
#' @return Loss(or derivative) value at x.
#' @examples
#' RSAVS_L2(1 : 10)
#' RSAVS_L1(1 : 10)
#' RSAVS_Huber(seq(from = -3, to = 3, by = 0.1), param = 1.345)
#' RSAVS_Huber(seq(from = -3, to = 3, by = 0.1), param = 1.345, derivative = TRUE)
#' @name loss function
NULL

#' @rdname loss function 
#' @export
RSAVS_L2 <- function(x, param){
  # The L2 loss function
  return(x ^ 2)
}

#' @rdname loss function 
#' @export
RSAVS_L1 <- function(x, param){
  # The L1 loss function
  return(abs(x))
}

#' @rdname loss function 
#' @export
RSAVS_Huber <- function(x, param, derivative = FALSE){
  # The huber loss function
  # Huber(x, c) = 0.5 * x ^ 2                 if abs(x) <= c
  #               c * abs(x) - 0.5 * c ^ 2    if abs(x) > c
  huber_c <- param[1]
  if(!derivative){
    res <- x
    id <- which(abs(x) <= huber_c)
    res[id] <- 0.5 * (x[id] ^ 2)
    res[-id] <- huber_c * abs(x[-id]) - 0.5 * huber_c ^ 2
  }else{
    res <- x
    id <- which(abs(x) <= huber_c)
    res[-id] <- huber_c * sign(x[-id])
  }
  return(res)
}

#' Generate the pair-wise different matrix
#'
#' This function generate the pairwise difference matrix(D matrix in the paper)
#' By default it returns a sparse matrix(matrix.csr) from the package SparseM
#'
#' @param n number of observations.
#' @param dense logical, whether the return type should be in dense matrix or not
#' @return a difference matrix with size \eqn{(n * (n - 1) / 2) \times n}
RSAVS_Generate_D_Matrix <- function(n, dense = FALSE){
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
  ra <- rep(c(1, -1), times = n * (n - 1) / 2)
  # ia: starting point of each row in ra and ja
  ia <- (1 : (n * (n - 1) / 2 + 1)) * 2 - 1
  # ja: column index, according to ra
  positive_column_indicator <- rep(1 : (n - 1), times = (n - 1) : 1)
  negative_column_indicator <- unlist(lapply(2 : n, function(x) x : n))
  ja_matrix <- t(cbind(positive_column_indicator, negative_column_indicator))
  ja <- as.vector(ja_matrix)
  res@ra <- ra
  res@ja <- ja
  res@ia <- as.integer(ia)
  # check whether the dense format is needed
  if(dense){
    res <- as.matrix(res)
  }
  return(res)
}


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


RSAVS_Determine_Mu <- function(mu_vec, group_res){
  # This function determines the final mu vector given the grouping results
  # Args: mu_vec: The given mu vector, length n, probability comes from the ADMM algorithm and not a very good grouping result
  #       group_res: A list, containing the grouping results. 
  #                  Each element of group_res is a list containing the indecs from the same subgroup
  #                  You can refer to RSAVS_S_to_Groups for the output structure
  # Returns: res: a new mu vector. 
  #               Current strategy is taking average value of mu_vec for those belong to the same subgroup in group_res

  if(missing(group_res)){    # new version, use cluster method
    n <- length(mu_vec)
    pamk_res <- try(pamk(mu_vec, krange = 1 : 7, usepam = (n <= 2000)), silent = T)
    if(!inherits(pamk_res, "try-error")){
      group_num <- pamk_res$nc
    }else{
      print("pamk fail!")
      group_num <- 1
    }
    
    if(group_num == 1){
      res <- rep(mean(mu_vec), n)
    }else{
      res <- pamk_res$pamobject$medoids[pamk_res$pamobject$clustering]
    }
    return(res)
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


RSAVS_Compute_BIC <- function(y_vec, x_mat, beta_vec, mu_vec, loss_type, loss_param, phi){
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
  #                              loss_type = "1", L1 loss, no actual parameter is needed
  #                              loss_type = "2", L2 loss, no actual parameter is needed
  #                              loss_type = "H", Huber loss, and c = loss_param[0]
  #       phi: a constant
  # Returns: bic.
  
  # prepare the data
  y_vec <- matrix(y_vec, ncol = 1)    # make sure y_vec is column vector
  n <- length(y_vec)
  x_mat <- matrix(x_mat, nrow = n)    # make sure x_mat is a matrix, even if it has only one column
  p <- ncol(x_mat)
  
  # compute PHI
  # phi <- phi * log(log(n + p)) * log(n + p) / n
  phi <- phi * log(log(n + p)) * log(n) / n    # another version
  
  # Find loss function
  loss_fun <- RSAVS_L1
  if(loss_type == "H"){
    loss_fun <- RSAVS_Huber
  } else{
    if(loss_type == "2"){
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
  bic_p1 <- log(1 / n * sum(loss_fun(y_vec - mu_vec - x_mat %*% beta_vec, loss_param)))
  bic_p2 <- (group_num + active_beta_num) * phi
  bic <- bic_p1 + bic_p2
  
  return(bic)
}

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

#' @export
RSAVS_Further_Improve <- function(y_vec, x_mat, l_type = "1", l_param = NULL, mu_vec, beta_vec){
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
      if(l_type == "1"){
          tmp <- rq(y_vec ~ newx_mat - 1, tau = 0.5)
          new_mu <- tmp$coefficients[1 : group_num]
          new_beta <- tmp$coefficients[(1 : length(active_id)) + group_num]
      } else{
        if(l_type == "2"){
            tmp <- lm(y_vec ~ newx_mat - 1)
            new_mu <- tmp$coefficients[1 : group_num]
            new_beta <- tmp$coefficients[(1 : length(active_id)) + group_num]
        } else{
          if(l_type == "H"){
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

#' @export
RSAVS_LargeN <- function(y_vec, x_mat, l_type = "L1", l_param = NULL, 
                         p1_type = "S", p1_param = c(2, 3.7), p2_type = "S", p2_param = c(2, 3.7), 
                         lam1_vec, lam2_vec, min_lam_ratio = 0.03, lam1_length, lam2_length, 
                         initial_vec, phi = 1.0, tol = 0.001, max_iter = 10, subgroup_benchmark = FALSE){
  # ADMM algorithm, for large n
  # Args: y_vec: response vector, length(y_vec) = n
  #       x_mat: covariate matrix, nrow(x_mat) = n, ncol(x_mat) = p
  #       l_type: character string, type of the loss function in the objective function
  #               "L2": l-2 loss, "L1": l-1 loss, "Huber": huber loss
  #       l_param: necessary parameters for the corresponding loss function.
  #                For l-2 and l-1 loss, this parameter is not necessary
  #                For Huber loss, c = l_param[1]
  #       p1_type, p2_type: character variable, type of the penalties applied over mu and beta
  #                         "L: Lasso, "S": SCAD, "M": MCP
  #       p1_param, p2_param: numeric vector, necessary parameters for the corresponding vector.
  #                           For lasso: lam = p_param[1]
  #                           For SCAD and MCP: lam = p_param[1], gamma = p_param[2]
  #       lam1_vec, lam2_vec:
  #       min_lam_ratio: the ration between the minimal and maximal lambda, equals to (minimal lambda) / (maximal lambda)
  #       lam1_length, lam2_length: the length of lam1_vec and lam2_vec
  #       initial_vec: list of vector, providing initial values for the algorithm
  #                    mu_initial = initial_vec$mu
  #                    beta_initial = initial_vec$beta
  #       phi: constant needed for computing BIC
  #       tol: tolerance
  #       max_iter: maximum number of iteration
  #       subgroup_benchmark: bool, default to FALSE. Whether this computation is only for subgroup benchmark
  #                           if true, then the lambda vector of penalty for covariate selection will be shrink to a small value.
  # Returns: A list containning:
  
  ### preparation ###
  # preparation for x and y #
  y_vec <- matrix(y_vec, ncol = 1)    # make sure y_vec is column vector
  n <- length(y_vec)
  x_mat <- matrix(x_mat, nrow = n)    # make sure x_mat is a matrix, even if it has only one column
  p <- ncol(x_mat)
  
  # preparation for loss function #
  if(l_type == "L2"){
    l_type = "2"
    l_param = 0
  } else{
    if(l_type == "L1"){
      l_type = "1"
      l_param = 0
    } else {
      if(l_type == "Huber"){
        l_type = "H"
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
  
  # preparation for lam_vec
  if(missing(lam1_vec)){
    # Old version, currently only designed for L-1 loss
    # lam1_max <- 1 / (n - 1) 
    # lam1_max <- 1 / (2 - 1)    #use 1 / 2 \sum\loss, not 1 / n \sum\loss
    
    # new version
    d_mat <- RSAVS_Generate_D_Matrix(n = n)
    if(l_type == "1"){    # L1 loss
      y_median <- median(y_vec)
      d_vec <- (-1) * sign(y_vec - y_median)
    }
    if(l_type == "2"){    # L2 loss
      y_mean <- mean(y_vec)
      d_vec <- (-2) * (y_vec - y_mean)
    }
    if(l_type == "H"){    # Huber loss
      tmp <- rlm(y_vec ~ 1, k = l_param[1])
      y_huber <- tmp$coefficients[1]
      d_vec <- (-1) * RSAVS_Huber(y_vec - y_huber, param = l_param[1], derivative = T)
      rm(tmp)
    }
    d_vec <- d_vec / n    # use 1 / n sum loss
    # d_vec <- d_vec / 2    # use 1 / 2 sum loss
    
    if(n <= 500){
      d_inv <- ginv(as.matrix(t(d_mat) %*% d_mat))
      lam1_max <- max(abs(as.matrix(d_mat %*% d_inv %*% d_vec)))
    }else{
      lam1_max <- 2 / n * max(abs(d_vec))
    }
    
    lam1_min <- lam1_max * min_lam_ratio
    lam1_vec <- exp(seq(from = log(lam1_max), to  = log(lam1_min), length.out = lam1_length))
    
  }else{
    
  }
  
  if(missing(lam2_vec)){
    if(l_type == "1"){    # L-1 loss
      d_vec <- (-1) * sign(y_vec - median(y_vec))
    }
    if(l_type == "2"){    # L-2 loss
      d_vec <- (-2) * (y_vec - mean(y_vec))
    }
    if(l_type == "H"){    # Huber loss
      tmp <- rlm(y_vec ~ 1, k = l_param[1])
      y_huber <- tmp$coefficients[1]
      d_vec <- (-1) * RSAVS_Huber(y_vec - y_huber, param = l_param[1], derivative = T)
      rm(tmp)
    }    
    d_vec <- d_vec / n    # use 1 / n \sum\loss
    # d_vec <- d_vec / 2    # use 1 / 2 \sum\loss
    lam2_max <- max(abs(t(x_mat) %*% d_vec))
    lam2_min <- lam2_max * min_lam_ratio
    lam2_vec <- exp(seq(from = log(lam2_max), to = log(lam2_min), length.out = lam2_length)) 
  }else{
    
  }
  if(subgroup_benchmark){
    # supress lambda for variable selection if this is subgroup benchmark
    lam2_vec <- lam2_vec * 0.05
  }
  
  # prepare other values
  if(l_type == "2"){
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
