RSAVS_Get_Mu0 <- function(y_vec, l_type = "L1", l_param = NULL){
  if(l_type == "L1"){
    mu0 <- median(y_vec)
  }else{
    if(l_type == "L2"){
      mu0 <- mean(y_vec)
    }else{
      if(l_type == "Huber"){
        tmp <- MASS::rlm(y_vec ~ 1, k = l_param[1])
        mu0 <- tmp$coefficients[1]
      }else{
        stop(paste("Unsupported type of loss function!, l_type = ", l_type, sep = ""))
      }
    }
  }
  return(mu0)
}

#' Compute psi_vec 
#' 
#' This function computes the psi_vec needed by the algorithm.
#' 
#' In the algorithm, when all \eqn{\mu} are shrinked to the same value, we have
#' \deqn{
#'   \mu_0 = \mathrm{argmin}_{\mu} \sum_i\rho(y_i - \mu)
#' }
#' and 
#' \deqn{
#'   psi\_vec = (\psi(y_1 - \mu_0), \cdots, \psi(y_n - \mu_0))
#' }
#' where \eqn{\psi = \partial \rho}.
#' 
#' @note there is a SIGN different between this \code{psi_vec} and 
#'   \eqn{\partial \rho} w.r.t \eqn{\mu}
#' @note for non-differential loss like L1 at 0, the derivative at 0+ is returned.
#' @param y_vec numeric vector
#' @param l_type character string, type of loss function.
#'   \itemize{
#'     \item "L1": l-1 loss(absolute value loss)
#'     \item "L2": l-2 loss(squared error loss)
#'     \item "Huber": Huber loss. Its parameter is given in l_param.
#'   }
#'   Default value is "L1".
#' @param l_param numeric vector containing necessary parameters of the corresponding loss function. 
#'   The default value is \code{NULL}.
#' @return the \code{psi_vec}
RSAVS_Get_Psi <- function(y_vec, l_type = "L1", l_param = NULL){
    # This function gets the default psi vec for the algorithm
    # Be default, we mean: 
    #     mu0 = argmin_{mu} \sum \rho(y_i - mu)
    # and 
    #     psi_vec = c(psi(y_1 - mu0), psi(y_2 - mu0), \cdots, psi(y_n - mu0))
    # where 
    #     psi = \partial \rho
    # NOTE: there is a SIGN different between this psi_vec and \partial \rho w.r.t \mu
    # NOTE: for non-differential loss like L1 at 0, see source code for its value
    # Args: y_vec: numerical vector of y
    #       l_type: type of loss function, character string
    #               "L1": median regression
    #               "L2": linear regression
    #               "Huber": Huber loss regression
    #       l_param: parameter vector for the loss function. 
    #                Not needed for "L1" and "L2" type of loss
    
    mu0 <- RSAVS_Get_Mu0(y_vec, l_type, l_param)
    if(l_type == "L1"){
        psi_vec <- sign(y_vec - mu0)
    }else{
        if(l_type == "L2"){
            psi_vec <- 2 * (y_vec - mu0)
        }else{
            if(l_type == "Huber"){
                psi_vec <- RSAVS_Huber(y_vec - mu0, param = l_param[1], derivative = T)
            }else{
                stop(paste("Unsupported type of loss function!, l_type = ", l_type, sep = ""))
            }
        }
    }
    return(psi_vec)
}

#' Compute lam_max
#' 
#' This function computes the values of \code{lam1_max} or \code{lam2_max} which 
#'   can shrink same type of covariates to the same value. 
#'   
#' For \code{beta}, they will be shrinked to 0. For \code{mu}, the final value 
#' is determined by the type of loss function. Please refer to \code{\link{RSAVS_Get_Psi}}
#' for more details. Basically, we have
#' \deqn{
#'   \mu_0 = \mathrm{argmin}_{\mu} \sum_i\rho(y_i - \mu)
#' }
#' and 
#' \deqn{
#'   psi\_vec = (\psi(y_1 - \mu_0), \cdots, \psi(y_n - \mu_0))
#' }
#' where \eqn{\psi = \partial \rho}.
#' 
#' @param y_vec numeric vector, the response
#' @param x_mat numeric matrix, the covariate matrix
#' @param l_type character string, type of loss function.
#'   \itemize{
#'     \item "L1": l-1 loss(absolute value loss)
#'     \item "L2": l-2 loss(squared error loss)
#'     \item "Huber": Huber loss. Its parameter is given in l_param.
#'   }
#'   Default value is "L1".
#' @param l_param numeric vector containing necessary parameters of the corresponding loss function. 
#'   The default value is \code{NULL}.
#' @param lam_ind integer, indicating computation for \code{mu} or \code{beta}.
#'   \itemize{
#'     \item 1: compute \code{lam_max} for \code{mu}.
#'     \item 2: compute \code{lam_max} for \code{beta}.
#'   }
#' @param const_abc a length-3 numeric vector, providing the scalars to adjust weight
#'   of regression function, penalty for subgroup identification and penalty for 
#'   variable selection in the overall objective function. Defaults to \code{c(1, 1, 1)}.
#' @param eps a samll safe guard constant
#' @return \code{lam_max}, a numerical variable
#' @seealso \code{\link{RSAVS_Get_Psi}}.
#' @note For convex penalty like lasso, this \code{lam_max} can guarantee the global
#'   condition for shrinking the corresponding variables to the same value. For non-convex
#'   penalties such as SCAD and MCP, this is just a local condition.
RSAVS_Get_Lam_Max <- function(y_vec, x_mat, l_type = "L1", l_param = NULL, lam_ind = 1, const_abc = rep(1, 3), eps = 10^(-6)){
    # lam_ind: 1: lam1
    #          2: lam2
    n <- length(y_vec)
    const_a <- const_abc[1]
    const_b <- const_abc[2]
    const_c <- const_abc[3]
    
    psi_vec <- RSAVS_Get_Psi(y_vec = y_vec, l_type = l_type, l_param = l_param)
    
    if(lam_ind == 1){
        lam_max <- 2 / n / const_a / const_b * max(abs(psi_vec)) * (1 + eps)
        lam_max <- lam_max * 1.1    # safe guard
    }else{
        if(lam_ind == 2){
            # `psi_vec` is in matrix form
            #   convert it to vector form
            # psi_vec <- as.vector(psi_vec)
            # tmp <- psi_vec * x_mat
            # tmp <- colSums(tmp)
            tmp <- t(x_mat) %*% psi_vec
            lam_max <- 1 / const_a / const_c * max(abs(tmp)) * (1 + eps)
        }else{
            stop(paste("Wrong index for lambda with lam_ind = ", lam_ind, sep = ""))
        }
    }
    
    return(lam_max)
}

#' Compute value of objective function
#' 
#' This function computes objective function's value for the ADMM algorithm.
#' 
#' The augmented lagrangian objective function for the ADMM algorithm contains
#'   \itemize{
#'     \item regression part, \code{1 / const_a * sum(rho(z_vec))}, where \code{rho} is
#'       the loss function. Refer to \code{\link{loss_function}} for more details.
#'     \item subgroup analysis part, \code{const_b * sum(P_1(s_vec))}.
#'     \item variable selection part, \code{const_c * sum(P_2(w_vec))}.
#'     \item augmented part1:
#'       \code{r_1 / 2 * norm(y_vec - mu_vec - x_mat * beta - z_vec) ^ 2 + inner_product(y_vec - mu_vec - x_mat * beta - z_vec, q1_vec)}
#'     \item augmented part1:
#'       \code{r_2 / 2 * norm(D_mat * mu_vec - s_vec) ^ 2 + inner_product(D_mat * mu_vec - s_vec, q2_vec)}
#'     \item augmented part1:
#'       \code{r_3 / 2 * norm(beta_vec - w_vec) ^ 2 + inner_product(beta_vec - w_vec, q3_vec)}              
#'   }
#' @note for "L2" loss, the \code{z_vec} could be eliminated, but currently this
#'   is not implemented. 
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
#'   The default value is \code{NULL}.
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
#' @param const_r123 a length-3 numerical vector, providing the scalars needed in the 
#'   augmented lagrangian part of the ADMM algorithm
#' @param const_abc a length-3 numeric vector, providing the scalars to adjust weight
#'   of regression function, penalty for subgroup identification and penalty for 
#'   variable selection in the overall objective function. Defaults to \code{c(1, 1, 1)}.
#' @param beta_vec,mu_vec,z_vec,s_vec,w_vec,q1_vec,q2_vec,q3_vec variables needed
#'   in the objective function
#' @return a list containing
#' \itemize{
#'     \item \code{loss}: overall loss value
#'     \item \code{loss_part1}: loss value from regression part
#'     \item \code{loss_part2}: loss value from subgroup analysis part
#'     \item \code{loss_part3}: loss value from variable selection part
#'     \item \code{loss_aug1}: loss value from augmented part1
#'     \item \code{loss_aug2}: loss value from augmented part2
#'     \item \code{loss_aug3}: loss value from augmented part3   
#'     \item \code{diff_z}: difference between \code{z_vec} and \code{y_vec - mu_vec - x_mat \%*\% beta_vec} 
#'     \item \code{diff_s}: difference between \code{s_vec} and \code{d_mat \%*\% mu_vec}
#'     \item \code{diff_w}: difference between \code{w_vec} and \code{beta_vec}
#'   }
#' @seealso \code{\link{loss_function}}
RSAVS_Compute_Loss_Value <- function(y_vec, x_mat, l_type = "L1", l_param = NULL, 
                                     p1_type = "S", p1_param = c(2, 3.7), p2_type = "S", p2_param = c(2, 3.7), 
                                     const_r123, const_abc, 
                                     beta_vec, mu_vec, z_vec, s_vec, w_vec, 
                                     q1_vec, q2_vec, q3_vec){
  # This functions computes the current value of loss function(augmented lagrangian form)
  # There are three parts in the main body of loss function:
  #   1 / a * sum(rho(z_vec))
  #   b * sum(P_1(s_vec))
  #   c * sum(P_2(w_vec))
  # There are three parts in the augmented part of loss function:
  #   r_1 / 2 * norm(y_vec - mu_vec - x_mat * beta - z_vec) ^ 2 + inner_product(y_vec - mu_vec - x_mat * beta - z_vec, q1_vec)
  #   r_2 / 2 * norm(D_mat * mu_vec - s_vec) ^ 2 + inner_product(D_mat * mu_vec - s_vec, q2_vec)
  #   r_3 / 2 * norm(beta_vec - w_vec) ^ 2 + inner_product(beta_vec - w_vec, q3_vec)
  # NOTE: for "L2" loss, there should not be a part about the `z_vec`, 
  #         but currently this is not implemented correctly.
  # ------ prepare some basic variables ------
  y_vec <- as.vector(y_vec)
  n <- length(y_vec)
  p <- nrow(x_mat)
  d_mat <- RSAVS_Generate_D_Matrix(n)    # pairwise difference matrix
  
  # --- check for loss function ---
  if(l_type == "L1"){
    loss_fun <- RSAVS_L1
  }else{
    if(l_type == "L2"){
      loss_fun <- RSAVS_L2
    }else{
      if(l_type == "Huber"){
        loss_fun <- RSAVS_Huber
      }else{
        stop("Unsupported type of loss function!")
      }
    }
  }
  
  # --- check for penalty 1 ---
  if(p1_type == "L"){
    p1_fun <- penalty_lasso
  }else{
    if(p1_type == "S"){
      p1_fun <- penalty_scad
    }else{
      if(p1_type == "M"){
        p1_fun <- penalty_mcp
      }else{
        stop("Wrong type of `p1_type`! Must be either `L`, `S` or `M`!")
      }
    }
  }
  
  # --- check for penalty 2 --- 
  if(p2_type == "L"){
    p2_fun <- penalty_lasso
  }else{
    if(p2_type == "S"){
      p2_fun <- penalty_scad
    }else{
      if(p2_type == "M"){
        p2_fun <- penalty_mcp
      }else{
        stop("Wrong type of `p2_type`! Must be either `L`, `S` or `M`!")
      }
    }
  }
  
  
  # ------ compute value of loss function -------
  # --- main loss function ---
  loss_part1 <- 1 / const_abc[1] * sum(loss_fun(z_vec, param = l_param))
  loss_part2 <- const_abc[2] * sum(p1_fun(s_vec, params = p1_param))
  loss_part3 <- const_abc[3] * sum(p2_fun(w_vec, params = p2_param))
  
  # --- augmented lagrangian part ---
  tmp <- y_vec - mu_vec - x_mat %*% beta_vec - z_vec
  loss_aug1 <- as.numeric(const_r123[1] / 2 * t(tmp) %*% tmp + t(tmp) %*% q1_vec)
  diff_z <- as.numeric(t(tmp) %*% tmp)
  
  tmp <- SparseM::as.matrix(d_mat %*% mu_vec) - s_vec
  loss_aug2 <- as.numeric(const_r123[2] / 2 * t(tmp) %*% tmp + t(tmp) %*% q2_vec)
  diff_s <- as.numeric(t(tmp) %*% tmp)
  
  tmp <- beta_vec - w_vec
  loss_aug3 <- as.numeric(const_r123[3] / 2 * t(tmp) %*% tmp + t(tmp) %*% q3_vec)
  diff_w <- as.numeric(t(tmp) %*% tmp)
  
  loss <- loss_part1 + loss_part2 + loss_part3 + loss_aug1 + loss_aug2 + loss_aug3
  
  res <- list(loss = loss, 
              loss_part1 = loss_part1, 
              loss_part2 = loss_part2, 
              loss_part3 = loss_part3, 
              loss_aug1 = loss_aug1, 
              loss_aug2 = loss_aug2, 
              loss_aug3 = loss_aug3, 
              diff_z = diff_z, 
              diff_s = diff_s, 
              diff_w = diff_w)
  return(res)
}


RSAVS_Compute_BIC_V2 <- function(id, rsavs_res, y_vec, x_mat, l_type, l_param, 
                                 update_mu = list(useS = FALSE, round_digits = NULL), 
                                 phi, const_a, doublle_log_lik = TRUE, 
                                 from_rsi = FALSE, pb = NULL){
  # report progress info if a progressor is provided
  if(!is.null(pb)){
    pb(message = sprintf("ID = %d", id))
  }
  
  # prepare mu_vec and beta_vec
  if(from_rsi){
    beta_vec <- rsavs_res$beta_mat[id, ]
    mu_vec <- rsavs_res$mu_updated_mat[id, ]
    if(!is.null(update_mu)){
      update_mu$useS <- FALSE
    }
  }else{
    beta_vec <- rsavs_res$w_mat[id, ]
    mu_vec <- rsavs_res$mu_updated_mat[id, ]
  }
  
  n <- length(mu_vec)
  # if mu_vec is already updated, which means it's in a meaningful subgroup strucutre
  #   then there's no need to perform RSAVS_Determine_Mu
  if(is.null(update_mu)){
    mu_improved <- mu_vec
  }else{
    useS <- update_mu$useS
    if(is.null(update_mu$klim)){
      klim <- c(2, 7, 4)
    }else{
      klim <- update_mu$klim
    }
    
    if(is.null(update_mu$usepam)){
      usepam <- length(res$mu_vec < 2000)
    }else{
      usepam <- update_mu$usepam
    }
    
    round_digits <- update_mu$round_digits
    
    if(useS){
      group_res <- RSAVS_S_to_Groups(round(rsavs_res$s_mat[id, ], 3), n)
      mu_improved <- RSAVS_Determine_Mu(mu_vec, group_res)
    }else{
      mu_improved <- RSAVS_Determine_Mu(mu_vec, klim = klim, usepam = usepam, round_digits = round_digits)
    }
  }
  
  # post-selection estimation
  post_est <- try(RSAVS_Further_Improve(y_vec, x_mat, l_type, l_param, mu_improved, beta_vec))
  if(!inherits(post_est, "try-error")){
    beta_vec <- post_est$beta_vec
    mu_improved <- post_est$mu_vec
  }else{
    message("further improve error at id = ", id)
  }
  
  
  bic <- RSAVS_Compute_BIC(y_vec, x_mat, beta_vec, mu_improved, l_type, l_param, phi, const_a, double_log_lik)
  active_num <- sum(beta_vec != 0)
  group_num <- length(unique(mu_improved))
  max_group_size <- max(table(mu_improved))
  
  res <- c(bic, group_num, active_num, max_group_size)
  names(res) <- c("bic", "group_num", "active_num", "max_group_size")
  
  return(list(
    id = id, 
    bic_info = res, 
    mu_vec = mu_improved, 
    beta_vec = beta_vec))
}
