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
    
    if(l_type == "L1"){
        mu0 <- median(y_vec)
        psi_vec <- sign(y_vec - mu0)
    }else{
        if(l_type == "L2"){
            mu0 <- mean(y_vec)
            psi_vec <- 2 * (y_vec - mu0)
        }else{
            if(l_type == "Huber"){
                tmp <- MASS::rlm(y_vec ~ 1, k = l_param[1])
                mu0 <- tmp$coefficients[1]
                psi_vec <- RSAVS_Huber(y_vec - mu0, param = l_param[1], derivative = T)
            }else{
                stop(paste("Unsupported type of loss function!, l_type = ", l_type, sep = ""))
            }
        }
    }
    return(psi_vec)
}

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
