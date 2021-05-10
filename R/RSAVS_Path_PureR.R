#' @export
RSAVS_Solver_PureR <- function(y_vec, x_mat, l_type = "L1", l_param = NULL, 
                         p1_type = "S", p1_param = c(2, 3.7), p2_type = "S", p2_param = c(2, 3.7), 
                         const_r123, const_abc = rep(1, 3), 
                         initial_values, additional, tol = 0.001, max_iter = 10, cd_max_iter = 1, cd_tol = 0.001, 
                         phi = 1.0, subgroup_benchmark = FALSE){
  # Core solver for small n situation
  
  # check the dataset
  y_vec <- matrix(y_vec, ncol = 1)    # make sure y_vec is column vector
  n <- length(y_vec)
  x_mat <- matrix(x_mat, nrow = n)    # make sure x_mat is a matrix, even if it has only one column
  p <- ncol(x_mat)
  
  # Under construction, no variable checking is implemented!
  r1 <- const_r123[1]
  r2 <- const_r123[2]
  r3 <- const_r123[3]
  
  const_a <- const_abc[1]
  const_b <- const_abc[2]
  const_c <- const_abc[3]
  
  # prepare update functions for z, s and w
  if(l_type == "L1"){
    UpdateZ <- RSAVS_UpdateZ_L1
  }else{
    if(l_type == "L2"){
      UpdateZ <- RSAVS_UpdateZ_L2
    }else{
      if(l_type == "Huber"){
        UpdateZ <- RSAVS_UpdateZ_Huber
      }else{
        stop("Unsupported type of loss function!")
      }
    }
  }
  
  if(p1_param[1] == 0){ 
    # No penalty for mu part
    UpdateS <- RSAVS_UpdateSW_Identity
    initial_values$q2_init <- rep(0, n * (n - 1) / 2)    # fix q2 at 0
    # r2 <- 0    # fix r2 to 0
  }else{
    if(p1_type == "L"){
      UpdateS <- RSAVS_UpdateSW_Lasso
    }else{
      if(p1_type == "S"){
        UpdateS <- RSAVS_UpdateSW_SCAD
      }else{
        if(p1_type == "M"){
          UpdateS <- RSAVS_UpdateSW_MCP
        }else{
          stop("Unsupported type of penalty for mu!")
        }
      }
    }
  }
  
  
  if(p2_param[1] == 0){
    # No penalty for beta part
    UpdateW <- RSAVS_UpdateSW_Identity
    initial_values$q3_init <- rep(0, p)    # fix q3 at 0
    # r3 <- 0    # fix r3 to 0
  }else{
    if(p2_type == "L"){
      UpdateW <- RSAVS_UpdateSW_Lasso
    }else{
      if(p2_type == "S"){
        UpdateW <- RSAVS_UpdateSW_SCAD
      }else{
        if(p2_type == "M"){
          UpdateW <- RSAVS_UpdateSW_MCP
        }else{
          stop("Unsupported type of penalty for beta!")
        }
      }
    }
  }
  
  
  #------ main algorithm ------
  #--- initial values ---
  beta_old = initial_values$beta_init
  mu_old = initial_values$mu_init
  z_old = initial_values$z_init
  s_old = initial_values$s_init
  w_old = initial_values$w_init
  q1_old = initial_values$q1_init
  q2_old = initial_values$q2_init
  q3_old = initial_values$q3_init
  
  
  loss_vec <- RSAVS_Compute_Loss_Value(y_vec = y_vec, x_mat = x_mat, l_type = l_type, l_param = l_param, 
                                       p1_type = p1_type, p1_param = p1_param, p2_type = p2_type, p2_param = p2_param, 
                                       const_r123 = const_r123, const_abc = const_abc, 
                                       beta_vec = beta_old, mu_vec = mu_old, 
                                       z_vec = z_old, s_vec = s_old, w_vec = w_old, 
                                       q1_vec = q1_old, q2_vec = q2_old, q3_vec = q3_old)
  loss_old <- loss_vec$loss
  loss_record = rep(loss_old, max_iter + 1)    # record the loss value at the start of each iteration, including the final value of last iteration
  #--- algorithm status ---
  current_step <- 1
  diff <- tol + 1
  
  
  while((current_step <= max_iter) && (diff > tol)){
    # message("Start the ADMM iteration with loss_value = ", loss_old)
    # update beta and mu
    if(cd_max_iter == 0){
      # update beta and mu together
      rhs <- matrix(0, nrow = n + p, ncol = 1)
      rhs[1 : n, 1] <- r1 * (y_vec - z_old) + SparseM::as.matrix(SparseM::t(additional$d_mat) %*% (r2 * s_old - q2_old)) * (p1_param[1] != 0) + q1_old
      rhs[n + (1 : p), 1] <- r1 * t(x_mat) %*% (y_vec - z_old) + r3 * w_old * (p2_param[1] != 0) + t(x_mat) %*% q1_old  - q3_old
      
      mu_beta <- additional$mu_beta_lhs %*% rhs
      
      mu_vec <- mu_beta[1 : n, 1, drop = FALSE]
      beta_vec <- mu_beta[n + (1 : p), 1, drop = FALSE]
      
    }else{
      cd_step <- 1
      cd_diff <- cd_tol + 1
      while((cd_step <= cd_max_iter) && (cd_diff > cd_tol)){
        # update beta
        beta_rhs <- r1 * t(x_mat) %*% (y_vec - z_old - mu_old) + r3 * w_old + t(x_mat) %*% q1_old - q3_old
        beta_vec <- additional$beta_lhs %*% beta_rhs
        
        # update mu
        mu_rhs <- r1 * (y_vec - x_mat %*% beta_vec - z_old) +  SparseM::as.matrix(SparseM::t(additional$d_mat) %*% (r2 * s_old - q2_old)) + q1_old
        mu_vec <- additional$mu_lhs %*% mu_rhs
        
        # check and update CD status
        cd_step <- cd_step + 1
        cd_diff <- max(sqrt(sum((beta_vec - beta_old) ^ 2)), 
                       sqrt(sum((mu_vec - mu_old) ^ 2)))
        
        beta_old <- beta_vec
        mu_old <- mu_old
      }
    }
    
    # update z
    invec <- y_vec - mu_vec - x_mat %*% beta_vec + q1_old / r1
    z_vec <- UpdateZ(invec, param = l_param, r1, const_a)
    
    # update s
    invec <- SparseM::as.matrix(additional$d_mat %*% mu_vec) + q2_old / r2
    s_vec <- UpdateS(invec, p1_param, r2, const_b)
    
    # update w
    invec <- beta_vec + q3_old / r3
    w_vec <- UpdateW(invec, p2_param, r3, const_c)
    
    # update q1, q2 and q3
    q1_vec <- q1_old + r1 * (y_vec - mu_vec - x_mat %*% beta_vec - z_vec)
    q2_vec <- q2_old + r2 * SparseM::as.matrix(additional$d_mat %*% mu_vec - s_vec)
    q3_vec <- q3_old + r3 * (beta_vec - w_vec)
    
    # possible post-selection estimation here?
    
    # check algorithm status
    diff <- max(sqrt(sum((y_vec - mu_vec - x_mat %*% beta_vec - z_vec) ^ 2)), 
                sqrt(sum(SparseM::as.matrix(additional$d_mat %*% mu_vec - s_vec) ^ 2)), 
                sqrt(sum((beta_vec - w_vec) ^ 2)))
    
    beta_old <- beta_vec
    mu_old <- mu_vec
    z_old <- z_vec
    s_old <- s_vec
    w_old <- w_vec
    q1_old <- q1_vec
    q2_old <- q2_vec
    q3_old <- q3_vec
    
    
    current_step <- current_step + 1
    
    loss_vec <- RSAVS_Compute_Loss_Value(y_vec = y_vec, x_mat = x_mat, l_type = l_type, l_param = l_param, 
                                         p1_type = p1_type, p1_param = p1_param, p2_type = p2_type, p2_param = p2_param, 
                                         const_r123 = const_r123, const_abc = const_abc, 
                                         beta_vec = beta_old, mu_vec = mu_old, 
                                         z_vec = z_old, s_vec = s_old, w_vec = w_old, 
                                         q1_vec = q1_old, q2_vec = q2_old, q3_vec = q3_old)
    loss_old <- loss_vec$loss
    loss_record[current_step] <- loss_old
  }
  # print(length(loss_record))
  res <- list(beta_vec = beta_vec, 
              mu_vec = mu_vec, 
              z_vec = z_vec, 
              s_vec = s_vec, 
              w_vec = w_vec, 
              q1_vec = q1_vec, 
              q2_vec = q2_vec, 
              q3_vec = q3_vec, 
              current_step = current_step - 1,
              diff = diff, 
              loss_vec = loss_record)
  return(res)
}

#' @export
RSAVS_Path_PureR <- function(y_vec, x_mat, l_type = "L1", l_param = NULL, 
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
    const_r123 <- c(r1, r2, r3)
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
  
  # preparation for lam_vec
  if(missing(lam1_vec)){
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
  
  # --- prepare initial values for the algorithm ---
  if(missing(initial_values)){
    message("`initial_values` is missing, use default settings!")
    # if(l_type == "L1"){
    #     mu0 <- median(y_vec)
    # }else{
    #     if(l_type == "L2"){
    #         mu0 <- mean(y_vec)
    #     }else{
    #         if(l_type == "Huber"){
    #             tmp <- MASS::rlm(y_vec ~ 1, k = l_param[1])
    #             mu0 <- tmp$coefficients[1]
    #         }
    #     }
    # }
    
    mu_init <- as.vector(y_vec)    # lam1 start from small values, hence we use y_vec as initial values for mu_vec
    
    initial_values <- list(beta_init = rep(0, p), 
                           mu_init = mu_init, 
                           z_init = y_vec - mu_init, 
                           s_init = rep(0, n * (n - 1) / 2), 
                           w_init = rep(0, p), 
                           q1_init = rep(0, n), 
                           q2_init = rep(0, n * (n - 1) / 2), 
                           q3_init = rep(0, p))
  }
  # --- prepare other values ---
  message("prepare intermediate variables needed by the algorithm")
  # intermediate variables needed for the algorithm
  message("generate pairwise difference matrix")
  d_mat <- RSAVS_Generate_D_Matrix(n)    # pairwise difference matrix
  
  message("compute `beta_lhs`")
  beta_lhs <- NA    # left part for updating beta
  if(n >= p){
    beta_lhs <- solve(r1 * t(x_mat) %*% x_mat + r3 * diag(nrow = p))
  }else{
    beta_lhs <- 1.0 / r3 * (diag(nrow = p) - r1 * t(x_mat) %*% solve(r1 * x_mat %*% t(x_mat) + r3 * diag(nrow = n)) %*% x_mat)
  }
  
  # left part for updating mu
  message("compute `mu_lhs`")
  mu_lhs <- solve(SparseM::as.matrix(r1 * diag(nrow = n) + r2 * SparseM::t(d_mat) %*% d_mat))
  
  
  mu_beta_lhs <- NA    # left part for updating beta and mu together
  if(cd_max_iter == 0){
    message("compute `mu_beta_lhs'")
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
  message("additional variables prepared!")
  
  # --- variables storing the results ---
  beta_mat <- matrix(0, nrow = lam1_len * lam2_len, ncol = p)
  mu_mat <- matrix(0, nrow = lam1_len * lam2_len, ncol = n)
  z_mat <- matrix(0, nrow = lam1_len *  lam2_len, ncol = n)
  s_mat <- matrix(0, nrow = lam1_len *  lam2_len, ncol = n * (n - 1) / 2)
  w_mat <- matrix(0, nrow = lam1_len *  lam2_len, ncol = p)
  q1_mat <- matrix(0, nrow = lam1_len *  lam2_len, ncol = n)
  q2_mat <- matrix(0, nrow = lam1_len *  lam2_len, ncol = n * (n - 1) / 2)
  q3_mat <- matrix(0, nrow = lam1_len *  lam2_len, ncol = p)
  step_vec <- rep(1, lam1_len * lam2_len)
  diff_vec <- rep(tol + 1, lam1_len * lam2_len)
  loss_mat <- matrix(0, nrow = lam1_len * lam2_len, ncol = max_iter + 1)
  
  # run algorithm over the path
  pb <- progressr::progressor(steps = lam1_len * lam2_len)
  for(i in 1 : lam1_len){
    for(j in 1 : lam2_len){
      # index for current solution and initial values
      ind <- (j - 1) * lam1_len + i
      
      #  grid points for lam1 and lam2 are designed in a matrix form:
      #
      #     +--- lam_2 ---+
      #     ||    -->     |
      #     ||    -->     |
      #     |+----------> |
      #   lam_1           |
      #     |             |
      #     |             |
      #     |             |
      #     +-------------+
      #
      # But the variables are stored in a col-major vector form
      #   with `ind <- (j - 1) * lam1_len + i` for coord `(i, j)`
      
      
      # if(j == 1){
      #     ind_old <- i - 1
      # }else{
      #     ind_old <- (j - 2) * lam1_len + i
      # }
      
      # prepare lambda's
      p1_param[1] <- lam1_vec[i]
      p2_param[1] <- lam2_vec[j]
      
      # carry out the algorithm
      res <- RSAVS_Solver_PureR(y_vec, x_mat, l_type = l_type, l_param = l_param, 
                                p1_type = p1_type, p1_param = p1_param, p2_type = p2_type, p2_param = p2_param, 
                                const_r123 = const_r123, const_abc = const_abc, 
                                initial_values = initial_values, additional = additional, 
                                tol = tol, max_iter = max_iter, cd_max_iter = cd_max_iter, cd_tol = cd_tol, 
                                phi = phi, subgroup_benchmark = subgroup_benchmark)
      
      
      # store the results
      # NOTE: we should store the result before prepare initial values for next iteration.
      #       Otherwise it causes the initial values to be all 0 when lam2_len = 1.
      beta_mat[ind, ] <- res$beta_vec
      mu_mat[ind, ] <- res$mu_vec
      z_mat[ind, ] <- res$z_vec
      s_mat[ind, ] <- res$s_vec
      w_mat[ind, ] <- res$w_vec
      q1_mat[ind, ] <- res$q1_vec
      q2_mat[ind, ] <- res$q2_vec
      q3_mat[ind, ] <- res$q3_vec
      step_vec[ind] <- res$current_step
      diff_vec[ind] <- res$diff
      loss_mat[ind, ] <- res$loss_vec
      
      # update inital values for next run
      if(j == lam2_len){
        ind_old <- i
        initial_values <- list(beta_init = beta_mat[ind_old, ], 
                               mu_init = mu_mat[ind_old, ], 
                               z_init = z_mat[ind_old, ], 
                               s_init = s_mat[ind_old, ], 
                               w_init = w_mat[ind_old, ], 
                               q1_init = q1_mat[ind_old, ], 
                               q2_init = q2_mat[ind_old, ], 
                               q3_init = q3_mat[ind_old, ])
      }else{
        initial_values <- list(beta_init = res$beta_vec, 
                               mu_init = res$mu_vec, 
                               z_init = res$z_vec, 
                               s_init = res$s_vec, 
                               w_init = res$w_vec, 
                               q1_init = res$q1_vec, 
                               q2_init = res$q2_vec, 
                               q3_init = res$q3_vec)
      }
      
      pb(message = paste("lam1: ", i, "/", lam1_len, ", lam2: ", j, "/", lam2_len, sep = ""))
    }
  }
  
  # --- Use mBIC to find the best solution in this solution plane ---
  
  res <- list(beta_mat = beta_mat, 
              mu_mat = mu_mat, 
              z_mat = z_mat, 
              s_mat = s_mat, 
              w_mat = w_mat, 
              q1_mat = q1_mat, 
              q2_mat = q2_mat, 
              q3_mat = q3_mat, 
              step_vec = step_vec, 
              diff_vec = diff_vec,
              lam1_vec = lam1_vec, 
              lam2_vec = lam2_vec, 
              loss_mat = loss_mat, 
              const_r123 = const_r123, 
              const_abc = const_abc, 
              initial_values = initial_values, 
              additional = additional)
  
  return(res)
  
  # if(l_type == "L2"){
  #     # beta_left_inv <- solve(t(x_mat) %*% x_mat / n + r3 / 2 * diag(nrow = p))
  #     # d_mat <- RSAVS_Generate_D_Matrix(n)
  #     # mu_left_inv <- solve(diag(nrow = n) / n + r2 / 2 * as.matrix(t(d_mat) %*% d_mat))
  #     res <- RSAVS_LargeN_L2_Rcpp(x_mat, y_vec, n, p, p1_type, p1_param, p2_type, p2_param, lam1_vec, lam2_vec, r2, r3, phi, tol, max_iter)
  # } else{
  #     res <- RSAVS_LargeN_Rcpp(x_mat, y_vec, n, p, l_type, l_param, p1_type, p1_param, p2_type, p2_param, lam1_vec, lam2_vec, r1, r2, r3, phi, tol, max_iter)
  # }
  
  # Idealy, the result of c(lam1[1], lam2[1]) should be mu being median and beta being 0
  # So the result is directly set to this, without actually computing using ADMM
  # But in actuality, our derivation of lam1[1] and lam2[1] is not that accurate, 
  # Hence this direct setting is not that accurate
  # Given this situation, when best_ind being 1
  # We re-choose the best result
  # if(res$best_ind == 1){
  #     min_bic_id <- which.min(res$bic_mat[-1])
  #     res$best_ind <- min_bic_id + 1
  #     res$best_i <- (min_bic_id + 1 - 1) %% length(lam1_vec) + 1
  #     res$best_j <- (min_bic_id + 1 - 1) %/% length(lam1_vec) + 1
  #     res$repick <- TRUE
  # }
}