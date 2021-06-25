RSI_Fit <- function(y_vec, x_mat, l_type = "L1", l_param = NULL, 
                    p1_type = "S", p1_param = c(2, 3.7), 
                    const_abc = rep(1, 3), initial_values, 
                    tol = 10 ^ {-3}, max_iter = 100, use_pamk = FALSE){
  ### preparation ###
  # preparation for x and y #
  y_vec <- matrix(y_vec, ncol = 1)    # make sure y_vec is column vector
  n <- length(y_vec)
  x_mat <- matrix(x_mat, nrow = n)    # make sure x_mat is a matrix, even if it has only one column
  p <- ncol(x_mat)
  
  # preparation for loss function #
  if(l_type != "L1"){
    warning("Currently, only `L1` is supported for `l_type`!")
    l_type <- "L1"
    l_param <- NULL
  }
  
  # preparation for penalties #
  if(p1_type == 'L'){
    p1_fun <- penalty_lasso
  } else{
    if(p1_type == 'S'){
      p1_fun <- penalty_scad
    } else {
      if(p1_type == 'M'){
        p1_fun <- penalty_mcp
      }
    }
  }
  
  # const_a <- const_abc[1]
  # const_b <- const_abc[2]
  # const_c <- const_abc[3]
  
  # prepare initial values
  if(missing(initial_values)){
    message("`initial_values` is missing. Use default settings!")
    initial_values <- list(beta_init = rep(0, p), 
                           mu_init = as.vector(y_vec))
  }
  
  # ------ main algorithm ------
  # --- initial values ---
  beta_old <- initial_values$beta_init
  mu_old <- initial_values$mu_init
  d_mat <- RSAVS_Generate_D_Matrix(n)
  loss_detail <- RSAVS_Compute_Loss_Value(y_vec = y_vec, x_mat = x_mat, l_type = l_type, l_param = l_param, 
                                          p1_type = p1_type, p1_param = p1_param, p2_type = "L", p2_param = 0, 
                                          const_r123 = rep(1, 3), const_abc = const_abc, 
                                          beta_vec = beta_old, mu_vec = mu_old, 
                                          z_vec = y_vec - mu_old - x_mat %*% beta_old, 
                                          s_vec = d_mat %*% mu_old, 
                                          w_vec = beta_old, 
                                          q1_vec = rep(0, n), q2_vec = rep(0, n * (n - 1) / 2), q3_vec = rep(0, p))
  
  loss_old <- loss_detail$loss_part1 + loss_detail$loss_part2 + loss_detail$loss_part3
  loss_vec <- rep(loss_old, max_iter + 1)
  
  # augment the part of x and y, which will not be changed during the iteration
  y_aug <- rbind(y_vec, matrix(0, nrow = n * (n - 1) / 2, ncol = 1))
  x_aug <- as.matrix.csr(0, nrow = n + n * (n - 1) / 2, ncol = n + p)
  x_aug[1 : n, 1 : n] <- diag(x = 1, nrow = n)
  x_aug[1 : n, (n + 1) : (n + p)] <- x_mat
  
  # --- algorithm status ---
  current_step <- 1
  diff <- tol + 1
  converge_status <- FALSE
  
  while((current_step <= max_iter) & (diff > tol)){
    # 1. augment the data
    # 1.1 current difference vector
    mu_diff_old <- d_mat %*% mu_old
    # 1.2 derivative vector
    w_vec <- p1_fun(x = abs(mu_diff_old), param = p1_param, derivative = TRUE)
    # 1.3 generate the w_mat
    w_mat <- RSAVS_Generate_D_Matrix(n, w_vec = as.vector(w_vec))
    # 1.4 update x_aug
    x_aug[(n + 1) : (n + n * (n - 1) / 2), 1 : n] <- const_abc[1] * const_abc[2] * w_mat
    
    # 2. fit the model
    tmp <- quantreg::rq.fit.sfn(a = x_aug, y = y_aug, tau = 0.5)
    mu_vec <- tmp$coefficients[1 : n]
    beta_vec <- tmp$coefficients[(n + 1) : (n + p)]
    
    # 3. check algorithm status
    mu_diff_vec <- d_mat %*% mu_vec
    w_vec_new <- p1_fun(x = abs(mu_diff_vec), param = p1_param, derivative = TRUE)
    diff <- sum(abs(w_vec_new - w_vec))
    if(diff < tol){
      converge_status <- TRUE
    }
    
    current_step <- current_step + 1
    
    beta_old <- beta_vec
    mu_old <- mu_vec
    
    loss_detail <- RSAVS_Compute_Loss_Value(y_vec = y_vec, x_mat = x_mat, l_type = l_type, l_param = l_param, 
                                            p1_type = p1_type, p1_param = p1_param, p2_type = "L", p2_param = 0, 
                                            const_r123 = rep(1, 3), const_abc = const_abc, 
                                            beta_vec = beta_old, mu_vec = mu_old, 
                                            z_vec = y_vec - mu_old - x_mat %*% beta_old, 
                                            s_vec = d_mat %*% mu_old, 
                                            w_vec = beta_old, 
                                            q1_vec = rep(0, n), q2_vec = rep(0, n * (n - 1) / 2), q3_vec = rep(0, p))
    
    loss_old <- loss_detail$loss_part1 + loss_detail$loss_part2 + loss_detail$loss_part3
    loss_vec[current_step] <- loss_old
  }
  
  # modify mu_vec into reasonable subgroups
  if(use_pamk){
    mu_updated <- RSAVS_Determine_Mu(mu_vec)
  }else{
    mu_updated <- RSAVS_Determine_Mu(mu_vec, round_digits = 2)
  }
  
  res <- list(beta_vec = beta_vec, 
              mu_vec = mu_vec, 
              mu_updated = mu_updated, 
              current_step = current_step - 1, 
              diff = diff, 
              loss_vec = loss_vec, 
              converge_status = converge_status)
  return(res)
}


RSI_Path <- function(y_vec, x_mat, l_type = "L1", l_param = NULL, p1_type = "S", p1_param = c(2, 3.7), 
                     lam1_vec, min_lam1_ratio = 0.03, lam1_len, 
                     initial_values, 
                     const_abc = rep(1, 3), use_pamk = FALSE, max_iter = 100, tol = 10 ^{-3}){
  # This function finds the solution path of Robust Subgroup Identification
  # Args: x: covariate matrix, NO intercept term
  #       y: response vector
  #       penalty: type of penalty function, character variable
  #       p_param: list variable, parameters for the penalty function
  #                'L': Lasso, lambda = p_param$lambda
  #                'S': SCAD, lambda = p_param$lambda, a = p_param$a
  #                'M': MCP, lambda = p_param$lambda, a = p_param$a
  #                Note: actually, the lambda will be provided from lam_seq
  #       lam_seq: lambda sequence, from big to small
  #       lam_len: length of the lambda sequence, 
  #                If lam_seq is missing, then use default method to get a sequence of lambda with length = lam_len
  #       min_lam_ratio: this ratio  = lam_min / lam_max
  #       use_pamk: boolen, whether pamk will be used for post-procedure clustering over \mu
  #       max_iter, tol: convergence and iteration related parameters
  # Returns:
  #       res: list of pathes of: \mu, \beta, k(number of iteration), and convergence indicator
  
  ### preparation ###
  # preparation for x and y #
  y_vec <- matrix(y_vec, ncol = 1)    # make sure y_vec is column vector
  n <- length(y_vec)
  x_mat <- matrix(x_mat, nrow = n)    # make sure x_mat is a matrix, even if it has only one column
  p <- ncol(x_mat)
  
  # preparation for loss function #
  if(l_type != "L1"){
    warning("Currently, only `L1` is supported for `l_type`!")
    l_type <- "L1"
    l_param <- NULL
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
  
  # check const_abc
  const_abc <- abs(const_abc)
  
  # preparation for lam_vec
  if(missing(lam1_vec)){
    # newer version
    message("lam1_vec is missing, use default values...")
    
    if(missing(lam1_len) | missing(min_lam1_ratio)){
      stop("Both `lam1_len` and `min_lam1_ratio` must be provided in order to generate default lambda vector!")
    }else{
      lam1_max <- RSAVS_Get_Lam_Max(y_vec = y_vec, l_type = l_type, l_param = l_param, 
                                    lam_ind = 1, 
                                    const_abc = const_abc, eps = 10^(-6))
      lam1_min <- lam1_max * min_lam1_ratio
      if(p1_type != "L"){
        lam1_max <- lam1_max * 100    # safe guard for non-convex penalties
      }
      
      # lam1_vec <- exp(seq(from = log(lam1_max), to = log(lam1_min), length.out = lam1_len))
      lam1_vec <- exp(seq(from = log(lam1_min), to = log(lam1_max), length.out = lam1_len))
    }
  }else{
    # lam1_vec <- sort(abs(lam1_vec), decreasing = T)
    # lam1_length = length(lam1_vec)
    lam1_len <- length(lam1_vec)
    lam1_vec <- sort(abs(lam1_vec), decreasing = FALSE)
  }
  
  # prepare initial values
  if(missing(initial_values)){
    message("`initial_values` is missing, use default settings!")
    initial_values <- list(beta_init = rep(0, p), 
                           mu_init = as.vector(y_vec))
  }
  
  # prepare variables to store the results
  beta_mat <- matrix(0, nrow = lam1_len, ncol = p)
  mu_mat <- matrix(0, nrow = lam1_len, ncol = n)
  mu_updated_mat <- matrix(0, nrow = lam1_len, ncol = n)
  step_vec <- rep(1, lam1_len)
  diff_vec <- rep(tol + 1, lam1_len)
  converge_vec <- rep(FALSE, lam1_len)
  loss_mat <- matrix(0, nrow = lam1_len, ncol = max_iter + 1)
  
  # main algorithm
  pb <- progressr::progressor(steps = lam1_len)
  for(i in 1 : lam1_len){
    # index for current solution
    ind <- i
    
    # prepare lambda 
    p1_param[1] <- lam1_vec[i]
    
    # carry out the algorithm
    res <- RSI_Fit(y_vec, x_mat, l_type = l_type, l_param = l_param, 
                   p1_type = p1_type, p1_param = p1_param, 
                   const_abc = const_abc, initial_values = initial_values, 
                   tol = tol, max_iter = max_iter, use_pamk = use_pamk)
    
    # store the results
    beta_mat[ind, ] <- res$beta_vec
    mu_mat[ind, ] <- res$mu_vec
    mu_updated_mat[ind, ] <- res$mu_updated
    step_vec[ind] <- res$current_step
    diff_vec[ind] <- res$diff
    converge_vec[ind] <- res$converge_status
    loss_mat[ind, ] <- res$loss_vec
    
    # prepare initial values for the next run
    initial_values <- list(beta_init = res$beta_vec, 
                           mu_init = res$mu_updated    # use updated mu vector as initial values
                                                       # it seems to provide better results than original mu_vec
                           )
    
    # progress information
    pb(message = paste("lam1: ", i, "/", lam1_len, sep = ""))
  }
  
  res <- list(beta_mat = beta_mat, 
              mu_mat = mu_mat, 
              mu_updated_mat = mu_updated_mat, 
              step_vec = step_vec, 
              diff_vec = diff_vec, 
              converge_vec = converge_vec, 
              loss_mat = loss_mat, 
              lam1_vec = lam1_vec, 
              const_abc = const_abc)
  return(res)
}