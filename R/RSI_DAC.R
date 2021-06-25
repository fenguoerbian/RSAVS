RSI_Subsample_ID <- function(n, sub_sample_size = 200){
  # Generate the subsampling id matrix.
  # Each row is a sub sampling result, with sub-sample size being sub_sample_size
  # And the number of sub-samples is the smallest one, but still guarantee every observation is presented at least once.
  # Args: n: number of observations
  #       sub_sample_size: size of each sub sample
  # Return: sample_res: a matrix, each row is one sub-sample
  reminder <- n %% sub_sample_size
  if(reminder == 0){
    sample_num <- n %/% sub_sample_size
    full_list <- 1 : n
  }else{
    sample_num <- n %/% sub_sample_size + 1
    full_list <- c(1 : n, sample(1 : n, size = sub_sample_size - reminder, replace = F))
  }
  
  sample_res <- sample(full_list, size = length(full_list), replace = F)
  sample_res <- matrix(sample_res, ncol = sub_sample_size)
  return(sample_res)
}

RSI_Mu_to_ID <- function(mu_vec, mu_unique){
  # This function converts the original mu vector to a resulting vector, with the same length.
  # But the entries is replaced with the index of this entry in unique(mu_vec), or the provided variable mu_unique
  # This functions is designed for the conquer part of the whole algorithm
  
  # mu_unique <- unique(mu_vec)
  if(missing(mu_unique)){
    mu_unique <- unique(mu_vec)
  }
  mu_vec <- sapply(mu_vec, FUN = function(x, unique_vec){
    return(which(unique_vec == x))
  }, unique_vec = mu_unique)
  return(mu_vec)
}

RSI_Conqure <- function(mu_response, p1_type, p1_param, const_abc = rep(1, 3), 
                        mu_initial, tol = 10^(-3), max_iter = 100, 
                        use_pamk = FALSE, round_digits = 3){
  # This is the conqure part of the algorithm. 
  # It tries to shrink the mu_vec from all subgroups into a more concentrated and reasonable results
  # This can be seen as RSI_Fit without the x_mat.
  # prepare the response
  mu_response <- matrix(mu_response, ncol = 1)
  n <- length(mu_response)
  
  # prepare penalty function
  if(p1_type == 'L'){
    p1_fun <- penalty_lasso
  } else{
    if(p1_type == 'S'){
      p1_fun <- penalty_scad
    } else {
      if(p1_type == 'M'){
        p1_fun <- penalty_mcp
      }else{
        stop("Wrong type of penalty! `p1_type` only supports `L`, `S` and `M`!")
      }
    }
  }
  
  if(missing(mu_initial)){
    message("`initial_values` is missing. Use default settings!")
    mu_initial <- as.vector(mu_response)
  }
  
  # ------ main algorithm ------
  # get initial values
  mu_old <- mu_initial
  d_mat <- RSAVS_Generate_D_Matrix(n)
  
  # compute initial loss
  tmpx <- matrix(rnorm(n * 2), ncol = 2)
  tmpbeta <- rep(0, 2)
  loss_detail <- RSAVS_Compute_Loss_Value(y_vec = mu_response, x_mat = tmpx, l_type = "L1", l_param = NULL, 
                                          p1_type = p1_type, p1_param = p1_param, p2_type = "L", p2_param = 0, 
                                          const_r123 = rep(1, 3), const_abc = const_abc, 
                                          beta_vec = tmpbeta, mu_vec = mu_old, 
                                          z_vec = mu_response - mu_old, 
                                          s_vec = d_mat %*% mu_old, 
                                          w_vec = tmpbeta, 
                                          q1_vec = rep(0, n), q2_vec = rep(0, n * (n - 1) / 2), q3_vec = rep(0, 2))
  
  loss_old <- loss_detail$loss_part1 + loss_detail$loss_part2
  loss_vec <- rep(loss_old, max_iter + 1)
  
  # augment the part of x and y, which will not be changed during the iteration
  x_aug <- as.matrix.csr(0, nrow = n + n * (n - 1) / 2, ncol = n)
  x_aug[1 : n, 1 : n] <- diag(x = 1, nrow = n)
  y_aug <- rbind(mu_response, matrix(0, nrow = n * (n - 1) / 2, ncol = 1))
  
  # algorithm status
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
    
    # 3. check algorithm status
    mu_diff_vec <- d_mat %*% mu_vec
    w_vec_new <- p1_fun(x = abs(mu_diff_vec), param = p1_param, derivative = TRUE)
    diff <- sum(abs(w_vec_new - w_vec))
    if(diff < tol){
      converge_status <- TRUE
    }
    
    current_step <- current_step + 1
    
    mu_old <- mu_vec
    
    loss_detail <- RSAVS_Compute_Loss_Value(y_vec = mu_response, x_mat = tmpx, l_type = "L1", l_param = NULL, 
                                            p1_type = p1_type, p1_param = p1_param, p2_type = "L", p2_param = 0, 
                                            const_r123 = rep(1, 3), const_abc = const_abc, 
                                            beta_vec = tmpbeta, mu_vec = mu_old, 
                                            z_vec = mu_response - mu_old, 
                                            s_vec = d_mat %*% mu_old, 
                                            w_vec = tmpbeta, 
                                            q1_vec = rep(0, n), q2_vec = rep(0, n * (n - 1) / 2), q3_vec = rep(0, 2))
    
    loss_old <- loss_detail$loss_part1 + loss_detail$loss_part2
    loss_vec[current_step] <- loss_old
  }
  
  # modify mu_vec into reasonable subgroups
  if(use_pamk){
    mu_updated <- RSAVS_Determine_Mu(mu_vec)
  }else{
    mu_updated <- RSAVS_Determine_Mu(mu_vec, round_digits = round_digits)
  }
  
  return(mu_updated)
}

RSI_DAC <- function(y_vec, x_mat, sub_sample_size = 200, 
                    l_type = "L1", l_param = NULL, 
                    p1_type = "S", p1_param = c(1, 3.7), 
                    lam1_vec, min_lam1_ratio = 0.03, lam1_len, 
                    initial_values, 
                    const_abc = rep(1, 3), 
                    divide_pamk = FALSE, divide_digit = 1, 
                    conquer_pamk = FALSE, conquer_digit = 1, 
                    max_iter = 100, tol = 10^(-3)){
  # ------ preparation ------
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
  
  # subsampling
  subsample_mat <- RSI_Subsample_ID(n, sub_sample_size)
  subsample_num <- nrow(subsample_mat)
  # variables to store divided results
  beta_mat <- matrix(0, nrow = subsample_num, ncol = p)
  mu_mat <- matrix(0, nrow = subsample_num, ncol = sub_sample_size)
  
  # ------ main algorithm ------
  # 1. divided analysis
  message("Start divided analysis!")
  for(sample_id in 1 : subsample_num){
    message(paste("Divided sub sample analysis at ", sample_id, "/", subsample_num, ".", sep = ""))
    y_subsample <- y_vec[subsample_mat[sample_id, ], , drop = F]
    x_subsample <- x_mat[subsample_mat[sample_id, ], , drop = F]
    # res_subsample <- RSI_Fit(y_vec = y_subsample, x_mat = x_subsample, l_type = l_type, l_param = l_param, 
    #                          p1_type = p1_type, p1_param = p1_param, 
    #                          const_abc = const_abc)
    # res_subsample <- RSI_Fit(x = x_subsample, y = y_subsample, penalty = p_type, p_param = p_param, 
    #                          min_lam_ratio = min_lam_ratio, use_pamk = divide_pamk, max_iter = max_iter, tol = tol)
    # beta_mat[sample_id, ] <- res_subsample$beta
    # mu_mat[sample_id, ] <- res_subsample$mu
  }
  # 2. conquer analysis
}