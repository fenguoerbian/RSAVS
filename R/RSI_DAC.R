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
                    initial_values, 
                    const_abc = rep(1, 3), 
                    divide_pamk = FALSE, divide_digit = 1, 
                    conquer_pamk = FALSE, conquer_digit = 1, conquer_use_updated_mu = TRUE, 
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
  mu_updated_mat <- matrix(0, nrow = subsample_num, ncol = sub_sample_size)
  
  # initial_values 
  if(missing(initial_values)){
    message("`initial_values` is missing. Use default settings!")
    initial_values <- list(beta_init = rep(0, p), 
                           mu_init = as.vector(y_vec))
  }
  # ------ main algorithm ------
  # 1. divided analysis
  message("Start divided analysis!")
  for(sample_id in 1 : subsample_num){
    message(paste("Divided sub sample analysis at ", sample_id, "/", subsample_num, ".", sep = ""))
    # prepare the subsample
    y_subsample <- y_vec[subsample_mat[sample_id, ], , drop = F]
    x_subsample <- x_mat[subsample_mat[sample_id, ], , drop = F]
    initial_subsample <- list(beta_init = initial_values$beta_init, 
                              mu_init = initial_values$mu_init[subsample_mat[sample_id, ]])
    
    # subsample analysis
    res_subsample <- RSI_Fit(y_vec = y_subsample, x_mat = x_subsample, l_type = l_type, l_param = l_param, 
                             p1_type = p1_type, p1_param = p1_param, 
                             const_abc = const_abc, initial_values = initial_subsample, 
                             tol = tol, max_iter = max_iter, 
                             use_pamk = divide_pamk, round_digits = divide_digit)
    
    # save the results
    beta_mat[sample_id, ] <- res_subsample$beta_vec
    mu_mat[sample_id, ] <- res_subsample$mu_vec
    mu_updated_mat[sample_id, ] <- res_subsample$mu_updated
  }
  
  # 2. conquer analysis
  message("Start conquered analysis!")
  # conquered analysis of mu
  # 1st, form the full mu vector from the mu_mat of Divided analysis
  mu_vec <- rep(0, n)
  for(i in 1 : n){
    id <- which(subsample_mat == i)
    if(conquer_use_updated_mu){
      mu_vec[i] <- mean(mu_updated_mat[id])
    }else{
      mu_vec[i] <- mean(mu_mat[id])
    }
  }
  
  # 2nd, find the ID form of this vector
  mu_unique <- unique(mu_vec)
  mu_ids <- RSI_Mu_to_ID(mu_vec = mu_vec, mu_unique = mu_unique)
  
  # 3rd, conquer this mu_unique vector
  message("Start conquering algorithm")
  mu_unique <- RSI_Conqure(mu_response = mu_unique, p1_type = p1_type, p1_param = p1_param, 
                           const_abc = c(length(mu_unique), 1, 1), tol = tol, max_iter = max_iter, 
                           use_pamk = conquer_pamk, round_digits = conquer_digit)
  
  # 4th, update the mu vector and its id form
  mu_conquered <- mu_unique[mu_ids]

  # conquered analysis of beta
  beta_vec <- apply(beta_mat, 2, mean)
  
  # return the result
  res <- list(mu_conquered_vec = mu_conquered, 
              beta_vec = beta_vec, 
              mu_origin_vec = mu_vec)
  
  return(res)
}

RSI_DAC_Path <- function(y_vec, x_mat, sub_sample_size, 
                         l_type = "L1", l_param = NULL, 
                         p1_type = "S", p1_param = c(1, 3.7), 
                         lam1_vec, min_lam1_ratio = 0.03, lam1_len, 
                         initial_values, 
                         const_abc = rep(1, 3), 
                         divide_pamk = FALSE, divide_digit = 1, 
                         conquer_pamk = FALSE, conquer_digit = 1, conquer_use_updated_mu = TRUE, 
                         max_iter = 100, tol = 10^(-3)){
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
  # step_vec <- rep(1, lam1_len)
  # diff_vec <- rep(tol + 1, lam1_len)
  # converge_vec <- rep(FALSE, lam1_len)
  # loss_mat <- matrix(0, nrow = lam1_len, ncol = max_iter + 1)
  
  # main algorithm
  pb <- progressr::progressor(steps = lam1_len)
  for(i in 1 : lam1_len){
    # index for current solution
    ind <- i
    
    # prepare lambda 
    p1_param[1] <- lam1_vec[i]
    
    # carry out the algorithm
    res <- RSI_DAC(y_vec = y_vec, x_mat = x_mat, l_type = l_type, l_param = l_param, 
                   p1_type = p1_type, p1_param = p1_param, 
                   initial_values = initial_values, const_abc = const_abc, 
                   divide_pamk = divide_pamk, divide_digit = divide_digit, 
                   conquer_pamk = conquer_pamk, conquer_digit = conquer_digit, conquer_use_updated_mu = conquer_use_updated_mu, 
                   tol = tol, max_iter = max_iter)
    
    # store the results
    beta_mat[ind, ] <- res$beta_vec
    mu_mat[ind, ] <- res$mu_origin_vec
    mu_updated_mat[ind, ] <- res$mu_conquered_vec

    # prepare initial values for the next run
    initial_values <- list(beta_init = res$beta_vec, 
                           mu_init = res$mu_conquered_vec    # use updated mu vector as initial values
                           # it seems to provide better results than original mu_vec
    )
    
    # progress information
    pb(message = paste("lam1: ", i, "/", lam1_len, sep = ""))
  }
  
  res <- list(beta_mat = beta_mat, 
              mu_mat = mu_mat, 
              mu_updated_mat = mu_updated_mat, 
              # step_vec = step_vec, 
              # diff_vec = diff_vec, 
              # converge_vec = converge_vec, 
              # loss_mat = loss_mat, 
              lam1_vec = lam1_vec, 
              const_abc = const_abc)
  return(res)  
}