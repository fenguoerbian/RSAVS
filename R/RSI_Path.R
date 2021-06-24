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
  
}