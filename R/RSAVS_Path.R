#' Robust subgroup analysis and variable selection simultaneously
#' 
#' This function utilize the cpp solver when carring out robust subgroup 
#'   analysis and variable selection simultaneously. And it is the core solver
#'   for \code{\link{RSAVS_Path}}.
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
#'   Note: This function searches the whole lam1_vec * lam2_vec grid for the best solution. 
#'   Hence the \code{lambda}s provided in these parameters serve only as placeholder 
#'   and will be ignored and overwritten in the actual computation.
#' @param lam1_vec,lam2_vec numerical vectors of customized lambda vectors. 
#'   For \code{lam1_vec}, it's preferred to be in the order from small to big.
#' @param min_lam_ratio the ratio between the minimal and maximal lambda, equals to (minimal lambda) / (maximal lambda).
#'   The default value is 0.03.
#' @param lam1_len,lam2_len integers, lengths of the auto-generated lambda vectors.
#' @param initial_values list of vector, providing initial values for the algorithm. 
#' @param additional a list providing additional variables needed during the algorithm.
#' @param phi numerical variable. A parameter needed for mBIC.
#' @param const_r123 a length-3 numerical vector, providing the scalars needed in the 
#'   augmented lagrangian part of the ADMM algorithm
#' @param const_abc a length-3 numeric vector, providing the scalars to adjust weight
#'   of regression function, penalty for subgroup identification and penalty for 
#'   variable selection in the overall objective function. Defaults to \code{c(1, 1, 1)}.
#' @param tol numerical, convergence tolerance for the algorithm.
#' @param cd_tol numerical, convergence tolerance for the coordinate descent part 
#'   when updating \code{mu} and \code{beta}.
#' @param max_iter integer, max number of iteration during the algorithm.
#' @param cd_max_iter integer, max number of iteration during the coordinate descent
#'   update of \code{mu} and \code{beta}. If set to 0, will use analytical solution(
#'   instead of coordinate descent algorithm) to update \code{mu} and \code{beta}.
#' @param subgroup_benchmark bool. Whether this call should be taken as a benchmark of subgroup identification. 
#'   If \code{TRUE}, then the penalty for variable selection will be surpressed to a minimal value.
#' @seealso \code{\link{RSAVS_Path}}, \code{\link{RSAVS_Path_PureR}}, \code{\link{RSAVS_LargeN}}
RSAVS_Solver <- function(y_vec, x_mat, l_type = "L1", l_param = NULL, 
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
  res <- RSAVS_Solver_Cpp(y_vec, x_mat, n, p, l_type, l_param, 
                          p1_type, p1_param, p2_type, p2_param, 
                          const_r123, const_abc, tol, max_iter, cd_tol, cd_max_iter, 
                          initial_values, additional, phi)
  
  return(res)
}


#' Robust subgroup analysis and variable selection simultaneously
#' 
#' This function carries out robust subgroup analysis and variable selection 
#' simultaneously. It utilizes the cpp core solver and supports different 
#'   types of loss functions and penalties.
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
#'   Note: This function searches the whole lam1_vec * lam2_vec grid for the best solution. 
#'   Hence the \code{lambda}s provided in these parameters serve only as placeholder 
#'   and will be ignored and overwritten in the actual computation.
#' @param lam1_vec,lam2_vec numerical vectors of customized lambda vectors. 
#'   For \code{lam1_vec}, it's preferred to be in the order from small to big.
#' @param min_lam_ratio the ratio between the minimal and maximal lambda, equals to (minimal lambda) / (maximal lambda).
#'   The default value is 0.03.
#' @param lam1_len,lam2_len integers, lengths of the auto-generated lambda vectors.
#' @param initial_vec list of vector, providing initial values for the algorithm. 
#' @param phi numerical variable. A parameter needed for mBIC.
#' @param const_r123 a length-3 numerical vector, providing the scalars needed in the 
#'   augmented lagrangian part of the ADMM algorithm
#' @param const_abc a length-3 numeric vector, providing the scalars to adjust weight
#'   of regression function, penalty for subgroup identification and penalty for 
#'   variable selection in the overall objective function. Defaults to \code{c(1, 1, 1)}.
#' @param tol numerical, convergence tolerance for the algorithm.
#' @param cd_tol numerical, convergence tolerance for the coordinate descent part 
#'   when updating \code{mu} and \code{beta}.
#' @param max_iter integer, max number of iteration during the algorithm.
#' @param cd_max_iter integer, max number of iteration during the coordinate descent
#'   update of \code{mu} and \code{beta}. If set to 0, will use analytical solution(
#'   instead of coordinate descent algorithm) to update \code{mu} and \code{beta}.
#' @param subgroup_benchmark bool. Whether this call should be taken as a benchmark of subgroup identification. 
#'   If \code{TRUE}, then the penalty for variable selection will be surpressed to a minimal value.
#' @seealso \code{\link{RSAVS_Solver}}, \code{\link{RSAVS_Path_PureR}}, \code{\link{RSAVS_LargeN}}
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
#' res <- RSAVS_Path(y_vec = y_vec, x_mat = x_mat, 
#'                  lam1_len = 50, lam2_len = 40, 
#'                  phi = 5)
#' @export
RSAVS_Path <- function(y_vec, x_mat, l_type = "L1", l_param = NULL, 
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
    d_mat <- RSAVS_Generate_D_Matrix(n)
    
    initial_values <- list(beta_init = rep(0, p), 
                           mu_init = mu_init, 
                           z_init = y_vec - mu_init, 
                           # s_init = rep(0, n * (n - 1) / 2), 
                           s_init = as.vector(d_mat %*% mu_init), 
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
  
  
  mu_beta_lhs <- matrix(0, nrow = 1, ncol = 1)    # left part for updating beta and mu together
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
      res <- RSAVS_Solver(y_vec, x_mat, l_type = l_type, l_param = l_param, 
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
              const_abc = const_abc)
  
  return(res)
}