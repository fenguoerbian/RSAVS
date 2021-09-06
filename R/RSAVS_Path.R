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
#' @param update_mu list of parameters for updating \code{mu_vec} in the algorithm into meaningful subgroup structure.
#'   Defaults to \code{NULL}, which means there is no update performed. The update of \code{mu_vec} is carried out through
#'   \code{RSAVS_Determine_Mu} and the necessary parameters in \code{update_mu} are:
#'   \itemize{
#'     \item \code{UseS}: a bool variable, whether the \code{s_vec} should be used to provide subgroup structure information.
#'     \item \code{klim}: a length-3 integer vector, given the range of number of cluster for considering.
#'     \item \code{usepam}: a bool variable, whether to use \code{pam} for clustering.
#'     \item \code{round_digits}: non-negative integer digits, indicating the rounding digits when merging \code{mu_vec}
#'   }
#'   Please refer to \code{RSAVS_Determine_Mu} to find out more details about how the algorithm works
#' @param loss_track boolen, whether to track the value of objective function(loss value) during each iteration.
#' @param diff_update boolen, whether to update the difference between each iteration. If set to \code{FALSE},
#'   the algorithm will still stop when it reaches \code{max_iter}.
#' @param omp_zsw a length-three integer vector, defaults to \code{c(1, 4, 1)}. It represents how many
#'   parallel threads to be used during the update of \code{z}, \code{s} and \code{w} respectively.
#' @param eigen_pnum integer number, representing the number of Eigen threads for matrix computation, 
#'   defaults to 4. 
#' @param s_v2 boolen, whether to use the updated and faster version during the computation of 
#'   \code{s} and \code{q2}
#' @seealso \code{\link{RSAVS_Path}}, \code{\link{RSAVS_Path_PureR}}, \code{\link{RSAVS_LargeN}}
RSAVS_Solver <- function(y_vec, x_mat, l_type = "L1", l_param = NULL, 
                         p1_type = "S", p1_param = c(2, 3.7), p2_type = "S", p2_param = c(2, 3.7), 
                         const_r123, const_abc = rep(1, 3), 
                         initial_values, additional, tol = 0.001, max_iter = 10, cd_max_iter = 1, cd_tol = 0.001, 
                         phi = 1.0, subgroup_benchmark = FALSE, update_mu = NULL, 
                         loss_track = FALSE, diff_update = TRUE, 
                         omp_zsw = c(1, 4, 1), eigen_pnum = 1, s_v2 = TRUE){
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
                          initial_values, additional, phi, loss_track, diff_update, 
                          omp_zsw, eigen_pnum, s_v2)
  # update the resulting mu vector into meaningfull subgroup results
  if(is.null(update_mu)){
    res$mu_updated_vec <- res$mu_vec
  }else{
    d_mat <- RSAVS_Generate_D_Matrix(n)
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
      group_res <- RSAVS_S_to_Groups(res$s_vec, n)
      mu_updated_vec <- RSAVS_Determine_Mu(res$mu_vec, group_res, klim = klim, usepam = usepam)
    }else{
      mu_updated_vec <- RSAVS_Determine_Mu(res$mu_vec, klim = klim, usepam = usepam, round_digits = round_digits)
    }
    res$mu_updated_vec <- mu_updated_vec
    res$s_vec <- SparseM::as.matrix(d_mat %*% mu_updated_vec)
  }
  
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
#' @param lam1_sort,lam2_sort boolen, whether to force sorting the provided \code{lam1_vec} and \code{lam2_vec}.
#'   By default, \code{lam1_vec} will sort in increasing order while \code{lam2_vec} in decreasing order.
#' @param min_lam1_ratio,min_lam2_ratio the ratio between the minimal and maximal 
#'   lambda, equals to (minimal lambda) / (maximal lambda). The default value is 0.03.
#' @param lam1_max_ncvguard a safe guard constant for \code{lam1_max} when the penalty is nonconvex such as SCAD and MCP.
#' @param lam1_len,lam2_len integers, lengths of the auto-generated lambda vectors.
#' @param initial_values list of vector, providing initial values for the algorithm. 
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
#' @param update_mu list of parameters for updating \code{mu_vec} in the algorithm into meaningful subgroup structure.
#'   Defaults to \code{NULL}, which means there is no update performed. The update of \code{mu_vec} is carried out through
#'   \code{RSAVS_Determine_Mu} and the necessary parameters in \code{update_mu} are:
#'   \itemize{
#'     \item \code{UseS}: a bool variable, whether the \code{s_vec} should be used to provide subgroup structure information.
#'     \item \code{klim}: a length-3 integer vector, given the range of number of cluster for considering.
#'     \item \code{usepam}: a bool variable, whether to use \code{pam} for clustering.
#'     \item \code{round_digits}: non-negative integer digits, indicating the rounding digits when merging \code{mu_vec}
#'   }
#'   Please refer to \code{RSAVS_Determine_Mu} to find out more details about how the algorithm works
#' @param loss_track boolen, whether to track the value of objective function(loss value) during each iteration.
#' @param diff_update boolen, whether to update the difference between each iteration. If set to \code{FALSE},
#'   the algorithm will still stop when it reaches \code{max_iter}.
#' @param omp_zsw a length-three integer vector, defaults to \code{c(1, 4, 1)}. It represents how many
#'   parallel threads to be used during the update of \code{z}, \code{s} and \code{w} respectively.
#' @param eigen_pnum integer number, representing the number of Eigen threads for matrix computation, 
#'   defaults to 4. 
#' @param s_v2 boolen, whether to use the updated and faster version during the computation of 
#'   \code{s} and \code{q2}
#' @param dry_run boolen, whether this is a so-called 'dry-run'. A dry-run will not carry out the real
#'   core solver, but only prepares and return the necessary initial values, additional values and 
#'   solution plane for the core solver.
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
                       lam1_vec = NULL, lam2_vec = NULL, lam1_sort = TRUE, lam2_sort = TRUE, 
                       min_lam1_ratio = 0.03, min_lam2_ratio = 0.03, lam1_max_ncvguard = 100,
                       lam1_len = 50, lam2_len = 40, 
                       const_r123 = NULL, const_abc = rep(1, 3), 
                       initial_values = NULL, phi = 1.0, tol = 0.001, max_iter = 10, 
                       cd_max_iter = 1, cd_tol = 0.001, 
                       subgroup_benchmark = FALSE, update_mu = NULL, loss_track = FALSE, diff_update = TRUE, 
                       omp_zsw = c(1, 4, 1), eigen_pnum = 4, s_v2 = TRUE, dry_run = FALSE){
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
  if(is.null(const_r123)){
    message("`const_r123` is missing. Use default settings!")
    r1 <- 2
    
    if(p1_type == "L"){
      r2 <- 2
    }
    if(p1_type == "S"){
      r2 <- max(1.3 * 1 / (p1_gamma - 1), 1 / (p1_gamma - 1) + 1, 2)
    }
    if(p1_type == "M"){
      r2 <- max(1.3 * 1 / p1_gamma, 1 / p1_gamma + 1, 2)
    }
    
    if(p2_type == "L"){
      r3 <- 2
    }
    if(p2_type == "S"){
      r3 <- max(1.3 * 1 / (p2_gamma - 1), 1 / (p2_gamma - 1) + 1, 2)
    }
    if(p2_type == "M"){
      r3 <- max(1.3 * 1 / p2_gamma, 1 / p2_gamma + 1, 2)
    }
    const_r123 <- c(r1, r2, r3)
  }else{
    const_r123 <- abs(const_r123)
    r1 <- const_r123[1]
    r2 <- const_r123[2]
    r3 <- const_r123[3]
    
    # check for r2
    r2min <- 2
    if(p1_type == "S"){
      r2min <- max(1.3 * 1 / (p1_param[2] - 1), 1 / (p1_param[2] - 1) + 1, 2)
    }
    if(p1_type == "M"){ 
      r2min <- max(1.3 * 1 / p1_param[2], 1 / p1_param[2] + 1, 2)
    }
    
    if(r2 < r2min){
      message("r2 = ", r2, "< r2min = ", r2min, ". Modify value of r2!")
      r2 <- r2min
      const_r123[2] <- r2
    }
    # check for r3
    r3min <- 2
    if(p2_type == "S"){
      r3min <- max(1.3 * 1 / (p2_param[2] - 1), 1 / (p2_param[2] - 1) + 1, 2)
    }
    if(p2_type == "M"){
      r3min <- max(1.3 * 1 / p2_param[2], 1 / p2_param[2] + 1, 2)
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
  if(is.null(lam1_vec)){
    # newer version
    message("lam1_vec is missing, use default values...")
    lam1_max <- RSAVS_Get_Lam_Max(y_vec = y_vec, l_type = l_type, l_param = l_param, 
                                  lam_ind = 1, 
                                  const_abc = const_abc, eps = 10^(-6))
    lam1_min <- lam1_max * min_lam1_ratio
    if(p1_type != "L"){
      lam1_max <- lam1_max * lam1_max_ncvguard    # safe guard for non-convex penalties
    }
    
    # lam1_vec <- exp(seq(from = log(lam1_max), to = log(lam1_min), length.out = lam1_len))
    lam1_vec <- exp(seq(from = log(lam1_min), to = log(lam1_max), length.out = lam1_len))
  }else{
    # lam1_vec <- sort(abs(lam1_vec), decreasing = T)
    # lam1_length = length(lam1_vec)
    lam1_len <- length(lam1_vec)
    if(lam1_sort){
      lam1_vec <- sort(abs(lam1_vec), decreasing = FALSE)
    }else{
      lam1_vec <- abs(lam1_vec)
    }
    
  }
  
  if(is.null(lam2_vec)){
    # newer version
    message("lam2_vec is missing, use default values...")
    lam2_max <- RSAVS_Get_Lam_Max(y_vec = y_vec, x_mat = x_mat, 
                                  l_type = l_type, l_param = l_param, 
                                  lam_ind = 2, 
                                  const_abc = const_abc, eps = 10^(-6))
    lam2_min <- lam2_max * min_lam2_ratio
    lam2_vec <- exp(seq(from = log(lam2_max), to = log(lam2_min), length.out = lam2_len))
  }else{
    lam2_len <- length(lam2_vec)
    if(lam2_sort){
      lam2_vec <- sort(abs(lam2_vec), decreasing = TRUE)
    }else{
      lam2_vec <- abs(lam2_vec)
    }
    
    
  }
  
  if(subgroup_benchmark){
    # supress lambda for variable selection if this is subgroup benchmark
    lam2_vec <- lam2_vec * 0.05
  }
  
  # --- prepare initial values for the algorithm ---
  if(is.null(initial_values)){
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
  
  if(dry_run){
    res <- list(
      lam1_vec = lam1_vec, 
      lam2_vec = lam2_vec, 
      const_r123 = const_r123, 
      const_abc = const_abc, 
      additional = additional
    )
    
    return(res)
  }
  
  # --- variables storing the results ---
  beta_mat <- matrix(0, nrow = lam1_len * lam2_len, ncol = p)
  mu_mat <- matrix(0, nrow = lam1_len * lam2_len, ncol = n)
  mu_updated_mat <- matrix(0, nrow = lam1_len * lam2_len, ncol = n)
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
                          phi = phi, subgroup_benchmark = subgroup_benchmark, update_mu = update_mu, 
                          loss_track = loss_track, diff_update = diff_update, 
                          omp_zsw = omp_zsw, eigen_pnum = eigen_pnum, s_v2 = s_v2)
      
      
      # store the results
      # NOTE: we should store the result before prepare initial values for next iteration.
      #       Otherwise it causes the initial values to be all 0 when lam2_len = 1.
      beta_mat[ind, ] <- res$beta_vec
      mu_mat[ind, ] <- res$mu_vec
      mu_updated_mat[ind, ] <- res$mu_updated_vec
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
                               mu_init = mu_updated_mat[ind_old, ], 
                               z_init = z_mat[ind_old, ], 
                               s_init = s_mat[ind_old, ], 
                               w_init = w_mat[ind_old, ], 
                               q1_init = q1_mat[ind_old, ], 
                               q2_init = q2_mat[ind_old, ], 
                               q3_init = q3_mat[ind_old, ])
      }else{
        initial_values <- list(beta_init = res$beta_vec, 
                               mu_init = res$mu_updated_vec, 
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
              mu_updated_mat = mu_updated_mat, 
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
              const_abc = const_abc
  )
  
  return(res)
}

RSAVS_Simple_Path <- function(y_vec, x_mat, l_type = "L1", l_param = NULL, 
                              p1_type = "S", p1_param = c(2, 3.7), p2_type = "S", p2_param = c(2, 3.7), 
                              lam1_vec = NULL, lam2_vec = NULL, 
                              min_lam1_ratio = 0.03, min_lam2_ratio = 0.03, lam1_max_ncvguard = 100,
                              lam1_len = 100, lam2_len = 80, 
                              const_r123_s1, const_r123_s2, 
                              const_abc_s1 = rep(1, 3), const_abc_s2 = rep(1, 3), 
                              bic_const_a_s1 = 1, bic_const_a_s2 = 1, 
                              bic_double_loglik_s1 = TRUE, bic_double_loglik_s2 = FALSE, 
                              initial_values_s1, initial_values_s2, 
                              phi_s1 = 5.0, phi_s2 = 5.0, 
                              tol_s1 = 0.001, tol_s2 = 0.001, 
                              max_iter_s1 = 100, max_iter_s2 = 100, 
                              cd_max_iter_s1 = 1, cd_max_iter_s2 = 1, 
                              cd_tol_s1 = 0.001, cd_tol_s2 = 0.001, 
                              subgroup_benchmark_s1 = FALSE, subgroup_benchmark_s2 = TRUE, 
                              update_mu_s11 = NULL, update_mu_s12 = list(useS = FALSE, round_digits = NULL),
                              update_mu_s21 = NULL, update_mu_s22 = list(useS = FALSE, round_digits = NULL),
                              complex_upper_bound = 15, 
                              loss_track_s1 = FALSE, loss_track_s2 = FALSE, 
                              diff_update_s1 = TRUE, diff_update_s2 = TRUE, 
                              omp_zsw = c(1, 4, 1), eigen_pnum = 4, s_v2 = TRUE){
  
  message("------ Analysis stage1: search over lam2_vec ------")
  message("--- a dry run to prepare some variables ---")
  dry_run_s1 <- RSAVS_Path(y_vec = y_vec, x_mat = x_mat, l_type = l_type, l_param = l_param, 
                           p1_type = p1_type, p1_param = p1_param, p2_type = p2_type, p2_param = p2_param, 
                           lam1_vec = lam1_vec, lam2_vec = lam2_vec, lam1_sort = FALSE, lam2_sort = FALSE, 
                           min_lam1_ratio = min_lam1_ratio, min_lam2_ratio = min_lam2_ratio, 
                           lam1_max_ncvguard = lam1_max_ncvguard, subgroup_benchmark = FALSE, 
                           lam1_len = lam1_len, lam2_len = lam2_len, 
                           const_r123 = const_r123_s1, const_abc = const_abc_s1, 
                           initial_values = initial_values_s1, 
                           cd_max_iter = cd_max_iter_s1, 
                           omp_zsw = omp_zsw, eigen_pnum = eigen_pnum, s_v2 = s_v2, 
                           dry_run = TRUE
                           )
  
  if(is.null(initial_values_s1)){
    mu0 <- RSAVS_Get_Mu0(y_vec, l_type, l_param)
    initial_values_s1 <- list(beta_init = rep(0, p), 
                              mu_init = rep(mu0, n * (n - 1) / 2), 
                              z_init = as.vector(y_vec - mu0), 
                              s_init = rep(0, n * (n - 1) / 2), 
                              w_init = rep(0, p), 
                              q1_init = rep(0, n), 
                              q2_init = rep(0, n * (n - 1) / 2), 
                              q3_init = rep(0, p))
  }
  
  lam1_chosen <- max(dry_run_s1$lam1_vec)
  lam2_vec_s1 <- dry_run_s1$lam2_vec
  
  message("--- real analysis ---")
  res_s1 <- RSAVS_Path(y_vec = y_vec, x_mat = x_mat, l_type = l_type, l_param = l_param, 
                       p1_type = p1_type, p1_param = p1_param, p2_type = p2_type, p2_param = p2_param, 
                       lam1_vec = lam1_chosen, lam2_vec = lam2_vec_s1, lam2_sort = FALSE, 
                       # min_lam2_ratio = min_lam2_ratio, lam2_len = lam2_len, 
                       subgroup_benchmark = subgroup_benchmark_s1, 
                       const_r123 = const_r123_s1, const_abc = const_abc_s1, 
                       initial_values = initial_values_s1, phi = phi_s1, 
                       tol = tol_s1, max_iter = max_iter_s1, 
                       cd_tol = cd_tol_s1, cd_max_iter = cd_max_iter_s1, 
                       update_mu = update_mu_s11, 
                       loss_track = loss_track_s1, diff_update = diff_update_s1, 
                       omp_zsw = omp_zsw, eigen_pnum = eigen_pnum, s_v2 = s_v2, dry_run = FALSE)
  
  message("------ mBIC check on results over lam2_vec ------")
  lam2_len <- length(res_s1$lam2_vec)
  mu_further_improve_mat <- matrix(0, nrow = 1 * lam2_len, ncol = n)
  beta_further_improve_mat <- matrix(0, nrow = 1 * lam2_len, ncol = p)
  bic_vec_s1 <- rep(0, 1 * lam2_len)
  group_num_vec_s1 <- rep(1, 1 * lam2_len)
  active_num_vec_s1 <- rep(0, 1 * lam2_len)
  max_bic <- Inf
  
  bic_res <- future.apply::future_lapply(1 : (1 * lam2_len), 
                                         FUN = RSAVS_Compute_BIC_V2, 
                                         rsavs_res = list(w_mat = res_s1$w_mat, 
                                                          mu_updated_mat = res_s1$mu_updated_mat), 
                                         y_vec = y_vec, x_mat = x_mat, 
                                         l_type = "L1", l_param = NULL, 
                                         phi = phi_s1, const_a = bic_const_a_s1, 
                                         update_mu = update_mu_s12, 
                                         double_log_lik = bic_double_loglik_s1, 
                                         from_rsi = FALSE)
  gc()
  
  for(j in 1 : (1 * lam2_len)){
    # print(j)
    tmp <- bic_res[[j]]
    mu_new <- as.vector(tmp$mu_vec)
    beta_new <- tmp$beta_vec
    mu_further_improve_mat[j, ] <- mu_new
    beta_further_improve_mat[j, ] <- beta_new
    bic_vec_s1[j] <- tmp$bic_info["bic"]
    
    # update bic according to complexsity upper bound
    current_group_num <- length(unique(mu_new))
    current_active_num <- sum(beta_new != 0)
    group_num_vec_s1[j] <- current_group_num
    active_num_vec_s1[j] <- current_active_num
    if((current_active_num + current_group_num) >= complex_upper_bound){
      bic_vec_s1[j] <- max_bic
    }
    
    # update bic according to active number
    if(current_active_num == 0){
      bic_vec_s1[j] <- max_bic
    }
  }
  
  # save the summarized information with the further improved estimation
  best_id_s1 <- which.min(bic_vec_s1)
  mu_s1 <- mu_further_improve_mat[best_id_s1, ]
  beta_s1 <- beta_further_improve_mat[best_id_s1, ]
  
  active_idx <- which(beta_s1 != 0)
  if(length(active_idx) == 0){
    active_idx <- sample(1 : p, size = 1)
  }
  if(length(active_idx) > (n / 2)){
    active_idx <- active_idx[1 : (n / 2)]
  }
  active_chosen_num <- length(active_idx)
  
  lam2_chosen <- res_s1$lam2_vec[best_id_s1]
  
  
  message("------ Analysis stage2: search over lam1_vec ------")
  message("--- a dry run to prepare some variables ---")
  dry_run_s2 <- RSAVS_Path(y_vec = y_vec, x_mat = x_mat[, active_idx, drop = FALSE],  
                           l_type = l_type, l_param = l_param, 
                           p1_type = p1_type, p1_param = p1_param, 
                           p2_type = p2_type, p2_param = p2_param, 
                           lam1_vec = lam1_vec, lam2_vec = lam2_vec, lam1_sort = FALSE, lam2_sort = FALSE, 
                           min_lam1_ratio = min_lam1_ratio, min_lam2_ratio = min_lam2_ratio, 
                           lam1_max_ncvguard = lam1_max_ncvguard, subgroup_benchmark = FALSE, 
                           lam1_len = lam1_len, lam2_len = lam2_len, 
                           const_r123 = const_r123_s2, const_abc = const_abc_s2, 
                           initial_values = initial_values_s2, 
                           cd_max_iter = cd_max_iter_s2, 
                           omp_zsw = omp_zsw, eigen_pnum = eigen_pnum, s_v2 = s_v2, 
                           dry_run = TRUE
  )
  
  if(is.null(initial_values_s2)){
    mu0 <- RSAVS_Get_Mu0(y_vec, l_type, l_param)
    initial_values_s2 <- list(beta_init = rep(0, active_chosen_num), 
                              mu_init = rep(mu0, n * (n - 1) / 2), 
                              z_init = as.vector(y_vec - mu0), 
                              s_init = rep(0, n * (n - 1) / 2), 
                              w_init = rep(0, active_chosen_num), 
                              q1_init = rep(0, n), 
                              q2_init = rep(0, n * (n - 1) / 2), 
                              q3_init = rep(0, active_chosen_num))
  }
  
  lam1_vec <- sort(dry_run_s2$lam1_vec, decreasing = TRUE)
  
  message("--- real analysis ---")
  res_s2 <- RSAVS_Path(y_vec = y_vec, x_mat = x_mat[, active_idx, drop = FALSE], 
                       l_type = l_type, l_param = l_param, 
                       p1_type = p1_type, p1_param = p1_param, 
                       p2_type = p2_type, p2_param = p2_param, 
                       lam1_vec = lam1_vec, lam2_vec = lam2_chosen, lam1_sort = FALSE, 
                       # min_lam2_ratio = min_lam2_ratio, lam2_len = lam2_len, 
                       subgroup_benchmark = subgroup_benchmark_s2, 
                       const_r123 = const_r123_s2, const_abc = const_abc_s2, 
                       initial_values = initial_values_s2, phi = phi_s2, 
                       tol = tol_s2, max_iter = max_iter_s2, 
                       cd_tol = cd_tol_s2, cd_max_iter = cd_max_iter_s2, 
                       update_mu = update_mu_s21, 
                       loss_track = loss_track_s2, diff_update = diff_update_s2, 
                       omp_zsw = omp_zsw, eigen_pnum = eigen_pnum, s_v2 = s_v2, dry_run = FALSE)
  
  message("------ mBIC check on results over lam1_vec ------")
  lam1_len <- length(lam1_vec)
  mu_further_improve_mat <- matrix(0, nrow = lam1_len * 1, ncol = n)
  beta_further_improve_mat <- matrix(0, nrow = lam1_len * 1, ncol = active_chosen_num)
  bic_vec_s2 <- rep(0, lam1_len * 1)
  group_num_vec_s2 <- rep(1, lam1_len * 1)
  active_num_vec_s2 <- rep(0, lam1_len * 1)
  max_bic <- Inf
  message("before bic")
  print(active_idx)
  bic_res <- future.apply::future_lapply(1 : (lam1_len * 1), 
                                         FUN = RSAVS_Compute_BIC_V2, 
                                         rsavs_res = list(w_mat = res_s2$w_mat, 
                                                          mu_updated_mat = res_s2$mu_updated_mat), 
                                         y_vec = y_vec, x_mat = x_mat[, active_idx, drop = FALSE], 
                                         l_type = l_type, l_param = l_param, 
                                         phi = phi_s2, const_a = bic_const_a_s2, 
                                         update_mu = update_mu_s22, 
                                         double_log_lik = bic_double_loglik_s2, 
                                         from_rsi = FALSE)
  print(active_idx)
  gc()
  print(active_idx)
  
  for(j in 1 : (lam1_len * 1)){
    # print(j)
    tmp <- bic_res[[j]]
    mu_new <- as.vector(tmp$mu_vec)
    beta_new <- tmp$beta_vec
    mu_further_improve_mat[j, ] <- mu_new
    beta_further_improve_mat[j, ] <- beta_new
    bic_vec_s2[j] <- tmp$bic_info["bic"]
    
    # update bic according to complexsity upper bound
    current_group_num <- length(unique(mu_new))
    current_active_num <- sum(beta_new != 0)
    group_num_vec_s2[j] <- current_group_num
    active_num_vec_s2[j] <- current_active_num
    if((current_active_num + current_group_num) >= complex_upper_bound){
      bic_vec_s2[j] <- max_bic
    }
    
    # update bic according to active number
    if(current_active_num == 0){
      bic_vec_s2[j] <- max_bic
    }
  }
  print(active_idx)
  # save the summarized information with the further improved estimation
  best_id_s2 <- which.min(bic_vec_s2)
  mu_s2 <- mu_further_improve_mat[best_id_s2, ]
  beta_s2 <- beta_further_improve_mat[best_id_s2, ]
  
  lam1_chosen <- res_s2$lam1_vec[best_id_s2]
  
  # combine the final result
  message("before best")
  print(active_idx)
  mu_best <- mu_s2
  beta_best <- rep(0, p)
  beta_best[active_idx] <- beta_s2
  
  res_full <- list(res_s1 = res_s1, 
                   res_s2 = res_s2, 
                   best_id_s1 = best_id_s1, 
                   best_id_s2 = best_id_s2, 
                   bic_vec_s1 = bic_vec_s1, 
                   bic_vec_s2 = bic_vec_s2, 
                   group_num_vec_s1 = group_num_vec_s1, 
                   group_num_vec_s2 = group_num_vec_s2, 
                   active_num_vec_s1 = active_num_vec_s1, 
                   active_num_vec_s2 = active_num_vec_s2, 
                   mu_best = mu_best, 
                   beta_best = beta_best
                   )
  
  return(res_full)
}