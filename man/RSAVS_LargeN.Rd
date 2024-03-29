% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pre_define_functions.r
\name{RSAVS_LargeN}
\alias{RSAVS_LargeN}
\title{Robust subgroup analysis and variable selection simultaneously for large scale dataset}
\usage{
RSAVS_LargeN(
  y_vec,
  x_mat,
  l_type = "L1",
  l_param = NULL,
  p1_type = "S",
  p1_param = c(2, 3.7),
  p2_type = "S",
  p2_param = c(2, 3.7),
  lam1_vec,
  lam2_vec,
  min_lam1_ratio = 0.03,
  min_lam2_ratio = 0.03,
  lam1_len,
  lam2_len,
  const_r123,
  const_abc = rep(1, 3),
  initial_values,
  phi = 1,
  tol = 0.001,
  max_iter = 10,
  cd_max_iter = 1,
  cd_tol = 0.001,
  subgroup_benchmark = FALSE
)
}
\arguments{
\item{y_vec}{numerical vector of response. n = length(y_vec) is the number of observations.}

\item{x_mat}{numerical matrix of covariates. Each row for one observation and 
\code{p = ncol(x_mat)} is the number of covariates.}

\item{l_type}{character string, type of loss function.
\itemize{
  \item "L1": l-1 loss(absolute value loss)
  \item "L2": l-2 loss(squared error loss)
  \item "Huber": Huber loss. Its parameter is given in l_param.
}
Default value is "L1".}

\item{l_param}{numeric vector containing necessary parameters of the corresponding loss function. 
The default value is NULL.}

\item{p1_type, p2_type}{a character indicating the penalty types for subgroup identification and variable selection.
\itemize{
  \item "S": SCAD
  \item "M": MCP
  \item "L": Lasso
}
Default values for both parameters are "S".}

\item{p1_param, p2_param}{numerical vectors providing necessary parameters for the corresponding penalties.
\itemize{
  \item For Lasso, lam = p_param[1]
  \item For SCAD and MCP, lam = p_param[1], gamma = p_param[2]
}
Default values for both parameters are \code{c(2, 3.7)}. 
Note: This function searches the whole lam1_vec * lam2_vec grid for the best solution. 
Hence the \code{lambda}s provided in these parameters will be ignored and overwritten.}

\item{lam1_vec, lam2_vec}{numerical vectors of customized lambda vectors. 
For \code{lam1_vec}, it's preferred to be in the order from small to big.}

\item{lam1_len, lam2_len}{integers, lengths of the auto-generated lambda vectors.}

\item{phi}{numerical variable. A parameter needed for mBIC.}

\item{tol}{numerical, convergence tolerance for the algorithm.}

\item{max_iter}{integer, max number of iteration during the algorithm.}

\item{subgroup_benchmark}{bool. Whether this call should be taken as a benchmark of subgroup identification. 
If \code{TRUE}, then the penalty for variable selection will be surpressed to a minimal value.}

\item{min_lam_ratio}{the ratio between the minimal and maximal lambda, equals to (minimal lambda) / (maximal lambda).
The default value is 0.03.}

\item{initial_vec}{list of vector, providing initial values for the algorithm. 
\code{mu_initial = initial_vec$mu}
and 
\code{beta_initial = initial_vec$beta}.
Note: Current the function uses auto-generated lambda vectors 
and the corresponding initial values are all just shrinked to 0.
Hence this parameter is ignored.}
}
\description{
This function carries out robust subgroup analysis and variable selection simultaneously. 
It's implemented in a parallel fashion. It supports different types of loss functions and penalties.
}
\examples{
# a toy example
# first we generate data
n <- 200    # number of observations
q <- 5    # number of active covariates
p <- 50    # number of total covariates
k <- 2    # number of subgroups

# k subgroup effect, centered at 0
group_center <- seq(from = 0, to = 2 * (k - 1), by = 2) - (k - 1)
# covariate effect vector
beta_true <- c(rep(1, q), rep(0, p - q))
# subgroup effect vector    
alpha_true <- sample(group_center, size = n, replace = TRUE)    
x_mat <- matrix(rnorm(n * p), nrow = n, ncol = p)    # covariate matrix
err_vec <- rnorm(n, sd = 0.1)    # error term
y_vec <- alpha_true + x_mat \%*\% beta_true + err_vec    # response vector

# a simple analysis using default loss and penalties
res <- RSAVS_LargeN(y_vec = y_vec, x_mat = x_mat, 
                    lam1_len = 50, lam2_len = 40, 
                    phi = 5)

# you can choose different loss and penalties
res_huber <- RSAVS_LargeN(y_vec = y_vec, x_mat = x_mat, 
                          l_type = "Huber", l_param = 1.345,
                          p1_type = "M", p2_type = "L",
                          lam1_len = 40, lam2_len = 30, 
                          phi = 5)
# you can do post-selection estimation by
ind <- res$best_ind    # pick an id of the solution
res2 <- RSAVS_Further_Improve(y_vec = y_vec, x_mat = x_mat, 
                              mu_vec = res$mu_improve_mat[ind, ], 
                              beta_vec = res$w_mat[ind, ])                           
}
