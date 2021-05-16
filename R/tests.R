RSAVS_Test_Cluster_Lasso <- function(y_vec, lam, a){
  # prepare data
  y_vec <- as.vector(y_vec)
  n <- length(y_vec)
  
  # prepare variables
  response <- rep(0, n + n * (n - 1) / 2)
  response[1 : n] <- y_vec / a
  
  d_mat <- RSAVS_Generate_D_Matrix(n, dense = FALSE)
  x_mat <- SparseM::as.matrix.csr(0, nrow = n + n * (n - 1) / 2, ncol = n)
  x_mat[1 : n, 1 : n] <- diag(x = 1, nrow = n)
  x_mat[(1 : (n * (n - 1) / 2)) + n, 1 : n] <- d_mat * lam
  
  res_fit <- quantreg::rq.fit.sfn(a = x_mat, y = response, tau = 0.5)
  mu_new <- res_fit$coefficients[1 : n]
  return(mu_new)
}