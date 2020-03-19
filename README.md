# RSAVS
This package offers **R**obust **S**ubgroup **A**nalysis and **V**ariable **S**election simultaneously. It's a companion `R` package for the paper **Robust Subgroup Analysis of High Dimensional Data**.

## Installation
You can use `devtools` to directly install from github
```r
#install.packages("devtools")
devtools::install_github("fenguoerbian/RSAVS")
```

## Example
We can generate a toy sample
```r
n <- 200    # number of observations
q <- 5    # number of active covariates
p <- 50    # number of total covariates
k <- 2    # number of subgroups
# k subgroup effect, centered at 0
group_center <- seq(from = 0, to = 2 * (k - 1), by = 2) - (k - 1)
beta_true <- c(rep(1, q), rep(0, p - q))    # covariate effect vector
alpha_true <- sample(group_center, size = n, replace = T)    # subgroup effect vector

x_mat <- matrix(rnorm(n * p), nrow = n, ncol = p)    # covariate matrix
err_vec <- rnorm(n, sd = 0.5)    # error term
y_vec <- alpha_true + x_mat %*% beta_true + err_vec    # response vector
```

Then we can conduct the analysis
```r
res <- RSAVS_LargeN(y_vec = y_vec, x_mat = x_mat, lam1_length = 50, lam2_length = 40, phi = 5)
```
where `phi` is the parameter needed by mBIC. By default, the function uses `L1` as regression loss function and `SCAD` penalty for both subgroup identification and variable selection. You can use other loss function or penalty term, e.g
```r
res_huber <- RSAVS_LargeN(y_vec = y_vec, x_mat = x_mat, l_type = "Huber", l_param = 1.345, 
                          lam1_length = 50, lam2_length = 40, p1_type = "M", p2_type = "L", 
                          phi = 5)
```
uses `Huber` loss with parameter 1.345, `MCP` penalty for subgroup identification and `Lasso` penalty for variable selection. More details about setting loss type, penalty type and their parameters can be found in the package documentation.

The function uses ADMM as a solver and the resulting variable `res` is a list containing all `lam1_length` \* `lam2_length` results. And `res$best_id` corresponds to the solution with the lowest mBIC.

You can do post-selection estimation by
```r
ind <- res$best_id    # pick an id of the solution
res2 <- RSAVS_Further_Improve(y_vec = y_vec, x_mat = x_mat, 
                              mu_vec = res$mu_improve_mat[ind, ], 
                              beta_vec = res$w_mat[ind, ])
```
This function carries out ordinary low dimensional estimation(without any penalty function) given the parameter structure indicated by `mu_vec` and `beta_vec`. And it supports all the loss function that `RSAVS_LargeN` dose.