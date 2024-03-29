# RSAVS
This package carries out the **R**obust **S**ubgroup **A**nalysis and **V**ariable **S**election simultaneously. It implements the computation in a parallel manner.

## Installation
You can use `devtools` to directly install the latest version from Github
```r
#install.packages("devtools")
devtools::install_github("fenguoerbian/RSAVS")
```
If you have trouble connecting to Github, this packages is also hosted on Gitlab.
```r
#install.packages("devtools")
devtools::install_gitlab("fenguoerbian/RSAVS")
```
___Note___: If you are interested in the original yet unoptimized version of this
package used in the simulation of our paper, you can check the `simulation_archive`
branch of this repo.

__Note__: It's recommended to read the vignettes on the package's website. 
But if you want to build the vignettes locally, you can run
```r
devtools::install_github("fenguoerbian/RSAVS", force = TRUE, build_vignettes = TRUE)
```
Just keep in mind that this will take some really __long__ time.

## Example
Here is a toy example:
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

Then we analyze the generated data with the function `RSAVS_LargeN`:
```r
res <- RSAVS_LargeN(y_vec = y_vec, x_mat = x_mat, lam1_len = 50, lam2_len = 40, phi = 5)
```
where `phi` is the parameter needed by mBIC. By default, this function uses `L1` as the loss function with the `SCAD` penalty for both subgroup identification and variable selection. You can use other losses or penalties, e.g
```r
res_huber <- RSAVS_LargeN(y_vec = y_vec, x_mat = x_mat, l_type = "Huber", l_param = 1.345, 
                          lam1_len = 50, lam2_len = 40, p1_type = "M", p2_type = "L", 
                          phi = 5)
```
uses `Huber` loss with parameter 1.345, `MCP` penalty for subgroup identification and `Lasso` penalty for variable selection. More details of options can be found in the package documentation.

The function uses the ADMM method to obtain the solution and the result stored in the variable `res` is a list containing all `lam1_len` \* `lam2_len` results. And `res$best_id` corresponds to the solution with the lowest mBIC.

You can do post-selection estimation by
```r
ind <- res$best_ind    # pick an id of the solution
res2 <- RSAVS_Further_Improve(y_vec = y_vec, x_mat = x_mat, 
                              mu_vec = res$mu_improve_mat[ind, ], 
                              beta_vec = res$w_mat[ind, ])
```
This function carries out ordinary low dimensional estimation(without any penalties) given the parameter structure indicated by `mu_vec` and `beta_vec`.

Besides `RSAVS_LargeN`, there are also `RSAVS_Path` and `RSAVS_Path_PureR` which can perform the same task. Furthermore, these two functions use the progress bar framework provided by package `progressr`.

```r
library(progressr)
if(as.numeric(version$major) >= 4){
    # For R >= 4.0, one can enable global handler so that
    #   there is NO need to enclose expression in `with_progress`
    handlers(global = TRUE)
}
handlers("progress")    # "progress" handlers requires the package `progress`
                        # defualt handler is `txtprogress` from base R

with_progress({
res3 <- RSAVS_Path(y_vec = y_vec, x_mat = x_mat, l_type = "Huber", l_param = 1.345, 
                   lam1_len = 50, lam2_len = 40, p1_type = "M", p2_type = "L", 
                   phi = 5)
})
```
