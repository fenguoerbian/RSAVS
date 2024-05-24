# Upcoming changes

## API changes

  - Make API about specifing loss function's type more consistent. See the changes of `RSAVS_Further_Improve` for more details.
  - Internally, add functions `Draw_Mu_Path` and `Draw_Beta_Path` to visualize the solution path regarding to subgroup identification and variable selection namely.
  - Add support to choose whether to use `double log likelihood` when computing the modified BIC.
  - Add a parameter `update_mu` to `RSAVS_Path` and `RSAVS_Path_PureR` so one can choose whether to update the mu vector into meaningful subgroup structure during the computation of solution plain.
  - Add a parameter `omp_zsw` to `RSAVS_Path` and `RSAVS_Solver` which controls number of threads during the update of `z`, `s` and `w` via OpenMP.
  - Add a parameter `s_v2` to `RSAVS_Path` and `RSAVS_Solver` which indicates whether to use the faster implementation to update `s_vec` and `q2_vec`.

## General changes
  - Tidy up some scirpt.
  - Add a vignette about applying the proposed method onto the `iris` dataset. It demonstrates the unsupervised nature of the proposed method.
  - Add `RSI` and `RSI_DAC` internally.
  - Add `RSAVS_Simple_Path` internally. This is a wrapper function for a simplified search over the solution plane.

# RSAVS 0.1.4

## Major changes

  - Tidy up cpp source files to make it compatible with newer version Eigen.

  - Update and add more detailed dependency info into DESCRIPTION.

# RSAVS 0.1.3

## API changes
  
  - Parameters `lam1_length` and `lam2_length` in function `RSAVS_LargeN` are renamed to `lam1_len` and `lam2_len`.
  - Possible choices for parameter `loss_type` are standardized to the same format. Currently they are `"L1"`, `"L2"` and `"Huber"`.

## Major changes

  - Add pure R implementation of RSAVS. Refer to `RSAVS_Path_PureR` and `RSAVS_Solver_PureR`for more details.
  - Add cpp implementation of RSAVS. The detailed algorithm uses the original design, which might be unstable when dealing with large scale dataset but is more flexible. Refer to `RSAVS_Path` and `RSAVS_Solver` for more details.
  - By default, `lam1_vec` is now in order from small to big.
  - The newly added `RSAVS_Path` and `RSAVS_Path_PureR` now supports progress bar via package `progressr`.

## Future plans

  - [done] Finish the vignette about the detail design of the internal ADMM algorithm
  - Clean up the code
  - [done] Add functionality for smaller dataset
  - Add RSI/LLA algorithm
  - Add a vignette about why lasso is not suitable for subgroup identificaion in this algorithm
  
# RSAVS 0.1.2

## Minor changes
  
  - Add student performance data into the package.
  - Add a vignette about using this package to analyze the student performance dataset
  - Update the website layout of this package
  
## Future plans

  - Add a vignette about the detail design of the internal ADMM algorithm
  - Clean up the code
  - Add functionality for smaller dataset
  - Add RSI/LLA algorithm


# RSAVS 0.1.1

## Initial release
   
  - Initial release of the package.
