## Upcoming changes

# RSAVS 0.1.3

## API changes
  
  - Parameters `lam1_length` and `lam2_length` in function `RSAVS_LargeN` are renamed to `lam1_len` and `lam2_len`.
  - Possible choices for parameter `loss_type` are standardized to the same format. Currently they are `"L1"`, `"L2"` and `"Huber"`.

## Major changes

  - Add pure R implementation of RSAVS. Refer to `RSAVS_Path_PureR` and `RSAVS_Solver_PureR`.
  - Add cpp implementation of RSAVS. The detailed algorithm uses the original design, which might be unstable when dealing with large scale dataset but is more flexible. Refer to `RSAVS_Path` and `RSAVS_Solver` for more details.
  - By default, `lam1_vec` is now in order from small to big.
  - The newly added `RSAVS_Path` and `RSAVS_Path_PureR` now supports progress bar via package `progressr`.

## Future plans

  - Finish the vignette about the detail design of the internal ADMM algorithm
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
