# this script is used to load all needed libraries
# This file is written to work with roxygen2.

# This package uses Cpp code.
#' @useDynLib RSAVS


# This packge needs Rcpp and RcppEigen.
#' @import Rcpp
#' @importFrom Rcpp evalCpp

#' @import RcppEigen

# Other packages needed for this package.
#' @import stats

# Package fpc for pamk. Use cluster method to better determine number of subgroups from the mu vector
#' @import fpc

# Package SparseM for generating D matrix in a sparse matrix form
#' @import SparseM

# Package quantreg for further improve result of L-1 loss
#' @import quantreg


# Package MASS for further improve result for Huber loss(rlm() function)
#' @import MASS

.onUnload <- function (libpath) {
  library.dynam.unload("RSAVS", libpath)
}

# Which functions of this package to be exported for users?



# --------------------------------------------
# output from previous devtools::check(document = F)
#  RSAVS_Determine_Mu: no visible global function definition for 'pamk'
#  RSAVS_Further_Improve: no visible global function definition for 'rq'
#  RSAVS_Further_Improve: no visible global function definition for 'lm'
#  RSAVS_Further_Improve: no visible global function definition for 'rlm'
#  RSAVS_Generate_D_Matrix: no visible global function definition for
#    'as.matrix.csr'
#  RSAVS_LargeN: no visible global function definition for 'median'
#  RSAVS_LargeN: no visible global function definition for 'rlm'
#  RSAVS_LargeN: no visible global function definition for 'ginv'
#  Undefined global functions or variables:
#    as.matrix.csr ginv lm median pamk rlm rq
#  Consider adding
#    importFrom("stats", "lm", "median")
#  to your NAMESPACE file.