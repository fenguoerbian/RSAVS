#' Soft thresh hold
#' 
#' This function performs the so-called `soft threshholding` on the input 
#' vector(\code{invec}) based on the given threshholding value(\code{thresh})
#' 
#' @param invec numerical vector, the original input value
#' @param thresh non-negative scalar, the threshholding point
#' @return \code{outvec}. For each input entry \code{invec[i]}, the output value
#'   \code{outvec[i]} is
#'   \itemize{
#'     \item \code{outvec[i] = invec[i] - thresh}: if \code{invec[i] > thresh}
#'     \item \code{outvec[i] = invec[i] + thresh}: if \code{invec[i] < -thresh}
#'     \item \code{outvec[i] = 0}: if otherwise
#'   }
#' @examples 
#' RSAVS:::RSAVS_Softthresh(seq(-2, 2, by = 0.1), thresh = 1)
RSAVS_Softthresh <- function(invec, thresh){
    thresh <- abs(thresh)
    
    outvec <- abs(invec) - thresh
    id0 <- which(outvec <= 0)
    outvec[id0] <- 0
    
    outvec <- sign(invec) * outvec
    return(outvec)
}

#' Update Z
#' 
#' Update \code{z_vec} in the ADMM algorithm
#' 
#' These functions correspond to different types of loss function defined in the
#'   original problem.
#' \itemize{
#'   \item \code{RSAVS_UpdateZ_L1} corresponds to the "L1" loss, which means \code{l_type = "L1"}.
#'   \item \code{RSAVS_UpdateZ_L2} corresponds to the "L2" loss, which means \code{l_type = "L2"}.
#'   \item \code{RSAVS_UpdateZ_Huber} corresponds to the "Huber" loss, which means \code{l_type = "Huber"}.
#' }
#' 
#' @param invec numerical vector, the orignal \code{z_vec}
#' @param param numerical vector, parameters needed for the loss function. For "L1"
#'   and "L2" loss, this will be ignored
#' @param r1 numerical scalar, parameter needed in the quadratic term in the augmented
#'   part for \code{z_vec}.
#' @param const_a numerical scalar, parameter controls the weight of regression part 
#'   in the orignal objective function
#' @name RSAVS_UpdateZ
#' @note When the loss type is "L2", it is possible to modify the algorithm to proceed 
#'   without \code{z_vec}. But currently for a more unified paradigm, this is not 
#'   implemented.
#' @note Please refer to the vignette about the algorithm detail design to find out
#'   more about how these parameters are defined.   
#' @seealso \code{\link{RSAVS_Path_PureR}}, \code{\link{RSAVS_Solver_PureR}}
NULL

#' Update S and W
#' 
#' These functions perform the updating step of \code{s_vec} and \code{w_vec} in 
#'   the ADMM algorithm for different types of penalty functions.
#'   
#' The updating steps of \code{s_vec} and \code{w_vec} share the same pattern in 
#'   the ADMM algorithm. These functions correspond to different penalty types
#'   (\code{p1_type} and \code{p2_type}) provided from the original solver function.
#'   \itemize{
#'     \item \code{UpdateSW_Identity}: identity update, this is needed when there
#'       is no shrinkage applied, possibly since the corresponding \code{lam} is 
#'       set to 0.
#'     \item \code{UpdateSW_Lasso}: Lasso update
#'     \item \code{UpdateSW_SCAD}: SCAD update
#'     \item \code{UpdateSW_MCP}: MCP update
#'   }
#' 
#' @param invec numerical vector, the orignal \code{s_vec} or \code{w_vec}.
#' @param param numerical vector, parameters needed for the penalty function.
#' @param const_r numerical scalar, parameter needed in the quadratic term in the 
#'   augmented part.
#' @param const_bc numerical scalar, parameter represents the weight of corresponding
#'   penalty part in the objective function
#' @note Please refer to the vignette about the algorithm detail design to find out
#'   more about how these parameters are defined.
#' @seealso \code{\link{RSAVS_Path_PureR}}, \code{\link{RSAVS_Solver_PureR}}
#' @name RSAVS_UpdateSW
NULL

#' @rdname RSAVS_UpdateZ
RSAVS_UpdateZ_L1 <- function(invec, param, r1, const_a){
    outvec <- RSAVS_Softthresh(invec, 1 / r1 / const_a)
    return(outvec)
}

#' @rdname RSAVS_UpdateZ
RSAVS_UpdateZ_L2 <- function(invec, param, r1, const_a){
    outvec <- invec / (1 + 2 / const_a / r1)
    return(outvec)
}

#' @rdname RSAVS_UpdateZ
RSAVS_UpdateZ_Huber <- function(invec, param, r1, const_a){
    outvec <- invec
    thresh <- param[1] / const_a / r1
    for(i in 1 : length(invec)){
        if(invec[i] > (thresh + param[1])){
            outvec[i] <- invec[i] - thresh
        }else{
            if(invec[i] < (-thresh - param[1])){
                outvec[i] <- invec[i] + thresh
            }else{
                outvec[i] <- invec[i] / (1 + 1 / const_a / r1)
            }
        }
        
    }
    return(outvec)
}

#' @rdname RSAVS_UpdateSW
RSAVS_UpdateSW_Identity <- function(invec, param, const_r, const_bc){
    return(invec)
}

#' @rdname RSAVS_UpdateSW
RSAVS_UpdateSW_Lasso <- function(invec, param, const_r, const_bc){
    outvec <- RSAVS_Softthresh(invec, thresh = const_bc * param[1] / const_r)
    return(outvec)
}

#' @rdname RSAVS_UpdateSW
RSAVS_UpdateSW_SCAD <- function(invec, param, const_r, const_bc){
    outvec <- invec
    lam <- param[1]
    gam <- param[2]
    
    for(i in 1 : length(invec)){
        if(abs(invec[i]) <= (1 + const_bc / const_r) * lam){
            outvec[i] <- RSAVS_Softthresh(invec[i], thresh = const_bc * lam / const_r)
        }else{
            if(abs(invec[i]) <= lam * gam){
                outvec[i] <- RSAVS_Softthresh(invec[i], thresh = const_bc * lam * gam / const_r / (gam - 1)) / (1 - const_bc / const_r / (gam - 1))
            }
        }
    }
    
    return(outvec)
}

#' @rdname RSAVS_UpdateSW
RSAVS_UpdateSW_MCP <- function(invec, param, const_r, const_bc){
    outvec <- invec
    lam <- param[1]
    gam <- param[2]
    
    for(i in 1 : length(invec)){
        if(abs(invec[i]) <= lam * gam){
            outvec[i] <- RSAVS_Softthresh(invec[i], thresh = const_bc * lam / const_r) / (1 - const_bc / const_r / gam)
        }
    }
    
    return(outvec)
}