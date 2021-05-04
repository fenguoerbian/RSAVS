#' Built-in loss functions
#'
#' These are built-in loss functions.
#'
#' @aliases RSAVS_L2 RSAVS_L1 RSAVS_Huber
#' @param x input numeric vector
#' @param param parameters needed for the function, takes the form of 
#' numeric vector. Unused for L1 and L2.
#' @param derivative logical, whether the return is the loss value or the
#' derivative.
#' @return Loss(or derivative) value at x.
#' @examples
#' RSAVS_L2(1 : 10)
#' RSAVS_L1(1 : 10)
#' RSAVS_Huber(seq(from = -3, to = 3, by = 0.1), param = 1.345)
#' RSAVS_Huber(seq(from = -3, to = 3, by = 0.1), param = 1.345, derivative = TRUE)
#' @name loss_function
NULL

#' @rdname loss_function 
#' @export
RSAVS_L2 <- function(x, param){
  # The L2 loss function
  return(x ^ 2)
}

#' @rdname loss_function 
#' @export
RSAVS_L1 <- function(x, param){
  # The L1 loss function
  return(abs(x))
}

#' @rdname loss_function 
#' @export
RSAVS_Huber <- function(x, param, derivative = FALSE){
  # The huber loss function
  # Huber(x, c) = 0.5 * x ^ 2                 if abs(x) <= c
  #               c * abs(x) - 0.5 * c ^ 2    if abs(x) > c
  huber_c <- param[1]
  if(!derivative){
    res <- x
    id <- which(abs(x) <= huber_c)
    res[id] <- 0.5 * (x[id] ^ 2)
    res[-id] <- huber_c * abs(x[-id]) - 0.5 * huber_c ^ 2
  }else{
    res <- x
    id <- which(abs(x) <= huber_c)
    res[-id] <- huber_c * sign(x[-id])
  }
  return(res)
}
