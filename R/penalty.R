penalty_lasso <- function(x, params, derivative = FALSE){
  # Lasso penalty function
  #     p(x) = lam * abs(x)
  # Args: x: inputs
  #       params: parameter vector, lam = params[1]
  #       derivative: bool, should the derivative at x be returned 
  #                   instead of the penalty value
  # Note: Strictly speaking, there is no derivative at 0 
  #       and the subgradient at 0 is [-1, 1].
  #       but in this function, 
  #       the derivative at 0 is taken to be the limit at 0+, hence lam.
  # Return: a vector, same length as x
  x <- as.vector(x)
  res <- x
  lam <- params[1]
  if(lam < 0){
    stop("lambda(`params[1]`) must not be negative!")
  }
  if(!derivative){
    res <- lam * abs(x)
  }else{    #compute derivative
    # NOTE: derivative at 0 is taken to be the same at 0+
    id0 <- which(x == 0)
    x[id0] <- 1
    
    pos_id <- which(x > 0)
    res[pos_id] <- lam
    res[-pos_id] <- -lam
  }
  return(res)
}

penalty_scad <- function(x, params, derivative = FALSE){
  # SCAD penalty function
  #     p(x) = lam * abs(x)
  # Args: x: inputs
  #       params: parameter vector, lam = params[1], gam = params[2]
  #       derivative: bool, should the derivative at x be returned 
  #                   instead of the penalty value
  # Note: Strictly speaking, there is no derivative at 0 
  #       and the subgradient at 0 is [-lam, lam].
  #       but in this function, 
  #       the derivative at 0 is taken to be the limit at 0+, hence lam.
  # Return: a vector, same length as x
  x <- as.vector(x)
  res <- abs(x)
  lam <- params[1]
  gam <- params[2]
  if(lam < 0){
    stop("Lambda(params[1]) must not be negative!")
  }
  if(gam <= 2){
    stop("Gamma(params[2]) must be greater than 2!")
  }
  if(lam > 0){
    if(!derivative){
      id1 <- which(res <= lam)
      id2 <- which((res > lam) & (res <= lam * gam))
      id3 <- which(res >= lam * gam)
      res[id1] <- lam * res[id1]
      res[id2] <- 1 / (gam - 1) * (- res[id2] ^ 2 / 2 + gam * lam * res[id2]) - lam ^ 2 / 2 / (gam - 1)
      res[id3] <- (gam + 1) * lam ^ 2 / 2
    }else{    # compute the derivative
      # NOTE: derivative at 0 is taken to be the same at 0+
      tmp <- x
      id0 <- which(tmp == 0)
      tmp[id0] <- lam / 2
      
      id1 <- which((0 < tmp) & (tmp <= lam))
      id2 <- which((-lam <= tmp) & (tmp < 0))
      id3 <- which((lam < tmp) & (tmp <= gam * lam))
      id4 <- which((-gam * lam <= tmp) & (tmp < -lam))
      id5 <- which(abs(tmp) >= gam * lam)
      
      res[id1] <- lam
      res[id2] <- -lam
      res[id3] <- (lam * gam - x[id3]) / (gam - 1)
      res[id4] <- (-lam * gam - x[id4]) / (gam - 1)
      res[id5] <- 0
      res[id0] <- lam
    }
  }else{    
    # when lam == 0, the penalty value and its derivative are all set to be 0
    res <- rep(0, length(x))
  }
  
  return(res)
}

penalty_mcp <- function(x, params, derivative = FALSE){
  # MCP penalty function
  #     
  # Args: x: inputs
  #       params: parameter vector, lam = params[1], gam = params[2]
  #       derivative: bool, should the derivative at x be returned 
  #                   instead of the penalty value
  # Note: Strictly speaking, there is no derivative at 0 
  #       and the subgradient at 0 is [-lam, lam].
  #       but in this function, 
  #       the derivative at 0 is taken to be the limit at 0+, hence lam.
  # Return: a vector, same length as x 
  x <- as.vector(x)
  res <- abs(x)
  lam <- params[1]
  gam <- params[2]
  if(lam < 0){
    stop("Lambda(params[1]) must not be negative!")
  }
  if(gam <= 1){
    stop("Gamma(params[2]) must be greater than 1!")
  }
  if(lam > 0){
    if(!derivative){
      id1 <- which(abs(x) <= gam * lam)
      id2 <- which(abs(x) > gam * lam)
      res[id1] <- lam * res[id1] - res[id1] ^ 2 / 2 / gam
      res[id2] <- gam * lam ^ 2 / 2
    }else{    # compute the derivative
      # NOTE: derivative at 0 is taken to be the same at 0+
      tmp <- x
      id0 <- which(tmp == 0)
      
      id1 <- which((0 < tmp) & (tmp <= gam * lam))
      id2 <- which((-gam * lam <= tmp) & (tmp < 0))
      id3 <- which(abs(tmp) > gam * lam)
      
      res[id0] <- lam
      res[id1] <- lam * (1 - x[id1] / gam / lam)
      res[id2] <- -lam * (1 + x[id2] / gam / lam)
      res[id3] <- 0
    }
  }else{
    res <- rep(0, length(x))
  }

  return(res)
}
