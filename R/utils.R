RSAVS_Get_Psi <- function(y_vec, l_type = "L1", l_param = NULL){
    # This function gets the default psi vec for the algorithm
    # Be default, we mean: 
    #     mu0 = argmin_{mu} \sum \rho(y_i - mu)
    # and 
    #     psi_vec = c(psi(y_1 - mu0), psi(y_2 - mu0), \cdots, psi(y_n - mu0))
    # where 
    #     psi = \partial \rho
    # NOTE: there is a SIGN different between this psi_vec and \partial \rho w.r.t \mu
    # NOTE: for non-differential loss like L1 at 0, see source code for its value
    # Args: y_vec: numerical vector of y
    #       l_type: type of loss function, character string
    #               "L1": median regression
    #               "L2": linear regression
    #               "Huber": Huber loss regression
    #       l_param: parameter vector for the loss function. 
    #                Not needed for "L1" and "L2" type of loss
    
    if(l_type == "L1"){
        mu0 <- median(y_vec)
        psi_vec <- sign(y_vec - mu0)
    }else{
        if(l_type == "L2"){
            mu0 <- mean(y_vec)
            psi_vec <- 2 * (y_vec - mu0)
        }else{
            if(l_type == "Huber"){
                tmp <- MASS::rlm(y_vec ~ 1, k = l_param[1])
                mu0 <- tmp$coefficients[1]
                psi_vec <- RSAVS_Huber(y_vec - mu0, param = l_param[1], derivative = T)
            }else{
                stop(paste("Unsupported type of loss function!, l_type = ", l_type, sep = ""))
            }
        }
    }
    return(psi_vec)
}

RSAVS_Get_Lam_Max <- function(y_vec, x_mat, l_type = "L1", l_param = NULL, lam_ind = 1, const_abc = rep(1, 3), eps = 10^(-6)){
    # lam_ind: 1: lam1
    #          2: lam2
    n <- length(y_vec)
    const_a <- const_abc[1]
    const_b <- const_abc[2]
    const_c <- const_abc[3]
    
    psi_vec <- RSAVS_Get_Psi(y_vec = y_vec, l_type = l_type, l_param = l_param)
    
    if(lam_ind == 1){
        lam_max <- 2 / n / const_a / const_b * max(abs(psi_vec)) * (1 + eps)
        lam_max <- lam_max * 1.1    # safe guard
    }else{
        if(lam_ind == 2){
            # `psi_vec` is in matrix form
            #   convert it to vector form
            # psi_vec <- as.vector(psi_vec)
            # tmp <- psi_vec * x_mat
            # tmp <- colSums(tmp)
            tmp <- t(x_mat) %*% psi_vec
            lam_max <- 1 / const_a / const_c * max(abs(tmp)) * (1 + eps)
        }else{
            stop(paste("Wrong index for lambda with lam_ind = ", lam_ind, sep = ""))
        }
    }
    
    return(lam_max)
}
