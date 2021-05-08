RSAVS_Softthresh <- function(invec, thresh){
    thresh <- abs(thresh)
    
    outvec <- abs(invec) - thresh
    id0 <- which(outvec <= 0)
    outvec[id0] <- 0
    
    outvec <- sign(invec) * outvec
    return(outvec)
}

RSAVS_UpdateZ_L1 <- function(invec, param, r1, const_a){
    outvec <- RSAVS_Softthresh(invec, 1 / r1 / const_a)
    return(outvec)
}

RSAVS_UpdateZ_L2 <- function(invec, param, r1, const_a){
    outvec <- invec / (1 + 2 / const_a / r1)
    return(outvec)
}

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

RSAVS_UpdateSW_Identity <- function(invec, param, const_r, const_bc){
    return(invec)
}

RSAVS_UpdateSW_Lasso <- function(invec, param, const_r, const_bc){
    outvec <- RSAVS_Softthresh(invec, thresh = const_bc * param[1] / const_r)
    return(outvec)
}

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