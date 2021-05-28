Draw_Mu_Path <- function(result, lam2_id, lam1_id, idx_vec = NULL, visual = TRUE){
  n <- ncol(result$mu_mat)
  lam1_len <- length(result$lam1_vec)
  lam2_len <- length(result$lam2_vec)
  
  if(!is.null(idx_vec)){
    muplot <- matrix(0, nrow = length(idx_vec), ncol = n)
    visual <- FALSE
  }else{
    if(!missing(lam2_id)){
      muplot <- matrix(0, nrow = lam1_len, ncol = n)
      idx_vec <- (1 : lam1_len) + lam1_len * (lam2_id - 1)
      lam_vec <- result$lam1_vec
    }else{
      if(!missing(lam1_id)){
        warning("`lam1_id` is supplied for constructing `mu_path`! Maybe you should double check your code!")
        muplot <- matrix(0, nrow = lam2_len, ncol = n)
        idx_vec <- (0 : (lam2_len - 1)) * lam1_len + lam1_id
        lam_vec <- result$lam2_vec
      }else{
        stop("idx_vec, lam2_id, lam1_id, at least one of them should be specified!")
      }
    }
  }
  
  for(i in 1 : length(idx_vec)){
    idx <- idx_vec[i]
    groups <- RSAVS_S_to_Groups(result$s_mat[idx, ], n)
    message("There are ", length(groups), " group(s).")
    for(j in 1 : length(groups)){
      group_idx <- groups[[j]]
      muplot[i, group_idx] <- mean(result$mu_mat[idx, group_idx])
    }
  }
  
  if(visual){
    ymin <- min(muplot)
    ymax <- max(muplot)
    
    plot(lam_vec, lam_vec, type = "n", ylim = c(ymin - 0.1 * (ymax - ymin), ymax + 0.1 * (ymax - ymin)))
    for(i in 1 : n){
      lines(lam_vec, muplot[, i])
      points(lam_vec, muplot[, i], cex = 0.5)
    }
  }
  return(muplot)
}

Draw_Beta_Path <- function(result, lam1_id, lam2_id, idx_vec = NULL, beta_true = NULL, visual = TRUE){
  n <- ncol(result$mu_mat)
  p <- ncol(result$beta_mat)
  lam1_len <- length(result$lam1_vec)
  lam2_len <- length(result$lam2_vec)
  
  
  if(!is.null(idx_vec)){
    # betaplot <- matrix(0, nrow = length(idx_vec), ncol = n)
    visual <- FALSE
  }else{
    if(!missing(lam1_id)){
      # betaplot <- matrix(0, nrow = lam2_len, ncol = n)
      idx_vec <- (0 : (lam2_len - 1)) * lam1_len + lam1_id
      lam_vec <- result$lam2_vec
    }else{
      if(!missing(lam2_id)){
        warning("`lam2_id` is supplied for constructing `beta_path`! Maybe you should double check your code!")
        # betaplot <- matrix(0, nrow = lam1_len, ncol = n)
        idx_vec <- (1 : lam1_len) + lam1_len * (lam2_id - 1)
        lam_vec <- result$lam1_vec
      }else{
        stop("idx_vec, lam2_id, lam1_id, at least one of them should be specified!")
      }
    }
  }
  
  betaplot <- result$w_mat[idx_vec, , drop = FALSE]
  
  active_num <- apply(betaplot, 1, function(x){
    return(sum(x != 0))
  })
  
  print("Number of active covariates:")
  print(active_num)
  
  if(visual){
    ymin <- min(betaplot)
    ymax <- max(betaplot)
    
    plot(lam_vec, lam_vec, type = "n", ylim = c(ymin - 0.1 * (ymax - ymin), ymax + 0.1 * (ymax - ymin)))
    for(i in 1 : p){
      lines(lam_vec, betaplot[, i])
      points(lam_vec, betaplot[, i], cex = 0.5)
    }
    
    # visualize the true value(`beta_true`)
    if(!is.null(beta_true)){
      idx <- which(beta_true != 0)
      for(i in 1 : length(idx)){
        lines(lam_vec, betaplot[, idx[i]], col = i + 1)
        points(lam_vec, betaplot[, idx[i]], col = i + 1, cex = 0.5)
        points(min(lam_vec), beta_true[idx[i]], col = i + 1, pch = 8)
        abline(h = beta_true[idx[i]], lty = 2, col = i + 1, lwd = 0.5)
      }
    }
  }
  
  return(betaplot)
}