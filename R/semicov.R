#' Minimum Semivariance Weights
#'
#' This function numerically computes the weights for the minimum semivariance portfolio.
#' Runs slowly with more than 3 assets. Cannot be used with more than 5 assets.
#'
#' @param mydata Matrix of returns. Each row is a period and each column is an asset.
#' @param ngrid Number of values in the grid of weights. Defaults to 201.
#' @param bench Benchmark for semivariance. Defaults to 0.
#' @param short Allowing short selling. Defaults to FALSE (no short positions allowed).
#' @return The minimum semivariance weights.
#' @export
numMinSvar <- function(mydata,ngrid=201,bench=0,short=FALSE){
  if (ncol(mydata)>5){stop("Maximum number of assets allowed is 5")}
  if ((ncol(mydata)>3)&(short==TRUE)){stop("Short selling only allowed with up to 3 assets")}
  if ((ncol(mydata)==3)&(short==TRUE)&(ngrid>201)){stop("Maximum ngrid allowed with 3 assets and short selling is 201")}
  if ((ncol(mydata)==4)&(ngrid>201)){stop("Maximum ngrid allowed with 4 assets is 201")}
  if ((ncol(mydata)==5)&(ngrid>101)){stop("Maximum ngrid allowed with 5 assets is 101")}
  if ((ncol(mydata)==5)&(nrow(mydata)>120)){stop("Longest estimation window allowed with 5 assets is 120")}
  
  comb <- seq(from=0,to=1,length.out=ngrid)
  step <- 1/(length(comb)-1)
  n_weights <- ncol(mydata)
  comb_mat <- as.matrix(t(partitions::restrictedparts(n = 1/step, m = n_weights) * step))
  perms <- gtools::permutations(n = ncol(mydata), r = ncol(mydata), v = 1:ncol(mydata))
  comb_perm <- list()
  for (i in 1:nrow(perms)){
    comb_perm[[i]] <- comb_mat[,perms[i,]]
  }
  comb_perm_mat <- do.call(rbind,comb_perm)  #matrix with all possible sets of weights

  if ((ncol(mydata)==2)&(short==TRUE)){
    comb_perm_mat1 <- cbind(comb_perm_mat[,1]-1,comb_perm_mat[,2]+1)
    comb_perm_mat2 <- cbind(comb_perm_mat[,1]+1,comb_perm_mat[,2]-1)
    comb_perm_mat <- rbind(comb_perm_mat,comb_perm_mat1,comb_perm_mat2)
    comb_perm_mat <- comb_perm_mat[comb_perm_mat[,1]>=-0.5,]
  }
  
  if ((ncol(mydata)==3)&(short==TRUE)){
    comb <- seq(from=-0.5,to=1.5,length.out=ngrid*2-1)
    comb_perm_mat <- list()
    for (x in 1:length(comb)) {
      combs <- list()
      for (y in 1:length(comb)) {
        combs[[y]] <- c(comb[x],comb[y])
      }
      combs <- do.call(rbind,combs)
      comb_perm_mat[[x]] <- combs
    }
    comb_perm_mat <- do.call(rbind,comb_perm_mat)
    comb_perm_list <- list()
    for (i in 1:length(comb)) {
      comb_i <- rep(comb[i],nrow(comb_perm_mat))
      comb_perm_list[[i]] <- cbind(comb_perm_mat,comb_i)
    }
    comb_perm_mat <- do.call(rbind,comb_perm_list)
    comb_perm_mat <- comb_perm_mat[rowSums(comb_perm_mat)==1,]  #only retain sets of weights that sum up to 1
    rm(comb_perm_list)
  }

  M_bench <- matrix(data=1, nrow=nrow(mydata)) %*% rbind(rep(bench,ncol(mydata)))  #matrix with benchmark (used to compute semicovariance matrix)
  r_combs <- comb_perm_mat%*%t(mydata)
  P_SVAR <- vector()
  for (z in 1:nrow(r_combs)) {
      w_comb <- comb_perm_mat[z,] #selected fixed portfolio weights
      r_comb <- r_combs[z,]  #returns of the portfolio with the selected weights
      pos_comb_ind <- which(r_comb>=bench)  #periods in which the portfolio has returns greater or equal to the benchmark
      est_dwn <- mydata
      est_dwn[pos_comb_ind,] <- bench  #set returns equal to the benchmark if above or equal to it

      D <- est_dwn - M_bench  #create a difference matrix
      SCN <- (nrow(mydata))^-1*t(D)%*%D  #compute semicovariance matrix with fixed weights
      P_SVAR[z] <- t(w_comb)%*%SCN%*%w_comb
    }
    p_svar_min <- which.min(P_SVAR)
    W <- comb_perm_mat[p_svar_min,]  #optimal weights that minimize portfolio semivariance
    names(W) <- NULL
    return(W)
}


#' Approximated semicovariance matrix
#'
#' This function computes an approximated semicovariance matrix as in Estrada (2008).
#'
#' @param mydata Matrix of returns. Each row is a period and each column is an asset.
#' @param bench Benchmark for semivariance. Defaults to 0.
#' @return The approximated semicovariance matrix as in Estrada (2008).
#' @export
semCov <- function(mydata,bench=0){
  D <- mydata - bench
  ZM <- matrix(data=0, nrow=nrow(mydata), ncol=ncol(mydata))
  DZ <- pmin(D,ZM)
  S <- (nrow(mydata))^(-1)*t(DZ)%*%DZ  #approximated semicovariance matrix
  return(S)
}