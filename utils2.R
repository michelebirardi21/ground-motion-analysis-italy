## Utilities for regularized nonparametric covariance reconstruction


library(Matrix)
library(lbfgs)
library(MASS)
library(psych)


obj_f <- function(C,Ro,Mo,rank,L,alpha,M_filter){
  
  K <- nrow(Mo)
  
  Co <- matrix(C,ncol=rank)
  temp <- Mo*(Ro-(Co%*%t(Co)))
  theta <- crossprod(t(Co))
  
  I <- diag(K)
  I_theta <- I*theta
  J <- matrix(1, ncol = K, nrow = K)
  
  # Ic <- 1-I
  Ic <- matrix(1, nrow = K, ncol = K)
  for (i in 1:K) {
    Ic[i, i] <- 0  # Diagonal
    if (i > 1)  Ic[i, i-1] <- 0  # Lower diagonal
    if (i < K)  Ic[i, i+1] <- 0  # Upper diagonal
  }
  # for (i in 1:K) {
  #   Ic[i, i] <- 1  # Main diagonal
  #   if (i > 1) Ic[i, i - 1] <- 0  # Lower diagonal
  #   if (i < K) Ic[i, i + 1] <- 0  # Upper diagonal
  #   if (i > 2) Ic[i, i - 2] <- 0  # Second lower diagonal
  #   if (i < K - 1) Ic[i, i + 2] <- 0  # Second upper diagonal
  # }
  
  M_filter_Ic <- M_filter*Ic
  M_filter_I <- M_filter*I
  
  temp2 <- M_filter_Ic*(L%*%theta+theta%*%L)
  temp3 <- M_filter_I*(L%*%I_theta%*%J)
  
  
  return(sum(temp*temp)+alpha*tr(t(temp2)%*%temp2)+alpha*tr(t(temp3)%*%temp3))
}

# gradient of the objective function
gradient_obj_f <- function(C,Ro,Mo,rank,L,alpha,M_filter){
  Co <- matrix(C,ncol=rank)
  theta <- crossprod(t(Co))
  
  K <- nrow(Mo)
  
  I <- diag(K)
  I_theta <- I*theta
  J <- matrix(1, ncol = K, nrow = K)
  
  # Ic <- 1-I
  Ic <- matrix(1, nrow = K, ncol = K)
  for (i in 1:K) {
    Ic[i, i] <- 0  # Diagonal
    if (i > 1)  Ic[i, i-1] <- 0  # Lower diagonal
    if (i < K)  Ic[i, i+1] <- 0  # Upper diagonal
  }
  # for (i in 1:K) {
  #   Ic[i, i] <- 1  # Main diagonal
  #   if (i > 1) Ic[i, i - 1] <- 0  # Lower diagonal
  #   if (i < K) Ic[i, i + 1] <- 0  # Upper diagonal
  #   if (i > 2) Ic[i, i - 2] <- 0  # Second lower diagonal
  #   if (i < K - 1) Ic[i, i + 2] <- 0  # Second upper diagonal
  # }
  
  M_filter_Ic <- M_filter*Ic
  M_filter_I <- M_filter*I
  
  T0 <- L%*%theta
  T1 <- theta%*%L
  T2 <- M_filter_Ic*(T0+T1)*M_filter_Ic
  
  T3 <- (M_filter_I*M_filter_I)*(J%*%I_theta%*%L)
  T4 <- (M_filter_I*M_filter_I)*(L%*%I_theta%*%J)
  
  temp <- -4*((Mo*(Ro-(Co%*%t(Co))))%*%Co)+
    alpha*(4*L%*%T2%*%Co+4*T2%*%L%*%Co)+
    alpha*2*(((J%*%T3%*%L)*I)%*%Co+((L%*%T4%*%J)*I)%*%Co)
  
  c(temp)
}


# covariance reconstruction
covar_reconstruction <- function(Rn,alpha,m,rk=NULL){
  
  K <- ncol(Rn)
  if(is.null(rk)){
    rk <- K
  }
  
  ## Build Laplacian
  # Penalization of second derivative --> Classic Laplacian
  L_cut <- 2*diag(K)
  diag_idxs <- (row(L_cut) == col(L_cut) + 1)
  L_cut[diag_idxs] <- -1
  diag_idxs <- (row(L_cut) == col(L_cut) - 1)
  L_cut[diag_idxs] <- -1
  L_cut[1,1] <- 1
  L_cut[K,K] <- 1
  
  # mask
  P_k = matrix(0, nrow = K, ncol = K)
  P_k[Rn!=0] = 1
  
  temp <- Rn
  old_diag <- diag(temp)
  old_diag[old_diag==0] <- NA
  new_diag <- approx(x=as.vector(which(!is.na(old_diag))), y=old_diag[as.vector(which(!is.na(old_diag)))], xout = 1:K, rule=2)
  diag(temp) <- new_diag$y
  
  temp <- interpolate_diag(mat=temp, offset=1)
  temp <- interpolate_diag(mat=temp, offset=-1)
  
  # Creation of mask entering the laplacian term
  mask_observed = matrix(0, nrow = nrow(covar), ncol = ncol(covar))
  mask_observed[covar!=0] = 1
  M_filter <- m*(1-mask_observed) + (1-m)*mask_observed
  
  # Creation of a matrix C (= gamma) of dim K x K as the initial starting value for the minimization problem
  svd_Rn <- eigen(temp*P_k)
  if(rk==1) {
    C_init <- svd_Rn$vectors*sqrt(svd_Rn$values)
  } else {
    C_init <- svd_Rn$vectors[,1:rk]%*%diag(sqrt(svd_Rn$values[1:rk]))
  }
  
  # fig1 <- plot_ly(z=C_init%*%t(C_init), type = "heatmap", coloraxis = "coloraxis")
  # fig1 <- fig1 %>% layout(title = 'Incomplete covariance',
  #                         coloraxis=list(colorscale='RdBu',cmin = -0.05, cmax = 0.15))
  # fig1
  
  # Solving the minimization problem 
  res_optim <- optim(par=c(C_init),
                     obj_f,
                     gradient_obj_f,
                     Ro=Rn,
                     Mo=P_k,
                     rank=rk,
                     L=L_cut,
                     alpha=alpha,
                     M_filter=M_filter,
                     method = "BFGS",
                     control = list(maxit = 1000))
  hat_C <- matrix(res_optim$par,ncol=rk)
  hat_R_k <- hat_C%*%t(hat_C)
  value <- res_optim$value
  
  return(list(cov = hat_R_k, value = value))
  
}


interpolate_diag <- function(mat, offset) {
  diag_idxs <- (row(mat) == col(mat) + offset)
  old_diag <- mat[diag_idxs]
  
  # Perform interpolation only if there are valid (non-NA) values
  if (sum(!is.na(old_diag)) > 1) {
    new_diag <- approx(x = as.vector(which(!is.na(old_diag))),
                       y = old_diag[as.vector(which(!is.na(old_diag)))],
                       xout = 1:(nrow(mat) - abs(offset)))
    mat[diag_idxs] <- new_diag$y
  }
  
  return(mat)
}



simulate_cov_mat_from_samples <- function(generated_samples, n_samples, q, indexes){
  
  n_samples <- dim(generated_samples)[1]
  K <- dim(generated_samples)[2]
  
  # Complete covariance estimate
  t.min <- 1
  t.max <- n_samples
  cov.mat.compl <- cov(generated_samples)
  cov.mat.compl <- cov.mat.compl*(t.max-t.min)/(t.max-t.min+1)
  
  
  # Incomplete covariance
  generated_samples[,indexes] <- NA
  
  # Incomplete covariance estimate
  t.min <- 1
  t.max <- n_samples
  cov.mat <- cov(generated_samples, use="pairwise.complete.obs")
  cov.mat <- cov.mat*(t.max-t.min)/(t.max-t.min+1)
  
  mask <- matrix(1, nrow = dim(cov.mat)[1], ncol = dim(cov.mat)[2])
  mask[is.na(cov.mat)] <- 0
  
  cov.mat.filled <- cov.mat
  
  #q <- 1
  if(q >= 0){
    old_diag <- diag(cov.mat)
    new_diag <- approx(x=as.vector(which(!is.na(old_diag))), y=old_diag[as.vector(which(!is.na(old_diag)))], xout = 1:K)
    diag(cov.mat.filled) <- new_diag$y
    if(q > 0){
      for (k in 1:q) {
        cov.mat.filled <- interpolate_diag(mat=cov.mat.filled, offset=k)
        cov.mat.filled <- interpolate_diag(mat=cov.mat.filled, offset=-k)
      }
    }
  }
  
  cov.mat <- cov.mat.filled
  
  return(list(cov=cov.mat, compl_cov=cov.mat.compl, mask=mask))
}



# funzione obiettivo per GA (1d)
f1 <- function(kappa, sigma2, l){
  
  phi = l/sqrt(2*kappa)
  pars <- c(sigma2, phi)
  # cov_stat <- cov.spatial(distance_matrix, cov.model = "matern", cov.pars = pars, kappa = kappa)
  # 
  cov_stat <- Matern(distance_matrix, range = l, phi = sigma2, smoothness = kappa)
  
  temp <- mask *(cov_stat-covar)
  sum(temp*temp)
}


# funzione obiettivo per GA (2d)
f <- function(length_scale, nu, sigma_f){
  
  # phi = l/sqrt(2*kappa)
  # pars <- c(sigma2, phi)
  # cov_stat <- cov.spatial(distance_matrix, cov.model = "matern", cov.pars = pars, kappa = kappa)
  cov_stat <- matern_covariance(distances, length_scale, nu, sigma_f)
  
  temp <- mask *(cov_stat-covar)
  sum(temp*temp)
}


# monitor function to print the progress
monitor_function <- function(obj) {
  cat("Iteration:", obj@iter, 
      "Best fitness:", max(obj@fitness), 
      "\n")
}



matern_covariance <- function(distances, length_scale = 1.0, nu = 1.5, sigma_f = 1.0) {
  # compute pairwise distances
  n <- nrow(distances)
  
  if (nu == 0.5) {
    # Exponential kernel (special case)
    K <- sigma_f^2 * exp(-distances / length_scale)
  } else if (nu == 1.5) {
    # mat?rn 3/2
    sqrt3_d_l <- sqrt(3) * distances / length_scale
    K <- sigma_f^2 * (1 + sqrt3_d_l) * exp(-sqrt3_d_l)
  } else if (nu == 2.5) {
    # mat?rn 5/2
    sqrt5_d_l <- sqrt(5) * distances / length_scale
    K <- sigma_f^2 * (1 + sqrt5_d_l + (5 * distances^2) / (3 * length_scale^2)) * exp(-sqrt5_d_l)
  } else {
    # general Mat?rn kernel using besselK
    sqrt_2nu_d_l <- sqrt(2 * nu) * distances / length_scale
    # avoid division by zero
    sqrt_2nu_d_l[sqrt_2nu_d_l == 0] <- 1e-10
    
    K <- sigma_f^2 * (2^(1-nu) / gamma(nu)) * 
      (sqrt_2nu_d_l)^nu * besselK(sqrt_2nu_d_l, nu)
    
    # handle the diagonal (distance = 0)
    diag(K) <- sigma_f^2
  }
  
  return(K)
}



# ============================================================
# Laplaciana di una griglia KxK con ordinamento a serpentina
# ============================================================


laplacian_2d <- function(n_rows, n_cols){
  # Input parameter:
  # n_rows      numero di righe della griglia del dato bi-dimensionale discretizzato
  # n_cols      numero di colonne della griglia del dato bi-dimensionale discretizzato
  
  n  <- n_rows * n_cols  # numero totale di celle
  
  # --- Funzione che mappa (riga, colonna) -> indice nodo ---
  # Righe pari  (0,2,4,...): sx -> dx  => indice = r*NC + c
  # Righe dispari (1,3,5,...): dx -> sx => indice = r*NC + (NC-1-c)
  node_index <- function(r, c, NC) {
    if (r %% 2 == 0) {
      r * NC + c
    } else {
      r * NC + (NC - 1 - c)
    }
  }
  
  # --- Costruzione della lista degli archi ---
  edges <- list()
  k <- 1
  
  for (r in 0:(n_rows - 1)) {
    for (c in 0:(n_cols - 1)) {
      i <- node_index(r, c, n_cols) + 1  # +1 per indici R (1-based)
      
      # Vicino a destra (stessa riga, colonna+1)
      if (c + 1 < n_cols) {
        j <- node_index(r, c + 1, n_cols) + 1
        edges[[k]] <- c(i, j)
        k <- k + 1
      }
      
      # Vicino in basso (riga+1, stessa colonna)
      if (r + 1 < n_rows) {
        j <- node_index(r + 1, c, n_cols) + 1
        edges[[k]] <- c(i, j)
        k <- k + 1
      }
    }
  }
  
  edge_mat <- do.call(rbind, edges)
  i_idx <- edge_mat[, 1]
  j_idx <- edge_mat[, 2]
  
  # --- Costruzione della Laplaciana sparsa ---
  # L = D - A  (D = matrice dei gradi, A = matrice di adiacenza)
  
  # Matrice di adiacenza (simmetrica)
  A <- sparseMatrix(
    i = c(i_idx, j_idx),
    j = c(j_idx, i_idx),
    x = -1,
    dims = c(n, n)
  )
  
  # Gradi (somma di ogni riga di A in valore assoluto)
  degrees <- rowSums(abs(A))  
  
  # Laplaciana
  L <- A + Diagonal(n, degrees)
  
  return(L)
}


# ============================================================
# COVARIANCE RECONSTRUCTION
# ============================================================


# covariance reconstruction
covar_reconstruction_2D <- function(Rn,alpha,m,rk=NULL){
  
  K <- ncol(Rn)
  if(is.null(rk)){
    rk <- min(K,1000)
  }
  
  ## Build Laplacian in 2D
  L_cut = laplacian_2d(n_rows=sqrt(K), n_cols=sqrt(K))
  
  # mask
  P_k = matrix(0, nrow = K, ncol = K)
  P_k[Rn!=0] = 1
  
  temp <- Rn
  old_diag <- diag(temp)
  old_diag[old_diag==0] <- NA
  new_diag <- approx(x=as.vector(which(!is.na(old_diag))), y=old_diag[as.vector(which(!is.na(old_diag)))], xout = 1:K, rule=2)
  diag(temp) <- new_diag$y
  
  temp <- interpolate_diag(mat=temp, offset=1)
  temp <- interpolate_diag(mat=temp, offset=-1)
  
  # Creation of mask entering the laplacian term
  mask_observed = matrix(0, nrow = nrow(covar), ncol = ncol(covar))
  mask_observed[covar!=0] = 1
  M_filter <- m*(1-mask_observed) + (1-m)*mask_observed
  
  # Creation of a matrix C (= gamma) of dim K x rk as the initial starting value for the minimization problem
  svd_Rn <- eigen(temp*P_k)
  if(rk==1) {
    C_init <- svd_Rn$vectors*sqrt(svd_Rn$values)
  } else {
    C_init <- svd_Rn$vectors[,1:rk]%*%diag(sqrt(svd_Rn$values[1:rk]))
  }
  
  
  # Solving the minimization problem 
  res_optim <- lbfgs(obj_f,
                     gradient_obj_f,
                     vars=c(C_init),
                     Ro=Rn,
                     Mo=P_k,
                     rank=rk,
                     L=L_cut,
                     alpha=alpha,
                     M_filter=M_filter,
                     invisible = 1,
                     max_iterations = 200)
  
  hat_C <- matrix(res_optim$par,ncol=rk)
  hat_R_k <- hat_C%*%t(hat_C)
  value <- res_optim$value
  
  return(list(cov = hat_R_k, value = value))
  
}





