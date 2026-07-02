#library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(viridis)
library(ggpubr)
library(rstan)
library(foreach)
library(doParallel)
library(radiant.data)
library(RGeostats)
library(plotly)
library(SpatialTools)
#options(mc.cores = parallel::detectCores(),Ncpus =  parallel::detectCores())
options(mc.cores = 6, Ncpus = 6)
rstan_options(auto_write = TRUE)


create_spatial_parameters = function(grid, funlist)
{
  temp = grid
  coordnames = names(grid)
  par_names = names(funlist)
  for ( name in par_names)
    temp = cbind(temp, funlist[[name]](grid$x_1,grid$x_2))
  temp = data.frame(temp)
  names(temp) = c(coordnames,par_names)
  return(temp)
}

plot_func_1 = function(df,name)
{
  t = ggplot(df,mapping = aes(x = x_1, y = x_2, fill = value)) +  
    geom_tile()  + coord_fixed() + theme_pubclean(base_size = 30) +
    scale_fill_viridis(option = "viridis", direction = -1, name = name) +
    theme(legend.text=element_text(size=30), legend.key.size = unit(1.8, 'cm'),
          legend.position = "top",axis.title = element_blank())
  return(t)
}

multiple_heatmaps = function(data, plot_func = plot_func_1)
{ 
  data$theta = NULL
  names = names(data)
  par_names = names[3:length(names)]
  pivoted = data %>% pivot_longer( cols = 3:length(names), names_to = "parameter" , values_to = "value")
  nested = pivoted %>% 
    group_by(parameter) %>% 
    nest() %>% 
    mutate(plots = map2(data,parameter,.f = plot_func)) 
  
  lambdas_lim = c(min(c(data$lambda_1,data$lambda_2)) , max(c(data$lambda_1,data$lambda_2)))
  Sigmas_lim = c(min(c(data$Sigma11,data$Sigma22)) , max(c(data$Sigma11,data$Sigma22)))
  
  #nested$plots[[3]] = nested$plots[[3]] + scale_fill_viridis(option = "viridis", direction = -1, name = nested$parameter[3], limits = lambdas_lim)
  #nested$plots[[4]] = nested$plots[[4]] + scale_fill_viridis(option = "viridis", direction = -1, name = nested$parameter[4], limits = lambdas_lim)
  #nested$plots[[7]] = nested$plots[[7]] + scale_fill_viridis(option = "viridis", direction = -1,name = nested$parameter[7], limits = Sigmas_lim)
  #nested$plots[[8]] = nested$plots[[8]] + scale_fill_viridis(option = "viridis", direction = -1, name = nested$parameter[8], limits = Sigmas_lim)
  
  return(gridExtra::grid.arrange(grobs = nested$plots,nrow = 1))
}

create_stan_list = function(covariate, spatial_data, nu, NN_ind, M, fit_mu = 1, fit_sigma = 1, fit_lambda = 1, fit_theta = 1)
{
  P_mu = P_sigma = P_lambda = P_theta = ncol(covariate)
  Y = spatial_data$process
  N = length(Y)
  X_mu = array(covariate,c(1,N,P_mu))
  X_sigma = array(covariate, c(1,N,P_sigma))
  X_lambda = array(covariate, c(1,N,P_lambda))
  X_theta = array(covariate, c(1,N,P_theta))
  uB_mu = array(0,c(1,P_mu))
  uB_sigma = array(0,c(1,P_sigma))
  uB_lambda = array(c(0,0),c(1,P_lambda))
  uB_lambda[1,1] = -5
  uB_theta = array(0,c(1,P_theta))
  VB_mu = array(diag(1,P_mu,P_mu),c(1,P_mu,P_mu))
  VB_sigma = array(diag(1,P_sigma,P_sigma),c(1,P_sigma,P_sigma) )
  VB_lambda = array(diag(1,P_lambda,P_lambda), c(1,P_lambda,P_lambda))
  VB_theta = array(diag(1,P_theta,P_theta), c(1, P_theta, P_theta))
  Locations = t(spatial_data[,1:2])
  
  
  data = list(
    fit_mu = fit_mu,
    fit_sigma = fit_sigma,
    fit_lambda = fit_lambda,
    fit_theta = fit_theta,
    N = N,
    P_mu = P_mu,
    P_sigma = P_sigma,
    P_lambda = P_lambda,
    P_theta = P_theta,
    nu = nu,
    Y = Y,
    X_mu = X_mu,
    X_sigma = X_sigma,
    X_lambda = X_lambda,
    X_theta = X_theta,
    Locations = Locations,
    uB_mu = uB_mu,
    uB_sigma = uB_sigma,
    uB_lambda = uB_lambda,
    uB_theta = uB_theta,
    VB_mu = VB_mu,
    VB_sigma = VB_sigma,
    VB_lambda = VB_lambda,
    VB_theta = VB_theta,
    M = M,
    NN_ind = NN_ind
  )
  return(data)
}


samples_from_chain = function(filename) {
  stanfit = rstan::read_stan_csv(filename)
  samples = data.frame(stanfit@sim$samples)
  samples %>% 
    filter(lp__ != 0) %>%  # Remove non-sampled iterations
    mutate(
      iter = 1:nrow(.),
      chain = gsub('.csv', '', last(strsplit(filename, '_')[[1]]))  # Get chain number [filename_1.csv]
    )
}


from_chains_to_pars_distr = function (chains, fitting_data)
{ 
  aux_list_mat = fitting_data[c("X_mu","X_sigma","X_lambda","X_lambda","X_theta")]
  aux_list_P = fitting_data[c("P_mu","P_sigma","P_lambda","P_lambda","P_theta")]
  names = c("mu","sigma","lambda_1","lambda_2","theta")
  to_ret = list()
  acc = 0
  for( i in 1:5)
  {
    to_ret[[names[i]]] = aux_list_mat[[1]][1,,]%*%t(chains[,(acc + 1):ifelse(i == 1, aux_list_P[[i]] , aux_list_P[[i]] + acc)])
    acc = acc + aux_list_P[[i]]
  }
  to_ret$sigma = exp(to_ret$sigma)
  to_ret$lambda_1 = exp(to_ret$lambda_1)
  to_ret$lambda_2 = exp(to_ret$lambda_2)
  to_ret$theta = pi/2 * rstanarm::invlogit(to_ret$theta)
  to_ret$theta_deg = 180/pi*to_ret$theta
  to_ret$phi = sqrt(to_ret$lambda_1/to_ret$lambda_2)
  
  return(to_ret)
}
  
from_chains_to_pars_est = function(chains,fitting_data,est_fun)
{ 
  names = c("mu","sigma","lambda_1","lambda_2","theta")
  estimates = as.numeric(apply(chains,2,est_fun))
  aux_list_mat = fitting_data[c("X_mu","X_sigma","X_lambda","X_lambda","X_theta")]
  aux_list_P = fitting_data[c("P_mu","P_sigma","P_lambda","P_lambda","P_theta")]
  to_ret = matrix(nrow = fitting_data$N,ncol = 8)
  to_ret = data.frame(to_ret)
  names(to_ret) = c(c("x_1","x_2"),names,"theta_deg")
  to_ret$x_1 = fitting_data$Locations[1,]
  to_ret$x_2 = fitting_data$Locations[2,]
  acc = 0
  for( i in 1:5)
  {
    to_ret[[names[i]]] = aux_list_mat[[1]][1,,]%*%estimates[(acc + 1):ifelse(i == 1, aux_list_P[[i]] , aux_list_P[[i]] + acc)]
    acc = acc + aux_list_P[[i]]
  }
  to_ret$sigma = exp(to_ret$sigma)
  to_ret$lambda_1 = exp(to_ret$lambda_1)
  to_ret$lambda_2 = exp(to_ret$lambda_2)
  to_ret$theta = pi/2 * rstanarm::invlogit(to_ret$theta)
  to_ret$theta_deg = 180/pi*to_ret$theta
  to_ret$phi = sqrt(to_ret$lambda_1/to_ret$lambda_2)
  return(to_ret)
}
 


compute_errors = function (data,real,error_funs, log_rescale = rep(FALSE,6)){
  data$theta = NULL
  par_names = names(data)[-c(1,2)]
  error_to_plot = data
  for ( i in 1:length(par_names))
  { 
    if(!log_rescale[i])
      error_to_plot[[par_names[i]]] = error_funs[[par_names[i]]](data[[par_names[i]]],real[[par_names[i]]])
    else
      error_to_plot[[par_names[i]]] = error_funs[[par_names[i]]](log(data[[par_names[i]]]),log(real[[par_names[i]]]))
    
  }
  return(error_to_plot)
} 

draw_pars_histograms = function (data){
  data$theta = NULL
  order = names(data)
  temp = data %>% pivot_longer(cols = 1:ncol(data),names_to = "parameter") 
  temp$parameter = factor(temp$parameter,levels = order)
  p = ggplot(temp, mapping = aes(x = value, fill = parameter)) +
      geom_histogram(aes(y = ..density..))+ facet_wrap(facets = "parameter",scale = "free") +
      theme_pubclean(base_size = 30) + rremove("legend")
  return(p)
} 


relative_err = Vectorize(function(est,real){return((est - real)/real)})

abs_err = Vectorize(function(est,real){return((est - real))})

rpd_err = Vectorize(function(est,real){return((est - real)/( abs(real) + abs(est) ))})
                         
assign_centers = function(data,centers,Ncenters,Ndata){
  distmat = matrix(0,nrow = Ndata, ncol = Ncenters)
  data = as.matrix(data)
  centers = as.matrix(centers)
  for ( i in 1:Ncenters)
    distmat[,i] = (data[,1] - centers[i,1])^2 + (data[,2] - centers[i,2])^2
  return(apply(distmat, 1, which.min))
}

extract_spatial_predictions_mean = function(array)
{ 
  N = dim(array)[1]
  npar = dim(array)[2]
  B = dim(array)[3]
  
  temp = matrix(nrow = N, ncol = npar)
  
  for(j in 1:N) {
    for(i in 1:npar) {
      if (i == 5) {
        # Trasformazione vettoriale circolare per theta_deg (dati assiali/ellissi)
        # Moltiplichiamo per 2 perché l'ellisse è simmetrica a 180 gradi
        theta_rad = array[j, 5, ] * (pi / 180)
        sin_2t = sin(2 * theta_rad)
        cos_2t = cos(2 * theta_rad)
        
        # Calcolo della media delle componenti direzionali
        mean_sin = mean(sin_2t, na.rm = TRUE)
        mean_cos = mean(cos_2t, na.rm = TRUE)
        
        # atan2 a quattro quadranti e dimezzamento per tornare all'angolo originale
        mean_theta_rad = 0.5 * atan2(mean_sin, mean_cos)
        temp[j, 5] = mean_theta_rad * (180 / pi)
      } else {
        # Media lineare standard per gli altri parametri (mu, sigma, lambda, nugget)
        temp[j, i] = mean(array[j, i, ], na.rm = TRUE)
      }
    }
  }
  
  temp = data.frame(temp)
  
  # Assegnazione dinamica dei nomi in base alla presenza del nugget
  if (npar == 5) {
    names(temp) = c("mu", "sigma", "lambda_1", "lambda_2", "theta_deg")
  } else {
    names(temp) = c("mu", "sigma", "lambda_1", "lambda_2", "theta_deg", "nugget")
  }
  
  # Aggiunta colonne calcolate
  temp$theta = (pi / 180) * temp$theta_deg
  temp$phi = sqrt(temp$lambda_1 / temp$lambda_2)
  
  return(temp)
}

extract_spatial_predictions_median = function(array, B) {
  N    = dim(array)[1]
  npar = dim(array)[2]
  B    = dim(array)[3]
  
  temp = matrix(NA_real_, nrow = N, ncol = npar)
  n_valid_per_cell = integer(N)   # diagnostica: quanti bootstrap hanno contribuito
  
  for (j in 1:N) {
    l1_all  = array[j, 3, ]
    l2_all  = array[j, 4, ]
    th_all  = array[j, 5, ] * (pi / 180)   # gradi -> radianti
    
    valid = !is.na(l1_all) & !is.na(l2_all) & !is.na(th_all)
    nv = sum(valid)
    n_valid_per_cell[j] = nv
    
    if (nv == 0) next   # resta NA
    
    l1 = l1_all[valid]; l2 = l2_all[valid]; th = th_all[valid]
    
    # Media tensoriale: media elemento per elemento delle matrici SPD
    A_bar = matrix(0, 2, 2)
    for (k in 1:nv) {
      ct = cos(th[k]); st = sin(th[k])
      R  = matrix(c(ct, st, -st, ct), 2, 2)         # rotazione
      A_bar = A_bar + R %*% diag(c(l1[k], l2[k])) %*% t(R)
    }
    A_bar = A_bar / nv
    
    eg = eigen(A_bar, symmetric = TRUE)             # autovalori decrescenti
    temp[j, 3] = eg$values[1]                        # lambda_1 (asse maggiore)
    temp[j, 4] = eg$values[2]                        # lambda_2 (asse minore)
    v = eg$vectors[, 1]                              # autovettore dell'asse maggiore
    temp[j, 5] = (atan2(v[2], v[1]) * 180 / pi) %% 180
    
    # mu, sigma, nugget: mediana sulle iterazioni valide
    temp[j, 1] = median(array[j, 1, ][valid], na.rm = TRUE)
    temp[j, 2] = median(array[j, 2, ][valid], na.rm = TRUE)
    if (npar >= 6) temp[j, 6] = median(array[j, 6, ][valid], na.rm = TRUE)
  }
  
  temp = data.frame(temp)
  if (npar == 5) {
    names(temp) = c("mu", "sigma", "lambda_1", "lambda_2", "theta_deg")
  } else {
    names(temp) = c("mu", "sigma", "lambda_1", "lambda_2", "theta_deg", "nugget")
  }
  temp$theta = (pi / 180) * temp$theta_deg
  temp$phi   = sqrt(temp$lambda_1 / temp$lambda_2)
  
  attr(temp, "n_valid_per_cell") = n_valid_per_cell  # leggilo per diagnosticare K alto
  return(temp)
}

extract_spatial_predictions_median_old = function(array, B) {
  N    = dim(array)[1]
  npar = dim(array)[2]
  B    = dim(array)[3]
  
  temp = matrix(nrow = N, ncol = npar)
  
  for (j in 1:N) {
    
    theta_rad_all = array[j, 5, ] * (pi / 180)
    l1_all        = array[j, 3, ]
    l2_all        = array[j, 4, ]
    
    # Maschera: teniamo solo le iterazioni con tutti e tre i valori validi
    valid = !is.na(theta_rad_all) & !is.na(l1_all) & !is.na(l2_all)
    
    if (sum(valid) == 0) {
      # Nessuna iterazione valida: lascia NA per questo punto
      temp[j, ] = NA
      next
    }
    
    theta_valid = theta_rad_all[valid]
    l1_valid    = l1_all[valid]
    l2_valid    = l2_all[valid]
    
    # Direzione di riferimento: mediana circolare sulle iterazioni valide
    ref_sin   = median(sin(2 * theta_valid), na.rm = TRUE)
    ref_cos   = median(cos(2 * theta_valid), na.rm = TRUE)
    ref_theta = 0.5 * atan2(ref_sin, ref_cos)
    
    # Allineamento angolare: per ogni bootstrap, decide se scambiare l1/l2
    n_valid     = sum(valid)
    l1_aligned  = numeric(n_valid)
    l2_aligned  = numeric(n_valid)
    
    for (k in 1:n_valid) {
      diff_angle = theta_valid[k] - ref_theta
      # Normalizza in [-pi/2, pi/2] (simmetria dell'ellisse a 180°, assi a 90°)
      diff_angle = ((diff_angle + pi/2) %% pi) - pi/2
      
      if (abs(diff_angle) > pi/4) {
        # Ellisse ruotata di ~90°: scambia gli assi
        l1_aligned[k] = l2_valid[k]
        l2_aligned[k] = l1_valid[k]
      } else {
        l1_aligned[k] = l1_valid[k]
        l2_aligned[k] = l2_valid[k]
      }
    }
    
    temp[j, 3] = median(l1_aligned, na.rm = TRUE)
    temp[j, 4] = median(l2_aligned, na.rm = TRUE)
    
    # mu (1) e sigma (2): mediana sulle iterazioni valide
    temp[j, 1] = median(array[j, 1, ][valid], na.rm = TRUE)
    temp[j, 2] = median(array[j, 2, ][valid], na.rm = TRUE)
    
    # theta_deg (5): mediana circolare
    sin_2t        = sin(2 * theta_valid)
    cos_2t        = cos(2 * theta_valid)
    med_theta_rad = 0.5 * atan2(median(sin_2t, na.rm = TRUE), median(cos_2t, na.rm = TRUE))
    temp[j, 5]    = med_theta_rad * (180 / pi)
    
    # nugget (6), se presente
    if (npar >= 6) {
      temp[j, 6] = median(array[j, 6, ], na.rm = TRUE)
    }
  }
  
  temp = data.frame(temp)
  if (npar == 5) {
    names(temp) = c("mu", "sigma", "lambda_1", "lambda_2", "theta_deg")
  } else {
    names(temp) = c("mu", "sigma", "lambda_1", "lambda_2", "theta_deg", "nugget")
  }
  temp$theta = (pi / 180) * temp$theta_deg
  temp$phi   = sqrt(temp$lambda_1 / temp$lambda_2)
  return(temp)
}


extract_spatial_predictions_median_old_old = function(array, B)
{ 
  N = dim(array)[1]
  npar = dim(array)[2]
  B = dim(array)[3]
  
  temp = matrix(nrow = N, ncol = npar)
  for( j in 1:N) {
    for( i in 1:npar) {
      if (i == 5) {
        # Trasformazione vettoriale circolare per theta_deg
        theta_rad = array[j, 5, ] * (pi / 180)
        sin_2t = sin(2 * theta_rad)
        cos_2t = cos(2 * theta_rad)
        
        med_sin = median(sin_2t, na.rm = TRUE)
        med_cos = median(cos_2t, na.rm = TRUE)
        
        # atan2 a quattro quadranti e dimezzamento
        med_theta_rad = 0.5 * atan2(med_sin, med_cos)
        temp[j, 5] = med_theta_rad * (180 / pi)
      } else {
        # Mediana lineare standard per gli altri parametri (mu, sigma, lambda, nugget)
        temp[j, i] = median(array[j, i, ], na.rm = TRUE)
      }
    }
  }
  
  temp = data.frame(temp)
  if (npar == 5) {
    names(temp) = c("mu", "sigma", "lambda_1", "lambda_2", "theta_deg")
  } else {
    names(temp) = c("mu", "sigma", "lambda_1", "lambda_2", "theta_deg", "nugget")
  }
  
  temp$theta = (pi / 180) * temp$theta_deg
  temp$phi = sqrt(temp$lambda_1 / temp$lambda_2)
  
  return(temp)
}


RDD_3d_plots = function(estimates, real)
{
  x_1 = unique(estimates[,1])
  x_2 = unique(estimates[,1])
  names = names(estimates)[3:length(names(estimates))]
  plist = list()
  for ( i in 1:length(names))
  { 
    z = matrix(data = real[[names[i]]], nrow = length(x_1),ncol = length(x_2),byrow = T)
    plist[[i]] = plot_ly(x = estimates[,1], y = estimates[,2], z = estimates[[names[i]]],
                         title = names[i], type = "scatter3d", size = 0.1,color = "red", 
                         showlegend = F, scene = paste0("scene",i)) %>%
                              add_surface(inherit = F, x = x_1, y = x_2, z = z,
                                          showscale=FALSE, scene = paste0("scene",i))  
  }   
  return(plist)
}



RDD_3d_plots_gp = function(parameter_names,dir_path, x_1,x_2, ord, real_pars, scaling_functions)
{
  nplots = length(parameter_names)
  plist = NULL
  for ( i in 1:nplots){
    if ( parameter_names[i] != "theta_deg")
      fit = rstan::extract(readRDS(paste0(dir_path,"/fit_",parameter_names[i],"_rdd_gp.rds")))$latent_functional
    else
      fit = rstan::extract(readRDS(paste0(dir_path,"/fit_","theta","_rdd_gp.rds")))$latent_functional
    fit = colMeans(data.frame(fit))
    fit = scaling_functions[[i]](fit)
    fit_rearrange = fit
    for ( j in 1:length(fit))
      fit_rearrange[ord[j]] = fit[j]
    z_real = matrix(data = real_pars[[parameter_names[i]]], nrow = length(x_1),ncol = length(x_2),byrow = T)
    fit_rearrange = matrix(data = fit_rearrange, nrow = length(x_1),ncol = length(x_2),byrow = T)
    plist[[i]] = plot_ly(x = x_1, y = x_2, z = fit_rearrange,scene = paste0("scene",i),
                         type = "surface",showscale = F) %>%
                               add_surface(x = x_1, y = x_2, z = z_real, scene = paste0("scene",i),
                                           color = "red",showscale = F) 
  }
  return(plist)
}

RDD_3d_plots_phi = function(dir_path, x_1,x_2, ord, real_pars)
{
  fit1 =  rstan::extract(readRDS(paste0(dir_path,"/fit_","lambda_1","_rdd_gp.rds")))$latent_functional
  fit2 =  rstan::extract(readRDS(paste0(dir_path,"/fit_","lambda_2","_rdd_gp.rds")))$latent_functional
  fit1 = exp(colMeans(fit1))
  fit2 = exp(colMeans(fit2))
  fit = sqrt(fit1/fit2)
  fit_rearrange = fit
  for ( j in 1:length(fit))
    fit_rearrange[ord[j]] = fit[j]
  z_real = sqrt(real_pars$lambda_1/real_pars$lambda_2)
  z_real = matrix(data = z_real, nrow = length(x_1),ncol = length(x_2),byrow = T)
  fit_rearrange = matrix(data = fit_rearrange, nrow = length(x_1),ncol = length(x_2),byrow = T)
  return(  plot_ly(x = x_1, y = x_2, z = fit_rearrange,scene = paste0("scene",6),type = "surface",showscale = F) %>%
             add_surface(x = x_1, y = x_2, z = z_real,scene = paste0("scene",6), color = "red",showscale = F) 
)
}

Sigma11 = function(lambda_1,lambda_2, theta){ return(sqrt(lambda_1 )* cos(theta)^2 + sqrt(lambda_2)*sin(theta)^2)}
Sigma22 = function(lambda_1,lambda_2, theta){ return(sqrt(lambda_2) * cos(theta)^2 + sqrt(lambda_1 )*sin(theta)^2)}
Sigma12 = function(lambda_1,lambda_2, theta){ return(sqrt(lambda_1 ) *sin(theta) * cos(theta) - sqrt(lambda_2)*sin(theta)*cos(theta))}

library(plotrix)

library(ggforce)
plot_ellipses = function(data, l_out = 10,scale = 1){
  seq_1 = unique(data$x_1)
  seq_2 = unique(data$x_2)
  seq_1 = seq_1[floor(seq(1,length(seq_1), length.out = l_out))]
  seq_2 = seq_2[floor(seq(1,length(seq_2), length.out = l_out))]
  data_t = data %>% dplyr::filter( x_1 %in% seq_1 , x_2 %in% seq_2)
  p = ggplot(data_t, aes(x0 = x_1, y0 = x_2, a = 3*scale*sqrt((lambda_1)), b = 3*scale*sqrt((lambda_2)),
             angle = theta_deg*pi/180)) + geom_ellipse() + geom_point(aes(x = x_1, y = x_2), size = 3) +
      theme_pubclean(base_size = 30) + coord_fixed() +
      ggtitle("Anisotropy Ellipses") +
      theme(axis.title = element_blank(), plot.title = element_text(margin = margin(60,0,0,0)))
  return(p)

}

add_ellipses = function(data, ggobj, color = "black",scale = 1){
  seq_1 = unique(data$x_1)
  seq_2 = unique(data$x_2)
  seq_1 = seq_1[floor(seq(1,length(seq_1), length.out = 6))]
  seq_2 = seq_2[floor(seq(1,length(seq_2), length.out = 6))]
  data_t = data %>% dplyr::filter( x_1 %in% seq_1 , x_2 %in% seq_2)
  p = ggobj +
      geom_ellipse(data = data_t, aes(x0 = x_1, y0 = x_2,
                   a = 3*scale*sqrt(lambda_1), b = 3*scale*sqrt(lambda_2),
                   angle = theta_deg*pi/180),inherit.aes = F, color = color, size = 0.01) +
      geom_point(data = data_t, aes(x = x_1, y = x_2), size = 3) + 
      coord_fixed() +theme(axis.title = element_blank())
  return(p)
  
}

#library(StanHeaders)
#stanFunction("modified_bessel_second_kind", v = 0.5, x = 10 )

Compute_Aniso_Old = function(lambda_1,lambda_2,theta)
{
  N = length(lambda_1)
  Anis = list()
  for (i in 1:N){
    rot = cbind(c(cos(theta[i]), -sin(theta[i])) , c(sin(theta[i]), cos(theta[i])))
    eig = cbind( c(lambda_1[i],0) , c(0,lambda_2[i]))
    Anis[[i]] = rot%*%eig%*%t(rot);
  }
  return(Anis);
}

Compute_Aniso = function(lambda_1, lambda_2, theta)
{
  N = length(lambda_1)
  Anis = list()
  for (i in 1:N){
    # Convenzione RGeostats (Antioraria da asse X)
    rot = matrix(c(cos(theta[i]), sin(theta[i]), 
                   -sin(theta[i]), cos(theta[i])), 
                 nrow=2, ncol=2) # Riempe per colonne [,1] e [,2]
    
    eig = matrix(c(lambda_1[i], 0, 
                   0, lambda_2[i]), 
                 nrow=2, ncol=2)
    
    Anis[[i]] = rot %*% eig %*% t(rot)
  }
  return(Anis)
}

matern_ns_corr = function(Locations, Aniso, dets, nu, sigma)
{ 
  N = ncol(Locations)
  sqrsqrdets = sqrt(sqrt(dets))
  norm_const = 2^(1-nu)/gamma(nu)
  C = matrix( nrow = N, ncol = N)
  cat(sprintf("\n[matern_ns_corr] Avvio calcolo matrice %d x %d...\n", N, N))
  for ( i in 1:N){
    
    if (i %% 100 == 0) {
      percentuale <- (i / N) * 100
      cat(sprintf("   -> Elaborata riga %d di %d (%.1f%% completato)\n", i, N, percentuale))
    }
    
    C[i,i] = sigma[i]^2;
    
    # La guardia logica evita la sequenza decrescente fatale in R
    if (i < N) {
      for ( j in (i+1):N)
      { 
        metric = 0.5 * (Aniso[[i]] + Aniso[[j]])
        normdet = sqrt(1/(metric[1,1]*metric[2,2] - (metric[1,2])^2));
        Q = sqrt( (Locations[,i] - Locations[,j]) %*% solve(metric, Locations[,i] - Locations[,j]));
        
        C[i,j] = sigma[i]*sigma[j]*sqrsqrdets[i]*sqrsqrdets[j] * normdet * norm_const * Q^nu * modified_bessel_second_kind(nu,Q);
        C[j,i] = C[i,j];
      }
    }
  }
  
  cat("[matern_ns_corr] Calcolo completato con successo!\n")
  return(C);
}

trim = Vectorize( function(x,lwr = -Inf, upr = +Inf)
{ if ( x < lwr) return(lwr)
  if(x > upr) return(upr)
  return(x)})

generate_centers = function (coords,center_grid, K , B, threshold = 10)
{
  N_grid = nrow(center_grid)
  N_obs = nrow(coords)
  center_list = list()
  for ( i in 1:B){
    centers = center_grid[sample(1:N_grid, size = K),]
    assignment = assign_centers(coords,centers,K,N_obs)
  while(min(table(assignment)) < threshold){
    #print(sum(table(assignment)))
    centers = center_grid[sample(1:N_grid, size = K),] 
    assignment = assign_centers(coords,centers,K,N_obs)
  }
    center_list[[i]] = centers
  }
  return(center_list)
}

RDD_fit_sequential = function (data, coords, center_grid,
                               K, B, clusters = 12, dir_path, struct = c("K-Bessel","Nugget Effect"),
                               dirs = seq(0,150, by = 30), param, lower, upper, threshold)
{ 
  # --- 1. SETUP ---
  # Rimuovo makePSOCKcluster e registerDoParallel
  
  # Assegnazioni originali
  assign_centers = assign_centers
  get_param = get_param
  
  # Importante: forziamo K numerico per evitare l'errore "matrice non-numerica"
  # che hai visto nei log precedenti, nel caso arrivi come testo.
  K = as.numeric(K) 
  
  N_obs = length(data)
  N_grid = nrow(center_grid)
  spatial_data = cbind(coords, data)
  
  # --- 2. GENERAZIONE CENTRI (Fuori dal ciclo, come in originale) ---
  # Se l'errore "matrix extent" persiste, è probabile che accada QUI dentro
  center_list = generate_centers(coords, center_grid, K, B, threshold)
  
  writeLines(c(""), paste0(dir_path,"/log.txt"))
  
  # Inizializzo la lista per i risultati
  ret = list()
  
  # --- 3. CICLO SEQUENZIALE (Sostituisce foreach) ---
  print("Avvio ciclo sequenziale...")
  
  for (i in 1:B) {
    
    # Visualizzazione progresso in console
    cat(paste("Iterazione:", i, "/", B, "\n"))
    
    # Logica interna originale
    centers = center_list[[i]]
    assignment_obs = assign_centers(coords, centers, K, N_obs)
    
    sink(paste0(dir_path,"/log.txt"), append = T)
    
    if(length(struct) == 1) 
      npar = 5 
    else 
      npar = 6 
    
    mat_param = t(matrix(NA, K, npar))
    
    for (j in 1:K)
    { 
      idxs_obs = which(assignment_obs == j)
      
      # Controllo di sicurezza: se il cluster è vuoto, saltiamo per evitare errori in db.create
      if (length(idxs_obs) == 0) next
      
      db = db.create(spatial_data[idxs_obs,], flag.grid = F, ndim = 2, autoname = F)
      
      # Nota: dirs = dirs è ridondante ma lo lascio per fedeltà
      vario = vario.calc(db, dir = dirs)
      
      # Blocco try() opzionale: se model.auto fallisce, non blocca tutto il ciclo
      fit = try(model.auto(vario, struct = struct,
                           param = param, lower = lower,
                           upper = upper, verbose = 0, flag.noreduce = T, flag.goulard = T, maxiter = 10000), silent=TRUE)
      
      if (!inherits(fit, "try-error")) {
        mat_param[2:npar, j] = unlist(get_param(fit))
        mat_param[1, j] = mean(spatial_data[idxs_obs, 3])
      }
    }
    
    print(i) # Stampa nel file di log
    sink()   # Chiude il log per questa iterazione
    
    # Salvo il risultato nella lista
    ret[[i]] = t(mat_param)
  }
  
  # Rimuovo stopCluster
  return(list(estimates = ret, center_list = center_list))
}

RDD_fit = function (data, coords, center_grid,
                    K, B, clusters = 12, dir_path, struct = c("K-Bessel","Nugget Effect"),
                    dirs = seq(0,150, by = 30), param, lower, upper, threshold)
{ 
  # --- 1. SETUP INIZIALE ---
  library(doParallel)
  library(foreach)
  
  # Forziamo K numerico (fix dalla versione sequenziale)
  K = as.numeric(K) 
  
  # Assicuriamoci che le funzioni helper siano visibili
  assign_centers_local = assign_centers
  get_param_local = get_param
  
  # Preparazione dati
  N_obs = length(data)
  N_grid = nrow(center_grid)
  spatial_data = cbind(coords, data)
  
  # Generazione centri (Fatta una volta sola sul master, come nella sequenziale)
  # Questo evita problemi di "matrix extent" all'interno dei worker
  center_list = generate_centers(coords, center_grid, K, B, threshold)
  
  # Reset del file di log (solo intestazione, niente scrittura parallela)
  writeLines(c("Start Parallel Process"), paste0(dir_path, "/log.txt"))
  
  # --- 2. CONFIGURAZIONE CLUSTER ---
  # Rileva core o usa il numero passato
  num_cores <- min(clusters, detectCores() - 1)
  cl <- makePSOCKcluster(num_cores)
  registerDoParallel(cl)
  
  print(paste("Avvio calcolo parallelo su", num_cores, "core. Attendere..."))
  
  # --- 3. CICLO PARALLELO ---
  # Nota: Rimosso 'sink' interno perché causa conflitti di scrittura (Race Conditions)
  # Aggiunto .export per sicurezza sulle funzioni custom
  ret = foreach(i = 1:B, 
                .packages = c("RGeostats", "radiant.data", "seqinr"),
                .export = c("assign_centers", "get_param"), 
                .errorhandling = "pass") %dopar% {
                  
                  # Recupera i centri per questa iterazione
                  centers = center_list[[i]]
                  
                  # Assegnazione punti ai centri
                  assignment_obs = assign_centers_local(coords, centers, K, N_obs)
                  
                  # Definizione numero parametri
                  if(length(struct) == 1) npar = 5 else npar = 6 
                  
                  mat_param = t(matrix(NA, K, npar))
                  
                  cluster_counts <- numeric(K)
                  
                  # Ciclo interno sui K cluster
                  for (j in 1:K) { 
                    idxs_obs = which(assignment_obs == j)
                    
                    # FIX 1: Controllo cluster vuoti (fondamentale per evitare crash)
                    n_punti_cluster <- length(idxs_obs)
                    cluster_counts[j] <- n_punti_cluster
                    
                    # Se ci sono meno di 15 stazioni, saltiamo per impossibilità matematica
                    if (n_punti_cluster < 15) next
                    
                    # Creazione DB RGeostats
                    db = db.create(spatial_data[idxs_obs,], flag.grid = F, ndim = 2, autoname = F)
                    
                    # Calcola estensione del cluster per adattare il variogramma
                    cluster_coords <- coords[idxs_obs, ]
                    max_dist <- max(dist(cluster_coords))
                    
                    # Sicurezza: Se max_dist è zero (o quasi), forza un valore minimo 
                    # per non far crashare l'algoritmo (es. 1 km)
                    if (max_dist < 1) max_dist <- 1
                    
                    # Usa al massimo metà della distanza massima del cluster come cutoff
                    cutoff_adattivo <- max_dist / 2
                    lag_adattivo  <- cutoff_adattivo / 10   # 10 lags
                    nlag_adattivo <- 10
                    
                    # FIX 2: Fallback dinamico sulle direzioni basato sulla densità del cluster
                    if (n_punti_cluster < 50) {
                      # Fallback: Variogramma Omnidirezionale (più robusto per pochi dati)
                      vario = vario.calc(db, lag = lag_adattivo, nlag = nlag_adattivo) 
                    } else {
                      # Variogramma Direzionale standard
                      vario = vario.calc(db, dir = dirs, lag = lag_adattivo, nlag = nlag_adattivo)
                    }
                    
                    # FIX 2: try() per gestire fallimenti di model.auto senza fermare il worker
                    fit = try(model.auto(vario, struct = struct,
                                         param = param, lower = lower,
                                         upper = upper, verbose = 0, 
                                         flag.noreduce = T, flag.goulard = T, maxiter = 10000), 
                              silent = TRUE)
                    
                    # Se il fit ha avuto successo, salva i parametri
                    if (!inherits(fit, "try-error")) {
                      mat_param[2:npar, j] = unlist(get_param_local(fit))
                      mat_param[1, j] = mean(spatial_data[idxs_obs, 3])
                    }
                  }
                  
                  # Ritorna la matrice trasposta per questa iterazione B
                  list(
                    params = t(mat_param),
                    counts = cluster_counts
                  )
                }
  
  # --- 4. CHIUSURA ---
  stopCluster(cl)
  print("Calcolo parallelo completato.")
  
  # --- 4. ESTRAZIONE DATI DAL FOREACH E CALCOLO STATISTICHE ---
  # Il foreach ora restituisce una lista di liste. Dobbiamo separare i parametri dai conteggi.
  
  final_estimates <- list()
  all_counts <- numeric() # Raccoglie tutti i conteggi di tutte le iterazioni per le statistiche globali
  
  for(i in 1:length(ret)) {
    if(!is.null(ret[[i]]$params)) {
      final_estimates[[i]] <- ret[[i]]$params
      all_counts <- c(all_counts, ret[[i]]$counts)
    } else {
      # Se il foreach ha restituito un errore puro per questa iterazione (molto raro con try)
      final_estimates[[i]] <- t(matrix(NA, K, ifelse(length(struct)==1, 5, 6)))
    }
  }
  
  # --- 5. STAMPA STATISTICHE IN CONSOLE ---
  cat("\n=========================================\n")
  cat(" STATISTICHE DISTRIBUZIONE STAZIONI NEI CLUSTER\n")
  cat(" (Media su", length(final_estimates), "iterazioni di Bootstrap)\n")
  cat("=========================================\n")
  
  # Quanti cluster in totale (K * B) erano troppo piccoli per essere calcolati? (< 5 stazioni)
  cluster_vuoti <- sum(all_counts < 5)
  perc_vuoti <- round((cluster_vuoti / length(all_counts)) * 100, 2)
  
  cat(sprintf("Totale cluster analizzati (K * B) : %d\n", length(all_counts)))
  cat(sprintf("Cluster 'vuoti' (< 5 stazioni)    : %d (%.2f%%)\n", cluster_vuoti, perc_vuoti))
  cat("Distribuzione del numero di stazioni per cluster:\n")
  print(summary(all_counts))
  cat("=========================================\n\n")
  
  # Restituiamo tutto: stime, centri, e i conteggi se volessi farci un grafico dopo
  return(list(
    estimates = final_estimates, 
    center_list = center_list,
    cluster_counts = all_counts
  ))
}

RDD_predict = function (RDD_output, grid)
{ 
  B = length(RDD_output$center_list)
  K = nrow(RDD_output$center_list[[1]])
  ret = list()
  for ( i in 1:B)
  {
    assign = assign_centers(grid,RDD_output$center_list[[i]],K,nrow(grid))
    ret[[i]] = RDD_output$estimates[[i]][assign,]
  }
  return(ret)
}
  
get_param = function(model, cost = 3.0) { 
  ret = list()
  ret$sigma = sqrt(model@basics[[1]]$sill)
  
  base_range = model@basics[[1]]$range
  coeffs = model@basics[[1]]$aniso.coeffs
  if (is.null(coeffs)) coeffs = c(1, 1)
  
  range_x = base_range * coeffs[1]
  range_y = base_range * coeffs[2]
  
  rotmat = model@basics[[1]]$aniso.rotmat
  if(!is.matrix(rotmat) || nrow(rotmat) == 1) {
    theta_rad = 0
  } else {
    theta_rad = atan2(rotmat[2,1], rotmat[1,1])
  }
  
  if (range_x >= range_y) {
    r1 = range_x; r2 = range_y
    final_theta_rad = theta_rad
  } else {
    r1 = range_y; r2 = range_x
    final_theta_rad = theta_rad + (pi / 2) 
  }
  
  # Hard cap coerente con upper bound * max aniso ratio atteso
  max_R = 400
  r1 = min(r1, max_R)
  r2 = min(r2, max_R)
  
  ret$lambda_1 = (r1 / cost)^2
  ret$lambda_2 = (r2 / cost)^2
  
  theta_deg = (final_theta_rad / pi * 180) %% 180
  if (theta_deg > 90) theta_deg = theta_deg - 180
  ret$theta_deg = theta_deg
  
  if(length(model@basics) >= 2) {
    raw_nugget = model@basics[[2]]$sill[1,1]
    # Hard cap: il nugget non può superare la varianza totale osservata
    ret$nugget = min(raw_nugget, 0.12)
  } else {
    ret$nugget = 0 
  }
  
  return(ret)
}

get_param_old = function(model, cost = 3.0) { 
  ret = list()
  ret$sigma = sqrt(model@basics[[1]]$sill)
  
  base_range = model@basics[[1]]$range
  coeffs = model@basics[[1]]$aniso.coeffs
  if (is.null(coeffs)) coeffs = c(1, 1) # Fallback sicuro
  
  # 1. Calcolo Diretto dei Veri Range
  range_x = base_range * coeffs[1]
  range_y = base_range * coeffs[2]
  
  # 2. Estrazione dell'angolo originale in radianti
  rotmat = model@basics[[1]]$aniso.rotmat
  if(!is.matrix(rotmat) || nrow(rotmat) == 1) {
    theta_rad = 0
  } else {
    theta_rad = atan2(rotmat[2,1], rotmat[1,1])
  }
  
  # 3. Ordinamento e Rotazione: Vogliamo SEMPRE r1 >= r2 (lambda_1 >= lambda_2)
  if (range_x >= range_y) {
    r1 = range_x
    r2 = range_y
    final_theta_rad = theta_rad
  } else {
    r1 = range_y
    r2 = range_x
    # Se invertiamo gli assi maggiore e minore, DOBBIAMO ruotare l'ellisse di 90 gradi
    final_theta_rad = theta_rad + (pi / 2) 
  }
  
  # 4. Hard Cap Anti-Esplosione (400 km)
  max_R = 400
  r1 = min(r1, max_R)
  r2 = min(r2, max_R)
  
  # 5. Assegnazione Parametri
  ret$lambda_1 = (r1 / cost)^2
  ret$lambda_2 = (r2 / cost)^2
  
  # 6. Conversione e normalizzazione angolo (tra -90° e 90°)
  theta_deg = (final_theta_rad / pi * 180) %% 180
  if (theta_deg > 90) theta_deg = theta_deg - 180
  ret$theta_deg = theta_deg
  
  # 7. Nugget
  if(length(model@basics) >= 2) {
    ret$nugget = model@basics[[2]]$sill[1,1]
  } else {
    ret$nugget = 0 
  }
  
  return(ret)
}

get_param_old_old = function(model, cost = 3.0)
{ 
  ret = list()
  
  # 1. Estrazione del Sill e del Range base
  ret$sigma = sqrt(model@basics[[1]]$sill)
  base_range = model@basics[[1]]$range
  
  # 2. Gestione Anisotropia e Autovalori (lambda)
  coeffs = model@basics[[1]]$aniso.coeffs
  if (is.null(coeffs)) coeffs = c(1, 1) # Fallback sicuro
  
  whichmin = which.min(coeffs)
  ratio = min(coeffs)
  
  if(whichmin == 1) {
    ret$lambda_1 = (base_range * ratio / cost)^2
    ret$lambda_2 = (base_range / cost)^2
  } else {
    ret$lambda_2 = (base_range * ratio / cost)^2
    ret$lambda_1 = (base_range / cost)^2
  }
  
  # 3. Estrazione SICURA della rotazione (Anti-Crash e Quadranti Corretti)
  rotmat = model@basics[[1]]$aniso.rotmat
  
  if(!is.matrix(rotmat) || nrow(rotmat) == 1) {
    # Modello Isotropico: Nessuna rotazione
    ret$theta_deg = 0
  } else {
    # Modello Anisotropico: Estrazione con atan2 (preserva i quadranti completi)
    # rotmat[1,1] è il coseno, rotmat[2,1] è il seno
    angolo_rad = atan2(rotmat[2,1], rotmat[1,1])
    ret$theta_deg = angolo_rad / pi * 180
  }
  
  # 4. Estrazione del Nugget (se presente nel modello)
  if(length(model@basics) >= 2) {
    ret$nugget = model@basics[[2]]$sill[1,1]
  } else {
    ret$nugget = 0 # Evita che ritorni NULL se il nugget non è fittato
  }
  
  return(ret)
}

get_param_old_old_old = function(model, cost = 4.9)
{ 
#  print(model)
  ret = list()
  ret$sigma = sqrt(model@basics[[1]]$sill)
  ret$lambda_1 = 0
  ret$lambda_2 = 0
  ret$theta_deg = 0
  whichmin = which.min(c(model@basics[[1]]$aniso.coeffs[1],model@basics[[1]]$aniso.coeffs[2]))
  ratio = min(c(model@basics[[1]]$aniso.coeffs[1],model@basics[[1]]$aniso.coeffs[2]))
  if(whichmin == 1){
    ret$lambda_1 = (model@basics[[1]]$range * ratio / cost)^2
    ret$lambda_2 = (model@basics[[1]]$range / cost)^2
    
    if(model@basics[[1]]$aniso.rotmat[2,1] < 0)
    {
        t = ret$lambda_1 
        ret$lambda_1 = ret$lambda_2
        ret$lambda_2 = t
    }

    ret$theta_deg = acos(model@basics[[1]]$aniso.rotmat[1,1])/pi*180
    if(ret$theta_deg > 90){
         ret$theta_deg = ret$theta_deg - 90
         t = ret$lambda_1 
         ret$lambda_1 = ret$lambda_2
         ret$lambda_2 = t
     }
  }
  
  if(whichmin == 2){
    ret$lambda_2 = (model@basics[[1]]$range * ratio / cost)^2
    ret$lambda_1 = (model@basics[[1]]$range / cost)^2
    if(model@basics[[1]]$aniso.rotmat[2,1] <0)
    {
      t = ret$lambda_1 
      ret$lambda_1 = ret$lambda_2
      ret$lambda_2 = t
    }
    ret$theta_deg = acos(model@basics[[1]]$aniso.rotmat[1,1])/pi*180
     if(ret$theta_deg > 90){
       ret$theta_deg = ret$theta_deg - 90
       t = ret$lambda_1 
       ret$lambda_1 = ret$lambda_2
       ret$lambda_2 = t
     }
  }
  if(length(model@basics) == 2)
    ret$nugget = model@basics[[2]]$sill[1,1]
#print(ret)
  return(ret)
}  
  

smooth_NSconvo = function(model, coords)
{
  weight_mat = as.matrix(proxy::dist(scale(coords),scale(model$mc.locations))^2)
  weight_mat = exp(-weight_mat/2*model$lambda.w)
  for ( i in 1:nrow(weight_mat)) weight_mat[i,] = weight_mat[i,]/sum(weight_mat[i,])
  weight_mat = as.matrix(weight_mat)
  knot_pars = model$MLEs.save
  knot_pars$n = NULL
  knot_pars$kappa = NULL
  knot_pars$theta = knot_pars$eta
  knot_pars$eta = NULL
  knot_pars$mu = model$beta.GLS[,1]
  knot_pars = as.matrix(knot_pars[,c(6,4,1,2,3,5)])
  par_est = weight_mat%*%knot_pars
  return(
    data.frame( 
    x_1 = coords[,1],
    x_2 = coords[,2],
    mu = par_est[,1],
    sigma = sqrt(par_est[,2]),
    lambda_1 = log(par_est[,3]),
    lambda_2 = log(par_est[,4]),
    nugget = par_est[,5],
    theta_deg = 180/pi * par_est[,6]
  )
  )
}

smooth_ellipses_locally <- function(grid_data, ell_data, raggio_km = 20) {
  ell_smoothed <- ell_data
  
  for (i in 1:nrow(ell_data)) {
    distanze_sq <- (grid_data$X_km - ell_data$X_km[i])^2 + 
      (grid_data$Y_km - ell_data$Y_km[i])^2
    
    vicini_idx <- which(distanze_sq <= raggio_km^2)
    vicini <- grid_data[vicini_idx, ]
    
    # Filtro di sicurezza fondamentale: scarta le righe con parametri NA
    vicini <- vicini[!is.na(vicini$lambda_1) & !is.na(vicini$lambda_2) & !is.na(vicini$theta_deg), ]
    peso_totale <- nrow(vicini)
    
    if (peso_totale > 3) {
      
      A_bar <- matrix(0, nrow = 2, ncol = 2)
      
      for (k in 1:peso_totale) {
        l1 <- vicini$lambda_1[k]
        l2 <- vicini$lambda_2[k]
        th_rad <- vicini$theta_deg[k] * (pi / 180)
        
        ct <- cos(th_rad)
        st <- sin(th_rad)
        R <- matrix(c(ct, st, -st, ct), 2, 2)
        
        A_bar <- A_bar + (R %*% diag(c(l1, l2)) %*% t(R))
      }
      
      A_bar <- A_bar / peso_totale
      eg <- eigen(A_bar, symmetric = TRUE)
      
      ell_smoothed$lambda_1[i] <- eg$values[1]
      ell_smoothed$lambda_2[i] <- eg$values[2]
      
      v <- eg$vectors[, 1]
      ell_smoothed$theta_deg[i] <- (atan2(v[2], v[1]) * 180 / pi) %% 180
    }
  }
  
  return(ell_smoothed)
}

smooth_ellipses_locally_old <- function(grid_data, ell_data, raggio_km = 20) {
  # Prepariamo un dataframe vuoto per i risultati
  ell_smoothed <- ell_data
  
  for (i in 1:nrow(ell_data)) {
    # 1. Trovo tutti i punti della griglia ORIGINALE che cadono nel raggio
    distanze_sq <- (grid_data$X_km - ell_data$X_km[i])^2 + 
      (grid_data$Y_km - ell_data$Y_km[i])^2
    
    vicini_idx <- which(distanze_sq <= raggio_km^2)
    vicini <- grid_data[vicini_idx, ]
    
    # Se ci sono abbastanza vicini, applico lo smoothing
    if (nrow(vicini) > 3) {
      
      # Mediana robusta per gli autovalori (scarta i K esplosi)
      ell_smoothed$lambda_1[i] <- median(vicini$lambda_1, na.rm = TRUE)
      ell_smoothed$lambda_2[i] <- median(vicini$lambda_2, na.rm = TRUE)
      
      # Mediana circolare assiale per l'angolo theta (la stessa che usi tu in Utilities)
      theta_rad <- vicini$theta_deg * (pi / 180)
      sin_2t <- sin(2 * theta_rad)
      cos_2t <- cos(2 * theta_rad)
      
      med_sin <- median(sin_2t, na.rm = TRUE)
      med_cos <- median(cos_2t, na.rm = TRUE)
      
      med_theta_rad <- 0.5 * atan2(med_sin, med_cos)
      ell_smoothed$theta_deg[i] <- med_theta_rad * (180 / pi)
    }
  }
  return(ell_smoothed)
}

plot_parameters = function(data, dir_path, K, borders_sf = NULL, n_grid = 10) {
  
  require(sf)
  require(rnaturalearth)
  require(dplyr)
  require(gridExtra)
  require(ggpubr)
  require(elevatr) 
  require(terra)  
  
  # ==============================================================================
  # 1. CONFINI ITALIA SENZA SARDEGNA + MASCHERATURA DATI
  # ==============================================================================
  ita_regioni <- ne_states(country = "Italy", returnclass = "sf")
  
  ita_no_sardegna <- ita_regioni %>%
    filter(!grepl("Sardegn", region, ignore.case = TRUE))
  
  italia_utm <- st_transform(st_union(ita_no_sardegna), 32633)
  
  borders_sf <- italia_utm
  
  df_temp <- data
  df_temp$X_m <- df_temp$X_km * 1000
  df_temp$Y_m <- df_temp$Y_km * 1000
  
  sf_all <- st_as_sf(df_temp, coords = c("X_m", "Y_m"), crs = 32633, remove = FALSE)
  
  inside <- st_intersects(sf_all, italia_utm, sparse = FALSE)
  data <- data[as.vector(inside), ] 
  
  # ==============================================================================
  # 2. HELPER PER TITOLI CON LETTERE GRECHE E VARIANZE
  # ==============================================================================
  get_plot_title <- function(par_name) {
    if (par_name == "mu") return(expression(mu))
    if (par_name == "sigma2") return(expression(sigma^2)) 
    if (par_name == "lambda_1") return(expression(lambda[1]))
    if (par_name == "lambda_2") return(expression(lambda[2]))
    if (par_name == "theta_azimuth") return(expression(theta~"(°)"))
    if (par_name == "nugget") return(expression(sigma[nugget]^2)) 
    return(par_name)
  }
  
  # --- PREPARAZIONE DATI ---
  if("sigma" %in% names(data)) {
    data$sigma2 <- data$sigma^2
  }
  
  # CONVERSIONE AL VOLO DA TRIGONOMETRICO AD AZIMUTALE (SOLO PER IL GRAFICO)
  if("theta_deg" %in% names(data)) {
    data$theta_azimuth <- (90 - data$theta_deg) %% 180
  }
  
  param_names = c("mu", "sigma2", "lambda_1", "lambda_2", "theta_azimuth", "nugget")
  plot_list = list()
  
  # --- 3. LIMITI ZOOM (METRI) ---
  limit_x = range(data$X_km * 1000, na.rm = TRUE)
  limit_y = range(data$Y_km * 1000, na.rm = TRUE)
  
  buffer_x = diff(limit_x) * 0.05
  buffer_y = diff(limit_y) * 0.05
  
  final_xlim = c(limit_x[1] - buffer_x, limit_x[2] + buffer_x)
  final_ylim = c(limit_y[1] - buffer_y, limit_y[2] + buffer_y)
  
  # --- 4. HEATMAPS ---
  for (par in param_names) {
    if (!par %in% names(data)) {
      p = ggplot() + theme_void() + ggtitle(paste("Missing:", par))
    } else {
      p = ggplot(data, aes(x = .data[["X_km"]] * 1000, y = .data[["Y_km"]] * 1000, fill = .data[[par]])) +
        geom_tile(color = NA) +
        {if(!is.null(borders_sf)) geom_sf(data = borders_sf, fill = NA, color = "black", size = 0.3, inherit.aes = FALSE)} +
        coord_sf(xlim = limit_x, ylim = limit_y, expand = FALSE, datum = st_crs(4326)) + 
        scale_fill_viridis_c(option = "plasma", name = NULL) +
        labs(title = get_plot_title(par), x = NULL, y = NULL) +
        theme_minimal(base_size = 12) +
        theme(
          legend.position = "right", 
          legend.key.width = unit(0.5, "cm"),
          legend.text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title = element_blank(), 
          axis.text = element_text(size = 8), 
          axis.text.x = element_text(angle = 45, hjust = 1), 
          panel.grid.major = element_line(color = "gray80", linetype = "dotted"), 
          panel.grid.minor = element_blank()
        )
    }
    plot_list[[par]] = p
  }
  
  combined_plot = ggpubr::ggarrange(plotlist = plot_list, nrow = 2, ncol = 3, align = "hv")
  print(combined_plot)
  
  # --- 5. PREPARAZIONE DEM (TOPOGRAFIA) ---
  cat("Download del DEM in corso (potrebbe richiedere qualche secondo)...\n")
  
  borders_sf_df <- st_as_sf(borders_sf)
  dem <- get_elev_raster(locations = borders_sf_df, z = 7)
  
  dem_terra <- terra::rast(dem)
  borders_vect <- terra::vect(borders_sf_df)
  
  dem_cropped <- terra::crop(dem_terra, borders_vect)
  dem_masked <- terra::mask(dem_cropped, borders_vect)
  
  dem_df <- as.data.frame(dem_masked, xy = TRUE)
  colnames(dem_df) <- c("x", "y", "elevation")
  dem_df <- dem_df[!is.na(dem_df$elevation), ]
  
  # --- 6. PLOT ELLISSI ---
  u_x = sort(unique(data$X_km))
  u_y = sort(unique(data$Y_km))
  
  n_x = min(length(u_x), n_grid)
  n_y = min(length(u_y), n_grid)
  
  idx_x = round(seq(1, length(u_x), length.out = n_x))
  idx_y = round(seq(1, length(u_y), length.out = n_y))
  
  data_ell = data %>% filter(X_km %in% u_x[idx_x] & Y_km %in% u_y[idx_y])
  data_ell <- smooth_ellipses_locally(grid_data = data, ell_data = data_ell, raggio_km = 25)
  
  dist_grid_km = (max(u_x[idx_x], na.rm=TRUE) - min(u_x[idx_x], na.rm=TRUE)) / (length(idx_x) - 1)
  
  data_ell$r_max <- pmax(sqrt(pmax(data_ell$lambda_1, 0, na.rm = TRUE)), sqrt(pmax(data_ell$lambda_2, 0, na.rm = TRUE)), na.rm = TRUE)
  
  raggio_rif <- quantile(data_ell$r_max[data_ell$r_max > 0], 0.80, na.rm = TRUE)
  if (is.na(raggio_rif) || raggio_rif == 0) raggio_rif <- 1
  
  scale_norm <- (dist_grid_km * 0.40) / raggio_rif
  
  data_ell$a_scaled <- sqrt(pmax(data_ell$lambda_1, 0, na.rm = TRUE)) * scale_norm
  data_ell$b_scaled <- sqrt(pmax(data_ell$lambda_2, 0, na.rm = TRUE)) * scale_norm
  
  max_draw_size <- dist_grid_km * 0.45
  data_ell$eccesso <- pmax(data_ell$a_scaled, data_ell$b_scaled) / max_draw_size
  data_ell$fattore_locale <- ifelse(data_ell$eccesso > 1, 1 / data_ell$eccesso, 1)
  
  data_ell$draw_a <- data_ell$a_scaled * data_ell$fattore_locale * 1000
  data_ell$draw_b <- data_ell$b_scaled * data_ell$fattore_locale * 1000
  data_ell$is_capped <- as.character(data_ell$eccesso > 1)
  data_ell <- data_ell %>% filter(!is.na(draw_a) & !is.na(draw_b))
  
  # Creazione del plot con DEM a colori e ellissi grigie
  p_ell = ggplot() +
    
    # Livello 1: Sfondo Topografico (A colori)
    # L'opacità è al 60% (alpha = 0.6) per renderlo gradevole e non far sparire i dati
    geom_raster(data = dem_df, aes(x = x, y = y, fill = elevation), alpha = 0.6) +
    # Utilizzo della palette standard topografica
    scale_fill_gradientn(colors = c("forestgreen", "green3",
                                    "chocolate", "saddlebrown", "white"),
                         guide = "none") +
    
    # Livello 2: Confini Italia
    {if(!is.null(borders_sf))
      geom_sf(data = borders_sf, fill = NA, color = "black", size = 0.4, inherit.aes = FALSE)} +
    
    # Livello 3: Ellissi (Tornate al colore originale gray20)
    geom_ellipse(data = data_ell, 
                 aes(x0 = X_km * 1000, y0 = Y_km * 1000, 
                     a = draw_a, 
                     b = draw_b, 
                     angle = theta_deg * pi / 180,
                     color = is_capped), 
                 fill = NA, 
                 size = 0.3, 
                 inherit.aes = FALSE) +
    
    # Livello 4: Centri
    geom_point(data = data_ell, 
               aes(x = X_km * 1000, y = Y_km * 1000, color = is_capped), 
               size = 0.6, inherit.aes = FALSE) +
    
    # Ripristinato gray20
    scale_color_manual(values = c("FALSE" = "gray20", "TRUE" = "gray20"), guide = "none") +    
    
    coord_sf(xlim = limit_x, ylim = limit_y, expand = FALSE, datum = st_crs(4326)) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      axis.title = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.grid.major = element_line(color = "gray60", linetype = "dotted") 
    )
  
  print(p_ell)
}

plot_parameters_sfondo_media = function(data, dir_path, K, borders_sf = NULL, n_grid = 10) {
  
  require(sf)
  require(rnaturalearth)
  require(dplyr)
  require(gridExtra)
  require(ggpubr) # Assicuriamoci che sia caricato per ggarrange
  
  # ==============================================================================
  # 1. CONFINI ITALIA SENZA SARDEGNA + MASCHERATURA DATI
  # ==============================================================================
  ita_regioni <- ne_states(country = "Italy", returnclass = "sf")
  
  ita_no_sardegna <- ita_regioni %>%
    filter(!grepl("Sardegn", region, ignore.case = TRUE))
  
  italia_utm <- st_transform(st_union(ita_no_sardegna), 32633)
  
  borders_sf <- italia_utm
  
  df_temp <- data
  df_temp$X_m <- df_temp$X_km * 1000
  df_temp$Y_m <- df_temp$Y_km * 1000
  
  sf_all <- st_as_sf(df_temp, coords = c("X_m", "Y_m"), crs = 32633, remove = FALSE)
  
  inside <- st_intersects(sf_all, italia_utm, sparse = FALSE)
  data <- data[as.vector(inside), ] 
  
  # ==============================================================================
  # 2. HELPER PER TITOLI CON LETTERE GRECHE, PEDICI E VARIANZE
  # ==============================================================================
  get_plot_title <- function(par_name) {
    if (par_name == "mu") return(expression(mu))
    if (par_name == "sigma2") return(expression(sigma^2)) 
    if (par_name == "lambda_1") return(expression(lambda[1]))
    if (par_name == "lambda_2") return(expression(lambda[2]))
    if (par_name == "theta_deg") return(expression(theta~"(°)"))
    if (par_name == "nugget") return(expression(sigma[nugget]^2)) 
    return(par_name)
  }
  # ==============================================================================
  
  # --- PREPARAZIONE DATI: CALCOLO VARIANZA (PARTIAL SILL) ---
  if("sigma" %in% names(data)) {
    data$sigma2 <- data$sigma^2
  }
  
  # MODIFICA PER CONVENZIONE AZIMUTALE
  
  param_names = c("mu", "sigma2", "lambda_1", "lambda_2", "theta_deg", "nugget")
  plot_list = list()
  
  # --- 3. CALCOLO LIMITI ZOOM (METRI) ---
  limit_x = range(data$X_km * 1000, na.rm = TRUE)
  limit_y = range(data$Y_km * 1000, na.rm = TRUE)
  
  buffer_x = diff(limit_x) * 0.05
  buffer_y = diff(limit_y) * 0.05
  
  final_xlim = c(limit_x[1] - buffer_x, limit_x[2] + buffer_x)
  final_ylim = c(limit_y[1] - buffer_y, limit_y[2] + buffer_y)
  
  # --- 4. HEATMAPS ---
  for (par in param_names) {
    if (!par %in% names(data)) {
      p = ggplot() + theme_void() + ggtitle(paste("Missing:", par))
    } else {
      p = ggplot(data, aes(x = .data[["X_km"]] * 1000, y = .data[["Y_km"]] * 1000, fill = .data[[par]])) +
        geom_tile(color = NA) +
        {if(!is.null(borders_sf)) geom_sf(data = borders_sf, fill = NA, color = "black", size = 0.3, inherit.aes = FALSE)} +
        
        coord_sf(xlim = limit_x, ylim = limit_y, expand = FALSE, datum = st_crs(4326)) + 
        
        scale_fill_viridis_c(option = "plasma", name = NULL) +
        labs(title = get_plot_title(par), x = NULL, y = NULL) +
        theme_minimal(base_size = 12) +
        theme(
          legend.position = "right", 
          legend.key.width = unit(0.5, "cm"),
          legend.text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title = element_blank(), 
          axis.text = element_text(size = 8), 
          axis.text.x = element_text(angle = 45, hjust = 1), 
          panel.grid.major = element_line(color = "gray80", linetype = "dotted"), 
          panel.grid.minor = element_blank()
        )
    }
    plot_list[[par]] = p
  }
  
  # MODIFICA QUI: Usiamo ggarrange con allineamento forzato orizzontale e verticale ('hv')
  combined_plot = ggpubr::ggarrange(plotlist = plot_list, nrow = 2, ncol = 3, align = "hv")
  
  # Mostra il grafico combinato ed allineato
  print(combined_plot)
  
  # --- 5. PLOT ELLISSI ---
  u_x = sort(unique(data$X_km))
  u_y = sort(unique(data$Y_km))
  
  n_x = min(length(u_x), n_grid)
  n_y = min(length(u_y), n_grid)
  
  idx_x = round(seq(1, length(u_x), length.out = n_x))
  idx_y = round(seq(1, length(u_y), length.out = n_y))
  
  data_ell = data %>% filter(X_km %in% u_x[idx_x] & Y_km %in% u_y[idx_y])
  data_ell <- smooth_ellipses_locally(grid_data = data, ell_data = data_ell, raggio_km = 25)
  
  dist_grid_km = (max(u_x[idx_x], na.rm=TRUE) - min(u_x[idx_x], na.rm=TRUE)) / (length(idx_x) - 1)
  
  data_ell$r_max <- pmax(sqrt(pmax(data_ell$lambda_1, 0, na.rm = TRUE)), sqrt(pmax(data_ell$lambda_2, 0, na.rm = TRUE)), na.rm = TRUE)
  
  raggio_rif <- quantile(data_ell$r_max[data_ell$r_max > 0], 0.80, na.rm = TRUE)
  if (is.na(raggio_rif) || raggio_rif == 0) raggio_rif <- 1
  
  scale_norm <- (dist_grid_km * 0.40) / raggio_rif
  
  data_ell$a_scaled <- sqrt(pmax(data_ell$lambda_1, 0, na.rm = TRUE)) * scale_norm
  data_ell$b_scaled <- sqrt(pmax(data_ell$lambda_2, 0, na.rm = TRUE)) * scale_norm
  
  max_draw_size <- dist_grid_km * 0.45
  data_ell$eccesso <- pmax(data_ell$a_scaled, data_ell$b_scaled) / max_draw_size
  
  data_ell$fattore_locale <- ifelse(data_ell$eccesso > 1, 1 / data_ell$eccesso, 1)
  
  data_ell$draw_a <- data_ell$a_scaled * data_ell$fattore_locale * 1000
  data_ell$draw_b <- data_ell$b_scaled * data_ell$fattore_locale * 1000
  
  data_ell$is_capped <- as.character(data_ell$eccesso > 1)
  
  data_ell <- data_ell %>% filter(!is.na(draw_a) & !is.na(draw_b))
  
  p_ell = ggplot(data, aes(x = X_km * 1000, y = Y_km * 1000)) +
    geom_raster(aes(fill = mu), alpha = 0.3) +
    scale_fill_viridis_c(option = "gray", guide = "none") +
    {if(!is.null(borders_sf))
      geom_sf(data = borders_sf, fill = NA, color = "gray20", size = 0.3, inherit.aes = FALSE)} +
    
    geom_ellipse(data = data_ell, 
                 aes(x0 = X_km * 1000, y0 = Y_km * 1000, 
                     a = draw_a, 
                     b = draw_b, 
                     angle = theta_deg * pi / 180,
                     color = is_capped), 
                 fill = NA, 
                 size = 0.1, 
                 inherit.aes = FALSE) +
    
    geom_point(data = data_ell, 
               aes(x = X_km * 1000, y = Y_km * 1000, color = is_capped), 
               size = 0.3) +
    
    scale_color_manual(values = c("FALSE" = "gray20", "TRUE" = "gray20"), guide = "none") +    
    coord_sf(xlim = limit_x, ylim = limit_y, expand = FALSE, datum = st_crs(4326)) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      axis.title = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.grid.major = element_line(color = "gray80", linetype = "dotted") 
    )
  
  print(p_ell)
}

plot_parameters_old = function(data, dir_path, K, borders_sf = NULL, n_grid = 10) {
  
  # --- SETUP ---
  param_names = c("mu", "sigma", "lambda_1", "lambda_2", "theta_deg", "nugget")
  
  plot_dir = file.path(dir_path, paste0("Plots_K", K))
  if(!dir.exists(plot_dir)) dir.create(plot_dir)
  
  plot_list = list()
  
  # --- 1. CALCOLO LIMITI ZOOM (METRI) ---
  # Calcoliamo i limiti esatti basati sui tuoi dati (la griglia colorata)
  limit_x = c(min(data$X_km * 1000), max(data$X_km * 1000))
  limit_y = c(min(data$Y_km * 1000), max(data$Y_km * 1000))
  
  # Buffer piccolo (5%) per non far toccare i dati ai bordi del grafico
  buffer_x = diff(limit_x) * 0.05
  buffer_y = diff(limit_y) * 0.05
  
  final_xlim = c(limit_x[1] - buffer_x, limit_x[2] + buffer_x)
  final_ylim = c(limit_y[1] - buffer_y, limit_y[2] + buffer_y)
  
  
  # --- 2. HEATMAPS ---
  for (par in param_names) {
    if (!par %in% names(data)) {
      p = ggplot() + theme_void() + ggtitle(paste("Missing:", par))
    } else {
        p = ggplot(data, aes(x = .data[["X_km"]] * 1000, y = .data[["Y_km"]] * 1000, fill = .data[[par]])) +
        
        # color = NA toglie la griglia grigia tra i pixel
        geom_tile(color = NA) +
        
        # Disegno i bordi ORIGINALI (senza tagliarli)
        # Il trucco è che coord_sf nasconderà le parti fuori dal grafico senza disegnare linee nere
        {if(!is.null(borders_sf)) geom_sf(data = borders_sf, fill = NA, color = "black", size = 0.3, inherit.aes = FALSE)} +
        
        # expand = FALSE è fondamentale: taglia esattamente dove diciamo noi
        coord_sf(xlim = limit_x, ylim = limit_y, expand = FALSE) + 
        
        scale_fill_viridis_c(option = "plasma", name = NULL) +
        labs(title = par, x = NULL, y = NULL) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "right", legend.key.width = unit(0.5, "cm"),
              legend.text = element_text(size = 12),
              plot.title = element_text(hjust = 0.5, face = "bold"),
              axis.text = element_blank(), panel.grid = element_blank())
    }
    plot_list[[par]] = p
  }
  
  # Modifica il nome file qui:
  heatmap_filename = file.path(plot_dir, paste0("heatmap_K", K, ".png"))
  combined_plot = arrangeGrob(grobs = plot_list, nrow = 2, ncol = 3)
  
  # Aggiungi questa riga per eseguire il salvataggio:
  ggsave(filename = heatmap_filename, plot = combined_plot, width = 14, height = 8, dpi = 300, bg = "white")
  
  grid.arrange(combined_plot)
  
  # --- 3. PLOT ELLISSI ---
  
  u_x = sort(unique(data$X_km))
  u_y = sort(unique(data$Y_km))
  
  n_x = min(length(u_x), n_grid)
  n_y = min(length(u_y), n_grid)
  
  idx_x = round(seq(1, length(u_x), length.out = n_x))
  idx_y = round(seq(1, length(u_y), length.out = n_y))
  
  data_ell = data %>% filter(X_km %in% u_x[idx_x] & Y_km %in% u_y[idx_y])
  data_ell <- smooth_ellipses_locally(grid_data = data, ell_data = data_ell, raggio_km = 25)
  
  # ----------------------------------------------------------------------------
  # NUOVA LOGICA DI SCALATURA (PERCENTILE ALTO + LOCAL CAPPING)
  # ----------------------------------------------------------------------------
  dist_grid_km = (max(u_x[idx_x]) - min(u_x[idx_x])) / (length(idx_x) - 1)
  
  # 1. Trovo il raggio maggiore per ogni cella della griglia
  data_ell$r_max <- pmax(sqrt(pmax(data_ell$lambda_1, 0, na.rm = TRUE)), sqrt(pmax(data_ell$lambda_2, 0, na.rm = TRUE)), na.rm = TRUE)
  
  # 2. LA SOLUZIONE: Non usiamo più la mediana! Usiamo l'80° percentile. 
  # In questo modo ignoriamo il "mare" (che ha valori piccoli) e calibriamo 
  # la scala sulle ellissi reali della terraferma.
  raggio_rif <- quantile(data_ell$r_max[data_ell$r_max > 0], 0.80, na.rm = TRUE)
  if (is.na(raggio_rif) || raggio_rif == 0) raggio_rif <- 1
  
  # Calcolo la scala globale affinché questa ellisse "vera" occupi l'80% della cella (raggio 40%)
  scale_norm <- (dist_grid_km * 0.40) / raggio_rif
  
  # 3. Applico la scala globale ai raggi
  data_ell$a_scaled <- sqrt(pmax(data_ell$lambda_1, 0, na.rm = TRUE)) * scale_norm
  data_ell$b_scaled <- sqrt(pmax(data_ell$lambda_2, 0, na.rm = TRUE)) * scale_norm
  
  # 4. CAP VISIVO LOCALE: Tetto massimo per le ellissi (semi-asse al 45% della cella)
  max_draw_size <- dist_grid_km * 0.45
  data_ell$eccesso <- pmax(data_ell$a_scaled, data_ell$b_scaled) / max_draw_size
  
  data_ell$fattore_locale <- ifelse(data_ell$eccesso > 1, 1 / data_ell$eccesso, 1)
  
  data_ell$draw_a <- data_ell$a_scaled * data_ell$fattore_locale * 1000
  data_ell$draw_b <- data_ell$b_scaled * data_ell$fattore_locale * 1000
  
  # Cast esplicito a stringa per massima compatibilità con scale_color_manual
  data_ell$is_capped <- as.character(data_ell$eccesso > 1)
  
  # Filtro di sicurezza anti-crash (rimuove NAs se il modello non ha convinto)
  data_ell <- data_ell %>% filter(!is.na(draw_a) & !is.na(draw_b))
  # ----------------------------------------------------------------------------
  
  # Creazione del plot con ggplot
  p_ell = ggplot(data, aes(x = X_km * 1000, y = Y_km * 1000)) +
    
    geom_raster(aes(fill = mu), alpha = 0.3) +
    scale_fill_viridis_c(option = "gray", guide = "none") +
    {if(!is.null(borders_sf))
      # ATTENZIONE: Tornato a "size = 0.3" al posto di linewidth
      geom_sf(data = borders_sf, fill = NA, color = "gray20", size = 0.3, inherit.aes = FALSE)} +
    
    # Plot delle ellissi 
    geom_ellipse(data = data_ell, 
                 aes(x0 = X_km * 1000, y0 = Y_km * 1000, 
                     a = draw_a, 
                     b = draw_b, 
                     angle = theta_deg * pi / 180,
                     color = is_capped), 
                 fill = NA, 
                 size = 0.1, # ATTENZIONE: Tornato a "size = 0.1" al posto di linewidth
                 inherit.aes = FALSE) +
    
    # Plot dei punti centrali
    geom_point(data = data_ell, 
               aes(x = X_km * 1000, y = Y_km * 1000, color = is_capped), 
               size = 0.3) +
    
    # La mappatura dei colori ora legge i testi "TRUE" e "FALSE"
    scale_color_manual(values = c("FALSE" = "gray20", "TRUE" = "gray20"), guide = "none") +    
    coord_sf(xlim = limit_x, ylim = limit_y, expand = FALSE) +
    labs( x = "Longitude (°)", y = "Latitude (°)" ) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(), # <-- AGGIUNTO QUI: spegne la griglia fitta
      panel.grid.major = element_line(color = "gray80")  # Decommenta questa se vuoi spegnere anche la griglia principale
    )
  
  ellipses_filename = file.path(plot_dir, paste0("ellipses_K", K, ".png"))
  ggsave(filename = ellipses_filename, plot = p_ell, width = 8, height = 8, dpi = 300, bg = "white")
  
  print(p_ell)
   
}

if (!exists("modified_bessel_second_kind")) {
  modified_bessel_second_kind <- function(nu, x) {
    # besselK è la funzione nativa di R identica a quella di Stan
    return(besselK(x, nu))
  }
}

get_induced_covariance_matrix <- function(preds_aggr, coords, nu = 0.5) {
  # --- 1. GESTIONE THETA (Gradi vs Radianti) ---
  
  # Caso 1: Esiste già la colonna "theta" (tipicamente radianti dal tuo codice)
  if("theta" %in% names(preds_aggr)) {
    theta_vals <- preds_aggr$theta
    
    # Controllo di sicurezza: se troviamo valori > 7 (2*pi è circa 6.28), 
    # è sospetto che siano gradi chiamati erroneamente "theta".
    if(max(abs(theta_vals), na.rm=T) > 7) {
      warning("ATTENZIONE: La colonna 'theta' ha valori molto alti (>7). 
              Sembrano gradi ma RGeostats/Trig richiedono radianti. 
              Verifica i dati!")
    }
    
  } else if ("theta_deg" %in% names(preds_aggr)) {
    # Caso 2: Abbiamo solo i gradi, li convertiamo noi
    theta_vals <- preds_aggr$theta_deg * (pi / 180)
  } else {
    stop("Errore: Non trovo né 'theta' né 'theta_deg' nel dataframe.")
  }
  
  # --- 2. Estrazione altri parametri ---
  lambda_1 <- preds_aggr$lambda_1
  lambda_2 <- preds_aggr$lambda_2
  sigma    <- preds_aggr$sigma
  
  # --- 3. Calcolo Anisotropia ---
  # Passiamo i theta validati
  Aniso_List <- Compute_Aniso(lambda_1, lambda_2, theta_vals)
  
  # Determinante = prodotto autovalori
  dets <- lambda_1 * lambda_2
  
  # --- 4. Coordinate ---
  if (is.data.frame(coords)) coords <- as.matrix(coords)
  Locations_t <- t(coords)
  
  # --- 5. Funzione Bessel di fallback ---
  if (!exists("modified_bessel_second_kind")) {
    modified_bessel_second_kind <- function(nu, x) {
      # Gestione x=0 (distanza nulla) -> correlazione 1 (limite)
      # besselK va a infinito a 0, ma la formula di Matern gestisce il limite.
      # Tuttavia, per sicurezza numerica su distanze piccolissime:
      val <- besselK(x, nu)
      return(val)
    }
  }
  
  # --- 6. Calcolo Covarianza ---
  Cov_Matrix <- matern_ns_corr(
    Locations = Locations_t, 
    Aniso = Aniso_List, 
    dets = dets, 
    nu = nu, 
    sigma = sigma
  )
  
  return(Cov_Matrix)
}
