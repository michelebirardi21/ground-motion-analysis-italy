
# SE SERVE, MODIFICARE "bin_width", "filter(n_pairs >= 10)" E "step_size_km"

# ------------------------------------------------------------------------------
# 1. librerie e dati
# ------------------------------------------------------------------------------

library(sf)
library(sp)
library(ggplot2)
library(viridis)
library(tidyr)
library(plotly)
library(reticulate)
library(psych)

library(foreach)
library(doParallel)

## run dei bambini:
#n_cores <- 3
#step <- 70

# run degli adulti:
n_cores <- 20
step <- 7

cl <- makePSOCKcluster(n_cores)
registerDoParallel(cl)

#setwd("C:/Users/miche/OneDrive/Desktop/prova3/INGV - pgv/Effetto stazione non tolto/Non Parametrico")
#source("utils2.R")
#setwd("C:/Users/miche/OneDrive/Desktop/prova3/INGV - pgv")
#source("ITA_18.R")
#data = read.csv("dataset_martina.csv")

# versione per cluster
setwd("/work/u10700026/nonparam/")
source("utils2.R")
source("ITA_18.R")
data = read.csv("dataset_martina.csv")

# --- INIZIO PULIZIA STAZIONI DUPLICATE/ANOMALE ---

# a) Rimuovo le stazioni da scartare
stazioni_da_buttare <- c("IT.GRG1", "E.ATQ", "E.ATT", "E.ATP", "XO.MN06")
data <- subset(data, !(NET_STA %in% stazioni_da_buttare))

# b) Aggrego le stazioni (rinomino i duplicati tenendo il nome principale)
data$NET_STA[data$NET_STA == "YI.SACS"] <- "IV.SACS"
data$NET_STA[data$NET_STA == "YI.ARCI"] <- "IV.ARCI"
data$NET_STA[data$NET_STA == "YI.MURB"] <- "IV.MURB"

# --- FINE PULIZIA STAZIONI ---

data$EVENT = as.factor(data$EVENT)
data$NET_STA = as.factor(data$NET_STA)
data$sf = as.factor(data$sf)

lat_min <- 41.5
lat_max <- 44
lon_min <- 10.5
lon_max <- 15

itacentraletemp <- subset(data, ev_lat >= lat_min & ev_lat <= lat_max &
                            ev_lon >= lon_min & ev_lon <= lon_max)

x1 = 14.4; y1 = 42.8; x2 = 12.6; y2 = 44
m = (y1 - y2)/(x1 - x2)
q = y1 - m*x1

itacentrale <- subset(itacentraletemp, (ev_lat < m*ev_lon + q) |
                        ( (ev_lon - 13.4)^2 + (ev_lat - 43.6)^2 < 0.05^2) )
itacentrale$EVENT <- droplevels(itacentrale$EVENT)

IM = -1
M = itacentrale$Mw_RCMT
R = itacentrale$RR
vs30 =itacentrale$vs30
sf = itacentrale$sf

ita18 = ITA18_RJB(IM, M, R, vs30, sf)

log_oss = log10(itacentrale$pgv)
log_pred = log10(ita18$pred)
l2b <- norm(matrix(log_oss - log_pred, ncol = 1), type = "2")

res = log_oss - log_pred
itacentrale$res = res

library(lme4)
library(lattice)

# 1. Fit del modello incrociato (stima simultanea degli effetti di evento e di stazione)
mod_completo <- lmer(res ~ 1 + (1 | EVENT) + (1 | NET_STA), data = itacentrale, REML = TRUE)

# 2. Calcoliamo i residui intra-evento (dWs = res - offset - dBe)
# Usiamo predict specificando di usare SOLO l'effetto random dell'evento
fitted_solo_evento <- predict(mod_completo, re.form = ~ (1 | EVENT))
itacentrale$dWs <- itacentrale$res - fitted_solo_evento

# ------------------------------------------------------------------------------
# 2. Aggregazione per stazione
# ------------------------------------------------------------------------------
library(dplyr)

obs <- itacentrale %>%
  filter(!is.na(dWs)) %>%
  group_by(NET_STA) %>%
  summarise(
    st_lon = first(st_lon),
    st_lat = first(st_lat),
    n_obs = n()
  ) %>%
  as.data.frame() 

# ------------------------------------------------------------------------------
# 3. proietta in coordinate planari
# ------------------------------------------------------------------------------

sf_pts <- st_as_sf(obs, coords = c("st_lon","st_lat"), crs = 4326)
sf_pts_proj <- st_transform(sf_pts, crs = 32633)
coords_proj <- st_coordinates(sf_pts_proj)
coords_proj_km <- coords_proj / 1000

# ------------------------------------------------------------------------------
# 4. Estrazione dS2S (BLUPs) e Stima Empirica con LOESS
# ------------------------------------------------------------------------------

# Estraiamo i BLUPs per le stazioni
dS2S_blup <- ranef(mod_completo)$NET_STA
dS2S_blup$NET_STA <- rownames(dS2S_blup)
colnames(dS2S_blup)[1] <- "Z_i"

# Uniamo i BLUPs alle coordinate delle stazioni (obs)
obs$X_km <- coords_proj_km[, 1]
obs$Y_km <- coords_proj_km[, 2]
obs <- merge(obs, dS2S_blup, by = "NET_STA") 

Z_vec <- obs$Z_i
n_sta <- nrow(obs)
mu_Z <- mean(Z_vec)

# Varianza h=0 esatta estratta dai componenti di varianza del modello misto
var_components <- as.data.frame(VarCorr(mod_completo))
var_h0 <- var_components$vcov[var_components$grp == "NET_STA"]

cat("\nVarianza di sito (h=0) stimata dal LMM:", var_h0, "\n")

# b) Calcolo delle coppie e Covarianza Empirica
dist_matrix <- as.matrix(dist(cbind(obs$X_km, obs$Y_km)))
pairs_df <- expand.grid(i = 1:n_sta, j = 1:n_sta)
pairs_df <- pairs_df[pairs_df$i < pairs_df$j, ]

pairs_df$dist <- dist_matrix[cbind(pairs_df$i, pairs_df$j)]
pairs_df$Z_i <- Z_vec[pairs_df$i]
pairs_df$Z_j <- Z_vec[pairs_df$j]
pairs_df$cov_product <- (pairs_df$Z_i - mu_Z) * (pairs_df$Z_j - mu_Z)

bin_width <- 10
breaks_dist <- seq(0, max(pairs_df$dist) + bin_width, by = bin_width)
pairs_df$bin_factor <- cut(pairs_df$dist, breaks = breaks_dist, include.lowest = TRUE)

empirical_cov <- pairs_df %>%
  group_by(bin_factor) %>%
  summarise(
    mean_dist = mean(dist),           
    cov_val = mean(cov_product),      
    n_pairs = n()                     
  ) %>%
  filter(n_pairs >= 10) 

# c) Smoothing Non-Parametrico (LOESS) pesato per il numero di coppie
# Questo risolve il problema dei "buchi" creando una curva continua
mod_loess <- loess(cov_val ~ mean_dist, data = empirical_cov, weights = n_pairs, span = 0.5)

# Aggiungiamo i valori fittati per il grafico
empirical_cov$loess_pred <- predict(mod_loess)

row_h0 <- data.frame(bin_factor = NA, mean_dist = 0, cov_val = var_h0, n_pairs = n_sta, loess_pred = var_h0)
empirical_cov <- bind_rows(row_h0, empirical_cov)

ggplot(empirical_cov, aes(x = mean_dist)) +
  geom_point(aes(y = cov_val, size = n_pairs), alpha = 0.5, color = "blue") +
  geom_line(aes(y = loess_pred), color = "red", linewidth = 1.2) + # La curva LOESS!
  geom_hline(yintercept = 0, color = "black") +
  labs(title = "Covarianza Empirica dS2S smussata con LOESS", x = "Distanza (km)", y = "Covarianza") +
  theme_minimal()

# ------------------------------------------------------------------------------
# 5. Creazione Griglia Spaziale (ALLINEATA CON RDD)
# ------------------------------------------------------------------------------
step_size_km <- step  # Usa 50 per test veloci, 7 per i risultati definitivi
buffer_km <- 20    # Stesso buffer dell'RDD

# 1. Troviamo i limiti massimi includendo il buffer
min_X <- min(obs$X_km) - buffer_km
max_X <- max(obs$X_km) + buffer_km
min_Y <- min(obs$Y_km) - buffer_km
max_Y <- max(obs$Y_km) + buffer_km

# 2. Forziamo un bounding box QUADRATO perfetto (Fondamentale per la serpentina)
span_X <- max_X - min_X
span_Y <- max_Y - min_Y
max_span <- max(span_X, span_Y)

mid_X <- (max_X + min_X) / 2
mid_Y <- (max_Y + min_Y) / 2

min_X <- mid_X - max_span / 2
max_X <- mid_X + max_span / 2
min_Y <- mid_Y - max_span / 2
max_Y <- mid_Y + max_span / 2

seq_X <- seq(min_X, max_X, by = step_size_km)
seq_Y <- seq(max_Y, min_Y, by = -step_size_km) 
N_side <- max(length(seq_X), length(seq_Y))

seq_X <- seq(min_X, by = step_size_km, length.out = N_side)
seq_Y <- seq(max_Y, by = -step_size_km, length.out = N_side)

# 3. Creazione griglia completa
grid_square <- expand.grid(col_idx = 0:(N_side-1), row_idx = 0:(N_side-1))
grid_square$X_km <- seq_X[grid_square$col_idx + 1]
grid_square$Y_km <- seq_Y[grid_square$row_idx + 1]

# 4. Indice a serpentina ESATTO per il Laplaciano
grid_square$serpentine_idx <- ifelse(
  grid_square$row_idx %% 2 == 0,
  grid_square$row_idx * N_side + grid_square$col_idx + 1,          
  grid_square$row_idx * N_side + (N_side - 1 - grid_square$col_idx) + 1 
)

grid_square <- grid_square[order(grid_square$serpentine_idx), ]
grid_square$cell_id <- grid_square$serpentine_idx
K_total <- nrow(grid_square)

cat("\nGriglia quadrata condivisa generata:", N_side, "x", N_side, "(", K_total, "celle totali )\n")

# Assegniamo le stazioni alla cella più vicina
trova_cella_vicina_km <- function(x_sta, y_sta, grid_df) {
  distanze_sq <- (grid_df$X_km - x_sta)^2 + (grid_df$Y_km - y_sta)^2
  return(grid_df$cell_id[which.min(distanze_sq)])
}
obs$cell_id <- mapply(trova_cella_vicina_km, obs$X_km, obs$Y_km, MoreArgs = list(grid_df = grid_square))

# NOTA BENE: Non filtriamo la terraferma qui! Lo faremo nel codice ponte.

# ------------------------------------------------------------------------------
# 6. Matrice di Covarianza sulle Celle Popolate e nearPD
# ------------------------------------------------------------------------------
library(Matrix)

# Isoliamo SOLO le celle attive
celle_attive_id <- unique(obs$cell_id)
grid_attiva <- grid_square[grid_square$cell_id %in% celle_attive_id, ]
n_attive <- nrow(grid_attiva)

dist_celle <- as.matrix(dist(cbind(grid_attiva$X_km, grid_attiva$Y_km)))
max_dist_loess <- max(empirical_cov$mean_dist[-1], na.rm = TRUE)
C_attiva_vec <- predict(mod_loess, newdata = data.frame(mean_dist = as.vector(dist_celle)))

# USIAMO L'HACK DELL' 1e-8!
C_attiva_vec[as.vector(dist_celle) > max_dist_loess] <- 1e-8
C_attiva_vec[is.na(C_attiva_vec)] <- 1e-8

C_attiva_matrix <- matrix(C_attiva_vec, nrow = n_attive, ncol = n_attive)
diag(C_attiva_matrix) <- var_h0

risultato_nearPD <- nearPD(C_attiva_matrix, corr = FALSE, base.matrix = TRUE)
C_attiva_PD <- as.matrix(risultato_nearPD$mat)

# ------------------------------------------------------------------------------
# 7. Espansione nella Matrice Sparsa Finale
# ------------------------------------------------------------------------------

coppie_celle <- expand.grid(idx_1 = 1:n_attive, idx_2 = 1:n_attive)
coppie_celle$cell_id_1 <- grid_attiva$cell_id[coppie_celle$idx_1]
coppie_celle$cell_id_2 <- grid_attiva$cell_id[coppie_celle$idx_2]
coppie_celle$cov_val <- as.vector(C_attiva_PD)

# Creazione Matrice Sparsa. I cell_id sono già gli indici a serpentina perfetti!
C_griglia_sparsa <- sparseMatrix(
  i = coppie_celle$cell_id_1,
  j = coppie_celle$cell_id_2,
  x = coppie_celle$cov_val,
  dims = c(K_total, K_total)
)

# ------------------------------------------------------------------------------
# 9. Ricostruzione
# ------------------------------------------------------------------------------

alpha_best <- 1
m_best <- 0.5

C_densa_input <- as.matrix(C_griglia_sparsa)
alpha_param <- alpha_best
m_param <- m_best
rank_param <- min(100, K_total)

cat("Avvio di covar_reconstruction_2D()...\n")

# --- INIZIO HACK R-SCOPING ---

# Hack 1: Variabile globale mancante in utils2.R
covar <- C_densa_input 

# Hack 2: Bypass del bug "m must be a square matrix"
tr <- function(m) {
  sum(diag(as.matrix(m)))
}

# Hack 3: Curiamo l'oggetto S4 per non far crashare LBFGS in C++
# Salviamo la funzione originale di utils2.R
laplacian_2d_originale <- laplacian_2d 

# La sovrascriviamo localmente forzando l'output a matrice base
laplacian_2d <- function(n_rows, n_cols) {
  as.matrix(laplacian_2d_originale(n_rows, n_cols))
}

eigen <- function(x, ...) {
  e <- base::eigen(x, ...);
  e$values[e$values < 0] <- 1e-12;
  return(e) }
svd <- function(x, ...) {
  s <- base::svd(x, ...);
  if(!is.null(s$d)) s$d[s$d < 0] <- 1e-12;
  if(!is.null(s$values)) s$values[s$values < 0] <- 1e-12;
  return(s) }

# --- FINE HACK ---

cat("Calcolo in corso (potrebbe richiedere qualche minuto)...\n")
t_inizio_ricostruzione <- Sys.time()

# Chiamata alla funzione standard
risultato_ricostruzione <- covar_reconstruction_2D(
  Rn = C_densa_input, 
  alpha = alpha_param, 
  m = m_param,
  rk = rank_param
)

t_fine_ricostruzione <- Sys.time()

C_globale_ricostruita <- risultato_ricostruzione$cov

minuti_totali_rec <- as.numeric(difftime(t_fine_ricostruzione, t_inizio_ricostruzione, units = "mins"))
ore_rec <- floor(minuti_totali_rec / 60)
minuti_restanti_rec <- round(minuti_totali_rec %% 60)

cat(sprintf("\n=> TEMPO DI ESECUZIONE RICOSTRUZIONE SEZ 9: %d ore e %d minuti\n", ore_rec, minuti_restanti_rec))

# Pulizia: Rimuoviamo i nostri hack
rm(covar, tr, laplacian_2d)
laplacian_2d <- laplacian_2d_originale # Ripristiniamo l'originale per sicurezza
rm(laplacian_2d_originale)

cat("\nRicostruzione completata con successo!\n")

# ------------------------------------------------------------------------------
# 10. Simple Kriging Data-Driven
# ------------------------------------------------------------------------------
library(MASS) # Necessario per ginv()

cat("\nAvvio del Kriging Spaziale...\n")

# 1. Prepariamo i dati osservati (Assicurati che l'ordine corrisponda!)
obs <- obs[order(obs$cell_id), ]
Z_obs <- obs$Z_i
idx_obs <- obs$cell_id

# 2. Estraiamo le sottomatrici dalla nostra matrice ricostruita
# Covarianza tra i punti osservati (Stazione vs Stazione)
Sigma_11 <- C_globale_ricostruita[idx_obs, idx_obs] 

# Covarianza tra TUTTA la griglia e i punti osservati
# Usiamo K_total per assicurarci di prendere tutte le righe
idx_all <- 1:K_total 
Sigma_all_1 <- C_globale_ricostruita[idx_all, idx_obs]

# 3. Calcolo dei pesi di Kriging
# Usiamo ginv() (Generalized Inverse) invece di solve() perché le stazioni 
# molto vicine potrebbero rendere la matrice Sigma_11 mal condizionata o singolare.

# Centratura preventiva (nel caso la media empirica non sia uno zero perfetto)
mu_global <- mean(Z_obs)

Sigma_11_inv <- ginv(Sigma_11)
pesi_kriging <- Sigma_all_1 %*% Sigma_11_inv

# 4. Predizione Simple Kriging con media esplicita
# Z_pred = mu + Pesi * (Z_obs - mu)
Z_pred_all <- mu_global + pesi_kriging %*% (Z_obs - mu_global)

# 5. Mappatura dei risultati sulla griglia geografica
grid_square$dS2S_kriging <- Z_pred_all[grid_square$cell_id]
# se partono warning: grid_square$dS2S_kriging <- as.vector(Z_pred_all)[grid_square$cell_id]

# --- PLOT DEL RISULTATO ---
fig_kriging <- ggplot(grid_square, aes(x = X_km, y = Y_km, fill = dS2S_kriging)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "dS2S (Kriging)") +
  # Aggiungiamo i puntini delle stazioni originali per vedere come fittano
  geom_point(data = obs, aes(x = X_km, y = Y_km), color = "black", size = 0.8, shape = 16, inherit.aes = FALSE) +
  labs(
    title = "Mappa Spaziale Predittiva del dS2S (Simple Kriging)",
    subtitle = "Interpolazione condizionata alla matrice di covarianza non parametrica",
    x = "Est-Ovest (km UTM)",
    y = "Nord-Sud (km UTM)"
  ) +
  theme_minimal() +
  coord_fixed()

print(fig_kriging)

cat("Calcolo dell'incertezza di Kriging in corso...\n")

# Estraiamo la varianza a priori di ogni cella (la diagonale della super-matrice)
# Dovrebbe essere molto simile al tuo var_h0
var_a_priori <- diag(C_globale_ricostruita)

# Calcolo EFFICIENTE della riduzione di varianza.
# Invece di moltiplicare matrici enormi per poi estrarre la diagonale, 
# facciamo il prodotto elemento per elemento (*) tra la matrice dei pesi 
# e la matrice di covarianza incrociata, e poi sommiamo per riga.
# Matematicamente equivale a estrarre la diagonale di (Pesi %*% Sigma_obs_all)
riduzione_varianza <- rowSums(pesi_kriging * Sigma_all_1)

# Varianza di Kriging = Varianza a priori - Riduzione garantita dai dati
var_kriging <- var_a_priori - riduzione_varianza

# Per sicurezza numerica (effetti di arrotondamento), forziamo eventuali valori < 0 a zero
var_kriging[var_kriging < 0] <- 0

# Convertiamo in Deviazione Standard (lo Standard Error della predizione)
grid_square$se_kriging <- sqrt(var_kriging)

# --- PLOT DELLA MAPPA DI INCERTEZZA ---
fig_incertezza <- ggplot(grid_square, aes(x = X_km, y = Y_km, fill = se_kriging)) +
  geom_tile() +
  scale_fill_viridis(option = "mako", direction = -1, name = "Std Error\n(dS2S)") +
  # Aggiungiamo le stazioni in rosso per vedere come l'incertezza crolli vicino ad esse
  geom_point(data = obs, aes(x = X_km, y = Y_km), color = "red", size = 1, shape = 16, inherit.aes = FALSE) +
  labs(
    title = "Mappa dell'Incertezza di Previsione (Standard Error)",
    subtitle = "I valori scuri (incertezza bassa) corrispondono alle aree vicine alle stazioni",
    x = "Est-Ovest (km UTM)",
    y = "Nord-Sud (km UTM)"
  ) +
  theme_minimal() +
  coord_fixed()

print(fig_incertezza)
# ------------------------------------------------------------------------------
# 11. Cross-Validation (Kriging RMSE)
# ------------------------------------------------------------------------------
cat("\n--- AVVIO CROSS-VALIDATION KRIGING RIGOROSA (NO LEAKAGE) ---\n")

n_folds <- 10
set.seed(123) 

# 1. Estraiamo gli ID univoci delle celle in cui cade almeno una stazione
celle_attive_univoche <- unique(obs$cell_id)

# 2. Assegniamo un fold (da 1 a 10) a ciascuna CELLA
fold_celle <- sample(rep(1:n_folds, length.out = length(celle_attive_univoche)))

# 3. Creiamo una mappa di associazione "Cella -> Fold"
mappa_fold <- data.frame(cell_id = celle_attive_univoche, fold = fold_celle)

# 4. Assegniamo il fold alle singole stazioni unendo la mappa (se avevi già una colonna 'fold', meglio rimuoverla prima)
if("fold" %in% names(obs)) obs$fold <- NULL 

obs$fold <- mappa_fold$fold[match(obs$cell_id, mappa_fold$cell_id)]

# Re-inizializziamo il vettore dei risultati (usando nrow(obs) perché le stazioni restano lo stesso numero)
Z_hat_cv <- numeric(nrow(obs))

# Re-inizializziamo gli hack necessari per la ricostruzione dentro il loop PARALLELO
risultati_cv_list <- foreach(k = 1:n_folds, .combine = 'rbind', .packages = c("Matrix", "lbfgs", "dplyr", "MASS")) %dopar% {
  
  # Re-inizializziamo le dipendenze per ogni worker "cieco"
  source("utils2.R", local = TRUE)
  
  # GLI HACK VENGONO CREATI SOLO QUI DENTRO
  tr <- function(m) sum(diag(as.matrix(m)))
  lap_orig <- laplacian_2d # Salviamo la funzione originale
  laplacian_2d <- function(n_rows, n_cols) as.matrix(lap_orig(n_rows, n_cols)) # La sovrascriviamo
  
  eigen <- function(x, ...) {
    e <- base::eigen(x, ...);
    e$values[e$values < 0] <- 1e-12;
    return(e) }
  svd <- function(x, ...) {
    s <- base::svd(x, ...);
    if(!is.null(s$d)) s$d[s$d < 0] <- 1e-12;
    if(!is.null(s$values)) s$values[s$values < 0] <- 1e-12;
    return(s) }
  
  idx_test_df <- which(obs$fold == k)
  idx_train_df <- which(obs$fold != k)
  
  if(length(idx_test_df) == 0) return(NULL)
  
  # --- 1. ISOLAMENTO DATI DI TRAINING ---
  obs_train <- obs[idx_train_df, ]
  n_sta_train <- nrow(obs_train)
  mu_Z_train <- mean(obs_train$Z_i)
  
  # --- 2. RICALCOLO COVARIANZA EMPIRICA E LOESS SUL TRAINING ---
  dist_matrix_tr <- as.matrix(dist(cbind(obs_train$X_km, obs_train$Y_km)))
  pairs_df_tr <- expand.grid(i = 1:n_sta_train, j = 1:n_sta_train)
  pairs_df_tr <- pairs_df_tr[pairs_df_tr$i < pairs_df_tr$j, ]
  
  pairs_df_tr$dist <- dist_matrix_tr[cbind(pairs_df_tr$i, pairs_df_tr$j)]
  pairs_df_tr$Z_i <- obs_train$Z_i[pairs_df_tr$i]
  pairs_df_tr$Z_j <- obs_train$Z_i[pairs_df_tr$j]
  pairs_df_tr$cov_product <- (pairs_df_tr$Z_i - mu_Z_train) * (pairs_df_tr$Z_j - mu_Z_train)
  
  pairs_df_tr$bin_factor <- cut(pairs_df_tr$dist, breaks = breaks_dist, include.lowest = TRUE)
  
  empirical_cov_tr <- pairs_df_tr %>%
    group_by(bin_factor) %>%
    summarise(
      mean_dist = mean(dist, na.rm=TRUE),           
      cov_val = mean(cov_product, na.rm=TRUE),      
      n_pairs = n()                     
    ) %>%
    filter(n_pairs >= 10)
  
  mod_loess_tr <- loess(cov_val ~ mean_dist, data = empirical_cov_tr, weights = n_pairs, span = 0.5)
  
  # --- 3. RICOSTRUZIONE MATRICE ATTIVA E NEAR-PD ---
  C_attiva_vec_tr <- predict(mod_loess_tr, newdata = data.frame(mean_dist = as.vector(dist_celle)))
  C_attiva_vec_tr[as.vector(dist_celle) > max_dist_loess] <- 1e-8
  C_attiva_vec_tr[is.na(C_attiva_vec_tr)] <- 1e-8
  
  C_attiva_mat_tr <- matrix(C_attiva_vec_tr, nrow = n_attive, ncol = n_attive)
  diag(C_attiva_mat_tr) <- var_h0
  
  C_attiva_PD_tr <- as.matrix(nearPD(C_attiva_mat_tr, corr = FALSE, base.matrix = TRUE)$mat)
  
  coppie_celle$cov_val_tr <- as.vector(C_attiva_PD_tr)
  C_griglia_sparsa_tr <- sparseMatrix(
    i = coppie_celle$cell_id_1,
    j = coppie_celle$cell_id_2,
    x = coppie_celle$cov_val_tr,
    dims = c(K_total, K_total)
  )
  
  # --- 4. RICOSTRUZIONE LBFGS SUL FOLD ---
  C_train_cv <- as.matrix(C_griglia_sparsa_tr)
  
  cell_test <- obs$cell_id[idx_test_df]
  C_train_cv[cell_test, ] <- 0
  C_train_cv[, cell_test] <- 0
  
  covar <- C_train_cv 
  
  res_cv_fold <- covar_reconstruction_2D(
    Rn = C_train_cv, 
    alpha = alpha_best,  
    m = m_best,
    rk = rank_param
  )
  C_fold_ricostruita <- res_cv_fold$cov
  
  # --- 5. PREDIZIONE KRIGING ---
  cell_train <- obs$cell_id[idx_train_df]
  cell_test  <- obs$cell_id[idx_test_df]
  Z_train    <- obs$Z_i[idx_train_df]
  
  Sigma_train_train <- C_fold_ricostruita[cell_train, cell_train]
  Sigma_test_train  <- C_fold_ricostruita[cell_test, cell_train]
  
  Sigma_train_inv <- ginv(Sigma_train_train)
  pesi_cv <- Sigma_test_train %*% Sigma_train_inv
  
  Z_hat_local <- mu_Z_train + pesi_cv %*% (Z_train - mu_Z_train)
  
  # Restituisce un dataframe con indici originali e predizioni
  data.frame(pos = idx_test_df, val = as.vector(Z_hat_local))
}

# Ricostruisce l'array Z_hat_cv aggregato fuori dal ciclo parallelo
for (i in 1:nrow(risultati_cv_list)) {
  Z_hat_cv[risultati_cv_list$pos[i]] <- risultati_cv_list$val[i]
}

stopCluster(cl)

# Inseriamo i risultati nel dataframe e calcoliamo l'RMSE
obs$Z_kriging_cv <- Z_hat_cv
rmse_kriging_loess <- sqrt(mean((obs$Z_i - obs$Z_kriging_cv)^2))

cat("\n========================================\n")
cat(" VERDETTO CROSS-VALIDATION RIGOROSA (RMSE) \n")
cat("========================================\n")
cat(sprintf(" RMSE Kriging (LOESS): %.5f\n", rmse_kriging_loess))
cat("========================================\n")

################################################################################
########################## SALVATAGGIO PER PONTE ###############################
################################################################################

# 1. DEFINIZIONE CARTELLA DI EXPORT COMUNE
#cartella_export <- "C:/Users/miche/OneDrive/Desktop/prova3/INGV - pgv"
cartella_export <- "/work/u10700026/export_pgv"

if(!dir.exists(cartella_export)) {
  dir.create(cartella_export, recursive = TRUE)
}

# 2. IMPACCHETTAMENTO
# --- SALVATAGGIO RISULTATI LOESS PER IL PONTE ---
risultati_nonparam <- list(
  # 1. Mappa con l'aggiunta vitale del cell_id
  mappa = grid_square[, c("cell_id", "X_km", "Y_km", "dS2S_kriging", "se_kriging")], 
  
  # 2. Variabili classiche
  cov = Sigma_11, 
  stazioni = obs$NET_STA, 
  rmse = rmse_kriging_loess,
  
  # 3. AGGIUNTE PER IL CONFRONTO SPAZIALE (FIGURA 2)
  cov_globale = C_globale_ricostruita, 
  stazioni_obs = obs[, c("NET_STA", "X_km", "Y_km", "cell_id", "st_lon", "st_lat")]
)

# 3. SALVATAGGIO
path_export <- file.path(cartella_export, "Risultati_NonParam.rds")
saveRDS(risultati_nonparam, file = path_export)
print(paste("Risultati NonParam (LOESS) esportati con successo in:", path_export))
