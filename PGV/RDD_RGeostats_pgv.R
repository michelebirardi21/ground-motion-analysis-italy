
# Modello RDD

# P = 0.5 -> cost = 3.0
# P = 1.0 -> cost = 4.0
# P = 1.5 -> cost = 4.74
# P = 2.0 -> cost = 5.37
# nel caso, cambiare: P in equal, cost in get_param e nu in matern_ns_corr

library(tidyr)
library(maps)
library(sp)
library(ggpubr)
library(reticulate)
library(ggplot2)
library(plotly)
library(ggforce)
library(viridis)
library(gridExtra)
library(rstan)
library(foreach)
library(doParallel)
library(SpatialTools)
#options(mc.cores = parallel::detectCores(), Ncpus = parallel::detectCores())
options(mc.cores = 6, Ncpus = 6)
rstan_options(auto_write = TRUE)

#setwd("C:/Users/miche/OneDrive/Desktop/prova3/Scripts")
#source("Utilities.r")
#setwd("C:/Users/miche/OneDrive/Desktop/prova3/INGV - pgv")
#source("ITA_18.R")
#data = read.csv("dataset_martina.csv")


# versione per cluster

setwd("/work/u10700026/")
source("Utilities.R")
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

# res = offset + dBe|evento + dWs
# quindi dWs contiene ancora degli effetti di stazione

library(lme4)
library(dplyr)
library(sf)

# ==============================================================================
# 0. PREPARAZIONE DATI: ESTRAZIONE BLUPs (dS2S) DAL MODELLO MISTO
# ==============================================================================

# 1. Fit del modello incrociato (stima simultanea degli effetti di evento e di stazione)
mod_completo <- lmer(res ~ 1 + (1 | EVENT) + (1 | NET_STA), data = itacentrale, REML = TRUE)

# 2. Estrazione dei BLUPs per le stazioni (Questo è il nostro dS2S target)
dS2S_blup <- ranef(mod_completo)$NET_STA
dS2S_blup$NET_STA <- rownames(dS2S_blup)
colnames(dS2S_blup)[1] <- "dS2S"

# 3. --- PROIEZIONE COORDINATE --- (Mancava questo pezzo!)
# Convertiamo prima le st_lon e st_lat in X_km e Y_km
itacentrale_sf <- st_as_sf(itacentrale, coords = c("st_lon", "st_lat"), crs = 4326, remove = FALSE)
itacentrale_utm <- st_transform(itacentrale_sf, 32633)
coords_mat <- st_coordinates(itacentrale_utm)
itacentrale$X_km <- coords_mat[,1] / 1000
itacentrale$Y_km <- coords_mat[,2] / 1000

# 4. Aggregazione delle coordinate uniche per stazione
stazioni_coords <- itacentrale %>%
  group_by(NET_STA) %>%
  summarize(
    X_km = mean(X_km), 
    Y_km = mean(Y_km)
  ) %>%
  ungroup()

# 5. Unione delle coordinate con i BLUPs
spatial_data_df <- merge(stazioni_coords, dS2S_blup, by = "NET_STA")

# 6. Preparazione input finali per RDD_fit
rdd_coords <- spatial_data_df[, c("X_km", "Y_km")]
rdd_data   <- as.matrix(spatial_data_df$dS2S) # Vettore colonna

cat("Dimensioni rdd_data:", nrow(rdd_data), "stazioni\n")
# ==============================================================================

# --- 2. PREPARAZIONE GRIGLIA PER RDD (ALLINEATA CON LOESS) ---
print("Generazione griglia spaziale (Quadrato perfetto)...")

step_size_km <- 7
buffer_km <- 20

min_X <- min(itacentrale$X_km) - buffer_km
max_X <- max(itacentrale$X_km) + buffer_km
min_Y <- min(itacentrale$Y_km) - buffer_km
max_Y <- max(itacentrale$Y_km) + buffer_km

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

grid_plain <- expand.grid(col_idx = 0:(N_side-1), row_idx = 0:(N_side-1))
grid_plain$X_km <- seq_X[grid_plain$col_idx + 1]
grid_plain$Y_km <- seq_Y[grid_plain$row_idx + 1]

# Ordinamento identico alla serpentina del LOESS per garantire lo stesso ordine di righe
grid_plain$serpentine_idx <- ifelse(
  grid_plain$row_idx %% 2 == 0,
  grid_plain$row_idx * N_side + grid_plain$col_idx + 1,          
  grid_plain$row_idx * N_side + (N_side - 1 - grid_plain$col_idx) + 1 
)
grid_plain <- grid_plain[order(grid_plain$serpentine_idx), ]

# Nomi attesi da RDD_fit
grid_plain$x_1 <- grid_plain$X_km
grid_plain$x_2 <- grid_plain$Y_km

print(paste("Punti griglia generati:", nrow(grid_plain)))
# NOTA BENE: Nessun filtro terraferma applicato. Lo faremo nei grafici finali.

# --- CONFIGURAZIONE RDD ---

B = 1
dir_path = "RDD_output_RGeostats_dS2S"
if(!dir.exists(dir_path)) dir.create(dir_path)
writeLines(c(""), paste0(dir_path, "/log.txt"))
rdd_data <- as.matrix(rdd_data)
Ks = c(1)
varianza_totale <- var(spatial_data_df$dS2S, na.rm = TRUE)

# ultima ipotesi: M2V = nugget, R = range, P = smoothing parameter, V = partial sill

param_init   = c("V=0.01", "R=110", "M2V=0.03")
lower_bounds = c("V=0.001", "R=5", "M2V=0.01")
upper_bounds = c("V=0.02", "R=400", "M2V=0.05")
#equal_bounds = c("P=0.5")

library(sf)
library(rnaturalearth)
ita_border <- ne_countries(scale = "medium", country = "Italy", returnclass = "sf")
ita_border_utm <- st_transform(ita_border, 32633)

# ------------------------------------------------------------------------------
# -------------- NUOVO FILTRO CENTRI: TERRAFERMA + LIMITI NORD/SUD -------------
# ------------------------------------------------------------------------------
print("Filtraggio dei centri candidati (Terraferma + Limiti dati)...")

# 1. Trova i limiti dalle stazioni (con 2 km di buffer di sicurezza)
limite_sud <- min(rdd_coords$Y_km) - 20
limite_nord <- max(rdd_coords$Y_km) + 20

# 2. Taglio Nord-Sud sulla griglia quadrata (in km)
grid_tagliata <- subset(grid_plain, Y_km >= limite_sud & Y_km <= limite_nord)

# 3. Taglio Terraferma: Costruzione sicura dell'oggetto spaziale in METRI
# Moltiplichiamo prima le colonne, SENZA toccare le geometrie
grid_tagliata_m <- grid_tagliata
grid_tagliata_m$X_m <- grid_tagliata_m$X_km * 1000
grid_tagliata_m$Y_m <- grid_tagliata_m$Y_km * 1000

# Creiamo l'oggetto sf già in metri e gli assegniamo esattamente il CRS dell'Italia
grid_sf <- st_as_sf(grid_tagliata_m, coords = c("X_m", "Y_m"), crs = st_crs(ita_border_utm))

# 4. Intersezione ultra-sicura con lengths()
# st_intersects restituisce una lista. lengths() conta quanti poligoni dell'Italia 
# toccano ogni punto. Se > 0, significa che il punto è su suolo italiano.
punti_in_italia_logico <- lengths(st_intersects(grid_sf, ita_border_utm)) > 0

# 5. Creiamo la griglia definitiva da cui pescare i centri
center_grid_ideale <- grid_tagliata[punti_in_italia_logico, ]

cat("Punti griglia quadrata totale (mantenuta per predizione):", nrow(grid_plain), "\n")
cat("Centri candidati validi (usati per pescare in RDD_fit):", nrow(center_grid_ideale), "\n")
# ------------------------------------------------------------------------------

# Ciclo principale
for (K in Ks) {
  
  print(paste(">>> Avvio RDD fit per K =", K))
  
  # 1. Fit del modello
  estimates = RDD_fit(
    data = rdd_data, 
    coords = rdd_coords,      
    center_grid = center_grid_ideale[, c("X_km", "Y_km")], 
    K = K, 
    B = B,
    clusters = 6,            
    dir_path = dir_path,
    struct = c("Exponential", "Nugget Effect"), 
    dirs = seq(0, 135, by = 45),
    #dirs = seq(0, 150, by = 30),
    #   param =  c("P=2"), 
    #   lower =  c("P=2"),  
    #   upper =  c("P=2"),
    param = param_init, 
    lower = lower_bounds,  
    upper = upper_bounds,
    #equal = equal_bounds,
    threshold = 10          
  )
  
  print("Fit completato. Calcolo predizioni...")
  
  # 2. Predizione
  ret = RDD_predict(estimates, grid_plain[, c("X_km", "Y_km")])
  ret = simplify2array(ret)
  
  # 3. Estrazione statistiche (Median)
  RDD_estimate = extract_spatial_predictions_median(ret, B = B)
  
  # Unione con la griglia originale (KM)
  plot_data = cbind(grid_plain, RDD_estimate)
  plot_data = as.data.frame(plot_data)
  
  # alias x_1/x_2 per le funzioni che usano questi nomi
  plot_data$x_1 <- plot_data$X_km
  plot_data$x_2 <- plot_data$Y_km
  
  # (opzionale ma utile) converti i campi in numeric per evitare tipi logical/factor
  cols_num <- c("mu","sigma","lambda_1","lambda_2","theta_deg","nugget")
  for(col in intersect(cols_num, names(plot_data))) {
    plot_data[[col]] <- as.numeric(as.character(plot_data[[col]]))
  }
  
  # =========================================================
  # CONTROLLO CONGELAMENTO PARAMETRI (FROZEN BOUNDS CHECK)
  # =========================================================
  cat("\n--- VERIFICA OTTIMIZZAZIONE PARAMETRI (K =", K, ") ---\n")
  
  # 1. Partial Sill (S). Nel dataframe abbiamo 'sigma' che è la radice quadrata di S
  S_stimato <- plot_data$sigma^2
  cat("Partial Sill (S) \n")
  print(summary(S_stimato))
  
  # 2. Nugget (V). Nel dataframe è salvato direttamente come 'nugget'
  if("nugget" %in% names(plot_data)) {
    cat("\nNugget (V) \n")
    print(summary(plot_data$nugget))
  }
  
  # 3. Range (R). Nel dataframe abbiamo lambda_1 e lambda_2. 
  # In Utilities.r, lambda = (Range / 3.0)^2. Invertiamo la formula per tornare al Range in km:
  Range_stimato <- sqrt(pmax(plot_data$lambda_1, plot_data$lambda_2)) * 3.0
  cat("\nRange (R) in km \n")
  print(summary(Range_stimato))
  
  
  # 1. Calcolo i tre parametri con il loro vero significato fisico
  S_stimato <- plot_data$sigma^2
  Range_stimato <- sqrt(pmax(plot_data$lambda_1, plot_data$lambda_2)) * 3.0
  Nugget_stimato <- if("nugget" %in% names(plot_data)) plot_data$nugget else rep(NA, nrow(plot_data))
  
  # 2. Creo un dataframe temporaneo per il plot
  check_df <- data.frame(
    Partial_Sill = plot_data$sigma^2,
    Range_MajorAxis_km = sqrt(plot_data$lambda_1) * 3.0,         # r1, asse maggiore
    Range_MinorAxis_km = sqrt(plot_data$lambda_2) * 3.0,         # r2, asse minore
    Aniso_Ratio = sqrt(plot_data$lambda_1 / plot_data$lambda_2), # phi = r1/r2
    Nugget = plot_data$nugget
  )
  
  # 3. Trasformo il dataframe nel formato "lungo" (ideale per facet_wrap)
  check_long <- gather(check_df, key = "Parametro", value = "Valore")
  
  # 4. Genero gli istogrammi
  p_check <- ggplot(check_long, aes(x = Valore)) +
    geom_histogram(bins = 40, fill = "steelblue", color = "black", alpha = 0.8) +
    facet_wrap(~ Parametro, scales = "free", ncol = 3) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste("Distribuzione Parametri Stimati (K =", K, ")"),
      subtitle = "Check visuale per l'effetto dei frozen bounds",
      x = "Valore stimato", 
      y = "Frequenza (N. celle griglia)"
    ) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  print(p_check)
  
  cat("---------------------------------------------------------\n\n")
  
  # --- CONVERSIONE COORDINATE (KM -> Metri -> Lat/Lon) ---
  # Fondamentale per far funzionare plot_parameter_heatmaps
  
  # Moltiplichiamo per 1000 per avere i Metri (UTM 33N richiede metri)
  plot_data$X_m <- plot_data$X_km * 1000
  plot_data$Y_m <- plot_data$Y_km * 1000
  
  # Creiamo oggetto spaziale temporaneo
  plot_data_sf <- st_as_sf(plot_data, coords = c("X_m", "Y_m"), crs = 32633) # UTM 33N
  
  # Trasformiamo in WGS84 (Lat/Lon)
  plot_data_wgs84 <- st_transform(plot_data_sf, 4326)
  
  # Estraiamo e assegniamo le nuove colonne
  coords_geo <- st_coordinates(plot_data_wgs84)
  plot_data$Longitude <- coords_geo[,1]
  plot_data$Latitude  <- coords_geo[,2]
  
  # --- GENERAZIONE GRAFICI ---
  print("Generazione mappe parametri...")
  
  # Chiamata alla funzione modificata
  # Try-catch per evitare che un errore grafico blocchi tutto il ciclo
  tryCatch({
    plot_parameters(plot_data, dir_path = dir_path, K = K, borders_sf = ita_border_utm)
  }, error = function(e) {
    print(paste("Errore durante il plotting per K =", K, ":", e$message))
  })
  
  print(paste("--- Completato K =", K, "---"))
}

################################################################################
###################### MATRICE COVARIANZA STAZIONI #############################
################################################################################

library(dplyr)

print(">>> Calcolo Matrice Covarianza Stazioni Uniche (Versione Safe)...")

# --- 1. Recupero e Pulizia Dati ---

# Cerchiamo il dataset sorgente (dallo script o dal bundle)
if(exists("spatial_data_df")) {
  source_df <- spatial_data_df
} else if(exists("bundle")) {
  source_df <- bundle$input_data$coords
  # Se mancano i nomi nel bundle, creiamo ID fittizi basati sulle coordinate
  if(!"NET_STA" %in% names(source_df)) {
    source_df$NET_STA <- paste0("ST_", round(source_df$X_km, 3), "_", round(source_df$Y_km, 3))
  }
} else {
  stop("Errore: Impossibile trovare i dati delle stazioni (spatial_data_df o bundle mancanti)")
}

# DEDUPLICAZIONE ROBUSTA:
# 1. Rimuove duplicati di nome (stessa stazione, eventi diversi)
# 2. Rimuove duplicati di coordinate (stazioni diverse ma sovrapposte, che causerebbero NaN)
unique_stations_df <- source_df %>%
  select(NET_STA, X_km, Y_km) %>%
  distinct(NET_STA, .keep_all = TRUE) %>%     # Unica riga per nome stazione
  distinct(X_km, Y_km, .keep_all = TRUE)      # Unica riga per coordinate (SAFETY CHECK)

print(paste("Numero stazioni fisiche uniche:", nrow(unique_stations_df)))

# --- 2. Predizione Parametri (Median) ---

# Input coordinate
stations_coords_input <- unique_stations_df[, c("X_km", "Y_km")]

# Selezioniamo il modello corretto (da memoria o bundle)
model_to_use <- if(exists("estimates")) estimates else bundle$estimates
B_val <- if(exists("B")) B else 140

# Predizione
ret_stations <- RDD_predict(model_to_use, stations_coords_input)
ret_stations_arr <- simplify2array(ret_stations)
preds_stations <- extract_spatial_predictions_median(ret_stations_arr, B = B_val)

# Assegniamo i nomi delle stazioni ai risultati per chiarezza
rownames(preds_stations) <- unique_stations_df$NET_STA

# --- 3. Preparazione Input Covarianza ---

# Theta (gestione gradi/radianti)
# Compute_Aniso in Utilities.r si aspetta RADIANTI (usa cos/sin)
theta_st <- if("theta" %in% names(preds_stations)) preds_stations$theta else preds_stations$theta_deg * (pi/180)

l1_st <- preds_stations$lambda_1
l2_st <- preds_stations$lambda_2
sigma_st <- preds_stations$sigma

# Nugget (se presente nel modello, altrimenti 0)
nugget_st <- if("nugget" %in% names(preds_stations)) preds_stations$nugget else rep(0, nrow(preds_stations))

# --- 4. Costruzione Matrice (Usando Utilities.r NATIVE) ---

# Importante: Usiamo Compute_Aniso definita in Utilities.r per coerenza con RGeostats
# (Non ridefiniamo funzioni custom che potrebbero invertire la rotazione)
if(!exists("Compute_Aniso")) source("Utilities.r") 

Aniso_List_Stations <- Compute_Aniso(l1_st, l2_st, theta_st)

# Calcolo Covarianza (Parte Spaziale)
coords_st_t <- t(as.matrix(stations_coords_input))
dets_st <- l1_st * l2_st

# Check funzione Bessel (se non caricata da Utilities)
if (!exists("modified_bessel_second_kind")) {
  modified_bessel_second_kind <- function(nu, x) besselK(x, nu)
}

Cov_Mat_Stations <- matern_ns_corr(
  Locations = coords_st_t, 
  Aniso     = Aniso_List_Stations, 
  dets      = dets_st, 
  nu        = 0.5, 
  sigma     = sigma_st
)

# --- 5. Aggiunta Nugget e Formattazione ---

# Aggiungiamo il Nugget SOLO sulla diagonale
diag(Cov_Mat_Stations) <- diag(Cov_Mat_Stations) + nugget_st

# Assegniamo i nomi a righe e colonne
rownames(Cov_Mat_Stations) <- unique_stations_df$NET_STA
colnames(Cov_Mat_Stations) <- unique_stations_df$NET_STA

print("Calcolo completato con successo.")
print(Cov_Mat_Stations[1:5, 1:5])

library(ggplot2)
library(reshape2)

# 1. DEFINIZIONE DELL'ORDINE GEOGRAFICO (Nord -> Sud)
#    Ordiniamo i nomi delle stazioni in base alla coordinata Y decrescente
ordine_geo <- unique_stations_df$NET_STA[order(unique_stations_df$Y_km, decreasing = TRUE)]

# 2. PREPARAZIONE DATI (Melting)
cov_melted <- melt(Cov_Mat_Stations)

# 3. APPLICAZIONE DELL'ORDINE
#    Trasformiamo le colonne Var1 e Var2 in Fattori con l'ordine specifico definito sopra
cov_melted$Var1 <- factor(cov_melted$Var1, levels = ordine_geo)
cov_melted$Var2 <- factor(cov_melted$Var2, levels = ordine_geo)

# 4. PLOT (Senza nomi)
p_cov_geo <- ggplot(cov_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  
  # Scala di colore per la COVARIANZA (Viridis è ottimo per valori assoluti/intensità)
  scale_fill_viridis_c(option = "magma", name = "Covarianza") +
  
  # Etichette generiche (dato che nascondiamo i nomi)
  labs(title = "Matrice di Covarianza (Ordinata N -> S)",
       x = "Stazioni (Nord -> Sud)", 
       y = "Stazioni (Nord -> Sud)") +
  
  theme_minimal() +
  coord_fixed() +
  theme(
    # QUESTI COMANDI NASCONDONO I NOMI:
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    
    # Pulizia generale
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face="bold")
  )

print(p_cov_geo)

################################################################################
############################# CROSS-VALIDATION #################################
################################################################################

# ==============================================================================
# 1. SETUP E SUBSETTING
# ==============================================================================
library(dplyr)
library(lme4)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(ggpubr)

print(">>> AVVIO SETUP CROSS-VALIDATION RIGOROSA (NO LEAKAGE) <<<")

# --- A. DEFINIZIONE FRAZIONE DATASET (Debug vs Full Run) ---
DATA_FRACTION <- 1  # 0.05 per test veloce, 1.0 per run completa

set.seed(123) 
if (DATA_FRACTION < 1.0) {
  n_sub <- round(nrow(itacentrale) * DATA_FRACTION)
  cat(paste("!!! ATTENZIONE: Esecuzione su subset di", n_sub, "righe (", DATA_FRACTION*100, "%) !!!\n"))
  idx_sub <- sample(1:nrow(itacentrale), n_sub)
  df_work <- itacentrale[idx_sub, ]
} else {
  cat(">>> Esecuzione su DATASET COMPLETO <<<\n")
  df_work <- itacentrale
}

if (!("res" %in% names(df_work))) {
  stop("Errore: Il dataframe deve contenere la colonna 'res' (residuo totale).")
}

# --- B. Gestione Gruppi (Stazioni) su df_work ---
station_ids <- as.character(df_work$NET_STA)
unique_stations <- unique(station_ids)
n_stations <- length(unique_stations)

print(paste("Numero stazioni uniche nel dataset di lavoro:", n_stations))

# --- C. Assegnazione Fold ---
n_folds <- 10
set.seed(123) 

folds_per_station <- sample(rep(1:n_folds, length.out = n_stations))
station_fold_map <- setNames(folds_per_station, unique_stations)
fold_indices <- station_fold_map[station_ids]

cv_results_global <- data.frame()

# Parametri CV
#Ks_CV = c(1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20)
#B_CV = 240   

Ks_CV = c(1)
B_CV = 1

build_standard_aniso <- function(l1, l2, th) {
  N_pts <- length(l1)
  Anis_List <- list()
  c_th <- cos(th); s_th <- sin(th)
  for(i in 1:N_pts) {
    R <- matrix(c(c_th[i], s_th[i], -s_th[i], c_th[i]), nrow=2, ncol=2)
    L <- diag(c(l1[i], l2[i]))
    Anis_List[[i]] <- R %*% L %*% t(R)
  }
  return(Anis_List)
}

if (!exists("modified_bessel_second_kind")) {
  modified_bessel_second_kind <- function(nu, x) besselK(x, nu)
}

# ==============================================================================
# 2. CICLO DI VALIDAZIONE
# ==============================================================================

for (K in Ks_CV) {
  
  print(paste("========================================"))
  print(paste(" VALUTAZIONE MODELLO CON K =", K))
  print(paste("========================================"))
  
  results_k <- data.frame()
  
  for (i in 1:n_folds) {
    
    test_idx <- which(fold_indices == i)
    if(length(test_idx) == 0) {
      cat(paste0("[K=", K, "] Fold ", i, " vuoto, salto...\n"))
      next
    }
    
    cat(paste0("[K=", K, "] Elaborazione Fold ", i, "/", n_folds, " (n_test=", length(test_idx), ")...\n"))
    
    # --- A. Splitting ---
    train_idx <- setdiff(1:nrow(df_work), test_idx)
    df_train <- df_work[train_idx, ]
    df_test  <- df_work[test_idx, ]
    
    # --- B. Ricalcolo LMM (ISOLAMENTO TOTALE TEST SET) ---
    skip_fold <- FALSE
    tryCatch({
      # 1. Modello LMM esclusivo per il training
      mod_cv <- lmer(res ~ 1 + (1 | EVENT) + (1 | NET_STA), data = df_train, REML = TRUE)
      
      # 2. Estrazione componenti apprese (Training)
      offset_train <- fixef(mod_cv)[1]
      
      blup_event_train <- ranef(mod_cv)$EVENT
      blup_event_train$EVENT <- rownames(blup_event_train)
      colnames(blup_event_train)[1] <- "dBe"
      
      blup_sta_train <- ranef(mod_cv)$NET_STA
      blup_sta_train$NET_STA <- rownames(blup_sta_train)
      colnames(blup_sta_train)[1] <- "dS2S"
      
      # 3. Preparazione Set di Addestramento RDD
      train_s2s <- df_train %>%
        group_by(NET_STA) %>%
        summarize(X_km = mean(X_km), Y_km = mean(Y_km)) %>%
        merge(blup_sta_train, by = "NET_STA")
      
      # 4. Preparazione Set di Test RDD (Target Empirico senza Leakage)
      
      # FIX: Identifichiamo gli eventi noti dal Training Set
      eventi_noti <- unique(df_train$EVENT)
      
      # Filtriamo il Test Set per tenere solo gli eventi noti, poi calcoliamo il residuo pulito
      df_test_clean <- df_test %>%
        filter(EVENT %in% eventi_noti) %>% # <-- FILTRO FONDAMENTALE
        left_join(blup_event_train, by = "EVENT") %>%
        mutate(
          empirical_dS2S = res - offset_train - dBe
        )
      
      # Se dopo il filtro il test set per le stazioni è vuoto, saltiamo
      if(nrow(df_test_clean) == 0) {
        skip_fold <<- TRUE
        next
      }
      test_s2s <- df_test_clean %>%
        group_by(NET_STA) %>%
        summarize(X_km = mean(X_km), 
                  Y_km = mean(Y_km), 
                  dS2S_target = mean(empirical_dS2S, na.rm = TRUE)) 
      
      # Assegnazione input RDD
      train_coords <- train_s2s[, c("X_km", "Y_km")]
      test_coords  <- test_s2s[, c("X_km", "Y_km")]
      train_data_vec <- as.numeric(train_s2s$dS2S)
      test_data_vec  <- as.numeric(test_s2s$dS2S_target) # Il "vero" osservato pulito
      
    }, error = function(e) {
      cat(paste("Errore LMM nel fold", i, ":", e$message, "- Salto Fold\n"))
      skip_fold <<- TRUE
    })
    
    if(skip_fold) next
    
    # --- C. Fitting RDD (Sul Training) ---
    fit_cv <- RDD_fit(
      data = as.matrix(train_data_vec), 
      coords = train_coords,      
      center_grid = center_grid_ideale[, c("X_km", "Y_km")], 
      K = K, 
      B = B_CV,
      clusters = 6,              
      dir_path = dir_path,
      struct = c("Exponential", "Nugget Effect"), 
      dirs = seq(0, 135, by = 45),
      param = param_init, lower = lower_bounds, upper = upper_bounds,
      threshold = 10           
    )
    
    # --- D. Kriging a Intorno Mobile (ROBUSTO E OTTIMIZZATO) ---
    n_max_neighbors <- 100
    
    all_coords_ordered <- rbind(train_coords, test_coords)
    preds_raw_all <- RDD_predict(fit_cv, all_coords_ordered)
    preds_all <- extract_spatial_predictions_median(simplify2array(preds_raw_all), B = B_CV)
    
    train_coords_mat <- as.matrix(train_coords)
    test_coords_mat  <- as.matrix(test_coords)
    
    n_tr <- nrow(train_coords)
    n_te <- nrow(test_coords)
    
    l1_tr <- preds_all$lambda_1[1:n_tr]
    l2_tr <- preds_all$lambda_2[1:n_tr]
    theta_tr <- if("theta" %in% names(preds_all)) preds_all$theta[1:n_tr] else preds_all$theta_deg[1:n_tr] * (pi/180)
    sigma_tr <- preds_all$sigma[1:n_tr]
    nugget_tr <- if(is.null(preds_all$nugget)) rep(0, n_tr) else preds_all$nugget[1:n_tr]
    mu_tr <- preds_all$mu[1:n_tr]
    res_tr <- train_data_vec - mu_tr 
    
    idx_te_start <- n_tr + 1
    idx_te_end   <- n_tr + n_te
    l1_te <- preds_all$lambda_1[idx_te_start:idx_te_end]
    l2_te <- preds_all$lambda_2[idx_te_start:idx_te_end]
    theta_te <- if("theta" %in% names(preds_all)) preds_all$theta[idx_te_start:idx_te_end] else preds_all$theta_deg[idx_te_start:idx_te_end] * (pi/180)
    sigma_te <- preds_all$sigma[idx_te_start:idx_te_end]
    mu_te <- preds_all$mu[idx_te_start:idx_te_end]
    
    nugget_te <- if(is.null(preds_all$nugget)) rep(0, n_te) else preds_all$nugget[idx_te_start:idx_te_end]
    Z_hat_kriging <- numeric(n_te)
    SE_kriging <- numeric(n_te)
    
    # 1. Identifichiamo i punti di training validi una volta sola
    valid_tr_idx <- which(!is.na(l1_tr) & !is.na(res_tr) & !is.na(sigma_tr))
    n_valid_tr <- length(valid_tr_idx)
    
    # Pre-filtriamo i vettori tenendo solo i dati validi per questo fold
    train_coords_valid <- train_coords_mat[valid_tr_idx, , drop = FALSE]
    l1_tr_valid <- l1_tr[valid_tr_idx]
    l2_tr_valid <- l2_tr[valid_tr_idx]
    th_tr_valid <- theta_tr[valid_tr_idx]
    sig_tr_valid <- sigma_tr[valid_tr_idx]
    nugget_tr_valid <- nugget_tr[valid_tr_idx]
    res_tr_valid <- res_tr[valid_tr_idx]
    
    # --- LOOP SUI PUNTI DI TEST ---
    for(j in 1:n_te) {
      if(j %% 10 == 0) cat("   -> Kriging stazione", j, "su", n_te, "\n")
      
      # Se mancano parametri nel test o ci sono troppi pochi dati di train, usiamo solo la media
      if(is.na(l1_te[j]) || is.na(sigma_te[j]) || n_valid_tr < 3) {
        Z_hat_kriging[j] <- mu_te[j]
        next
      }
      
      risultato_kriging <- tryCatch({
        
        # --- 1. TROVA I VICINI PIÙ PROSSIMI (INTORNO MOBILE) ---
        dist_sq <- (train_coords_valid[,1] - test_coords_mat[j,1])^2 + 
          (train_coords_valid[,2] - test_coords_mat[j,2])^2
        
        n_neighbors <- min(n_max_neighbors, n_valid_tr)
        idx_vicini <- order(dist_sq)[1:n_neighbors]
        
        # --- 2. CREA IL SOTTOINSIEME LOCALE ---
        loc_sub <- rbind(train_coords_valid[idx_vicini, , drop=FALSE], test_coords_mat[j, , drop=FALSE])
        
        l1_sub  <- c(l1_tr_valid[idx_vicini], l1_te[j])
        l2_sub  <- c(l2_tr_valid[idx_vicini], l2_te[j])
        th_sub  <- c(th_tr_valid[idx_vicini], theta_te[j])
        sig_sub <- c(sig_tr_valid[idx_vicini], sigma_te[j])
        nugget_sub <- c(nugget_tr_valid[idx_vicini])
        res_sub <- res_tr_valid[idx_vicini]
        
        # --- 3. CALCOLA LA COVARIANZA SOLO PER QUESTO PICCOLO INTORNO ---
        Aniso_List_Sub <- Compute_Aniso(l1_sub, l2_sub, th_sub)
        
        Cov_Sub <- matern_ns_corr(t(loc_sub), Aniso_List_Sub, l1_sub * l2_sub, 0.5, sig_sub)
        
        # --- 4. ESTRAI MATRICE DI TRAINING LOCALE E CROSS-COVARIANZA ---
        C_nn_locale <- Cov_Sub[1:n_neighbors, 1:n_neighbors, drop = FALSE]
        diag(C_nn_locale) <- diag(C_nn_locale) + nugget_sub + 1e-6*max(diag(C_nn_locale)) # Aggiunta nugget locale
        
        c_nt <- Cov_Sub[1:n_neighbors, n_neighbors + 1]
        
        if(any(is.na(c_nt))) stop("NaN in Cross-Covariance")
        
        # --- 5. RISOLVI IL SISTEMA (Kriging Weights) ---
        U_locale <- chol(C_nn_locale)
        weights_loc <- backsolve(U_locale, forwardsolve(t(U_locale), c_nt))
        
        correction <- sum(weights_loc * res_sub)
        
        if(is.na(correction)) stop("Correction is NA")
        
        # --- 6. CALCOLO DELLA VARIANZA DI KRIGING ---
        var_a_priori <- sigma_te[j]^2 + nugget_te[j]
        riduzione_var <- sum(weights_loc * c_nt)
        se_val <- sqrt(max(0, var_a_priori - riduzione_var))
        
        # Restituisce una lista con predizione e incertezza
        list(pred = mu_te[j] + correction, se = se_val)
        
      }, error = function(e) {
        
        var_a_priori <- sigma_te[j]^2 + nugget_te[j]
        list(pred = mu_te[j], se = sqrt(var_a_priori))
      })
      
      Z_hat_kriging[j] <- risultato_kriging$pred
      SE_kriging[j]    <- risultato_kriging$se
      
    }
    
    # --- E. Salvataggio ---
    fold_res <- data.frame(
      K_poly = K,
      Fold = i,
      Observed = test_data_vec,
      Predicted_MeanOnly = mu_te,
      Predicted_Kriging = Z_hat_kriging,
      SE_Kriging = SE_kriging
    )
    
    results_k <- rbind(results_k, fold_res)
  }
  
  if(nrow(results_k) > 0) {
    cv_results_global <- rbind(cv_results_global, results_k)
    rmse_kr <- sqrt(mean((results_k$Observed - results_k$Predicted_Kriging)^2, na.rm=T))
    print(paste(">>> FINE K =", K, "| RMSE Kriging:", round(rmse_kr, 4)))
  }
}

# ==============================================================================
# 3. ANALISI FINALE
# ==============================================================================

if(nrow(cv_results_global) > 0) {
  # 1. Calcolo metriche standard
  # Step 1: Calcolo l'RMSE per ogni K e per ogni Fold
  fold_summary <- cv_results_global %>%
    group_by(K_poly, Fold) %>%
    summarise(
      Fold_RMSE_MeanOnly = sqrt(mean((Observed - Predicted_MeanOnly)^2, na.rm=T)),
      Fold_RMSE_Kriging  = sqrt(mean((Observed - Predicted_Kriging)^2, na.rm=T)),
      .groups = "drop"
    )
  
  # Step 2: Calcolo media e deviazione standard degli RMSE tra i fold
  summary_table <- fold_summary %>%
    group_by(K_poly) %>%
    summarise(
      RMSE_MeanOnly    = mean(Fold_RMSE_MeanOnly, na.rm=T),
      RMSE_MeanOnly_sd = sd(Fold_RMSE_MeanOnly, na.rm=T),
      RMSE_Kriging     = mean(Fold_RMSE_Kriging, na.rm=T),
      RMSE_Kriging_sd  = sd(Fold_RMSE_Kriging, na.rm=T)
    )
  
  # 2. Trasformazione nel formato "long"
  summary_long <- summary_table %>%
    pivot_longer(
      cols = c(RMSE_MeanOnly, RMSE_Kriging),
      names_to = "Metodo",
      values_to = "Valore"
    ) %>%
    mutate(
      # Assegniamo la SD corretta in base al metodo
      Valore_sd = case_when(
        Metodo == "RMSE_MeanOnly" ~ RMSE_MeanOnly_sd,
        Metodo == "RMSE_Kriging"  ~ RMSE_Kriging_sd
      ),
      # Rinominiamo il metodo per il plot finale
      Metodo = case_when(
        Metodo == "RMSE_MeanOnly" ~ "RMSE Mean",
        Metodo == "RMSE_Kriging"  ~ "RMSE Kriging"
      )
    )
  
  # 3. Generazione del plot
  p_res <- ggplot(summary_long, aes(x = K_poly, y = Valore, color = Metodo)) +
    
    # 1. Linea
    geom_line(linewidth = 0.6) +
    
    # 2. Barre d'errore (Nere e trasparenti, sotto il pallino)
    geom_errorbar(aes(ymin = Valore - Valore_sd, ymax = Valore + Valore_sd),
                  width = 0.2, color = "black", alpha = 0.6) +
    
    # 3. Pallini (Disegnati per ultimi così restano in primo piano!)
    geom_point(size = 1) +
    
    scale_color_manual(values = c(
      "RMSE Mean"    = "red", 
      "RMSE Kriging" = "blue"
    )) +
    
    scale_x_continuous(breaks = unique(summary_long$K_poly)) +
    
    labs(
      x = "K",
      y = "RMSE", 
      color = "Method"
    ) +
    
    theme_minimal() +
    theme(
      legend.position = "bottom",
      
      # Font maggiorati come da tua richiesta precedente
      legend.title = element_text(size = 14),
      legend.text  = element_text(size = 12),
      
      axis.title.x = element_text(size = 14),
      axis.text.x  = element_text(size = 12),
      
      axis.title.y = element_text(size = 14),
      axis.text.y  = element_text(size = 12),
      
      panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
    )
  
  # 4. SALVATAGGIO
  ggsave(filename = "CV_RMSE_K_plot_pgv.png", 
         plot = p_res, 
         width = 8, 
         height = 5, 
         dpi = 300, 
         bg = "white")
  saveRDS(p_res, file = "plot_cv_pgv.rds")
} else {
  print("Nessun risultato generato (forse subset troppo piccolo?)")
}

################################################################################
########################## VARIOGRAMMA EMPIRICO ################################
################################################################################

# --- 1. RIEPILOGO VARIANZA ---
# Usiamo spatial_data_df che ora contiene le coordinate e i BLUPs (dS2S)
var_empirica <- var(spatial_data_df$dS2S, na.rm = TRUE)
cat(paste("Varianza Totale dS2S (Sill teorico max):", var_empirica, "\n"))

# --- 2. VARIOGRAMMA EMPIRICO ---
# Assicuriamoci di passarlo come data.frame pulito
db_temp <- db.create(as.data.frame(spatial_data_df[, c("X_km", "Y_km", "dS2S")]), flag.grid=FALSE)
db_temp <- db.locate(db_temp, names=c("X_km", "Y_km"), loctype="x")
db_temp <- db.locate(db_temp, names="dS2S", loctype="z")

# Parametri variogramma (proviamo a spingerci fino a 150km per vedere il plateau)
lag_step <- 50 
n_lags <- 10

# Calcolo
vario_exp <- vario.calc(db_temp, lag=lag_step, nlag=n_lags)

# Plot
plot(vario_exp, main="Variogramma Empirico Globale", pch=19, col="blue", type="b")
abline(h=var_empirica, col="red", lty=2) # Linea varianza totale

# Variogramma semplice senza RGeostats

library(gstat)
library(sf)
library(ggplot2)

# 1. Creazione dell'oggetto spaziale (sf)
# Usiamo il dataframe delle stazioni uniche con i dS2S (BLUPs)
df_vario <- spatial_data_df %>% 
  filter(!is.na(dS2S)) %>%
  st_as_sf(coords = c("X_km", "Y_km"))

# 2. Calcolo del variogramma empirico
v_empirico <- variogram(dS2S ~ 1, data = df_vario, width = 20, cutoff = 300)

# 3. Plot Standard (Semplice e pulito)
ggplot(v_empirico, aes(x = dist, y = gamma)) +
  geom_point(color = "black", size = 2) +        # Punti del variogramma
  geom_line(color = "gray", linetype = "dotted") + # Guida visiva
  geom_hline(yintercept = var(spatial_data_df$dS2S), 
             color = "red", linetype = "dashed") + # Linea della varianza totale (Sill teorico)
  labs(x = "Distanza (km)", 
       y = "Semivarianza", 
       title = "Variogramma Empirico dS2S") +
  theme_bw() # Tema bianco e nero, molto standard

# --- 1. PREPARAZIONE DATI ---
# (Assumiamo che tu abbia già il tuo spatial_data_df)
df_vario <- spatial_data_df %>% 
  filter(!is.na(dS2S)) %>%
  st_as_sf(coords = c("X_km", "Y_km"))

# --- 2. CALCOLO VARIOGRAMMA DIREZIONALE (gstat) ---
# Specificando il parametro alpha, gstat calcola il variogramma nelle direzioni indicate
# Includiamo tol.hor (tolleranza orizzontale), tipicamente si usa 22.5 per 4 direzioni (45/2)
v_empirico_dir <- variogram(dS2S ~ 1, 
                            data = df_vario, 
                            width = 20,         # Equivalente al tuo lag
                            cutoff = 300,       # Equivalente al tuo limite massimo
                            alpha = c(0, 45, 90, 135), 
                            tol.hor = 22.5)

# --- 3. PLOT DIREZIONALE (stile ggplot2) ---
# Usiamo facet_wrap per creare 4 riquadri, uno per ogni direzione, 
# per facilitare il confronto visuale
ggplot(v_empirico_dir, aes(x = dist, y = gamma)) +
  geom_point(color = "black", size = 2) +        
  geom_line(color = "gray", linetype = "dotted") + 
  
  # Aggiungiamo la linea della varianza totale in ogni riquadro
  geom_hline(yintercept = var(spatial_data_df$dS2S, na.rm = TRUE), 
             color = "red", linetype = "dashed") + 
  
  # Suddivide il plot in base alla direzione (alpha)
  facet_wrap(~ dir.hor, 
             labeller = labeller(dir.hor = function(x) paste0("Direzione: ", x, "°"))) +
  
  labs(x = "Distanza (km)", 
       y = "Semivarianza", 
       title = "Variogramma Empirico Direzionale dS2S") +
  theme_bw() +
  
  # (Opzionale) Personalizziamo i titoli dei riquadri
  theme(strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(face = "bold", size = 11))

# più bello

ggplot(v_empirico_dir, aes(x = dist, y = gamma)) +
  geom_point(color = "black", size = 2) +        
  geom_line(color = "gray", linetype = "dotted") + 
  
  # Aggiungiamo la linea della varianza totale 
  geom_hline(yintercept = var(spatial_data_df$dS2S, na.rm = TRUE), 
             color = "skyblue", linetype = "dashed") + 
  
  # Suddivide il plot in base alla direzione
  facet_wrap(~ dir.hor, 
             labeller = labeller(dir.hor = function(x) paste0("Direction: ", x, "°"))) +
  
  labs(x = "Distance (km)", 
       y = "Semivariance") +
  theme_bw() +
  
  # Personalizzazione completa dei font
  theme(
    # 1. Riquadri "Direction" (sfondo azzurro e font ingrandito)
    strip.background = element_rect(fill = "skyblue", color = "black"),
    strip.text = element_text(face = "bold", size = 14, color = "white"),
    
    # 2. Titoli degli assi "Distance (km)" e "Semivariance"
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.text.x  = element_text(size = 12, angle = 45, hjust = 1),
    
    # 3. Valori numerici sulle tacchette degli assi
    axis.text = element_text(size = 12, color = "black")
  )

################################################################################
######## GLOBAL KRIGING E MAPPA DELL'INCERTEZZA SU TUTTA LA GRIGLIA ############
################################################################################

print(">>> Avvio Global Kriging (Predizione e Incertezza su mappa)...")

# 1. PREPARAZIONE DATI DI TRAINING (Le stazioni originali)
train_df <- spatial_data_df %>% distinct(NET_STA, .keep_all = TRUE)
train_coords_mat <- as.matrix(train_df[, c("X_km", "Y_km")])
train_data_vec <- train_df$dS2S
n_tr <- nrow(train_coords_mat)

# 2. ESTRAZIONE PARAMETRI SULLE STAZIONI (Servono per la matrice di covarianza)
model_to_use <- if(exists("estimates")) estimates else bundle$estimates
B_val <- if(exists("B")) B else 140

preds_tr_raw <- RDD_predict(model_to_use, train_coords_mat)
preds_tr <- extract_spatial_predictions_median(simplify2array(preds_tr_raw), B = B_val)

l1_tr <- preds_tr$lambda_1
l2_tr <- preds_tr$lambda_2
theta_tr <- if("theta" %in% names(preds_tr)) preds_tr$theta else preds_tr$theta_deg * (pi/180)
sigma_tr <- preds_tr$sigma
nugget_tr <- if("nugget" %in% names(preds_tr)) preds_tr$nugget else rep(0, n_tr)
mu_tr <- preds_tr$mu

# Residuo spaziale sui punti di training (Z_obs - mu_locale)
res_tr <- train_data_vec - mu_tr 

# --- (Il blocco 3. Matrice Covarianza Globale di Train è stato rimosso per salvare RAM) ---

# 4. PREPARAZIONE GRIGLIA DI TEST
grid_coords_mat <- as.matrix(plot_data[, c("X_km", "Y_km")])
n_te <- nrow(grid_coords_mat)

l1_te <- plot_data$lambda_1
l2_te <- plot_data$lambda_2
theta_te <- if("theta" %in% names(plot_data)) plot_data$theta else plot_data$theta_deg * (pi/180)
sigma_te <- plot_data$sigma
mu_te <- plot_data$mu
nugget_te <- if("nugget" %in% names(plot_data)) plot_data$nugget else rep(0, n_te)

Z_hat_kriging_global <- numeric(n_te)
SE_kriging_global <- numeric(n_te)

# 5. LOOP SULLA GRIGLIA (Kriging e Calcolo Incertezza)
cat(sprintf("Calcolo Kriging su %d pixel della griglia...\n", n_te))

for(j in 1:n_te) {
  if(j %% 100 == 0) cat("   -> Pixel", j, "su", n_te, "\n")
  
  if(is.na(l1_te[j]) || is.na(sigma_te[j])) {
    Z_hat_kriging_global[j] <- mu_te[j]
    SE_kriging_global[j] <- NA
    next
  }
  
  # Eseguiamo il calcolo e salviamo il risultato in una lista temporanea
  risultato_pixel <- tryCatch({
    
    # 1. TROVA I VICINI PIÙ PROSSIMI (INTORNO MOBILE)
    # Calcolo distanza al quadrato dal pixel j a tutte le stazioni di training
    dist_sq <- (train_coords_mat[,1] - grid_coords_mat[j,1])^2 + 
      (train_coords_mat[,2] - grid_coords_mat[j,2])^2
    
    # Seleziono i 100 vicini più prossimi
    n_neighbors <- min(100, n_tr) 
    idx_vicini <- order(dist_sq)[1:n_neighbors]
    
    # 2. CREA IL SOTTOINSIEME LOCALE
    loc_tr_sub <- train_coords_mat[idx_vicini, , drop=FALSE]
    loc_sub    <- rbind(loc_tr_sub, grid_coords_mat[j, , drop=FALSE])
    
    l1_sub     <- c(l1_tr[idx_vicini], l1_te[j])
    l2_sub     <- c(l2_tr[idx_vicini], l2_te[j])
    th_sub     <- c(theta_tr[idx_vicini], theta_te[j])
    sig_sub    <- c(sigma_tr[idx_vicini], sigma_te[j])
    nugget_sub <- c(nugget_tr[idx_vicini])
    res_tr_sub <- res_tr[idx_vicini]
    
    # 3. CALCOLA LA COVARIANZA SOLO PER QUESTO PICCOLO INTORNO
    Aniso_List_Sub <- Compute_Aniso(l1_sub, l2_sub, th_sub)
    Cov_Sub <- matern_ns_corr(t(loc_sub), Aniso_List_Sub, l1_sub * l2_sub, 0.5, sig_sub)
    
    # 4. ESTRAI MATRICE DI TRAINING LOCALE E CROSS-COVARIANZA
    C_nn_locale <- Cov_Sub[1:n_neighbors, 1:n_neighbors, drop = FALSE]
    diag(C_nn_locale) <- diag(C_nn_locale) + nugget_sub + 1e-6 # Stabilità numerica
    
    c_nt <- Cov_Sub[1:n_neighbors, n_neighbors + 1]
    
    # 5. RISOLVI IL SISTEMA (Kriging Weights)
    # Su una matrice 40x40, solve() o chol() sono istantanei
    U_locale <- chol(C_nn_locale)
    weights_loc <- backsolve(U_locale, forwardsolve(t(U_locale), c_nt))
    
    # 6. CALCOLO PREDIZIONE E INCERTEZZA
    pred_val <- mu_te[j] + sum(weights_loc * res_tr_sub)
    
    var_a_priori  <- sigma_te[j]^2 + nugget_te[j]
    riduzione_var <- sum(weights_loc * c_nt)
    se_val        <- sqrt(max(0, var_a_priori - riduzione_var))
    
    # Restituiamo i due valori come lista
    list(pred = pred_val, se = se_val)
    
  }, error = function(e) {
    # Fallback in caso di errore (es. Cholesky fallisce)
    list(pred = mu_te[j], se = sqrt(sigma_te[j]^2 + nugget_te[j]))
  })
  
  # Assegnazione finale ai vettori globali (senza <<- )
  Z_hat_kriging_global[j] <- risultato_pixel$pred
  SE_kriging_global[j]    <- risultato_pixel$se
}

# Salviamo i risultati nel dataframe di plotting
plot_data$Pred_Kriging <- Z_hat_kriging_global
plot_data$SE_Kriging <- SE_kriging_global

print("Kriging completato con successo!")


################################################################################
########################## SALVATAGGIO PER PONTE ###############################
################################################################################

print(">>> Preparazione esportazione per il Ponte...")

# --- AGGIUNTA FONDAMENTALE: CREAZIONE cell_id ---
# Assegniamo un ID riga alla griglia di predizione
plot_data$cell_id <- 1:nrow(plot_data)

# --- AGGIUNTA PER IL PONTE: CALCOLO MATRICE COVARIANZA GLOBALE RDD ---
print("Calcolo Matrice di Covarianza Globale per il plot spaziale (Stile Fig.2)...")

C_globale_rdd <- get_induced_covariance_matrix(
  preds_aggr = plot_data, 
  coords = plot_data[, c("X_km", "Y_km")], 
  nu = 0.5
)

# Aggiunta del nugget (se presente) sulla diagonale principale
nugget_grid <- if("nugget" %in% names(plot_data)) plot_data$nugget else rep(0, nrow(plot_data))
diag(C_globale_rdd) <- diag(C_globale_rdd) + nugget_grid
# --------------------------------------------------------------------

# 1. DEFINIZIONE CARTELLA DI EXPORT COMUNE
#cartella_export <- "C:/Users/miche/OneDrive/Desktop/prova3/INGV - pgv"
cartella_export <- "/work/u10700026/export_pgv"

if(!dir.exists(cartella_export)) {
  dir.create(cartella_export, recursive = TRUE)
}

# 2. IDENTIFICAZIONE DINAMICA DEL NOME FILE BASATA SU K
k_corrente <- if(exists("K_val")) K_val else K

if (k_corrente == 1) {
  nome_file_export <- "Risultati_RDD_K1.rds"
} else {
  nome_file_export <- "Risultati_RDD_Kopt.rds"
}

path_export <- file.path(cartella_export, nome_file_export)

# 3. PREPARAZIONE STAZIONI OBS (Uniformata a quella del LOESS)
stazioni_per_ponte <- spatial_data_df %>% 
  distinct(NET_STA, .keep_all = TRUE) %>% 
  select(NET_STA, X_km, Y_km)

# Assegniamo a ogni stazione il cell_id della griglia più vicina
# Questo è vitale per il "Confronto 9" del Ponte
stazioni_per_ponte$cell_id <- sapply(1:nrow(stazioni_per_ponte), function(i) {
  dist_sq <- (plot_data$X_km - stazioni_per_ponte$X_km[i])^2 + 
    (plot_data$Y_km - stazioni_per_ponte$Y_km[i])^2
  return(plot_data$cell_id[which.min(dist_sq)])
})

# Aggiungiamo latitudine e longitudine dal dataframe originale itacentrale
info_geo <- itacentrale %>% distinct(NET_STA, .keep_all = TRUE) %>% select(NET_STA, st_lon, st_lat)
stazioni_per_ponte <- merge(stazioni_per_ponte, info_geo, by = "NET_STA", all.x = TRUE)

# 4. IMPACCHETTAMENTO RISULTATI
risultati_rdd <- list(
  mappa = plot_data,
  cov = Cov_Mat_Stations, 
  stazioni = rownames(Cov_Mat_Stations), 
  rmse = NA,   
  
  # Aggiunte per abilitare il Confronto 9 nel Ponte:
  cov_globale = NULL,
  stazioni_obs = stazioni_per_ponte
)

# 5. SALVATAGGIO
saveRDS(risultati_rdd, file = path_export)
print(paste("Risultati RDD esportati con successo in:", path_export))

################################################################################

# 1. Carichi il file salvato
setwd("C:/Users/miche/OneDrive/Desktop/prova3/INGV - pgv")
res_rdd <- readRDS("Risultati_RDD_Kopt.rds")

# 2. Carichi le tue funzioni
setwd("C:/Users/miche/OneDrive/Desktop/prova3/Scripts")
source("Utilities.R")

# 3. Calcoli la covarianza globale in post-produzione
C_globale_rdd <- get_induced_covariance_matrix(
  preds_aggr = res_rdd$mappa, 
  coords = res_rdd$mappa[, c("X_km", "Y_km")], 
  nu = 0.5
)

# 4. Aggiungi il nugget sulla diagonale (come nel tuo script principale)
nugget_grid <- if("nugget" %in% names(res_rdd$mappa)) res_rdd$mappa$nugget else rep(0, nrow(res_rdd$mappa))
diag(C_globale_rdd) <- diag(C_globale_rdd) + nugget_grid

# 5. INSERISCI la matrice calcolata all'interno della lista
res_rdd$cov_globale <- C_globale_rdd

# 6. SALVA il file aggiornato sul disco (sovrascrivendo il precedente)
setwd("C:/Users/miche/OneDrive/Desktop/prova3/INGV - pgv")
saveRDS(res_rdd, file = "Risultati_RDD_Kopt.rds")

cat("File aggiornato con la cov_globale!")
