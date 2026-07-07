
################################################################################
# SCRIPT 4: PONTE ESTESO - CONFRONTO 4 MODELLI                                 #
# (LOESS, RDD_K_Opt, RDD_K_1, Stazionario)                                     #
################################################################################

# SOMMARIO
#
# 1. VERDETTO NUMERICO (RMSE):
#    Classifica delle performance predittive basata sull'errore stimato in Cross-Validation.
# 2. MAPPE PREDITTIVE (dS2S):
#    Visualizzazione spaziale affiancata delle superfici dell'effetto di sito stimate dai 4 modelli.
# 3. MAPPE INCERTEZZA (SE):
#    Distribuzione spaziale dell'errore standard di stima per valutare dove i modelli sono più incerti.
# 4. STRUTTURA MATRICI COVARIANZA:
#    Plot delle matrici teoriche per le stazioni in comune, ordinate geograficamente (Nord -> Sud).
# 5. DIFFERENZA COVARIANZE (TUTTI vs LOESS):
#    Mappe di calore della differenza numerica punto-punto tra le matrici parametriche e il riferimento LOESS.
# 6. DIFFERENZE PREDITTIVE E SCATTERPLOT:
#    Mappe dei Delta spaziali (Modello - LOESS) e grafici di correlazione per valutare l'accordo puntuale.
# 7. VARIANZA DI SITO (DIAGONALE):
#    Confronto della varianza pura catturata ai siti (diagonale delle matrici) rispetto al non-parametrico.
# 8. DISTRIBUZIONE DEGLI SCARTI:
#    Istogrammi delle differenze spaziali per individuare eventuali bias sistematici rispetto al LOESS.


library(ggplot2)
library(viridis)
library(dplyr)
library(tidyr)
library(sf)
library(rnaturalearth)

print(">>> Avvio Ponte Esteso: Confronto 4 Modelli...")

# ==============================================================================
# 1. CARICAMENTO DATI
# ==============================================================================

setwd("C:/Users/miche/OneDrive/Desktop/prova3/INGV - pgv") 
files_necessari <- c("Risultati_NonParam.rds", "Risultati_RDD_Kopt.rds",
                     "Risultati_RDD_K1.rds", "Risultati_Stazionario.rds")

for(f in files_necessari) {
  if(!file.exists(f)) stop(paste("Errore: Manca il file", f))
}

res_loess <- readRDS("Risultati_NonParam.rds")
res_rdd   <- readRDS("Risultati_RDD_Kopt.rds")
res_rdd_k1<- readRDS("Risultati_RDD_K1.rds")
res_staz  <- readRDS("Risultati_Stazionario.rds")

# ------------------------------------------
# CONFRONTO 1: VERDETTO NUMERICO (CV RMSE)
# ------------------------------------------
df_rmse <- data.frame(
  Modello = c("1. Non-Parametrico (LOESS)", "2. RDD Multi-Cluster", "3. RDD K=1",
              paste("4.", res_staz$modello_scelto)),
              RMSE = c(res_loess$rmse, res_rdd$rmse, res_rdd_k1$rmse, res_staz$rmse)
)
df_rmse <- df_rmse %>% arrange(RMSE)

cat("\n--------------------------------------\n")
cat("      CLASSIFICA CROSS-VALIDATION       \n")
cat("\n--------------------------------------\n")
print(df_rmse, row.names = FALSE)
cat("----------------------------------------\n\n")

# ==============================================================================
# 2. ALLINEAMENTO GRIGLIE SPAZIALI
# ==============================================================================
# Dato che tutte e 4 le griglie sono ora state generate con la stessa logica,
# i valori di X_km e Y_km combaceranno esattamente. Li arrotondiamo al 3° decimale 
# (metro) solo per evitare noie con la precisione in virgola mobile di R.

res_loess$mappa$X_km <- round(res_loess$mappa$X_km, 3)
res_loess$mappa$Y_km <- round(res_loess$mappa$Y_km, 3)

res_rdd$mappa$X_km <- round(res_rdd$mappa$X_km, 3)
res_rdd$mappa$Y_km <- round(res_rdd$mappa$Y_km, 3)

res_rdd_k1$mappa$X_km <- round(res_rdd_k1$mappa$X_km, 3)
res_rdd_k1$mappa$Y_km <- round(res_rdd_k1$mappa$Y_km, 3)

res_staz$mappa$X_km <- round(res_staz$mappa$X_km, 3)
res_staz$mappa$Y_km <- round(res_staz$mappa$Y_km, 3)

# Fusione a catena
mappa_base <- merge(res_loess$mappa, res_rdd$mappa, by = c("X_km", "Y_km"))
mappa_base <- merge(mappa_base, res_rdd_k1$mappa, by = c("X_km", "Y_km"), suffixes = c("_RDD", "_RDD_K1"))
mappa_base <- merge(mappa_base, res_staz$mappa, by = c("X_km", "Y_km"))

# Rinominiamo le colonne per chiarezza e compatibilità con il resto dello script
colnames(mappa_base)[colnames(mappa_base) == "dS2S_kriging"] <- "Pred_LOESS"
colnames(mappa_base)[colnames(mappa_base) == "se_kriging"]   <- "SE_LOESS"
colnames(mappa_base)[colnames(mappa_base) == "Pred_Kriging_RDD"] <- "Pred_RDD"
colnames(mappa_base)[colnames(mappa_base) == "SE_Kriging_RDD"]   <- "SE_RDD"
colnames(mappa_base)[colnames(mappa_base) == "Pred_Kriging_RDD_K1"] <- "Pred_RDD_K1"
colnames(mappa_base)[colnames(mappa_base) == "SE_Kriging_RDD_K1"]   <- "SE_RDD_K1"

# Filtro Terraferma per pulire le mappe (lo facciamo alla fine così i dati sono completi)
sf_base <- st_as_sf(mappa_base, coords=c("X_km", "Y_km"), crs=32633, remove=FALSE)
ita_border <- ne_countries(scale = "medium", country = "Italy", returnclass = "sf")
ita_border_utm <- st_transform(ita_border, 32633)

st_geometry(sf_base) <- st_geometry(sf_base) * 1000 # Sf in metri per la mascheratura
st_crs(sf_base) <- 32633

inside <- st_intersects(sf_base, ita_border_utm, sparse = FALSE)
mappa_plot <- mappa_base[as.vector(inside), ]

# ------------------------------------------
# CONFRONTO 2: MAPPE PREDITTIVE
# ------------------------------------------
mappa_long_pred <- mappa_plot %>%
  select(X_km, Y_km, Pred_LOESS, Pred_RDD, Pred_RDD_K1, Pred_Staz) %>%
  pivot_longer(cols = starts_with("Pred_"), names_to = "Modello", values_to = "dS2S") %>%
  mutate(Modello = factor(Modello, levels = c("Pred_LOESS", "Pred_RDD", "Pred_RDD_K1", "Pred_Staz"),
                          labels = c("1. LOESS", "2. RDD Multi", "3. RDD K=1", "4. Stazionario")))

lim_val <- max(abs(mappa_long_pred$dS2S), na.rm=TRUE)

p_pred <- ggplot(mappa_long_pred, aes(x = X_km, y = Y_km, fill = dS2S)) +
  geom_tile() +
  facet_wrap(~ Modello, nrow = 1) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       limits = c(-lim_val, lim_val)) +
  coord_fixed() + theme_minimal() + 
  labs(title = "Confronto Mappe Predittive dS2S") +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

print(p_pred)

# ------------------------------------------
# CONFRONTO 3: MAPPE INCERTEZZA
# ------------------------------------------
mappa_long_se <- mappa_plot %>%
  select(X_km, Y_km, SE_LOESS, SE_RDD, SE_RDD_K1, SE_Staz) %>%
  pivot_longer(cols = starts_with("SE_"), names_to = "Modello", values_to = "SE") %>%
  mutate(Modello = factor(Modello, levels = c("SE_LOESS", "SE_RDD", "SE_RDD_K1", "SE_Staz"),
                          labels = c("1. LOESS", "2. RDD Multi", "3. RDD K=1", "4. Stazionario")))

lim_se <- max(mappa_long_se$SE, na.rm=TRUE)

p_se <- ggplot(mappa_long_se, aes(x = X_km, y = Y_km, fill = SE)) +
  geom_tile() +
  facet_wrap(~ Modello, nrow = 1) +
  scale_fill_viridis_c(option = "mako", direction = -1, limits = c(0, lim_se)) +
  coord_fixed() + theme_minimal() + 
  labs(title = "Confronto Mappe Incertezza (Standard Error)") +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

print(p_se)

# ------------------------------------------
# CONFRONTO 4: STRUTTURA MATRICI COVARIANZA
# ------------------------------------------
# Identifichiamo le stazioni in comune a tutti e 4 i modelli
nomi_comuni <- Reduce(intersect, list(res_loess$stazioni, res_rdd$stazioni,
                                      res_rdd_k1$stazioni, res_staz$stazioni))

# --- FIX: ORDINAMENTO NORD-SUD REALE ---
# Recuperiamo le coordinate dal dataset originale per essere sicuri
dataset_raw <- read.csv("dataset_martina.csv")
stazioni_coords <- dataset_raw %>%
  filter(NET_STA %in% nomi_comuni) %>%
  group_by(NET_STA) %>%
  summarise(st_lat = first(st_lat)) %>%
  arrange(desc(st_lat)) # Ordine decrescente (Nord -> Sud)

nomi_ordinati_NS <- stazioni_coords$NET_STA

# Estrazione e Melting
extract_cov <- function(res_obj, nome_modello) {
  mat <- res_obj$cov
  rownames(mat) <- res_obj$stazioni
  colnames(mat) <- res_obj$stazioni
  mat_sub <- mat[nomi_ordinati_NS, nomi_ordinati_NS]
  df <- reshape2::melt(mat_sub)
  df$Modello <- nome_modello
  return(df)
}

df_cov_all <- rbind(
  extract_cov(res_loess, "1. LOESS"),
  extract_cov(res_rdd, "2. RDD Multi"),
  extract_cov(res_rdd_k1, "3. RDD K=1"),
  extract_cov(res_staz, "4. Stazionario")
)

df_cov_all$Var1 <- factor(df_cov_all$Var1, levels = nomi_ordinati_NS)
df_cov_all$Var2 <- factor(df_cov_all$Var2, levels = nomi_ordinati_NS)

p_cov_all <- ggplot(df_cov_all, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  facet_wrap(~ Modello, nrow = 2) +
  scale_fill_viridis_c(option = "magma", name = "Covarianza") +
  theme_minimal() + coord_fixed() +
  labs(title = "Struttura delle Matrici di Covarianza Teoriche (Ordinate N -> S)") +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank())

print(p_cov_all)

library(gridExtra)

# ------------------------------------------
# CONFRONTO 5: DIFFERENZA TRA MATRICI DI COVARIANZA (TUTTI vs LOESS)
# ------------------------------------------
cat("\nCalcolo delle differenze numeriche tra le Matrici di Covarianza...\n")

# Funzione per calcolare la differenza rispetto al LOESS (Arricchita con distanze numeriche)
calc_diff_cov <- function(res_obj, ref_obj, nome_modello) {
  
  # Sottrazione delle matrici (Ordinate N-S)
  mat_diff <- res_obj$cov[nomi_ordinati_NS, nomi_ordinati_NS] - 
    ref_obj$cov[nomi_ordinati_NS, nomi_ordinati_NS]
  
  # CALCOLO AGGREGATO DELLE DIFFERENZE
  radice_somma_quadratica <- sqrt(sum(mat_diff^2)) # Norma di Frobenius
  rmse_matrice <- sqrt(mean(mat_diff^2))           # RMSE della matrice
  
  cat(sprintf("\n=> Distanza Matrici (%s vs LOESS):\n", nome_modello))
  cat(sprintf("   - Radice Somma Quadrati : %.4f\n", radice_somma_quadratica))
  cat(sprintf("   - RMSE Matrice          : %.6f\n", rmse_matrice))
  # =======================================================
  
  df <- reshape2::melt(as.matrix(mat_diff))
  df$Confronto <- paste(nome_modello, "vs LOESS")
  return(df)
}

# Uniamo le 3 differenze in un unico dataframe
df_diff_cov_all <- rbind(
  calc_diff_cov(res_rdd,    res_loess, "2. RDD Multi"),
  calc_diff_cov(res_rdd_k1, res_loess, "3. RDD K=1"),
  calc_diff_cov(res_staz,   res_loess, "4. Stazionario")
)

# Plot delle matrici di differenza
lim_cov_diff <- max(abs(df_diff_cov_all$value), na.rm=TRUE)

p_diff_mat_all <- ggplot(df_diff_cov_all, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  facet_wrap(~ Confronto, nrow = 1) +
  scale_fill_gradient2(low = "purple", mid = "white", high = "orange", 
                       midpoint = 0, limits = c(-lim_cov_diff, lim_cov_diff),
                       name = "Diff Cov") +
  labs(title = "Differenza Punto-Punto tra le Matrici di Covarianza (vs LOESS)",
       subtitle = "Arancione: Covarianza maggiore del LOESS | Viola: Covarianza minore del LOESS") +
  theme_minimal() + coord_fixed() +
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())

print(p_diff_mat_all)

# ------------------------------------------
# CONFRONTO 6: SCATTERPLOT, DISTRIBUZIONE E MAPPA DELTA
# ------------------------------------------
# 1. Calcolo dei Delta rispetto al LOESS
mappa_delta <- mappa_plot %>%
  mutate(
    Delta_RDD = Pred_RDD - Pred_LOESS,
    Delta_RDD_K1 = Pred_RDD_K1 - Pred_LOESS,
    Delta_Staz = Pred_Staz - Pred_LOESS
  ) %>%
  select(X_km, Y_km, Delta_RDD, Delta_RDD_K1, Delta_Staz) %>%
  pivot_longer(cols = starts_with("Delta_"), names_to = "Modello", values_to = "Differenza")

# 2. Mappe dei Delta (Dove i modelli dissentono dal LOESS?)
lim_diff <- max(abs(mappa_delta$Differenza), na.rm=TRUE)

p_delta_all <- ggplot(mappa_delta, aes(x = X_km, y = Y_km, fill = Differenza)) +
  geom_tile() +
  facet_wrap(~ Modello, nrow = 1) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       limits = c(-lim_diff, lim_diff), name = "Delta") +
  coord_fixed() + theme_minimal() +
  labs(title = "Differenza Spaziale dei Modelli rispetto al LOESS",
       subtitle = "Valori Rossi/Blu indicano dove il modello si allontana dal Non-Parametrico") +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

print(p_delta_all)

# 3. Scatterplot Combinato (LOESS vs ALTRI)
mappa_scatter <- mappa_plot %>%
  select(Pred_LOESS, Pred_RDD, Pred_RDD_K1, Pred_Staz) %>%
  pivot_longer(cols = c(Pred_RDD, Pred_RDD_K1, Pred_Staz), 
               names_to = "Modello", values_to = "Valore_Pred")

p_scatter_all <- ggplot(mappa_scatter, aes(x = Pred_LOESS, y = Valore_Pred)) +
  geom_point(alpha = 0.1, color = "darkblue", size = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  facet_wrap(~ Modello) +
  labs(title = "Correlazione delle Predizioni", 
       x = "dS2S (LOESS)", y = "dS2S (Altri Modelli)") +
  theme_minimal() + coord_fixed()

print(p_scatter_all)

# ------------------------------------------
# CONFRONTO 7: VARIANZA DI SITO (DIAGONALE MATRICI)
# ------------------------------------------
# Estraiamo le diagonali da tutti e 4 i modelli
cov_loess_ord  <- res_loess$cov[nomi_ordinati_NS, nomi_ordinati_NS]
cov_rdd_ord    <- res_rdd$cov[nomi_ordinati_NS, nomi_ordinati_NS]
cov_rdd_k1_ord <- res_rdd_k1$cov[nomi_ordinati_NS, nomi_ordinati_NS] # <--- AGGIUNTO K=1
cov_staz_ord   <- res_staz$cov[nomi_ordinati_NS, nomi_ordinati_NS]

var_df <- data.frame(
  Stazione = nomi_ordinati_NS,
  Var_LOESS = diag(cov_loess_ord),
  Var_RDD = diag(cov_rdd_ord),
  Var_RDD_K1 = diag(cov_rdd_k1_ord), # <--- AGGIUNTO K=1
  Var_Staz = diag(cov_staz_ord)
)

# Plot comparativo globale
p_var <- ggplot(var_df) +
  geom_point(aes(x = Var_LOESS, y = Var_RDD, color = "2. RDD Multi"), size = 2, alpha = 0.6) +
  geom_point(aes(x = Var_LOESS, y = Var_RDD_K1, color = "3. RDD K=1"), size = 2, alpha = 0.6) +
  geom_point(aes(x = Var_LOESS, y = Var_Staz, color = "4. Stazionario"), size = 2, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  scale_color_manual(values = c("2. RDD Multi" = "forestgreen", 
                                "3. RDD K=1" = "orange", 
                                "4. Stazionario" = "red")) +
  labs(title = "Confronto Varianze di Sito (Diagonale)",
       subtitle = "Ogni punto è una stazione. Il riferimento (X) è il modello LOESS.",
       x = "Varianza (LOESS)", y = "Varianza Modelli", color = "Modello") +
  theme_minimal() + coord_fixed(ratio = 1)

print(p_var)

# ------------------------------------------
# CONFRONTO 8: DISTRIBUZIONE DEGLI SCARTI PREDITTIVI (ISTOGRAMMI)
# ------------------------------------------
# Usiamo il dataframe mappa_delta creato nella Sezione 5
p_hist_all <- ggplot(mappa_delta, aes(x = Differenza, fill = Modello)) +
  geom_histogram(bins = 50, color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  facet_wrap(~ Modello, nrow = 1) +
  scale_fill_viridis_d(option = "turbo") +
  labs(title = "Distribuzione delle Differenze Spaziali (Modello - LOESS)", 
       x = "Delta dS2S", y = "Conteggio Celle") +
  theme_minimal() +
  theme(legend.position = "none")

print(p_hist_all)

cat("\n>>> Analisi completata. Tutti i confronti sono stati generati rispetto al LOESS.\n")

# ==============================================================================
# CONFRONTO 9: MAPPE DI COVARIANZA SPAZIALE (STILE FIGURA 2 PAPER)
# ==============================================================================
cat("\nGenerazione del confronto spaziale della covarianza (Stile Figura 2)...\n")

# Assicuriamoci che la mappa geometrica di base abbia la colonna cell_id
mappa_fig2_base <- res_loess$mappa %>% select(X_km, Y_km, cell_id)

# ------------------------------------------------------------------------------
# CASO A: STAZIONE DI ALTO INTERESSE SISMICO (es. Zona Appenninica)
# ------------------------------------------------------------------------------
# Inserisci una stazione nota
stazione_target <- "IT.AQK" 

if(stazione_target %in% res_loess$stazioni_obs$NET_STA) {
  cella_rif_A <- res_loess$stazioni_obs$cell_id[res_loess$stazioni_obs$NET_STA == stazione_target]
  pt_x_A      <- res_loess$stazioni_obs$X_km[res_loess$stazioni_obs$NET_STA == stazione_target]
  pt_y_A      <- res_loess$stazioni_obs$Y_km[res_loess$stazioni_obs$NET_STA == stazione_target]
  
  mappa_A <- mappa_fig2_base
  mappa_A$Cov_LOESS <- as.numeric(res_loess$cov_globale[cella_rif_A, mappa_A$cell_id])
  
  # Scommenta le prossime righe se hai esportato la cov_globale anche per gli altri modelli
  # mappa_A$Cov_Staz <- as.numeric(res_staz$cov_globale[cella_rif_A, mappa_A$cell_id])
  # mappa_A$Cov_RDD  <- as.numeric(res_rdd$cov_globale[cella_rif_A, mappa_A$cell_id])
  
  df_A_long <- mappa_A %>%
    pivot_longer(cols = starts_with("Cov_"), names_to = "Modello", values_to = "Covarianza")
  
  lim_A <- max(abs(df_A_long$Covarianza), na.rm = TRUE)
  
  plot_A <- ggplot(df_A_long, aes(x = X_km, y = Y_km, fill = Covarianza)) +
    geom_tile() +
    facet_wrap(~ Modello, nrow = 1) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-lim_A, lim_A)) +
    geom_point(aes(x = pt_x_A, y = pt_y_A), color = "black", size = 3, shape = 16) +
    coord_fixed() + theme_minimal() +
    labs(title = paste("CASO A: Covarianza spaziale - Stazione", stazione_target),
         subtitle = "Propagazione della covarianza da un sito osservato ad alta densità", fill = "Cov") +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
  
  print(plot_A)
} else {
  cat("Attenzione: La stazione target non è nel dataset.\n")
}

# ------------------------------------------------------------------------------
# CASO B: LUOGO DI INTERESSE NAZIONALE SENZA STAZIONE
# ------------------------------------------------------------------------------
# Inserisci le coordinate WGS84 (Gradi Decimali) della località che ti interessa
lat_target_B <- 42  
lon_target_B <- 13.75

# Trasformazione automatica in UTM 33N (km)
pt_sf <- st_sfc(st_point(c(lon_target_B, lat_target_B)), crs = 4326)
pt_utm <- st_transform(pt_sf, crs = 32633)
coords_km <- st_coordinates(pt_utm) / 1000

x_no_staz <- coords_km[1]
y_no_staz <- coords_km[2]

# Troviamo la cella geometrica più vicina a queste coordinate
distanze_sq <- (mappa_fig2_base$X_km - x_no_staz)^2 + (mappa_fig2_base$Y_km - y_no_staz)^2
cella_rif_B <- mappa_fig2_base$cell_id[which.min(distanze_sq)]

pt_x_B <- mappa_fig2_base$X_km[mappa_fig2_base$cell_id == cella_rif_B]
pt_y_B <- mappa_fig2_base$Y_km[mappa_fig2_base$cell_id == cella_rif_B]

mappa_B <- mappa_fig2_base
mappa_B$Cov_LOESS <- as.numeric(res_loess$cov_globale[cella_rif_B, mappa_B$cell_id])

# Scommenta le prossime righe se hai esportato la cov_globale anche per gli altri modelli
# mappa_B$Cov_Staz <- as.numeric(res_staz$cov_globale[cella_rif_B, mappa_B$cell_id])
# mappa_B$Cov_RDD  <- as.numeric(res_rdd$cov_globale[cella_rif_B, mappa_B$cell_id])

df_B_long <- mappa_B %>%
  pivot_longer(cols = starts_with("Cov_"), names_to = "Modello", values_to = "Covarianza")

lim_B <- max(abs(df_B_long$Covarianza), na.rm = TRUE)

plot_B <- ggplot(df_B_long, aes(x = X_km, y = Y_km, fill = Covarianza)) +
  geom_tile() +
  facet_wrap(~ Modello, nrow = 1) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-lim_B, lim_B)) +
  geom_point(aes(x = pt_x_B, y = pt_y_B), color = "black", size = 3, shape = 16) +
  coord_fixed() + theme_minimal() +
  labs(title = sprintf("CASO B: Covarianza in un gap strumentale (es. Roma) [Lat: %.2f, Lon: %.2f]", lat_target_B, lon_target_B),
       subtitle = "Capacità del modello di estrapolare la struttura spaziale bidimensionale", fill = "Cov") +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

print(plot_B)
