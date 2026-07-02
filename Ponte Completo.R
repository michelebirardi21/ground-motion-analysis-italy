
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
# 4. STRUTTURA MATRICI COVARIANZA: (ni, nell'appendice)
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

setwd("C:/Users/miche/OneDrive/Desktop/prova3/INGV") 
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
# NUOVO CONFRONTO 9: MAPPE DI COVARIANZA SPAZIALE RICOSTRUITA 2x3 (PGA e PGV)
# ==============================================================================
cat("\nGenerazione del confronto spaziale della covarianza ricostruita (2x3)...\n")

library(rnaturalearth)
library(sf)
library(ggplot2)
library(dplyr)

# 1. Carichiamo i confini dell'Italia e proiettiamoli in UTM 33N (metri)
ita_border <- ne_countries(scale = "medium", country = "Italy", returnclass = "sf")
ita_border_utm <- st_transform(ita_border, 32633)

# 2. Caricamento dei 6 file RDS (Modifica i percorsi se necessario)
path_pga_staz <- "C:/Users/miche/OneDrive/Desktop/prova3/INGV/Risultati_Stazionario.rds"
path_pga_rdd  <- "C:/Users/miche/OneDrive/Desktop/prova3/INGV/Risultati_RDD_Kopt.rds"
path_pga_np   <- "C:/Users/miche/OneDrive/Desktop/prova3/INGV/Risultati_NonParam.rds"

path_pgv_staz <- "C:/Users/miche/OneDrive/Desktop/prova3/INGV - pgv/Risultati_Stazionario.rds"
path_pgv_rdd  <- "C:/Users/miche/OneDrive/Desktop/prova3/INGV - pgv/Risultati_RDD_Kopt.rds"
path_pgv_np   <- "C:/Users/miche/OneDrive/Desktop/prova3/INGV - pgv/Risultati_NonParam.rds"

pga_staz <- readRDS(path_pga_staz)
pga_rdd  <- readRDS(path_pga_rdd)
pga_np   <- readRDS(path_pga_np)

pgv_staz <- readRDS(path_pgv_staz)
pgv_rdd  <- readRDS(path_pgv_rdd)
pgv_np   <- readRDS(path_pgv_np)

# 3. Prepariamo la griglia globale di base della terraferma (arrotondamento salvavita)
mappa_fig2_base <- pga_np$mappa %>% 
  mutate(X_km = round(X_km, 3), Y_km = round(Y_km, 3)) %>% 
  select(X_km, Y_km, cell_id)

griglia_sf <- st_as_sf(mappa_fig2_base, coords = c("X_km", "Y_km"), crs = 32633, remove = FALSE)
st_geometry(griglia_sf) <- st_geometry(griglia_sf) * 1000
st_crs(griglia_sf) <- 32633

inside_italy <- st_intersects(griglia_sf, ita_border_utm, sparse = FALSE)
griglia_terraferma <- mappa_fig2_base[as.vector(inside_italy), ]

# 4. Configurazione della stazione target
stazione_target <- "IT.NRC"

# Identifichiamo dinamicamente il cell_id di riferimento nei rispettivi cataloghi
cella_pga <- pga_np$stazioni_obs$cell_id[grep(stazione_target, pga_np$stazioni_obs$NET_STA)[1]]
cella_pgv <- pgv_np$stazioni_obs$cell_id[grep(stazione_target, pgv_np$stazioni_obs$NET_STA)[1]]

# Coordinate reali della stazione per posizionare il punto nero al centro dello zoom
# Aggiunto [1] per garantire che estragga un singolo numero
pt_x <- pga_np$stazioni_obs$X_km[pga_np$stazioni_obs$cell_id == cella_pga][1]
pt_y <- pga_np$stazioni_obs$Y_km[pga_np$stazioni_obs$cell_id == cella_pga][1]
# 5. Funzione interna per estrarre la riga di covarianza globale ricostruita (Modello Liscio)
estrai_cov_spaziale <- function(res_obj, cella_rif, nome_modello, nome_im) {
  df <- griglia_terraferma
  # Estraiamo la riga corrispondente alla cella della stazione
  df$Covarianza <- as.numeric(res_obj$cov_globale[cella_rif, df$cell_id])
  df$Modello <- nome_modello
  df$IM <- nome_im
  return(df)
}

# Estrazione coordinata dai 6 modelli
df_pga_staz_cov <- estrai_cov_spaziale(pga_staz, cella_pga, "Stationary", "PGA")
df_pga_rdd_cov  <- estrai_cov_spaziale(pga_rdd,  cella_pga, "RDD", "PGA")
df_pga_np_cov   <- estrai_cov_spaziale(pga_np,   cella_pga, "Non-Parametric", "PGA")

df_pgv_staz_cov <- estrai_cov_spaziale(pgv_staz, cella_pgv, "Stationary", "PGV")
df_pgv_rdd_cov  <- estrai_cov_spaziale(pgv_rdd,  cella_pgv, "RDD", "PGV")
df_pgv_np_cov   <- estrai_cov_spaziale(pgv_np,   cella_pgv, "Non-Parametric", "PGV")

# Unione in un solo grande dataset
mappa_all_cov <- bind_rows(df_pga_staz_cov, df_pga_rdd_cov, df_pga_np_cov,
                           df_pgv_staz_cov, df_pgv_rdd_cov, df_pgv_np_cov)

# Forziamo i livelli dei fattori per impaginare correttamente la griglia
mappa_all_cov$Modello <- factor(mappa_all_cov$Modello, levels = c("Stationary", "RDD", "Non-Parametric"))
mappa_all_cov$IM <- factor(mappa_all_cov$IM, levels = c("PGA", "PGV"))

library(patchwork)

# 6. Calcolo dei limiti geometrici dello zoom (70 km attorno alla stazione)
raggio_zoom_km <- 200
buffer_zoom_m <- raggio_zoom_km * 1000
xlim_m <- c((pt_x * 1000) - buffer_zoom_m, (pt_x * 1000) + buffer_zoom_m)
ylim_m <- c((pt_y * 1000) - buffer_zoom_m, (pt_y * 1000) + buffer_zoom_m)

# Calcolo dei limiti geometrici SEPARATI per le due scale colore
lim_pga <- max(abs(mappa_all_cov$Covarianza[mappa_all_cov$IM == "PGA"]), na.rm = TRUE)
lim_pgv <- max(abs(mappa_all_cov$Covarianza[mappa_all_cov$IM == "PGV"]), na.rm = TRUE)


# 7. Generazione dei due plot separati e unione con Patchwork

library(ggpubr)

# --- RIGA SUPERIORE: PGA ---
p_pga <- ggplot(subset(mappa_all_cov, IM == "PGA")) +
  geom_tile(aes(x = X_km * 1000, y = Y_km * 1000, fill = Covarianza)) +
  geom_sf(data = ita_border_utm, fill = NA, color = "black", linewidth = 0.4, inherit.aes = FALSE) +
  annotate("point", x = pt_x * 1000, y = pt_y * 1000, fill = "white", color = "black", size = 2, shape = 21, stroke = 1) +
  
  facet_wrap(~ Modello, nrow = 1) + 
  
  scale_fill_viridis_c(option = "C", 
                       limits = c(-lim_pga, lim_pga), direction = -1,
                       oob = scales::squish, na.value = "grey85", 
                       name = "Cov") +
  coord_sf(xlim = xlim_m, ylim = ylim_m, expand = FALSE, datum = st_crs(4326)) + 
  theme_minimal() +
  labs(x = NULL, y = NULL) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.title.x = element_blank(),  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3, linetype = "dashed"),
    panel.grid.minor = element_blank()
  )

# --- RIGA INFERIORE: PGV ---
p_pgv <- ggplot(subset(mappa_all_cov, IM == "PGV")) +
  geom_tile(aes(x = X_km * 1000, y = Y_km * 1000, fill = Covarianza)) +
  geom_sf(data = ita_border_utm, fill = NA, color = "black", linewidth = 0.4, inherit.aes = FALSE) +
  annotate("point", x = pt_x * 1000, y = pt_y * 1000, fill = "white", color = "black", size = 2, shape = 21, stroke = 1) +
  
  facet_wrap(~ Modello, nrow = 1) +
  
  scale_fill_viridis_c(option = "C", 
                       limits = c(-lim_pgv, lim_pgv), direction = -1,
                       oob = scales::squish, na.value = "grey85", 
                       name = "Cov") +
  coord_sf(xlim = xlim_m, ylim = ylim_m, expand = FALSE, datum = st_crs(4326)) + 
  theme_minimal() +
  labs(x = "Longitude (°)", y = NULL) +
  theme(
    strip.text = element_blank(), 
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3, linetype = "dashed"),
    panel.grid.minor = element_blank()
  )

# --- UNIONE VERTICALE CON GGPUBR (CORRETTA) ---
# ggarrange allinea perfettamente i grafici, rispettando gli assi
p_cov_2x3_base <- ggarrange(p_pga, p_pgv, ncol = 1, nrow = 2, align = "v")

# --- AGGIUNTA ASSE Y GLOBALE ---
p_cov_2x3 <- annotate_figure(p_cov_2x3_base, 
                             left = text_grob("Latitude (°)", rot = 90, size = 12))

# Mostra il grafico a schermo
print(p_cov_2x3)

################################################################################
################### PLOT PREV E INCERT KRIGING PGA E PGV #######################
################################################################################

# Caricamento librerie necessarie
library(ggplot2)
library(viridis)
library(dplyr)
library(sf)
library(rnaturalearth)
library(scales)

print(">>> Avvio generazione Figure: PGA vs PGV (Senza Sardegna, Assi Lat/Lon)...")

# ==============================================================================
# 1. CARICAMENTO E UNIONE DEI DATI
# ==============================================================================
path_pga <- "C:/Users/miche/OneDrive/Desktop/prova3/INGV/Risultati_RDD_Kopt.rds"
path_pgv <- "C:/Users/miche/OneDrive/Desktop/prova3/INGV - pgv/Risultati_RDD_Kopt.rds"

mappa_pga <- readRDS(path_pga)$mappa %>% mutate(IM = "PGA")
mappa_pgv <- readRDS(path_pgv)$mappa %>% mutate(IM = "PGV")

mappa_all <- bind_rows(mappa_pga, mappa_pgv)

if("pred" %in% names(mappa_all)) mappa_all <- rename(mappa_all, Pred_Kriging = pred)
if("std_dev" %in% names(mappa_all)) mappa_all <- rename(mappa_all, SE_Kriging = std_dev)

# ==============================================================================
# 2. CONFINI ITALIA SENZA SARDEGNA + MASCHERATURA
# ==============================================================================
# Carichiamo l'Italia a livello di regioni/province (qui c'è la colonna 'region')
ita_regioni <- ne_states(country = "Italy", returnclass = "sf")

# Teniamo tutto TRANNE la Sardegna. grepl copre eventuali varianti di nome.
ita_no_sardegna <- ita_regioni %>%
  filter(!grepl("Sardegn", region, ignore.case = TRUE))

# Unifichiamo le province in un unico poligono nazionale e proiettiamo in UTM 33N
ita_border  <- st_union(ita_no_sardegna)
italia_utm  <- st_transform(ita_border, 32633)

# Griglia dati -> oggetto spaziale in metri
sf_all <- st_as_sf(mappa_all, coords = c("X_km", "Y_km"), crs = 32633, remove = FALSE)
st_geometry(sf_all) <- st_geometry(sf_all) * 1000
st_crs(sf_all) <- 32633

# Teniamo solo i punti che cadono sulla terraferma (Sardegna ora esclusa)
inside <- st_intersects(sf_all, italia_utm, sparse = FALSE)
mappa_plot <- mappa_all[as.vector(inside), ]

# Limiti assi (in metri UTM)
limiti_x <- range(mappa_plot$X_km * 1000, na.rm = TRUE)
limiti_y <- range(mappa_plot$Y_km * 1000, na.rm = TRUE)

# Risoluzione del pixel (in metri)
RISOLUZIONE <- (sort(unique(mappa_plot$X_km))[2] - sort(unique(mappa_plot$X_km))[1]) * 1000

# ==============================================================================
# FIGURE 1: PREDICTIONS (PGA & PGV)
# ==============================================================================
lim_pred <- range(mappa_plot$Pred, na.rm = TRUE)
#lim_pred <- quantile(mappa_plot$Pred_Kriging, probs = c(0.01, 0.99), na.rm = TRUE)

p_pred <- ggplot(mappa_plot) +
  # Sfondo Italia senza Sardegna
  geom_sf(data = italia_utm, fill = "white", color = "grey80") +
  
  # Rasterizzazione della griglia
  geom_tile(aes(x = X_km * 1000, y = Y_km * 1000, fill = Pred_Kriging), 
            width = RISOLUZIONE, height = RISOLUZIONE) +
  
  # Disegno dei bordi
  geom_sf(data = italia_utm, fill = NA, color = "black", linewidth = 0.3) +
  facet_wrap(~ IM, ncol = 2) +
  
  # Scala viridis con espressione matematica per stampare delta*S2S
  scale_fill_viridis_c(option = "viridis", 
                       name = expression(delta*"S2S"~Prediction), 
                       limits = lim_pred, 
                       oob = scales::squish) +
  
  # coord_sf proietta i dati UTM, ma imposta assi e reticolo in Lat/Lon (4326)
  coord_sf(xlim = limiti_x, ylim = limiti_y, expand = FALSE, datum = st_crs(4326)) +
  scale_x_continuous(breaks = seq(10, 17, by = 2)) +
  scale_y_continuous(breaks = seq(40, 46, by = 2)) +
  
  labs(title = NULL, x = "Longitude (°)", y = "Latitude (°)") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 14),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3, linetype = "dashed"), 
    panel.grid.minor = element_blank()
  )

print(p_pred)
ggsave("Figure1_Predictions.png", plot = p_pred, width = 12, height = 6, dpi = 300, bg = "white")

# ==============================================================================
# FIGURE 2: UNCERTAINTIES (PGA & PGV)
# ==============================================================================
lim_se <- range(mappa_plot$SE_Kriging, na.rm = TRUE)

p_se <- ggplot(mappa_plot) +
  # Sfondo Italia senza Sardegna
  geom_sf(data = italia_utm, fill = "white", color = "grey80") +
  
  # Rasterizzazione della griglia
  geom_tile(aes(x = X_km * 1000, y = Y_km * 1000, fill = SE_Kriging), 
            width = RISOLUZIONE, height = RISOLUZIONE) +
  
  # Disegno dei bordi
  geom_sf(data = italia_utm, fill = NA, color = "black", linewidth = 0.3) +
  facet_wrap(~ IM, ncol = 2) +
  
  # Palette magma invertita e nome corretto "Standard Error"
  scale_fill_viridis_c(option = "magma", 
                       direction = -1, 
                       name = "Standard Error", 
                       limits = lim_se, 
                       oob = scales::squish) +
  
  # coord_sf proietta i dati UTM, ma imposta assi e reticolo in Lat/Lon (4326)
  coord_sf(xlim = limiti_x, ylim = limiti_y, expand = FALSE, datum = st_crs(4326)) +
  scale_x_continuous(breaks = seq(10, 17, by = 2)) +
  scale_y_continuous(breaks = seq(40, 46, by = 2)) +
  
  labs(title = NULL, x = "Longitude (°)", y = "Latitude (°)") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 14),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3, linetype = "dashed"),
    panel.grid.minor = element_blank()
  )

print(p_se)
ggsave("Figure2_Uncertainties.png", plot = p_se, width = 12, height = 6, dpi = 300, bg = "white")

print(">>> Figure generate e salvate con successo!")

################################################################################
# SCRIPT FINALE: PLOT 2x2 PREVISIONI E INCERTEZZE (PGA e PGV | RDD e Non-Param)
################################################################################

library(ggplot2)
library(viridis)
library(dplyr)
library(sf)
library(rnaturalearth)
library(scales)

print(">>> Avvio generazione Figure 2x2: RDD vs Non-Parametrico...")

# ==============================================================================
# 1. FUNZIONE ESTRATTRICE ROBUSTA
# ==============================================================================
estrai_dati <- function(path, nome_modello, nome_im) {
  df <- readRDS(path)$mappa
  
  # LA CORREZIONE: Arrotondiamo a 3 decimali per allineare la griglia del modello Stazionario
  df$X_km <- round(df$X_km, 3)
  df$Y_km <- round(df$Y_km, 3)
  
  nomi_col <- names(df)
  
  # 1. Trova dinamicamente la colonna di PREVISIONE
  if ("dS2S_kriging" %in% nomi_col) {
    df$Pred <- df$dS2S_kriging
  } else if ("Pred_Kriging" %in% nomi_col) {
    df$Pred <- df$Pred_Kriging
  } else if ("pred" %in% nomi_col) {
    df$Pred <- df$pred
  } else {
    idx <- grep("pred|ds2s", nomi_col, ignore.case = TRUE)[1]
    df$Pred <- df[[idx]]
  }
  
  # 2. Trova dinamicamente la colonna di INCERTEZZA
  if ("se_kriging" %in% nomi_col) {
    df$SE <- df$se_kriging
  } else if ("SE_Kriging" %in% nomi_col) {
    df$SE <- df$SE_Kriging
  } else if ("std_dev" %in% nomi_col) {
    df$SE <- df$std_dev
  } else {
    idx <- grep("se_|std", nomi_col, ignore.case = TRUE)[1]
    df$SE <- df[[idx]]
  }
  
  # 3. Seleziona le colonne e aggiunge le etichette
  df_final <- df %>% 
    select(X_km, Y_km, Pred, SE) %>% 
    mutate(Modello = nome_modello, IM = nome_im)
  
  return(df_final)
}

# ==============================================================================
# 2. CARICAMENTO DATI (Modifica i path se necessario)
# ==============================================================================

path_pga_rdd <- "C:/Users/miche/OneDrive/Desktop/prova3/INGV/Risultati_RDD_Kopt.rds"
path_pga_np  <- "C:/Users/miche/OneDrive/Desktop/prova3/INGV/Risultati_NonParam.rds"

path_pgv_rdd <- "C:/Users/miche/OneDrive/Desktop/prova3/INGV - pgv/Risultati_RDD_Kopt.rds"
path_pgv_np  <- "C:/Users/miche/OneDrive/Desktop/prova3/INGV - pgv/Risultati_NonParam.rds"

df_pga_rdd <- estrai_dati(path_pga_rdd, "RDD", "PGA")
df_pga_np  <- estrai_dati(path_pga_np,  "Non-Parametric", "PGA")

df_pgv_rdd <- estrai_dati(path_pgv_rdd, "RDD", "PGV")
df_pgv_np  <- estrai_dati(path_pgv_np,  "Non-Parametric", "PGV")

# Unione in un unico grande dataset
mappa_all <- bind_rows(df_pga_rdd, df_pga_np, df_pgv_rdd, df_pgv_np)

# Forziamo l'ordine dei fattori per la griglia del plot
mappa_all$Modello <- factor(mappa_all$Modello, levels = c("RDD", "Non-Parametric"))
mappa_all$IM <- factor(mappa_all$IM, levels = c("PGA", "PGV"))

# ==============================================================================
# 3. MASCHERATURA ITALIA SENZA SARDEGNA
# ==============================================================================
ita_regioni <- ne_states(country = "Italy", returnclass = "sf")
ita_no_sardegna <- ita_regioni %>% filter(!grepl("Sardegn", region, ignore.case = TRUE))

ita_border  <- st_union(ita_no_sardegna)
italia_utm  <- st_transform(ita_border, 32633)

# Creazione oggetto spaziale per il filtro
sf_all <- st_as_sf(mappa_all, coords = c("X_km", "Y_km"), crs = 32633, remove = FALSE)
st_geometry(sf_all) <- st_geometry(sf_all) * 1000
st_crs(sf_all) <- 32633

inside <- st_intersects(sf_all, italia_utm, sparse = FALSE)
mappa_plot <- mappa_all[as.vector(inside), ]

# Calcolo limiti spaziali e risoluzione
limiti_x <- range(mappa_plot$X_km * 1000, na.rm = TRUE)
limiti_y <- range(mappa_plot$Y_km * 1000, na.rm = TRUE)
RISOLUZIONE <- (sort(unique(mappa_plot$X_km))[2] - sort(unique(mappa_plot$X_km))[1]) * 1000

# ==============================================================================
# 4. PLOT 1: PREVISIONI 2x2 (facet_grid)
# ==============================================================================
lim_pred <- quantile(mappa_plot$Pred, probs = c(0.01, 0.99), na.rm = TRUE)

p_pred <- ggplot(mappa_plot) +
  geom_sf(data = italia_utm, fill = "white", color = "grey80") +
  geom_tile(aes(x = X_km * 1000, y = Y_km * 1000, fill = Pred), 
            width = RISOLUZIONE, height = RISOLUZIONE) +
  geom_sf(data = italia_utm, fill = NA, color = "black", linewidth = 0.3) +
  
  # LA MAGIA: facet_grid divide le righe per IM e le colonne per Modello
  facet_grid(IM ~ Modello) +
  scale_fill_viridis_c(
    option = "viridis", 
    name = expression(delta*"S2S"~Prediction), 
    limits = lim_pred, 
    oob = scales::squish,
    # Aggiungi questa riga per accorciare la barra delle previsioni
    guide = guide_colorbar(barheight = unit(3, "cm"), title.vjust = 1) ) +

  coord_sf(xlim = limiti_x, ylim = limiti_y, expand = FALSE, datum = st_crs(4326)) +
  scale_x_continuous(breaks = seq(10, 17, by = 2)) +
  scale_y_continuous(breaks = seq(40, 46, by = 2)) +
  
  labs(x = "Longitude (°)", y = "Latitude (°)") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 14, face = "bold"), # Testo Intestazioni
    axis.title = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 16),
    legend.key.height = unit(1.5, "cm"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3, linetype = "dashed"), 
    panel.grid.minor = element_blank()
  )

print(p_pred)
#ggsave("Final_Predictions_2x2.pdf", plot = p_pred, width = 12, height = 10, dpi = 300, bg = "white")

# ==============================================================================
# 5. PLOT 2: INCERTEZZE 2x2 (facet_grid)
# ==============================================================================
library(patchwork)

print(">>> Generazione Mappe di Incertezza (Legende unificate per Modello)...")

# ==============================================================================
# 1. Calcolo dei limiti GLOBALI divisi per MODELLO (Colonne)
# ==============================================================================
lim_rdd <- range(mappa_plot$SE[mappa_plot$Modello == "RDD"], na.rm = TRUE)
#lim_np  <- range(mappa_plot$SE[mappa_plot$Modello == "Non-Parametric"], na.rm = TRUE)
min_np <- min(mappa_plot$SE[mappa_plot$Modello == "Non-Parametric"], na.rm = TRUE)
max_np <- quantile(mappa_plot$SE[mappa_plot$Modello == "Non-Parametric"], probs = 0.99, na.rm = TRUE)
lim_np <- c(min_np, max_np)

# ==============================================================================
# 2. Definizione delle due Scale (una per RDD, una per Non-Parametrico)
# ==============================================================================
altezza_barra <- unit(3, "cm")

scala_rdd <- scale_fill_viridis_c(option = "magma", direction = -1, 
                                  name = "Standard Error", 
                                  limits = lim_rdd, oob = scales::squish,
                                  guide = guide_colorbar(barheight = altezza_barra, title.vjust = 1))

scala_np  <- scale_fill_viridis_c(option = "magma", direction = -1, 
                                  name = "Standard Error", 
                                  limits = lim_np, oob = scales::squish,
                                  guide = guide_colorbar(barheight = altezza_barra, title.vjust = 1))

# ==============================================================================
# 3. Funzione base del Plot (Modificata per accomodare l'asse destro)
# ==============================================================================
crea_plot_base <- function(dati) {
  ggplot(dati) +
    geom_sf(data = italia_utm, fill = "white", color = "grey80") +
    geom_tile(aes(x = X_km * 1000, y = Y_km * 1000, fill = SE), 
              width = RISOLUZIONE, height = RISOLUZIONE) +
    geom_sf(data = italia_utm, fill = NA, color = "black", linewidth = 0.3) +
    coord_sf(xlim = limiti_x, ylim = limiti_y, expand = FALSE, datum = st_crs(4326)) +
    scale_x_continuous(breaks = seq(10, 17, by = 2)) +
    # Rimuoviamo scale_y_continuous da qui, lo definiamo nei singoli plot per gestire il lato destro
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5), # Stile per "RDD" e "Non-Parametric"
      axis.title.x = element_text(size = 12),
      axis.title.y.left = element_text(size = 12),
      
      # --- TRUCCO PER "PGA" e "PGV" A DESTRA ---
      axis.title.y.right = element_text(size = 15, face = "bold", angle = 270, margin = margin(l = 10)),
      axis.text.y.right = element_blank(),  # Nasconde i numeri a destra
      axis.ticks.y.right = element_blank(), # Nasconde le tacchette a destra
      # ----------------------------------------
      
      axis.text = element_text(size = 10),
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.key.height = unit(1.2, "cm"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3, linetype = "dashed"),
      panel.grid.minor = element_blank()
    )
}

# ==============================================================================
# 4. Assemblaggio dei 4 Pannelli (Simulando l'estetica di facet_grid)
# ==============================================================================

# Filtriamo i dataset
d_pga_rdd <- subset(mappa_plot, IM == "PGA" & Modello == "RDD")
d_pga_np  <- subset(mappa_plot, IM == "PGA" & Modello == "Non-Parametric")
d_pgv_rdd <- subset(mappa_plot, IM == "PGV" & Modello == "RDD")
d_pgv_np  <- subset(mappa_plot, IM == "PGV" & Modello == "Non-Parametric")

# --- COLONNA SINISTRA: RDD ---
# Titolo in alto, asse Y standard
p1 <- crea_plot_base(d_pga_rdd) + scala_rdd + 
  scale_y_continuous(breaks = seq(40, 46, by = 2)) +
  labs(title = "RDD", y = "Latitude (°)", x = "")

p3 <- crea_plot_base(d_pgv_rdd) + scala_rdd + 
  scale_y_continuous(breaks = seq(40, 46, by = 2)) +
  labs(title = "", y = "Latitude (°)", x = "Longitude (°)")

# --- COLONNA DESTRA: Non-Parametric ---
# Titolo in alto, Asse Y secondario (destro) per scrivere PGA/PGV
p2 <- crea_plot_base(d_pga_np) + scala_np + 
  scale_y_continuous(breaks = seq(40, 46, by = 2), sec.axis = sec_axis(~ ., name = "PGA")) +
  labs(title = "Non-Parametric", y = "", x = "")

p4 <- crea_plot_base(d_pgv_np) + scala_np + 
  scale_y_continuous(breaks = seq(40, 46, by = 2), sec.axis = sec_axis(~ ., name = "PGV")) +
  labs(title = "", y = "", x = "Longitude (°)")


# ==============================================================================
# 5. Costruzione della griglia con patchwork e unione legende
# ==============================================================================

# Uniamo in verticale e raccogliamo la legenda (una sola legenda per tutta la colonna RDD)
colonna_rdd <- (p1 / p3) + plot_layout(guides = "collect")

# Uniamo in verticale e raccogliamo la legenda (una sola legenda per tutta la colonna NP)
colonna_np  <- (p2 / p4) + plot_layout(guides = "collect")

# Uniamo le due colonne in orizzontale
layout_se <- (colonna_rdd | colonna_np)

print(layout_se)


################################################################################
# MAPPA DELLE DIFFERENZE SPAZIALI (RDD vs Non-Parametrico)
################################################################################

library(ggplot2)
library(dplyr)
library(sf)
library(rnaturalearth)
library(scales)

print(">>> Avvio generazione Mappa Differenze (RDD - Non-Parametrico)...")

# ==============================================================================
# 1. FUNZIONE ESTRATTRICE ROBUSTA (Eredita dai tuoi ultimi script)
# ==============================================================================
estrai_pred <- function(path, nome_modello) {
  df <- readRDS(path)$mappa
  
  # Allineamento griglia
  df$X_km <- round(df$X_km, 3)
  df$Y_km <- round(df$Y_km, 3)
  
  nomi_col <- names(df)
  
  # Trova dinamicamente la colonna di PREVISIONE
  if ("dS2S_kriging" %in% nomi_col) {
    df$Pred <- df$dS2S_kriging
  } else if ("Pred_Kriging" %in% nomi_col) {
    df$Pred <- df$Pred_Kriging
  } else if ("pred" %in% nomi_col) {
    df$Pred <- df$pred
  } else {
    idx <- grep("pred|ds2s", nomi_col, ignore.case = TRUE)[1]
    df$Pred <- df[[idx]]
  }
  
  df_final <- df %>% 
    select(X_km, Y_km, Pred) %>% 
    rename_with(~ paste0("Pred_", nome_modello), Pred)
  
  return(df_final)
}

# ==============================================================================
# 2. CARICAMENTO DATI E CALCOLO DELTA
# ==============================================================================
# Sostituisci i percorsi con quelli del tuo PC (qui sto usando PGA come esempio)
path_rdd <- "C:/Users/miche/OneDrive/Desktop/prova3/INGV/Risultati_RDD_Kopt.rds"
path_np  <- "C:/Users/miche/OneDrive/Desktop/prova3/INGV/Risultati_NonParam.rds"

df_rdd <- estrai_pred(path_rdd, "RDD")
df_np  <- estrai_pred(path_np, "NP")

# Unione dei due dataframe tramite coordinate e calcolo della differenza
mappa_diff <- merge(df_rdd, df_np, by = c("X_km", "Y_km")) %>%
  mutate(Delta = Pred_RDD - Pred_NP) # Differenza: RDD - Non-Parametrico

# ==============================================================================
# 3. MASCHERATURA ITALIA SENZA SARDEGNA
# ==============================================================================
ita_regioni <- ne_states(country = "Italy", returnclass = "sf")
ita_no_sardegna <- ita_regioni %>% filter(!grepl("Sardegn", region, ignore.case = TRUE))

ita_border  <- st_union(ita_no_sardegna)
italia_utm  <- st_transform(ita_border, 32633)

# Creazione oggetto spaziale per il filtro in metri
sf_diff <- st_as_sf(mappa_diff, coords = c("X_km", "Y_km"), crs = 32633, remove = FALSE)
st_geometry(sf_diff) <- st_geometry(sf_diff) * 1000
st_crs(sf_diff) <- 32633

# Teniamo solo le celle interne alla penisola
inside <- st_intersects(sf_diff, italia_utm, sparse = FALSE)
mappa_plot <- mappa_diff[as.vector(inside), ]

# Calcolo limiti spaziali e risoluzione
limiti_x <- range(mappa_plot$X_km * 1000, na.rm = TRUE)
limiti_y <- range(mappa_plot$Y_km * 1000, na.rm = TRUE)
RISOLUZIONE <- (sort(unique(mappa_plot$X_km))[2] - sort(unique(mappa_plot$X_km))[1]) * 1000

# Limite massimo assoluto per centrare lo zero nella scala colori
lim_max <- max(abs(mappa_plot$Delta), na.rm = TRUE)

# ==============================================================================
# 4. PLOT DELLA MAPPA DELLE DIFFERENZE
# ==============================================================================
p_delta <- ggplot(mappa_plot) +
  # Sfondo Italia
  geom_sf(data = italia_utm, fill = "white", color = "grey80") +
  
  # Rasterizzazione della griglia delle differenze
  geom_tile(aes(x = X_km * 1000, y = Y_km * 1000, fill = Delta), 
            width = RISOLUZIONE, height = RISOLUZIONE) +
  
  # Disegno dei bordi esterni
  geom_sf(data = italia_utm, fill = NA, color = "black", linewidth = 0.3) +
  
  # Scala colori: freddi (negativi), bianco (zero), caldi (positivi)
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limits = c(-lim_max, lim_max), 
    name = expression(Delta~"dS2S"~"(RDD - NP)"),
    oob = scales::squish,
    guide = guide_colorbar(barheight = unit(4, "cm"))
  ) +
  
  # Gestione coordinate e reticolo Lat/Lon
  coord_sf(xlim = limiti_x, ylim = limiti_y, expand = FALSE, datum = st_crs(4326)) +
  scale_x_continuous(breaks = seq(10, 17, by = 2)) +
  scale_y_continuous(breaks = seq(40, 46, by = 2)) +
  
  labs(
    title = "Differenza Predittiva Spaziale",
    subtitle = "Colori Caldi: RDD stima valori > Non-Parametrico\nColori Freddi: RDD stima valori < Non-Parametrico",
    x = "Longitude (°)", y = "Latitude (°)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 12),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3, linetype = "dashed"), 
    panel.grid.minor = element_blank()
  )

print(p_delta)

# Se vuoi salvarla de-commenta la riga sotto
# ggsave("Mappa_Differenza_RDD_vs_NP.png", plot = p_delta, width = 8, height = 8, dpi = 300, bg = "white")

################################################################################
# MAPPA DELLE DIFFERENZE SPAZIALI (Stazionario vs RDD)
################################################################################

print(">>> Avvio generazione Mappa Differenze (Stazionario - RDD)...")

# 1. Caricamento ed estrazione dati Stazionario
path_staz <- "C:/Users/miche/OneDrive/Desktop/prova3/INGV/Risultati_Stazionario.rds"
df_staz <- estrai_pred(path_staz, "Staz")

# 2. Calcolo Delta: Stazionario - RDD
mappa_diff_staz_rdd <- merge(df_staz, df_rdd, by = c("X_km", "Y_km")) %>%
  mutate(Delta = Pred_Staz - Pred_RDD) 

# 3. Mascheratura spaziale
sf_diff_staz_rdd <- st_as_sf(mappa_diff_staz_rdd, coords = c("X_km", "Y_km"), crs = 32633, remove = FALSE)
st_geometry(sf_diff_staz_rdd) <- st_geometry(sf_diff_staz_rdd) * 1000
st_crs(sf_diff_staz_rdd) <- 32633

inside_staz_rdd <- st_intersects(sf_diff_staz_rdd, italia_utm, sparse = FALSE)
mappa_plot_staz_rdd <- mappa_diff_staz_rdd[as.vector(inside_staz_rdd), ]

# Limite massimo per la scala colori
lim_max_staz_rdd <- max(abs(mappa_plot_staz_rdd$Delta), na.rm = TRUE)

# 4. Plot Stazionario vs RDD
p_delta_staz_rdd <- ggplot(mappa_plot_staz_rdd) +
  geom_sf(data = italia_utm, fill = "white", color = "grey80") +
  geom_tile(aes(x = X_km * 1000, y = Y_km * 1000, fill = Delta), 
            width = RISOLUZIONE, height = RISOLUZIONE) +
  geom_sf(data = italia_utm, fill = NA, color = "black", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limits = c(-lim_max_staz_rdd, lim_max_staz_rdd), 
    name = expression(Delta~"dS2S"~"(Staz - RDD)"),
    oob = scales::squish, guide = guide_colorbar(barheight = unit(4, "cm"))
  ) +
  coord_sf(xlim = limiti_x, ylim = limiti_y, expand = FALSE, datum = st_crs(4326)) +
  scale_x_continuous(breaks = seq(10, 17, by = 2)) +
  scale_y_continuous(breaks = seq(40, 46, by = 2)) +
  labs(
    title = "Differenza Predittiva Spaziale (Stazionario vs RDD)",
    subtitle = "Colori Caldi: Stazionario > RDD | Colori Freddi: Stazionario < RDD",
    x = "Longitude (°)", y = "Latitude (°)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 12),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3, linetype = "dashed"), 
    panel.grid.minor = element_blank()
  )

print(p_delta_staz_rdd)
# ggsave("Mappa_Differenza_Staz_vs_RDD.png", plot = p_delta_staz_rdd, width = 8, height = 8, dpi = 300, bg = "white")


################################################################################
# MAPPA DELLE DIFFERENZE SPAZIALI (Stazionario vs Non-Parametrico)
################################################################################

print(">>> Avvio generazione Mappa Differenze (Stazionario - Non-Parametrico)...")

# 1. Calcolo Delta: Stazionario - Non-Parametrico (df_np e df_staz sono già caricati)
mappa_diff_staz_np <- merge(df_staz, df_np, by = c("X_km", "Y_km")) %>%
  mutate(Delta = Pred_Staz - Pred_NP)

# 2. Mascheratura spaziale
sf_diff_staz_np <- st_as_sf(mappa_diff_staz_np, coords = c("X_km", "Y_km"), crs = 32633, remove = FALSE)
st_geometry(sf_diff_staz_np) <- st_geometry(sf_diff_staz_np) * 1000
st_crs(sf_diff_staz_np) <- 32633

inside_staz_np <- st_intersects(sf_diff_staz_np, italia_utm, sparse = FALSE)
mappa_plot_staz_np <- mappa_diff_staz_np[as.vector(inside_staz_np), ]

# Limite massimo per la scala colori
lim_max_staz_np <- max(abs(mappa_plot_staz_np$Delta), na.rm = TRUE)

# 3. Plot Stazionario vs Non-Parametrico
p_delta_staz_np <- ggplot(mappa_plot_staz_np) +
  geom_sf(data = italia_utm, fill = "white", color = "grey80") +
  geom_tile(aes(x = X_km * 1000, y = Y_km * 1000, fill = Delta), 
            width = RISOLUZIONE, height = RISOLUZIONE) +
  geom_sf(data = italia_utm, fill = NA, color = "black", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limits = c(-lim_max_staz_np, lim_max_staz_np), 
    name = expression(Delta~"dS2S"~"(Staz - NP)"),
    oob = scales::squish, guide = guide_colorbar(barheight = unit(4, "cm"))
  ) +
  coord_sf(xlim = limiti_x, ylim = limiti_y, expand = FALSE, datum = st_crs(4326)) +
  scale_x_continuous(breaks = seq(10, 17, by = 2)) +
  scale_y_continuous(breaks = seq(40, 46, by = 2)) +
  labs(
    title = "Differenza Predittiva Spaziale (Stazionario vs NP)",
    subtitle = "Colori Caldi: Stazionario > NP | Colori Freddi: Stazionario < NP",
    x = "Longitude (°)", y = "Latitude (°)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 12),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3, linetype = "dashed"), 
    panel.grid.minor = element_blank()
  )

print(p_delta_staz_np)
# ggsave("Mappa_Differenza_Staz_vs_NP.png", plot = p_delta_staz_np, width = 8, height = 8, dpi = 300, bg = "white")
