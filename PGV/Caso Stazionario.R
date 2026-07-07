
# Modello Stazionario

setwd("C:/Users/miche/OneDrive/Desktop/prova3/INGV - pgv")
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

mod <- lmer(res ~ 1 + (1 | EVENT), data = itacentrale, REML = TRUE)
dBe <- ranef(mod)$EVENT[, "(Intercept)"]
fitted_vals <- fitted(mod)
offset <- fixef(mod)[1]
dWs <- res - fitted_vals
itacentrale$dWs = dWs

# NUOVO Modello ausiliario (serve SOLO per stimare phi_0, il rumore puro)
# Aggiungiamo (1 | NET_STA) per catturare la varianza di sito separatamente
mod_sep <- lmer(res ~ 1 + (1 | EVENT) + (1 | NET_STA), data = itacentrale, REML = TRUE)

# Estraiamo la varianza dei residui di questo modello, ovvero varianza del rumore
phi_0_sq <- sigma(mod_sep)^2 

# Sostituiamo la variabile var_intra con questo nuovo valore
var_intra <- phi_0_sq 

# Controllo (opzionale): Stampa la differenza
cat("Varianza Totale (S2S + Rumore):", sigma(mod)^2, "\n")
cat("Varianza Rumore Puro (phi_0^2):", phi_0_sq, "\n")

# Al Atik: var(dWs) = var(dS2S) + var(eps)

#---------------------------------------------------------------------
# 0. SETUP: CARICARE TUTTE LE LIBRERIE NECESSARIE
#---------------------------------------------------------------------
# Assicurati di averle installate: 
# install.packages(c("sf", "gstat", "dplyr", "sp", "ggplot2"))

library(sf)      # Per gestire i dati spaziali (Proiezione)
library(dplyr)   # Per la manipolazione dei dati (Pairing)
library(sp)      # Per la funzione spDists (Calcolo distanze)
library(gstat)   # Per fit.variogram, vgm
library(ggplot2) # Per il plot finale (opzionale)

#---------------------------------------------------------------------
# 1. PREPARAZIONE DATI (PRO: Rigore Geospaziale con SF e UTM)
#---------------------------------------------------------------------

# 1. Converti in 'sf' usando Lat/Lon (CRS 4326)
itacentrale_sf <- st_as_sf(itacentrale, 
                           coords = c("st_lon", "st_lat"), 
                           crs = 4326,
                           remove = FALSE) # Tieni le colonne originali

# 2. Trasforma in un sistema proiettato (es. UTM 33N per l'Italia)
#    Questo è FONDAMENTALE per calcolare distanze Euclidee corrette (in metri)
cat("Proiezione dei dati in UTM 33N (EPSG:32633)...\n")
itacentrale_utm <- st_transform(itacentrale_sf, crs = 32633)

# 3. Estrai le coordinate X/Y proiettate e prepara per il join
# Sezione 1: Preparazione Dati
itacentrale_df <- itacentrale_utm %>%
  mutate(
    x_utm = st_coordinates(.)[,1],
    y_utm = st_coordinates(.)[,2]
  ) %>%
  st_drop_geometry() %>%
  mutate(
    NET_STA = as.character(NET_STA),
    EVENT   = as.character(EVENT)
  ) %>%
  filter(!is.na(x_utm) & !is.na(dWs) & !is.na(res)) %>% # <-- Aggiunto res
  select(EVENT, NET_STA, dWs, res, x_utm, y_utm) # <-- Aggiunto res

cat("Preparazione dati completata.\n")

#---------------------------------------------------------------------
# 2. CREAZIONE COPPIE (PRO: Efficienza con dplyr::inner_join)
#---------------------------------------------------------------------

print("Inizio creazione coppie (metodo inner_join)...")
df_pairs <- inner_join(
  itacentrale_df, 
  itacentrale_df, 
  by = "EVENT",
  suffix = c("_i", "_j"),
  relationship = "many-to-many"
) %>%
  # Filtro per coppie uniche
  filter(NET_STA_i < NET_STA_j)

cat("Coppie create:", nrow(df_pairs), "\n")

#---------------------------------------------------------------------
# 3. CALCOLO VARIOGRAM CLOUD (Distanze Euclidee in Metri)
#---------------------------------------------------------------------

print("Inizio calcolo distanze (h) e semivarianza (gamma)...")

# Estrai le coordinate UTM in matrici
coords_i <- as.matrix(df_pairs[, c("x_utm_i", "y_utm_i")])
coords_j <- as.matrix(df_pairs[, c("x_utm_j", "y_utm_j")])

# Calcola le distanze Euclidee (in METRI) usando le coordinate proiettate
# spDists con longlat=FALSE è il modo corretto per farlo su UTM
df_for_vario <- data.frame(
  h_m = spDists(coords_i, coords_j, longlat = FALSE, diagonal = TRUE),
  gamma = 0.5 * (df_pairs$dWs_i - df_pairs$dWs_j)^2
) %>%
  na.omit()

cat("Calcolo Cloud completato.\n")
summary(df_for_vario$h_m)
quantile(df_for_vario$h_m, probs = c(0.90, 0.95, 0.99))

#---------------------------------------------------------------------
# 4. BINNING (PRO: Robusto con Quantile e Mediana)
#---------------------------------------------------------------------
print("Inizio binning (metodo robusto: quantili + mediana)...")

# Definisci un cutoff in METRI (es. 230 km = 230000 m)
MAX_DIST_METRI <- 300000

df_cut <- df_for_vario %>% filter(h_m > 0 & h_m <= MAX_DIST_METRI)

# 1. Binning basato su quantili
boundaries <- quantile(df_cut$h_m, probs = seq(0, 1, length.out = 11), na.rm = TRUE)
# Assicurati che i boundaries siano unici
boundaries <- unique(boundaries)
boundaries[1] <- 0 # Forza l'inizio a 0

exp_vario_df <- df_cut %>%
  mutate(dist_bin = cut(h_m, breaks = boundaries, include.lowest = TRUE)) %>%
  filter(!is.na(dist_bin)) %>%
  group_by(dist_bin) %>%
  summarize(
    # Distanza media del bin per l'asse X
    dist = mean(h_m, na.rm = TRUE),
    # Mediana della semivarianza (Pro da Script 1)
    gamma = median(gamma, na.rm = TRUE), 
    np = as.numeric(n())
  ) %>%
  ungroup() %>%
  filter(np > 0) # Rimuovi bin vuoti

#---------------------------------------------------------------------
# 5. FORMATTAZIONE PER GSTAT (Il "Trucco")
#---------------------------------------------------------------------
print("Applico il 'trucco' per gstat...")

exp_vario <- as.data.frame(exp_vario_df)
exp_vario$dist_bin <- NULL
class(exp_vario) <- c("gstatVariogram", "data.frame")
attr(exp_vario, "model") <- vgm(model = "Nug")

print("Variogramma sperimentale (binned):")
print(exp_vario)

#---------------------------------------------------------------------
# 6. FIT E PLOT
#---------------------------------------------------------------------
print("Inizio fit (metodo gstat)...")

# Stime iniziali
nugget_start <- 0.02
sill_start <- 0.04
range_start_metri <- 250000 # (200 km)

vgm_sph <- vgm(psill = sill_start-nugget_start, model = "Sph", range = range_start_metri, nugget = nugget_start)
vgm_exp <- vgm(psill = sill_start-nugget_start, model = "Exp", range = range_start_metri/3, nugget = nugget_start)
vgm_mat <- vgm(psill = sill_start-nugget_start, model = "Mat", range = range_start_metri/4, nugget = nugget_start, kappa = 1)

# Esegui il fit
fit_sph <- try(fit.variogram(exp_vario, model = vgm_sph, fit.method = 7), silent = TRUE)
fit_exp <- try(fit.variogram(exp_vario, model = vgm_exp, fit.method = 7), silent = TRUE)
fit_mat <- try(fit.variogram(exp_vario, model = vgm_mat, fit.method = 7), silent = TRUE)

print("--- Modello Sferico Fittato ---")
print(fit_sph)
print("--- Modello Esponenziale Fittato ---")
print(fit_exp)
print("--- Modello Matérn Fittato ---")
print(fit_mat)

# --- PLOT (Base R, come Script 1) ---

maxdist <- max(exp_vario$dist, na.rm = TRUE)

plot(
  exp_vario$dist,
  exp_vario$gamma,
  main = "Variogramma Ibrido (UTM + Mediana)",
  xlab = "Distanza h (metri)", # Nota: l'unità è cambiata!
  ylab = "Semivarianza (gamma)",
  pch = 19,
  col = "black",
  cex = 1.5)
grid()

# Aggiungi linee fittate
if(!inherits(fit_sph, "try-error")) {
  lines(variogramLine(fit_sph, maxdist = maxdist, n = 200), col = "red", lwd = 2)}

if(!inherits(fit_exp, "try-error")) {
  lines(variogramLine(fit_exp, maxdist = maxdist, n = 200), col = "blue", lwd = 2)}

if(!inherits(fit_mat, "try-error")) {
  lines(variogramLine(fit_mat, maxdist = maxdist, n = 200), col = "forestgreen", lwd = 2)}

legend(
  "bottomright",
  legend = c("Sperimentale (Mediana)", "Sferico Fittato", "Esponenziale Fittato", "Matérn Fittato"),
  pch = c(19, NA, NA, NA),
  lty = c(NA, 1, 1, 1),
  lwd = c(NA, 2, 2, 2),
  col = c("black", "red", "blue", "forestgreen"),
  cex = 1)

#---------------------------------------------------------------------
# 9. CALCOLO MATRICE DI COVARIANZA (METODO UNIVERSALE)
#---------------------------------------------------------------------

stazioni_uniche <- itacentrale_utm %>%
  select(NET_STA, geometry) %>%     # Seleziona solo ID e Coordinate
  distinct(NET_STA, .keep_all = TRUE) # Tieni solo la prima occorrenza per ogni stazione

cat("Numero totale di stazioni uniche trovate:", nrow(stazioni_uniche), "\n")

# 1. Calcola la Matrice delle Distanze (in metri)
#    Usa la colonna 'geometry' che hai nel dataset
dist_matrix_sf <- st_distance(stazioni_uniche)

# Converti in matrice numerica pura (rimuove le unità di misura)
dist_matrix <- matrix(as.numeric(dist_matrix_sf), 
                      nrow = nrow(dist_matrix_sf), 
                      ncol = ncol(dist_matrix_sf))

# 2. Assegna i nomi delle stazioni (NET_STA)
#    Questo rende la matrice leggibile
rownames(dist_matrix) <- stazioni_uniche$NET_STA
colnames(dist_matrix) <- stazioni_uniche$NET_STA

# 3. Calcola la Covarianza usando gstat::variogramLine
#    Questo è il trucco potente: 'covariance = TRUE' trasforma il variogramma in covarianza.
#    Funziona con fit_exp, fit_sph, fit_mat, ecc.

# Scegli il modello che vuoi usare (es. fit_exp)
MODELLO_SCELTO <- fit_exp 

# Calcola i valori
cov_values <- variogramLine(
  MODELLO_SCELTO, 
  dist_vector = as.vector(dist_matrix), # Appiattisce la matrice in un vettore lungo
  covariance = TRUE                     # <--- CHIEDE LA COVARIANZA, NON IL VARIOGRAMMA
)

# 4. Ricostruisci la Matrice di Covarianza
#    'cov_values' restituisce un dataframe con colonna 'gamma' (che qui contiene la covarianza)
cov_matrix <- matrix(cov_values$gamma, 
                     nrow = nrow(dist_matrix), 
                     ncol = ncol(dist_matrix))

# Riassegna i nomi
rownames(cov_matrix) <- rownames(dist_matrix)
colnames(cov_matrix) <- colnames(dist_matrix)

# Stampa le prime 5x5 stazioni
print(cov_matrix[1:5, 1:5])

mat_for_plot <- cov_matrix[nrow(cov_matrix):1, ] 

# Plot (PESANTISSIMO)
print(
  levelplot(
    mat_for_plot,
    main = "Matrice di Covarianza Spaziale",
    xlab = "Stazione i",
    ylab = "Stazione j",
    col.regions = hcl.colors(50, "viridis"), 
    scales = list( 
      x = list(lab = NULL),  # <-- CORREZIONE: Rimuove solo il testo
      y = list(lab = NULL)   # <-- CORREZIONE: Rimuove solo il testo
    ),
    colorkey = list( 
      space = "right")))

# ---------------------------------------------------------------------
# 9-BIS. MATRICE DI COVARIANZA TEORICA DEL SEGNALE PURO (dS2S)
# ---------------------------------------------------------------------

cat("Generazione Matrice Covarianza Segnale (Senza Nugget)...\n")

# 1. Definizione del Modello per dS2S (Solo Segnale)
#    Partiamo dal modello fittato sui dWs (fit_sph)
fit_site_pure <- fit_exp

#    AZZERIAMO IL NUGGET.
#    Il nugget rappresenta il rumore intra-evento (errore record-to-record).
#    Il termine di sito dS2S è la componente stabile, quindi ha nugget ~ 0.
fit_site_pure$psill[1] <- 0  # Azzera la componente Nugget (assumendo sia la prima riga)

#    Nota: Se il modello ha una struttura complessa, assicurati che la riga 1 sia il Nugget.
#    Di solito gstat mette il Nugget alla riga 1 e le strutture spaziali alla riga 2.
print("Modello per il Segnale Puro (dS2S):")
print(fit_site_pure)

# 2. Calcolo della Covarianza
#    Usiamo le stesse distanze calcolate nella Sezione 9 (dist_matrix)
if(!exists("dist_matrix")) stop("Esegui prima la Sezione 9 per avere 'dist_matrix'")

cov_values_signal <- variogramLine(
  fit_site_pure, 
  dist_vector = as.vector(dist_matrix), 
  covariance = TRUE
)

# 3. Ricostruzione Matrice
Sigma_signal_pure <- matrix(cov_values_signal$gamma, 
                            nrow = nrow(dist_matrix), 
                            ncol = ncol(dist_matrix))

rownames(Sigma_signal_pure) <- rownames(dist_matrix)
colnames(Sigma_signal_pure) <- colnames(dist_matrix)

# 4. Plot (Levelplot) (PESANTISSIMO)
#    Usiamo gli stessi colori della Sezione 9 per confronto diretto
mat_for_plot_sig <- Sigma_signal_pure[nrow(Sigma_signal_pure):1, ] 

print(
  levelplot(
    mat_for_plot_sig,
    main = "Matrice Covarianza Segnale Puro (dS2S)",
    xlab = "Stazione i",
    ylab = "Stazione j",
    col.regions = hcl.colors(50, "viridis"), 
    scales = list(x = list(draw = FALSE), y = list(draw = FALSE)), # Nasconde etichette
    colorkey = list(space = "right")
  )
)

# 5. CONFRONTO NUMERICO (Diagonale)
var_totale <- fit_exp$psill[1] + fit_exp$psill[2] # Nugget + Sill Parziale
var_segnale <- fit_site_pure$psill[1] + fit_site_pure$psill[2] # Solo Sill Parziale

cat("\n--- CONFRONTO VARIANZE ---\n")
cat("Varianza Totale (dWs) sulla diagonale: ", round(var_totale, 4), "\n")
cat("Varianza Segnale (dS2S) sulla diagonale:", round(var_segnale, 4), "\n")
cat("La differenza è il rumore rimosso.\n")

# ---------------------------------------------------------------------
# 10. CALCOLO COVARIANZA EMPIRICA (DATA-DRIVEN)
# ---------------------------------------------------------------------
print("Calcolo Matrice di Covarianza Empirica dai residui grezzi...")

library(tidyr)

# 1. RISTRUTTURAZIONE DEI DATI (PIVOTING)
# Dobbiamo creare una matrice dove:
# - Le RIGHE sono gli EVENTI
# - Le COLONNE sono le STAZIONI
# - Le CELLE contengono il valore dWs

matrix_data <- itacentrale_df %>%
  select(EVENT, NET_STA, dWs) %>%
  # Trasforma in formato "largo" (wide)
  pivot_wider(names_from = NET_STA, values_from = dWs) %>%
  # Rimuoviamo la colonna EVENT perché serve solo per allineare le righe
  select(-EVENT)

# Convertiamo in matrice pura per efficienza
mat_dWs <- as.matrix(matrix_data)

cat("Dimensioni matrice (Eventi x Stazioni):", dim(mat_dWs), "\n")

# 2. CALCOLO DELLA COVARIANZA
# Il trucco è l'argomento 'use = "pairwise.complete.obs"'.
# Poiché non tutte le stazioni registrano tutti gli eventi, la matrice è piena di NA.
# Questo comando dice a R: "Per calcolare la covarianza tra Stazione A e Stazione B,
# usa SOLO gli eventi che ENTRAMBE hanno registrato".

emp_cov_matrix <- cov(mat_dWs, use = "pairwise.complete.obs")

# 3. GESTIONE DEI VALORI MANCANTI (NA)
# Se due stazioni non hanno MAI registrato lo stesso evento in comune,
# la loro covarianza sarà NA. Possiamo lasciarla NA o metterla a 0.
# (Per visualizzazione mettiamo 0, per calcoli matematici NA è più onesto)

cat("Numero di celle NA nella matrice di covarianza:", sum(is.na(emp_cov_matrix)), "\n")
emp_cov_matrix_plot <- emp_cov_matrix
emp_cov_matrix_plot[is.na(emp_cov_matrix_plot)] <- 0 

# 4. PLOT DELLA MATRICE EMPIRICA
# Usiamo levelplot come prima per confrontare visivamente

library(lattice)

# Ruotiamo per il plot (come nel codice precedente)
mat_for_plot_emp <- emp_cov_matrix_plot[nrow(emp_cov_matrix_plot):1, ]

# Trasforma la matrice in un vettore per vedere la distribuzione
valori <- as.vector(mat_for_plot_emp)

print("Riassunto statistico dei valori nella matrice:")
summary(valori)

print("99esimo percentile (per vedere se ci sono outlier estremi):")
print(quantile(valori, 0.99, na.rm=TRUE))

limiti_colore <- quantile(as.vector(mat_for_plot_emp), probs = c(0.01, 0.99), na.rm = TRUE)
if(diff(limiti_colore) == 0) limiti_colore <- c(-0.01, 0.01)
breaks_seq <- seq(from = limiti_colore[1], to = limiti_colore[2], length.out = 51)

# Plot (PESANTISSIMO)
print(
  levelplot(
    mat_for_plot_emp,
    main = "Matrice di Covarianza Empirica (Scala Saturata)",
    xlab = "Stazioni",
    ylab = "Stazioni",
    col.regions = hcl.colors(50, "viridis"), # O prova "Spectral" o "RdBu"
    at = breaks_seq,                         # <--- QUESTO È IL TRUCCO
    scales = list(draw = FALSE),
    colorkey = list(space = "right", 
                    labels = list(at = pretty(breaks_seq))) # Etichette leggibili
  )
)

# 5. CONFRONTO DIAGONALE (VARIANZA DELLE SINGOLE STAZIONI)
# La diagonale di questa matrice rappresenta la varianza reale osservata 
# di ogni singola stazione su tutti i suoi eventi.

diag_var <- diag(emp_cov_matrix)
hist(diag_var, main = "Istogramma delle Varianze per Stazione", 
     xlab = "Varianza (dWs)", breaks = 20, col = "grey")
abline(v = mean(diag_var, na.rm=TRUE), col="red", lwd=2)

#---------------------------------------------------------------------
# 11. CONFRONTO TRA LE DUE MATRICI
#---------------------------------------------------------------------

# 1. Estrai i valori dalla matrice EMPIRICA (Sez 10)
#    Nota: emp_cov_matrix ha NA, li togliamo dopo
cov_emp_vec <- as.vector(emp_cov_matrix)

# 2. Estrai i valori dalla matrice TEORICA (Sez 9)
#    ATTENZIONE: Assicurati che l'ordine delle stazioni sia identico!
#    Per sicurezza, ricostruiamo la teorica usando GLI STESSI NOMI dell'empirica
stazioni_nomi <- colnames(emp_cov_matrix) 

# Coordinate delle stazioni nell'ordine della matrice empirica
coords_match <- itacentrale_utm %>%
  filter(NET_STA %in% stazioni_nomi) %>%
  distinct(NET_STA, .keep_all=TRUE) %>%
  arrange(factor(NET_STA, levels = stazioni_nomi)) # Forza l'ordine

# Calcola distanze per queste stazioni ordinate
dists_check <- st_distance(st_as_sf(coords_match))

# Calcola Covarianza Teorica per queste distanze
cov_theo_vals <- variogramLine(fit_exp, dist_vector = as.numeric(dists_check), covariance = TRUE)
cov_theo_vec <- cov_theo_vals$gamma

# Ricostruiamo la Matrice Teorica Quadrata (per il plot)
Sigma_theo_check <- matrix(cov_theo_vec, 
                           nrow = length(stazioni_nomi), 
                           ncol = length(stazioni_nomi))
rownames(Sigma_theo_check) <- stazioni_nomi
colnames(Sigma_theo_check) <- stazioni_nomi

# 3. PLOT DI CONFRONTO: Empirico vs Distanza (con sopra il Modello)
df_check <- data.frame(
  dist = as.numeric(dists_check),
  cov_emp = cov_emp_vec,
  cov_theo = cov_theo_vec
) %>% filter(!is.na(cov_emp)) # Rimuovi i NA dell'empirica

# 4. HEATMAP PER CONFRONTO VISIVO

cat("Generazione Heatmap di confronto...\n")

# Funzione per convertire matrice in dataframe per ggplot
matrix_to_df <- function(mat, label) {
  # reshape2::melt o tidyr::pivot_longer manuale
  df <- as.data.frame(as.table(mat))
  colnames(df) <- c("Sta1", "Sta2", "Covarianza")
  df$Tipo <- label
  return(df)
}

# Creiamo i dataframe
df_mat_emp <- matrix_to_df(emp_cov_matrix, "1. Matrice Empirica")
df_mat_theo <- matrix_to_df(Sigma_theo_check, "2. Matrice Teorica (fit_exp)")

# Uniamo tutto
df_mat_combo <- rbind(df_mat_emp, df_mat_theo)

# CALCOLO LIMITI COMUNI (Per avere la stessa scala colori)
# Ignoriamo NA e diagonali troppo alte se necessario
limiti_cov <- quantile(df_mat_combo$Covarianza, probs = c(0.01, 0.99), na.rm = TRUE)

print(
  ggplot(df_mat_combo, aes(x = Sta1, y = Sta2, fill = Covarianza)) +
    # Usiamo geom_raster o geom_tile. 
    # raster è più veloce per matrici grandi.
    geom_raster() + 
    
    # Dividiamo i due grafici
    facet_wrap(~ Tipo) +
    
    # Scala Colori Unificata
    scale_fill_viridis_c(
      option = "plasma",  # Plasma è ottimo per evidenziare le correlazioni
      name = "Covarianza",
      limits = limiti_cov,
      oob = scales::squish # Schiaccia i valori fuori scala (es. diagonale)
    ) +
    
    # Cosmetica
    labs(title = "Confronto Struttura di Covarianza",
         subtitle = "Sinistra: Dati Grezzi (Rumorosi) | Destra: Modello (Liscio)",
         x = "Stazioni", y = "Stazioni") +
    
    theme_minimal() +
    theme(
      axis.text = element_blank(), # Nascondiamo i nomi stazioni (illeggibili se tanti)
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      aspect.ratio = 1 # Forza formato quadrato
    )
)

cat("Confronto matrici completato.\n")

# Plot (PESANTISSIMO)
plot(df_check$dist, df_check$cov_emp,
     main = "Coerenza: Covarianza Empirica vs Modello",
     xlab = "Distanza (m)", ylab = "Covarianza",
     pch = 19, col = rgb(0,0,0,0.1), cex = 0.5,
     ylim = c(-0.05, max(df_check$cov_emp, na.rm=TRUE)))

# Aggiungi la linea del modello teorico (dovrebbe passare in mezzo)
points(df_check$dist, df_check$cov_theo, col="red", pch=".", cex=2)

legend("topright", legend=c("Dati Empirici (Matrice Sez 10)", "Modello Teorico (Matrice Sez 9)"),
       col=c("black", "red"), pch=c(19, 19))
grid()

# PLOT (PESANTISSIMO) "BISETTRICE": Covarianza Empirica vs Teorica
plot(df_check$cov_theo, df_check$cov_emp,
     main = "Confronto Diretto: Teorico vs Empirico",
     xlab = "Covarianza Teorica (Modello)",
     ylab = "Covarianza Empirica (Dati)",
     pch = 19, col = rgb(0,0,0,0.1), cex = 0.5,
     xlim = c(-0.02, max(df_check$cov_theo)),
     ylim = c(-0.02, max(df_check$cov_emp, na.rm=TRUE)))

# Qui aggiungi la bisettrice
abline(0, 1, col="red", lwd=2) 
grid()

# SCATTERPLOT BINNATO (CON BARRE DI ERRORE)

library(dplyr)
library(ggplot2)

cat("Generazione Scatterplot Binnato...\n")

# A. DEFINIZIONE DEI BIN (Intervalli)
# Raggruppiamo le distanze ogni 20 km (20000 metri).
# Puoi cambiare questo valore: 10000 per più dettagli, 30000 per più pulizia.
BIN_WIDTH <- 20000 

# Creiamo il dataframe per il binning
df_binned <- df_check %>%
  filter(!is.na(cov_emp)) %>%
  mutate(
    # Calcola il centro del bin (es. distanze 0-20km -> bin_center = 10km)
    bin_center = (floor(dist / BIN_WIDTH) * BIN_WIDTH + BIN_WIDTH/2) / 1000 # convertito in km
  ) %>%
  group_by(bin_center) %>%
  summarize(
    cov_mean = mean(cov_emp, na.rm = TRUE),       # Media della covarianza nel bin
    cov_sd   = sd(cov_emp, na.rm = TRUE),         # Deviazione standard (variabilità)
    n_pairs  = n(),                               # Numero di coppie nel bin
    .groups = 'drop'
  )

# B. GENERAZIONE LINEA DEL MODELLO TEORICO
# Creiamo una linea liscia per il modello (fit_exp) da sovrapporre
dummy_dist <- seq(0, max(df_check$dist), length.out = 200)
dummy_cov  <- variogramLine(fit_exp, dist_vector = dummy_dist, covariance = TRUE)
df_model_line <- data.frame(dist_km = dummy_dist / 1000, gamma = dummy_cov$gamma)

# C. PLOT
print(
  ggplot() +
    # 1. BARRE DI ERRORE (Grigie)
    # Rappresentano la variabilità dei dati (Media +/- 1 Deviazione Standard)
    # Se la linea rossa passa dentro le barre, il modello è statisticamente buono.
    geom_errorbar(data = df_binned, 
                  aes(x = bin_center, ymin = cov_mean - cov_sd, ymax = cov_mean + cov_sd),
                  width = 2, color = "grey60", alpha = 0.6) +
    
    # 2. PUNTI MEDI (Neri)
    # Il "cuore" dei dati empirici.
    geom_point(data = df_binned, 
               aes(x = bin_center, y = cov_mean, size = n_pairs), 
               color = "black", fill = "black", shape = 21) +
    
    # 3. LINEA MODELLO (Rossa)
    geom_line(data = df_model_line, 
              aes(x = dist_km, y = gamma), 
              color = "red", size = 1.2) +
    
    # Cosmetica
    labs(title = "Validazione Modello: Covarianza Binnata",
         subtitle = sprintf("Binning ogni %d km | Barre = ±1 Dev.Std.", BIN_WIDTH/1000),
         x = "Distanza (km)",
         y = "Covarianza",
         size = "N. Coppie") +
    
    scale_size_continuous(range = c(1, 6)) + # Dimensione pallini basata su quante coppie ci sono
    theme_minimal() +
    theme(legend.position = "bottom")
)

# ---------------------------------------------------------------------
# 13. KRIGING AVANZATO: PESATURA SU NUMERO DI EVENTI (N_i)
# Chilès, J. P., & Delfiner, P. (2012). Geostatistics: Modeling Spatial Uncertainty.
# Approccio: Kriging of the mean with measurement error.
# ---------------------------------------------------------------------
# Obiettivo: Non buttare via le stazioni con pochi eventi, ma usarle
#            con un "peso" minore (permettendo al modello di non passare per forza su di loro).
# Metodo:    Kriging con termine di Errore di Misura (Nugget Effect locale).

print("Inizio Smoothing Kriging (pesato su N eventi)...")

# 1. STIMA DELLA VARIANZA DEL "RUMORE" INTRA-EVENTO
#    Calcoliamo la varianza globale dei residui.
#    Assumiamo che questa varianza contenga sia il segnale che il rumore.
#    vecchia varianza: var_intra <- var(itacentrale_df$dWs, na.rm = TRUE)

#stima REML (Restricted Maximum Likelihood) della varianza residua del modello a effetti misti

cat("Varianza Rumore Intra-Evento:", round(var_intra, 4), "\n")

# 2. PREPARAZIONE DEI DATI DI SITO
site_terms_all <- itacentrale_df %>%
  group_by(NET_STA) %>%
  summarize(
    dS2S = mean(dWs, na.rm = TRUE),
    n_eventi = n(), 
    x_utm = mean(x_utm),
    y_utm = mean(y_utm)
  ) %>%
  ungroup()

# Calcoliamo l'Errore di Misura (Varianza dell'errore) per ogni stazione:
# Var_err = Var_intra / N
site_terms_all$meas_error_var <- var_intra / site_terms_all$n_eventi

# Non usiamo scaling empirico. Usiamo la decomposizione del variogramma 'fit_exp'.
# fit_exp = Nugget (Rumore) + Partial Sill (Segnale S2S).
# Per il Kriging del segnale, dobbiamo usare un modello con Nugget = 0.
# Creiamo una copia del modello globale
# Usiamo la covarianza fittata e non empirica perché:
# Definita positiva
# Dobbiamo creare il termine noto (b) su tutta la griglia e quindi dovremmo interpolare (instabile)

fit_site <- fit_exp

# Impostiamo il Nugget a 0 per isolare solo la continuità spaziale.
fit_site$psill[1] <- 0 

cat("Modello 'fit_site' per il segnale (Nugget imposto a 0):\n")
print(fit_site)

# 3. CREAZIONE DEL SISTEMA DI KRIGING MANUALE (Weighted)

# Coordinate delle stazioni
coords_staz <- site_terms_all[, c("x_utm", "y_utm")]
dist_matrix_sf <- st_distance(st_as_sf(coords_staz, coords=c("x_utm","y_utm"), crs=32633))
dist_matrix <- matrix(as.numeric(dist_matrix_sf), nrow=nrow(coords_staz))

# A. Calcolo Covarianza Spaziale (SOLO SEGNALE)
# Usiamo fit_site che ha nugget=0. Quindi C(0) = Partial Sill.
# C(h) = Sill - Gamma(h)
cov_values <- variogramLine(fit_site, dist_vector = as.vector(dist_matrix), covariance = TRUE)
Sigma_signal <- matrix(cov_values$gamma, nrow = nrow(dist_matrix))

# B. AGGIUNTA DELL'ERRORE DI MISURA (SOLO RUMORE)
# Aggiungiamo alla diagonale la varianza d'errore specifica (Var_intra / N).
# Ora la diagonale contiene: Segnale + Rumore specifico.
# Diag_totale = C(0) + Var_errore(i)
Sigma_total <- Sigma_signal
diag(Sigma_total) <- diag(Sigma_total) + site_terms_all$meas_error_var

library(maps) 

# 1. PREPARAZIONE SHAPEFILE ITALIA (Ci serve solo per gli sfondi dei plot)
italia_sf <- st_as_sf(maps::map("world", regions = "Italy", plot = FALSE, fill = TRUE))
st_crs(italia_sf) <- 4326             
italia_utm <- st_transform(italia_sf, 32633) 

# 2. CREAZIONE GRIGLIA (PERFETTAMENTE ALLINEATA A RDD E LOESS)
# Lavoriamo in metri (UTM), quindi step 5 km = 5000 m, buffer 20 km = 20000 m
RISOLUZIONE <- 7000
buffer_m <- 20000

min_X <- min(itacentrale_df$x_utm) - buffer_m
max_X <- max(itacentrale_df$x_utm) + buffer_m
min_Y <- min(itacentrale_df$y_utm) - buffer_m
max_Y <- max(itacentrale_df$y_utm) + buffer_m

span_X <- max_X - min_X
span_Y <- max_Y - min_Y
max_span <- max(span_X, span_Y)

mid_X <- (max_X + min_X) / 2
mid_Y <- (max_Y + min_Y) / 2

min_X <- mid_X - max_span / 2
max_X <- mid_X + max_span / 2
min_Y <- mid_Y - max_span / 2
max_Y <- mid_Y + max_span / 2

seq_X <- seq(min_X, max_X, by = RISOLUZIONE)
seq_Y <- seq(max_Y, min_Y, by = -RISOLUZIONE)
N_side <- max(length(seq_X), length(seq_Y))

seq_X <- seq(min_X, by = RISOLUZIONE, length.out = N_side)
seq_Y <- seq(max_Y, by = -RISOLUZIONE, length.out = N_side)

grid_plain <- expand.grid(col_idx = 0:(N_side-1), row_idx = 0:(N_side-1))
grid_plain$x_utm <- seq_X[grid_plain$col_idx + 1]
grid_plain$y_utm <- seq_Y[grid_plain$row_idx + 1]

# Assicuriamo lo stesso ordine a serpentina
grid_plain$serpentine_idx <- ifelse(
  grid_plain$row_idx %% 2 == 0,
  grid_plain$row_idx * N_side + grid_plain$col_idx + 1,          
  grid_plain$row_idx * N_side + (N_side - 1 - grid_plain$col_idx) + 1 
)
grid_plain <- grid_plain[order(grid_plain$serpentine_idx), ]

grid_coords <- as.matrix(grid_plain[, c("x_utm", "y_utm")])
cat("Griglia allineata con RDD generata:", nrow(grid_coords), "punti totali\n")

# 3. CICLO DI KRIGING
preds_weighted <- numeric(nrow(grid_coords))
preds_var      <- numeric(nrow(grid_coords))

# Il Sill totale del segnale (Partial Sill originale)
sill_signal_value <- sum(fit_site$psill) 

cat("Inizio Kriging (Previsione + Varianza)... \n")

# Pre-inversione della matrice A
N <- nrow(site_terms_all)
A_mat <- rbind(cbind(Sigma_total, 1), 1)
A_mat[N+1, N+1] <- 0
A_inv <- solve(A_mat) 

# Loop sui punti filtrati
for(i in 1:nrow(grid_coords)) {
  
  pt_x <- grid_coords[i, 1]
  pt_y <- grid_coords[i, 2]
  
  # Distanze punto-stazioni
  dists_pt <- sqrt((site_terms_all$x_utm - pt_x)^2 + (site_terms_all$y_utm - pt_y)^2)
  
  # Vettore b (Covarianza Segnale tra punto ignoto e stazioni)
  # Usiamo fit_site (nugget=0)
  cov_target_vals <- variogramLine(fit_site, dist_vector = dists_pt, covariance = TRUE)
  b_vec <- c(cov_target_vals$gamma, 1)
  
  # Risoluzione pesi lambda
  weights_mu <- A_inv %*% b_vec
  weights <- weights_mu[1:N]      
  mu_lagrange <- weights_mu[N+1]  
  
  # Calcolo Previsione
  preds_weighted[i] <- sum(weights * site_terms_all$dS2S)
  
  # Calcolo Varianza Kriging
  # Sigma^2 = C_segnale(0) - Sum(lambda * C_segnale(h)) - mu
  kriging_var <- sill_signal_value - sum(weights * b_vec[1:N]) - mu_lagrange
  
  preds_var[i] <- max(0, kriging_var)
  
  if(i %% 100 == 0) cat(sprintf("\rProgresso: %.1f%%", i/nrow(grid_coords)*100))
}
cat("\nCalcolo completato.\n")

# Calcola i limiti per il plot
padding <- 10000 
min_x <- min(site_terms_all$x_utm) - padding
max_x <- max(site_terms_all$x_utm) + padding
min_y <- min(site_terms_all$y_utm) - padding
max_y <- max(site_terms_all$y_utm) + padding
limiti_x <- c(min_x, max_x)
limiti_y <- c(min_y, max_y)

# --- PLOT PREVISIONE ---

df_map_weighted <- data.frame(
  x = grid_coords[,1], 
  y = grid_coords[,2], 
  pred = preds_weighted)

print(
  ggplot() +
    geom_sf(data = italia_utm, fill = "grey95", color = "grey50") +
    geom_tile(data = df_map_weighted, 
              aes(x=x, y=y, fill=pred), 
              width = RISOLUZIONE, height = RISOLUZIONE) +
    geom_point(data = site_terms_all, aes(x=x_utm, y=y_utm, size=n_eventi), 
               shape = 21, fill = "white", stroke = 0.5) +
    scale_fill_viridis_c(option = "viridis", name = "dWs medio") +
    scale_size_continuous(range = c(1, 5), name = "Numero di Eventi") +
    coord_sf(xlim = limiti_x, ylim = limiti_y, expand = FALSE) +
    labs(title = "Mappa di Previsione di dWs (Site Term)") +
    theme_minimal() +
    theme(axis.title = element_blank())
)

# --- PLOT INCERTEZZA ---

df_map_error <- data.frame(
  x = grid_coords[,1], 
  y = grid_coords[,2], 
  std_dev = sqrt(preds_var)
)

print(
  ggplot() +
    geom_sf(data = italia_utm, fill = "white", color = "grey80") +
    geom_tile(data = df_map_error, 
              aes(x = x, y = y, fill = std_dev), 
              width = RISOLUZIONE, height = RISOLUZIONE) +
    geom_point(data = site_terms_all, aes(x=x_utm, y=y_utm), 
               shape = 21, size = 0.5, fill = "black", alpha = 0.5) +
    scale_fill_viridis_c(option = "magma", 
                         name = "Dev. Std. (sigma)",
                         direction = -1) + 
    coord_sf(xlim = limiti_x, ylim = limiti_y, expand = FALSE) +
    labs(title = "Mappa dell'Incertezza di Stima") +
    theme_minimal() +
    theme(axis.title = element_blank())
  #    theme(panel.grid = element_blank(),
  #          axis.title = element_blank(),
  #          axis.text = element_blank(),
  #          axis.ticks = element_blank())
)

# ---------------------------------------------------------------------
# 14. SMOOTHING
# Wiener, N. (1949)
# Rasmussen & Williams 2006, Gaussian Processes for Machine Learning
# Extrapolation, Interpolation, and Smoothing of Stationary Time Series
# Stafford, 2014; Chilès & Delfiner, 2012
# ---------------------------------------------------------------------
# Obiettivo: Raffinare i termini di sito (dS2S) usando il modello
#            variografico teorico (fit_exp) per correggere le stazioni
#            poco affidabili basandosi sui vicini.

print("Inizio Smoothing basato su Modello Teorico (Isotropo)...")

# --- 1. PREPARAZIONE DATI ---

# Calcoliamo le statistiche per ogni stazione
site_terms_temp <- itacentrale_df %>%
  group_by(NET_STA) %>%
  summarize(
    dS2S_obs = mean(dWs, na.rm = TRUE), 
    n_eventi = n(),
    meas_error = var_intra / n()
  ) %>%
  ungroup()

# Aggiungiamo le coordinate e rimuoviamo chi non le ha
site_terms_smooth <- site_terms_temp %>%
  inner_join(site_terms_all %>% select(NET_STA, x_utm, y_utm), by = "NET_STA") %>%
  na.omit() 

# --- 2. COSTRUZIONE MATRICI (con fit_site) ---

# A. Matrice delle Distanze tra stazioni
coords_sf <- st_as_sf(site_terms_smooth, coords = c("x_utm", "y_utm"), crs = 32633)
dist_mat_smooth <- st_distance(coords_sf)

# B. Matrice di Covarianza Segnale (Sigma_signal) basata su fit_exp
cov_vals_smooth <- variogramLine(fit_site, 
                                 dist_vector = as.numeric(dist_mat_smooth), 
                                 covariance = TRUE)

Sigma_signal <- matrix(cov_vals_smooth$gamma, 
                       nrow = nrow(site_terms_smooth), 
                       ncol = nrow(site_terms_smooth))

# C. Matrice del Sistema (Sigma_sys = Segnale + Errore Misura)
Sigma_sys <- Sigma_signal
diag(Sigma_sys) <- diag(Sigma_sys) + site_terms_smooth$meas_error

# D. Inversione (Stabile)
Sigma_inv <- solve(Sigma_sys)

# --- 3. CALCOLO SMOOTHING (Valore e Varianza) ---

# A. Calcolo del Valore Smoothed (Wiener Filter)
# Formula: Smoothed = C * (C + R)^-1 * Z_obs
weights_vector <- Sigma_inv %*% site_terms_smooth$dS2S_obs
dS2S_smooth_vals <- Sigma_signal %*% weights_vector

# B. [NUOVO] Calcolo della Varianza Posteriore (Incertezza residua)
# Formula: Var_post = C - C * (C + R)^-1 * C
# Rappresenta quanto l'incertezza è scesa grazie ai dati vicini
reduction_term <- Sigma_signal %*% Sigma_inv %*% Sigma_signal
var_post_matrix <- Sigma_signal - reduction_term
sd_smooth_vals <- sqrt(pmax(0, diag(var_post_matrix))) # pmax evita radici negative per errori numerici

# --- 4. SALVATAGGIO RISULTATI ---
df_results <- data.frame(
  NET_STA = site_terms_smooth$NET_STA,
  N = site_terms_smooth$n_eventi,
  Osservato = site_terms_smooth$dS2S_obs,
  Smoothed = as.vector(dS2S_smooth_vals),
  SD_Smooth_Iso = as.vector(sd_smooth_vals) # <--- NUOVA COLONNA
)
df_results$Delta <- df_results$Smoothed - df_results$Osservato

print("Prime 10 stazioni corrette:")
print(head(df_results, 10))

# SALVATAGGIO NEL DATAFRAME PRINCIPALE (site_terms_all)
# Pulizia preventiva
if("Smoothed" %in% names(site_terms_all)) site_terms_all$Smoothed <- NULL
if("SD_Smooth_Iso" %in% names(site_terms_all)) site_terms_all$SD_Smooth_Iso <- NULL

# Join
site_terms_all <- site_terms_all %>%
  left_join(df_results %>% select(NET_STA, Smoothed, SD_Smooth_Iso), by = "NET_STA")

cat("Salvataggio completato. Aggiunte colonne 'Smoothed' e 'SD_Smooth_Iso'.\n")

# --- 5. PLOT VALIDAZIONE (Scatterplot) ---
library(ggplot2)
print(
  ggplot(df_results, aes(x = Osservato, y = Smoothed)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
    geom_point(aes(size = N, color = abs(Delta)), alpha = 0.7) +
    scale_color_viridis_c(option = "magma", name = "|Correzione|") +
    scale_size_continuous(range = c(1, 6)) +
    labs(
      title = "Validazione Smoothing (Modello Teorico)",
      subtitle = "Correzione Bayesiana basata su fit_exp",
      x = "Media Osservata (dS2S grezzo)",
      y = "Media Pesata (dS2S smoothed)"
    ) +
    theme_minimal()
)

# --- 6. PLOT MAPPA COMBINATA PREVISIONE ---
# (Il codice del plot resta identico al tuo, usa i nuovi dati)
df_punti_plot <- site_terms_all %>% filter(!is.na(Smoothed))
tutti_i_valori <- c(df_map_weighted$pred, df_punti_plot$Smoothed)
limiti_robusti <- quantile(tutti_i_valori, probs = c(0.01, 0.99), na.rm = TRUE)
if(diff(limiti_robusti) < 0.1) limiti_robusti <- c(-0.5, 0.5)

bbox_stazioni <- st_bbox(st_as_sf(site_terms_all, coords = c("x_utm", "y_utm"), crs = 32633))
padding <- 10000 
xlim_plot <- c(bbox_stazioni["xmin"] - padding, bbox_stazioni["xmax"] + padding)
ylim_plot <- c(bbox_stazioni["ymin"] - padding, bbox_stazioni["ymax"] + padding)

print(
  ggplot() +
    geom_sf(data = italia_utm, fill = "white", color = "grey80") +
    geom_tile(data = df_map_weighted, 
              aes(x = x, y = y, fill = pred), 
              width = RISOLUZIONE, height = RISOLUZIONE) +
    geom_point(data = df_punti_plot,
               aes(x = x_utm, y = y_utm, 
                   fill = Smoothed, 
                   size = n_eventi),
               shape = 21, color = "black", stroke = 0.5) +
    scale_fill_viridis_c(option = "viridis", name = "dWs", limits = limiti_robusti, oob = scales::squish) +
    scale_size_continuous(range = c(1, 5), name = "N. Eventi") +
    coord_sf(xlim = xlim_plot, ylim = ylim_plot, expand = FALSE) +
    labs(title = "Mappa Combinata Isotropa", subtitle = "Sfondo: Kriging of the Mean with Error | Pallini: Smoothing") +
    theme_minimal() + theme(panel.grid = element_blank(), axis.title = element_blank())
)

# --- PLOT MAPPA COMBINATA INCERTEZZA ---
# Uniamo i valori della mappa e dei punti per avere una scala colori coerente
valori_iso <- c(df_map_error$std_dev, site_terms_all$SD_Smooth_Iso)
limiti_sigma_iso <- range(valori_iso, na.rm = TRUE)

print(
  ggplot() +
    # Sfondo Italia
    geom_sf(data = italia_utm, fill = "white", color = "grey80") +
    
    # Mappa Raster dell'Errore
    geom_tile(data = df_map_error, 
              aes(x = x, y = y, fill = std_dev), 
              width = RISOLUZIONE, height = RISOLUZIONE) +
    
    # Punti Stazioni (ORA COLORATI CON LA LORO SIGMA)
    geom_point(data = site_terms_all, 
               aes(x = x_utm, y = y_utm, 
                   fill = SD_Smooth_Iso), # <--- Qui usiamo la Sigma smoothed
               shape = 21, 
               size = 3,          # Grandi abbastanza da vedere il colore
               color = "black",   # Bordino nero per staccare
               stroke = 0.5) +
    
    # Scala Colori Magma (Invertita e BLOCCATA sui limiti comuni)
    scale_fill_viridis_c(option = "magma", 
                         name = "Dev. Std.",
                         direction = -1,
                         limits = limiti_sigma_iso, # <--- FONDAMENTALE
                         oob = scales::squish) +    
    
    coord_sf(xlim = limiti_x, ylim = limiti_y, expand = FALSE) +
    
    labs(title = "Mappa dell'Incertezza Isotropo",
         subtitle = "Sfondo: Errore Kriging | Pallini: Smoothing") +
    
    theme_minimal() +
    theme(panel.grid = element_blank(), axis.title = element_blank())
)

# ---------------------------------------------------------------------
# 15. ANALISI DI ANISOTROPIA (VARIOGRAMMI DIREZIONALI)
# ---------------------------------------------------------------------

cat("\nInizio analisi anisotropia (in km)...\n")

# Se le coordinate sono in km, usiamo cutoff e width in km.
# Cutoff = 230 km
# Width = 15 km (bin)

vario_dir <- variogram(dWs ~ 1, 
                       itacentrale_sf, 
                       alpha = c(0, 45, 90, 135), 
                       cutoff = 300,  
                       width = 30)    

# Plot diagnostico
print(
  plot(vario_dir, 
       type = "b", pch = 16, col = "blue",
       main = "Variogrammi Direzionali (0=N, 90=E)",
       xlab = "Distanza (km)",         # <--- Etichetta corretta
       ylab = "Semivarianza")
)

cat("...Analisi completata.\n")

# ---------------------------------------------------------------------
# 16. FITTING DEL MODELLO ANISOTROPO (ISPEZIONE VISIVA + GRID SEARCH 3D)
# ---------------------------------------------------------------------

get_ani_coords <- function(x, y, alpha, ratio) {
  # A. Rotazione (Allinea l'asse maggiore all'asse X)
  x_rot <- x * cos(alpha) - y * sin(alpha)
  y_rot <- x * sin(alpha) + y * cos(alpha)
  
  # B. Scaling (Stira l'asse Y minore per renderlo lungo quanto X)
  # Dividendo per il ratio (<1), la coordinata aumenta -> distanza aumenta -> correlazione scende
  y_scaled <- y_rot / ratio
  
  return(data.frame(x_ani = x_rot, y_ani = y_scaled))
}

# Parametri stimati dal plot della Sezione 15
NUGGET_STIMATO <- 0.05   # Punto di partenza sull'asse Y
SILL_PARZIALE  <- 0.05   # Quanto sale la curva (Sill Totale - Nugget = 0.175 - 0.075)
RANGE_MAJOR    <- 250    # Km (Direzione 135°)
ANGLE_MAJOR    <- 135    # Gradi (Direzione Appenninica)
RATIO          <- 0.5    # Rapporto (Range 45° / Range 135°)

# Definizione del modello anisotropo
# psill totale sarà ~0.175
model_ani <- vgm(psill  = SILL_PARZIALE, 
                 model  = "Sph",       
                 range  = RANGE_MAJOR, 
                 nugget = NUGGET_STIMATO,
                 anis   = c(ANGLE_MAJOR, RATIO))

cat("Modello Anisotropo Definito:\n")
print(model_ani)

# Plot di verifica: Sovrapponiamo il modello ai dati sperimentali
# Questo plot ti mostrerà come il modello si adatta alle 4 direzioni
print(plot(vario_dir, model_ani, main = "Fit Modello Anisotropo (Curve continue)"))

cat("\n--- INIZIO GRID SEARCH ANISOTROPIA (3 PARAMETRI) ---\n")

# 1. Definisci la griglia di ricerca (Occhio a non esagerare con i valori!)
angoli_test <- seq(135, 135, by = 10) 
#angoli_test <- seq(110, 160, by = 10)             # 6 valori
ratio_test  <- seq(0.3, 0.7, by = 0.1)            # 5 valori
range_test  <- c(100000, 125000, 150000, 175000, 200000, 225000, 250000, 275000, 300000)
#range_test  <- c(150000, 200000, 250000, 300000)  # 4 valori (in metri!)


grid_search <- expand.grid(Angle = angoli_test, Ratio = ratio_test, Range = range_test)
grid_search$RMSE <- NA

cat(sprintf("Totale combinazioni da testare: %d\n", nrow(grid_search)))
cat("Questa operazione richiederà tempo. Mettiti comodo...\n")

F_SILL   <- SILL_PARZIALE # Lo teniamo fisso per non far esplodere i calcoli
F_NUGGET <- 0             

run_kfold_cv <- function(data_df, model_signal, coords_cols, meas_err_col, k = 10) {
  
  data_df <- as.data.frame(data_df)
  n_stazioni <- nrow(data_df)
  preds <- numeric(n_stazioni) # Vettore per salvare le predizioni (allineato ai dati originali)
  z_vals <- data_df$dS2S 
  
  # 1. Pre-calcolo della Matrice di Covarianza Segnale COMPLETA
  #    (Esattamente come nel LOOCV, per efficienza)
  coords_mat <- as.matrix(data_df[, coords_cols])
  dist_full <- as.matrix(dist(coords_mat))
  cov_full_list <- variogramLine(model_signal, dist_vector = as.vector(dist_full), covariance = TRUE)
  Sigma_signal_full <- matrix(cov_full_list$gamma, nrow = n_stazioni, ncol = n_stazioni)
  
  # 2. Creazione dei Folds (Casuale ma riproducibile)
  set.seed(123) # Importante per confrontare Isotropo vs Anisotropo sugli stessi gruppi
  folds <- sample(rep(1:k, length.out = n_stazioni))
  
  # 3. Loop sui K Folds
  for(f in 1:k) {
    # cat(sprintf("Fold %d/%d... ", f, k))
    
    # Identifica indici di Training e Test
    idx_test <- which(folds == f)
    idx_train <- setdiff(1:n_stazioni, idx_test)
    
    # --- A. COSTRUZIONE SISTEMA DI TRAINING (Una volta per fold!) ---
    
    # Matrice Covarianza Segnale (Training vs Training)
    Sigma_train_signal <- Sigma_signal_full[idx_train, idx_train]
    
    # Aggiunta Errore di Misura (Solo sulla diagonale del training)
    error_vals_train <- data_df[[meas_err_col]][idx_train]
    R_err_train <- diag(error_vals_train, nrow = length(error_vals_train))
    
    Sigma_train_tot <- Sigma_train_signal + R_err_train
    
    # Matrice Kriging Completa (LHS)
    N_train <- length(idx_train)
    A_mat <- rbind(cbind(Sigma_train_tot, 1), 1)
    A_mat[N_train+1, N_train+1] <- 0
    
    # INVERSIONE MATRICE (L'ottimizzazione chiave del K-Fold)
    # Calcoliamo l'inversa una volta sola e la usiamo per tutti i punti del test set
    A_inv <- try(solve(A_mat), silent = TRUE)
    
    if(inherits(A_inv, "try-error")) {
      cat("ERRORE: Inversione matrice fallita nel fold", f, "\n")
      preds[idx_test] <- NA
      next
    }
    
    # --- B. PREDIZIONE SUI PUNTI DI TEST ---
    for(i_test in idx_test) {
      
      # Vettore b (RHS): Covarianza Segnale tra Punto Test e Punti Training
      # Estraiamo la colonna corrispondente dalla matrice globale pre-calcolata
      cov_target <- Sigma_signal_full[idx_train, i_test]
      b_vec <- c(cov_target, 1)
      
      # Risoluzione Pesi: w = A^-1 * b
      weights_mu <- A_inv %*% b_vec
      weights <- weights_mu[1:N_train]
      
      # Calcolo Predizione
      preds[i_test] <- sum(weights * z_vals[idx_train])
    }
  }
  cat("\nCalcolo K-Fold completato.\n")
  return(preds)
}

# 4. Loop di ottimizzazione
for(i in 1:nrow(grid_search)) {
  
  test_angle <- grid_search$Angle[i]
  test_ratio <- grid_search$Ratio[i]
  test_range <- grid_search$Range[i] # <--- Nuovo parametro estratto!
  
  cat(sprintf("\rTest %d/%d: Angolo=%d°, Ratio=%.1f, Range=%.0f m ... ", 
              i, nrow(grid_search), test_angle, test_ratio, test_range))
  
  # A. Definisci il modello con il Range variabile
  fit_base_ani <- vgm(psill = F_SILL, model = "Sph", range = test_range, nugget = F_NUGGET)
  
  # B. Deformazione Spazio temporanea
  alpha_rad_test <- (test_angle - 90) * pi / 180
  coords_test <- get_ani_coords(site_terms_all$x_utm, site_terms_all$y_utm, alpha_rad_test, test_ratio)
  
  df_test <- site_terms_all
  df_test$x_ani_test <- coords_test$x_ani
  df_test$y_ani_test <- coords_test$y_ani
  
  # C. Esegui 10-Fold CV
  preds_test <- suppressMessages(
    run_kfold_cv(data_df = df_test, model_signal = fit_base_ani, 
                 coords_cols = c("x_ani_test", "y_ani_test"), 
                 meas_err_col = "meas_error_var", k = 10)
  )
  
  # D. Calcolo RMSE
  rmse_corrente <- sqrt(mean((df_test$dS2S - preds_test)^2, na.rm=TRUE))
  grid_search$RMSE[i] <- rmse_corrente
}

cat("\nRicerca completata!\n")

# 5. Estrazione dei risultati migliori
migliore <- grid_search[which.min(grid_search$RMSE), ]

cat("\n--- RISULTATO OTTIMIZZAZIONE 3D ---\n")
cat(sprintf("L'RMSE MINIMO (%.4f) si ottiene con:\n", migliore$RMSE))
cat(sprintf("Angolo Ottimale: %d°\n", migliore$Angle))
cat(sprintf("Ratio Ottimale:  %.2f\n", migliore$Ratio))
cat(sprintf("Range Ottimale:  %.0f metri\n", migliore$Range))

# 6. Sovrascriviamo i parametri globali
ANGLE_MAJOR     <- migliore$Angle
RATIO           <- migliore$Ratio
RANGE_MAJOR     <- migliore$Range 

# ---------------------------------------------------------------------
# 17. PREPARAZIONE DATASET ANISOTROPO
# ---------------------------------------------------------------------

# Chilès, J. P., & Delfiner, P. (2012). Geostatistics: Modeling Spatial Uncertainty. Wiley.
# Obiettivo: Ruotare e "stirare" l'Italia affinché l'ellisse di anisotropia
# diventi un cerchio. In questo spazio deformato, possiamo usare formule isotrope.

# 1. Recupero Parametri dal Modello Anisotropo (definito in Sez 6.B)
# Se non li hai salvati in variabili, li ridefiniamo qui per sicurezza:
# (Questi devono coincidere con quelli che hai visto nel variogramma direzionale)
ANI_RANGE_MAJOR <- 250000  # Range Lungo (Metri) -> Direzione Appenninica
ANI_ANGLE_MAJOR <- 135     # Angolo Azimuth (Gradi) -> 135° SE
ANI_RATIO       <- 0.5     # Rapporto (Range Corto / Range Lungo)
ANI_SILL_SIGNAL <- 0.05
ANI_NUGGET      <- 0
#ANI_NUGGET      <- 0.075   # Nugget (da Sez 6.B)
#ANI_SILL        <- 0.10    # Sill parziale (da Sez 6.B)

cat("--- Parametri Anisotropia ---\n")
cat(sprintf("Direzione: %d° | Ratio: %.2f | Range Max: %.0f m\n", 
            ANI_ANGLE_MAJOR, ANI_RATIO, ANI_RANGE_MAJOR))

# 2. Funzione di Trasformazione (Rotazione + Scaling)
# Converte Azimuth (0=N, 90=E) in Angolo Matematico (0=E, 90=N)
alpha_rad <- (ANI_ANGLE_MAJOR - 90) * pi / 180

# 3. Applicazione alle STAZIONI
coords_ani_staz <- get_ani_coords(site_terms_all$x_utm, site_terms_all$y_utm, alpha_rad, ANI_RATIO)
site_terms_all$x_ani <- coords_ani_staz$x_ani
site_terms_all$y_ani <- coords_ani_staz$y_ani

# 4. Applicazione alla GRIGLIA DI PREVISIONE
coords_ani_grid <- get_ani_coords(grid_coords[,1], grid_coords[,2], alpha_rad, ANI_RATIO)
grid_coords_ani <- as.matrix(coords_ani_grid) # Matrice per velocità nel loop

# 5. Definizione del "Modello Isotropo Equivalente"
# Nello spazio trasformato, usiamo SOLO il Range Maggiore, perché abbiamo "allungato"
# le distanze corte per farle sembrare lunghe.
fit_site_ani_equiv <- vgm(psill  = ANI_SILL_SIGNAL, 
                          model  = "Sph", 
                          range  = ANI_RANGE_MAJOR, 
                          nugget = ANI_NUGGET)

cat("Trasformazione completata. Coordinate '_ani' pronte.\n")

# ---------------------------------------------------------------------
# 18. KRIGING AVANZATO ANISOTROPO: PESATURA SU NUMERO DI EVENTI (N_i)
# ---------------------------------------------------------------------

# Distanze tra stazioni nello spazio deformato
dist_matrix_ani <- as.matrix(dist(site_terms_all[, c("x_ani", "y_ani")]))

# Calcoliamo la covarianza usando il modello equivalente
cov_matrix_val <- variogramLine(fit_site_ani_equiv, dist_vector = as.vector(dist_matrix_ani), covariance = TRUE)
Sigma_emp_ani <- matrix(cov_matrix_val$gamma, nrow = nrow(site_terms_all), ncol = nrow(site_terms_all))

# Aggiungiamo vincoli di Lagrange (Ordinary Kriging)
N <- nrow(site_terms_all)
meas_error_vec <- site_terms_all$meas_error_var
Sigma_total_ani <- Sigma_emp_ani
diag(Sigma_total_ani) <- diag(Sigma_total_ani) + meas_error_vec

A_mat_ani <- rbind(cbind(Sigma_total_ani, 1), 1)
A_mat_ani[N+1, N+1] <- 0

# Inversione della matrice (Fatta una volta sola!)
A_inv_ani <- solve(A_mat_ani)

# 2. Ciclo di Previsione sulla Griglia
preds_ani <- numeric(nrow(grid_coords))
vars_ani  <- numeric(nrow(grid_coords))
sill_total <- sum(fit_site_ani_equiv$psill)

cat("Calcolo Kriging Anisotropo in corso...\n")

for(i in 1:nrow(grid_coords)) {
  
  # PRENDIAMO IL PUNTO TRASFORMATO
  pt_x <- grid_coords_ani[i, 1]
  pt_y <- grid_coords_ani[i, 2]
  
  # Calcolo distanze dal punto a tutte le stazioni (nello spazio deformato)
  dists_pt <- sqrt((site_terms_all$x_ani - pt_x)^2 + (site_terms_all$y_ani - pt_y)^2)
  
  # Vettore b (Covarianza Segnale)
  cov_vals <- variogramLine(fit_site_ani_equiv, dist_vector = dists_pt, covariance = TRUE)
  b_vec <- c(cov_vals$gamma, 1)
  
  # Risoluzione pesi
  weights_mu <- A_inv_ani %*% b_vec
  weights <- weights_mu[1:N]
  mu <- weights_mu[N+1]
  
  # Previsione e Varianza
  preds_ani[i] <- sum(weights * site_terms_all$dS2S)
  vars_ani[i]  <- sill_total - sum(weights * cov_vals$gamma) - mu
  vars_ani[i] <- max(0, vars_ani[i])
  
  if(i %% 100 == 0) cat(sprintf("\rProgresso: %.1f%%", i/nrow(grid_coords)*100))
}
cat("\nCalcolo completato.\n")

# PLOT PREVISIONE

# 3. Creazione Dataframe per Plotting (Usiamo coordinate ORIGINALI per il plot)
df_map_ani <- data.frame(
  x = grid_coords[,1], # Coordinate reali UTM (per disegnare l'Italia dritta)
  y = grid_coords[,2],
  pred = preds_ani,
  var = vars_ani
)

# 4. Visualizzazione Mappa Anisotropa
print(
  ggplot() +
    geom_sf(data = italia_utm, fill = "grey95", color = "grey50") +
    geom_tile(data = df_map_ani, aes(x=x, y=y, fill=pred), width=RISOLUZIONE, height=RISOLUZIONE) +
    scale_fill_viridis_c(option = "viridis", name = "dS2S") +
    labs(title = "Mappa di Previsione Anisotropa (135°)") +
    coord_sf(xlim = limiti_x, ylim = limiti_y, expand = FALSE) +
    theme_minimal() +
    theme(axis.title = element_blank() ) )

# PLOT INCERTEZZA

# 1. Creiamo il dataframe usando i risultati anisotropi (vars_ani)
# Importante: Usiamo x, y ORIGINALI (grid_coords) per disegnare l'Italia dritta,
# ma coloriamo con la varianza calcolata nello spazio deformato.

df_map_error_ani <- data.frame(
  x = grid_coords[,1], 
  y = grid_coords[,2], 
  std_dev = sqrt(vars_ani) # Radice della varianza anisotropa
)

print(
  ggplot() +
    # Sfondo Italia
    geom_sf(data = italia_utm, fill = "white", color = "grey80") +
    
    # Mappa Raster dell'Errore
    geom_tile(data = df_map_error_ani, 
              aes(x = x, y = y, fill = std_dev), 
              width = RISOLUZIONE, height = RISOLUZIONE) +
    
    # Punti Stazioni (Coordinate UTM originali)
    geom_point(data = site_terms_all, aes(x=x_utm, y=y_utm), 
               shape = 21, size = 0.5, fill = "black", alpha = 0.5) +
    
    # Scala Colori Magma (Invertita: Scuro = Alto Errore / Pochi dati)
    scale_fill_viridis_c(option = "magma", 
                         name = "Dev. Std.",
                         direction = -1,
                         oob = scales::squish) + 
    
    # Inquadratura (Usa gli stessi limiti della mappa isotropa)
    coord_sf(xlim = limiti_x, ylim = limiti_y, expand = FALSE) +
    
    labs(title = "Mappa dell'Incertezza Anisotropa") +
    
    theme_minimal() +
    theme(axis.title = element_blank() ) )

# ---------------------------------------------------------------------
# 19. SMOOTHING ANISOTROPO
# --------------------------------------------------------------------- 

cat("\n--- INIZIO SMOOTHING ANISOTROPO ---\n")

# --- 1. PREPARAZIONE MATRICI ---

# A. Distanze Anisotrope
dist_matrix_smooth_ani <- as.matrix(dist(site_terms_all[, c("x_ani", "y_ani")]))

# B. Covarianza Segnale (Modello Anisotropo)
cov_vals_smooth_ani <- variogramLine(fit_site_ani_equiv, 
                                     dist_vector = as.vector(dist_matrix_smooth_ani), 
                                     covariance = TRUE)

Sigma_signal_ani <- matrix(cov_vals_smooth_ani$gamma, 
                           nrow = nrow(site_terms_all), 
                           ncol = nrow(site_terms_all))

# C. Matrice Totale (Segnale + Rumore)
R_err <- diag(site_terms_all$meas_error_var) 
Sigma_sys_ani <- Sigma_signal_ani + R_err

# --- 2. CALCOLO SMOOTHING E VARIANZA (Wiener Filter) ---

cat("Inversione matrice per smoothing anisotropico...\n")
Sigma_inv_ani <- solve(Sigma_sys_ani)

# A. Calcolo Valore Smoothed
weights_ani <- Sigma_inv_ani %*% site_terms_all$dS2S
dS2S_smooth_ani_vals <- Sigma_signal_ani %*% weights_ani

# B. [NUOVO] Calcolo Varianza Posteriore Anisotropa
# Var_post = C - C * (C+R)^-1 * C (tutto con matrici anisotrope)
reduction_term_ani <- Sigma_signal_ani %*% Sigma_inv_ani %*% Sigma_signal_ani
var_post_matrix_ani <- Sigma_signal_ani - reduction_term_ani
sd_smooth_ani_vals <- sqrt(pmax(0, diag(var_post_matrix_ani)))

# --- 3. SALVATAGGIO RISULTATI ---

site_terms_all$Smoothed_Ani <- as.vector(dS2S_smooth_ani_vals)
site_terms_all$SD_Smooth_Ani <- as.vector(sd_smooth_ani_vals) # <--- NUOVA COLONNA

# --- 4. PLOT MAPPA COMBINATA PREVISIONE ---
tutti_i_valori_ani <- c(df_map_ani$pred, site_terms_all$Smoothed_Ani)
limiti_robusti_ani <- quantile(tutti_i_valori_ani, probs = c(0.01, 0.99), na.rm = TRUE)
if(diff(limiti_robusti_ani) < 0.1) limiti_robusti_ani <- c(-0.5, 0.5)

if(!exists("xlim_plot")) {
  bb <- st_bbox(st_as_sf(site_terms_all, coords=c("x_utm","y_utm"), crs=32633))
  xlim_plot <- c(bb["xmin"]-10000, bb["xmax"]+10000)
  ylim_plot <- c(bb["ymin"]-10000, bb["ymax"]+10000)
}

print(
  ggplot() +
    geom_sf(data = italia_utm, fill = "white", color = "grey80") +
    geom_tile(data = df_map_ani, aes(x = x, y = y, fill = pred), width = RISOLUZIONE, height = RISOLUZIONE) +
    geom_point(data = site_terms_all,
               aes(x = x_utm, y = y_utm, fill = Smoothed_Ani, size = n_eventi),
               shape = 21, color = "black", stroke = 0.5, alpha = 0.9) +
    scale_fill_viridis_c(option = "viridis", name = "dS2S Aniso", limits = limiti_robusti_ani, oob = scales::squish) +
    scale_size_continuous(range = c(1, 5), name = "N. Eventi") +
    coord_sf(xlim = xlim_plot, ylim = ylim_plot, expand = FALSE) +
    labs(title = "Mappa Combinata Anisotropa", subtitle = "Sfondo: Kriging of the Mean with Error | Pallini: Smoothing") +
    theme_minimal() + theme(panel.grid = element_blank(), axis.title = element_blank())
)


# --- 5. PLOT MAPPA COMBINATA INCERTEZZA ---
valori_ani <- c(df_map_error_ani$std_dev, site_terms_all$SD_Smooth_Ani)
limiti_sigma_ani <- range(valori_ani, na.rm = TRUE)

print(
  ggplot() +
    # Sfondo Italia
    geom_sf(data = italia_utm, fill = "white", color = "grey80") +
    
    # Mappa Raster dell'Errore Anisotropo
    geom_tile(data = df_map_error_ani, 
              aes(x = x, y = y, fill = std_dev), 
              width = RISOLUZIONE, height = RISOLUZIONE) +
    
    # Punti Stazioni (COLORATI CON LA LORO SIGMA ANISOTROPA)
    geom_point(data = site_terms_all, 
               aes(x = x_utm, y = y_utm, 
                   fill = SD_Smooth_Ani), # <--- Sigma Anisotropa
               shape = 21, 
               size = 3, 
               color = "black", 
               stroke = 0.5) +
    
    # Scala Colori Magma (Invertita e BLOCCATA)
    scale_fill_viridis_c(option = "magma", 
                         name = "Dev. Std. (sigma)",
                         direction = -1,
                         limits = limiti_sigma_ani, # <--- FONDAMENTALE
                         oob = scales::squish) +    
    
    coord_sf(xlim = limiti_x, ylim = limiti_y, expand = FALSE) +
    
    labs(title = "Mappa dell'Incertezza Anisotropa (135°)",
         subtitle = "Sfondo: Errore Kriging | Pallini: Smoothing") +
    
    theme_minimal() +
    theme(panel.grid = element_blank(), axis.title = element_blank())
)

# ---------------------------------------------------------------------
# 20. CONFRONTO FINALE VISIVO: ISOTROPO vs ANISOTROPO
# ---------------------------------------------------------------------
# Affianchiamo le mappe forzando la stessa scala colori.

library(ggplot2)
library(ggpubr) 
library(scales) 

# 1. MAPPA DELLE DIFFERENZE

if("Smoothed" %in% names(site_terms_all)) {
  site_terms_all$Delta_Ani_Iso <- site_terms_all$Smoothed_Ani - site_terms_all$Smoothed
  cat("Delta (Aniso - Iso) calcolato.\n")
} else {
  site_terms_all$Delta_Ani_Iso <- NA
  cat("AVVISO: Colonna 'Smoothed' (Isotropa) assente. Delta non calcolato.\n")
}

cat("Smoothing Anisotropo completato. Risultati in site_terms_all$Smoothed_Ani e $SD_Smooth_Ani\n")

cat("--- Generazione Plot di Confronto (Sezione 20) ---\n")

print(
  ggplot() +
    geom_sf(data = italia_utm, fill = "grey95", color = "grey80") +
    geom_point(data = site_terms_all, 
               aes(x = x_utm, y = y_utm, color = Delta_Ani_Iso, size = n_eventi), alpha = 0.8) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Diff (Aniso - Iso)") +
    scale_size_continuous(range = c(2, 6), name = "N. Eventi") +
    coord_sf(xlim = xlim_plot, ylim = ylim_plot, expand = FALSE) +
    labs(title = "Impatto dell'Anisotropia sulle Stazioni", subtitle = "Differenza tra valore Smoothed Anisotropo e Isotropo") +
    theme_minimal() + theme(panel.grid = element_blank(), axis.title = element_blank())
)

# 2. PREPARAZIONE DATI E CONTROLLI

# Verifichiamo che i calcoli delle sezioni precedenti esistano
if(!exists("df_map_error")) stop("ERRORE: Manca 'df_map_error' (Sez 13).")
if(!exists("df_map_error_ani")) stop("ERRORE: Manca 'df_map_error_ani' (Sez 18).")
if(is.null(site_terms_all$SD_Smooth_Iso)) stop("ERRORE: Manca SD_Smooth_Iso (Sez 14).")
if(is.null(site_terms_all$SD_Smooth_Ani)) stop("ERRORE: Manca SD_Smooth_Ani (Sez 19).")

# Fix nomi colonne previsione se necessario
if(is.null(df_map_weighted$pred)) df_map_weighted$pred <- df_map_weighted$var1.pred
if(is.null(df_map_ani$pred))      df_map_ani$pred      <- df_map_ani$var1.pred

# 3. CALCOLO LIMITI COMUNI (COMMON SCALES)

# A. LIMITI PREVISIONE (dS2S)
# Uniamo tutto: Mappa Iso, Mappa Ani, Punti Iso, Punti Ani
all_preds <- c(df_map_weighted$pred, df_map_ani$pred, 
               site_terms_all$Smoothed, site_terms_all$Smoothed_Ani)

lim_pred_common <- quantile(all_preds, probs = c(0.01, 0.99), na.rm = TRUE)
if(diff(lim_pred_common) < 0.1) lim_pred_common <- c(-0.5, 0.5)

# B. LIMITI INCERTEZZA (Sigma)
# Uniamo tutto: Mappa Iso, Mappa Ani, Punti Iso (Residui), Punti Ani (Residui)
all_sigmas <- c(df_map_error$std_dev, df_map_error_ani$std_dev,
                site_terms_all$SD_Smooth_Iso, site_terms_all$SD_Smooth_Ani)

lim_sigma_common <- range(all_sigmas, na.rm = TRUE)

cat("Scala Colori Previsione:", round(lim_pred_common, 2), "\n")
cat("Scala Colori Incertezza (Min/Max):", round(lim_sigma_common, 3), "\n")

# Recupero Zoom
if(!exists("xlim_plot")) {
  bb <- st_bbox(st_as_sf(site_terms_all, coords=c("x_utm","y_utm"), crs=32633))
  pad <- 10000
  xlim_plot <- c(bb["xmin"]-pad, bb["xmax"]+pad)
  ylim_plot <- c(bb["ymin"]-pad, bb["ymax"]+pad)
}

# 4. GENERAZIONE GRAFICI DI PREVISIONE (dS2S)

# Previsione Isotropo
p_pred_iso <- ggplot() +
  geom_sf(data = italia_utm, fill = "white", color = "grey80") +
  geom_tile(data = df_map_weighted, aes(x=x, y=y, fill=pred), width=RISOLUZIONE, height=RISOLUZIONE) +
  geom_point(data = site_terms_all, aes(x=x_utm, y=y_utm, fill=Smoothed, size=n_eventi),
             shape=21, color="black", stroke=0.5, alpha=0.9) +
  scale_fill_viridis_c(option="viridis", name="dS2S", limits=lim_pred_common, oob=scales::squish) +
  scale_size_continuous(range=c(1,5), guide="none") +
  coord_sf(xlim=xlim_plot, ylim=ylim_plot, expand=FALSE) +
  labs(subtitle="ISOTROPO (Sferico)") +
  theme_minimal() + theme(axis.title=element_blank(), legend.position="bottom")

# Previsione Anisotropo
p_pred_ani <- ggplot() +
  geom_sf(data = italia_utm, fill = "white", color = "grey80") +
  geom_tile(data = df_map_ani, aes(x=x, y=y, fill=pred), width=RISOLUZIONE, height=RISOLUZIONE) +
  geom_point(data = site_terms_all, aes(x=x_utm, y=y_utm, fill=Smoothed_Ani, size=n_eventi),
             shape=21, color="black", stroke=0.5, alpha=0.9) +
  scale_fill_viridis_c(option="viridis", name="dS2S", limits=lim_pred_common, oob=scales::squish) +
  scale_size_continuous(range=c(1,5), guide="none") +
  coord_sf(xlim=xlim_plot, ylim=ylim_plot, expand=FALSE) +
  labs(subtitle=sprintf("ANISOTROPO (Dir: %d°)", ANI_ANGLE_MAJOR)) +
  theme_minimal() + theme(axis.title=element_blank(), legend.position="bottom")

# 5. GENERAZIONE GRAFICI DI INCERTEZZA (Sigma) - COMBINATI

# Incertezza Isotropo (Mappa + Punti colorati)
p_err_iso <- ggplot() +
  geom_sf(data = italia_utm, fill = "white", color = "grey80") +
  
  # Raster
  geom_tile(data = df_map_error, aes(x=x, y=y, fill=std_dev), width=RISOLUZIONE, height=RISOLUZIONE) +
  
  # Punti (Colorati con SD_Smooth_Iso)
  geom_point(data = site_terms_all, aes(x=x_utm, y=y_utm, fill=SD_Smooth_Iso),
             shape=21, color="black", stroke=0.5, size=2.5) +
  
  # Scala Magma Invertita
  scale_fill_viridis_c(option="magma", direction=-1, name="Dev.Std", 
                       limits=lim_sigma_common, oob=scales::squish) +
  
  coord_sf(xlim=xlim_plot, ylim=ylim_plot, expand=FALSE) +
  labs(subtitle="Incertezza ISOTROPA") +
  theme_minimal() + theme(axis.title=element_blank(), legend.position="bottom")

# Incertezza Anisotropo (Mappa + Punti colorati)
p_err_ani <- ggplot() +
  geom_sf(data = italia_utm, fill = "white", color = "grey80") +
  
  # Raster
  geom_tile(data = df_map_error_ani, aes(x=x, y=y, fill=std_dev), width=RISOLUZIONE, height=RISOLUZIONE) +
  
  # Punti (Colorati con SD_Smooth_Ani)
  geom_point(data = site_terms_all, aes(x=x_utm, y=y_utm, fill=SD_Smooth_Ani),
             shape=21, color="black", stroke=0.5, size=2.5) +
  
  # Stessa Scala
  scale_fill_viridis_c(option="magma", direction=-1, name="Dev.Std", 
                       limits=lim_sigma_common, oob=scales::squish) +
  
  coord_sf(xlim=xlim_plot, ylim=ylim_plot, expand=FALSE) +
  labs(subtitle="Incertezza ANISOTROPA") +
  theme_minimal() + theme(axis.title=element_blank(), legend.position="bottom")

# 6. OUTPUT FINALE (AFFIANCAMENTO)

cat("Creazione pannelli affiancati...\n")

# Pannello Previsioni
fig_preds <- ggarrange(p_pred_iso, p_pred_ani, 
                       ncol=2, nrow=1, common.legend=TRUE, legend="right")
fig_preds <- annotate_figure(fig_preds, 
                             top = text_grob("Confronto Previsioni", face = "bold", size=14))

# Pannello Incertezze
fig_errs <- ggarrange(p_err_iso, p_err_ani, 
                      ncol=2, nrow=1, common.legend=TRUE, legend="right")
fig_errs <- annotate_figure(fig_errs, 
                            top = text_grob("Confronto Incertezze", face = "bold", size=14))

# Stampa
print(fig_preds)
print(fig_errs)
ggarrange(fig_preds, fig_errs, ncol=1, nrow=2, common.legend=FALSE, legend="right")

cat("Sezione 20 completata. Grafici generati con successo.\n")

# ---------------------------------------------------------------------
# 21. CONFRONTO FINALE LOOCV: ISOTROPO vs ANISOTROPO
# ---------------------------------------------------------------------

cat("\n--- INIZIO CROSS-VALIDATION (Leave-One-Out) ---\n")

# Funzione generica per il Kriging LOO con Errore di Misura
run_loo_cv <- function(data_df, model_signal, coords_cols, meas_err_col) {
  
  # Assicuriamoci che sia un data.frame classico per evitare problemi con tibble/sf
  data_df <- as.data.frame(data_df)
  
  n_stazioni <- nrow(data_df)
  preds <- numeric(n_stazioni)
  z_vals <- data_df$dS2S # Valori osservati (medie)
  
  # Matrice delle distanze completa
  coords_mat <- as.matrix(data_df[, coords_cols])
  dist_full <- as.matrix(dist(coords_mat))
  
  # Pre-calcolo Covarianza Segnale Completa
  cov_full_list <- variogramLine(model_signal, dist_vector = as.vector(dist_full), covariance = TRUE)
  Sigma_signal_full <- matrix(cov_full_list$gamma, nrow = n_stazioni, ncol = n_stazioni)
  
  # Loop su ogni stazione
  for(i in 1:n_stazioni) {
    
    idx_obs <- setdiff(1:n_stazioni, i)
    
    # 1. Costruzione Matrice LHS
    Sigma_obs <- Sigma_signal_full[idx_obs, idx_obs]
    
    # --- CORREZIONE QUI ---
    # Estraiamo i valori come VETTORE puro usando [[ ]]
    error_vals <- data_df[[meas_err_col]][idx_obs]
    
    # Creiamo la matrice diagonale (specificando nrow per sicurezza se fosse 1 solo numero)
    R_err <- diag(error_vals, nrow = length(error_vals))
    
    Sigma_tot <- Sigma_obs + R_err
    
    # 2. Costruzione Vettore RHS
    cov_target <- Sigma_signal_full[idx_obs, i]
    
    # 3. Sistema di Kriging
    N_obs <- length(idx_obs)
    A_mat <- rbind(cbind(Sigma_tot, 1), 1)
    A_mat[N_obs+1, N_obs+1] <- 0
    
    b_vec <- c(cov_target, 1)
    
    # 4. Soluzione
    weights_res <- try(solve(A_mat, b_vec), silent = TRUE)
    
    if(!inherits(weights_res, "try-error")) {
      weights <- weights_res[1:N_obs]
      preds[i] <- sum(weights * z_vals[idx_obs])
    } else {
      preds[i] <- NA
    }
    
    if(i %% 50 == 0) cat("iterazione", i, "\n")
  }
  cat("\n")
  return(preds)
}

# --- A. CV MODELLO ISOTROPO ---
cat("Esecuzione CV Isotropa (coordinate UTM)...")
cv_iso <- run_loo_cv(
  data_df = site_terms_all,
  model_signal = fit_site,          # Modello da Sez. 13 (Nugget=0)
  coords_cols = c("x_utm", "y_utm"),
  meas_err_col = "meas_error_var"   # Colonna con phi_0^2 / N
)

# --- B. CV MODELLO ANISOTROPO ---
cat("Esecuzione CV Anisotropa (coordinate Trasformate)...")
cv_ani <- run_loo_cv(
  data_df = site_terms_all,
  model_signal = fit_site_ani_equiv, # Modello da Sez. 17 (Equivalente)
  coords_cols = c("x_ani", "y_ani"), # Coordinate stirate/ruotate
  meas_err_col = "meas_error_var"
)

# --- C. CALCOLO METRICHE E CONFRONTO ---

# Creazione dataframe risultati
df_cv_results <- data.frame(
  NET_STA = site_terms_all$NET_STA,
  Osservato = site_terms_all$dS2S,
  Pred_Iso = cv_iso,
  Pred_Ani = cv_ani
)

# Residui (Osservato - Predetto)
df_cv_results$Res_Iso <- df_cv_results$Osservato - df_cv_results$Pred_Iso
df_cv_results$Res_Ani <- df_cv_results$Osservato - df_cv_results$Pred_Ani

# Calcolo RMSE (Root Mean Square Error)
rmse_iso <- sqrt(mean(df_cv_results$Res_Iso^2, na.rm=TRUE))
rmse_ani <- sqrt(mean(df_cv_results$Res_Ani^2, na.rm=TRUE))

cat("\n--- RISULTATI CROSS-VALIDATION ---\n")
cat(sprintf("RMSE Isotropo:   %.4f\n", rmse_iso))
cat(sprintf("RMSE Anisotropo: %.4f\n", rmse_ani))

pct_improvement <- (rmse_iso - rmse_ani) / rmse_iso * 100
cat(sprintf("Miglioramento con Anisotropia: %.2f%%\n", pct_improvement))

# --- D. PLOT DI CONFRONTO ---
# Scatterplot Osservato vs Predetto per i due modelli

lims <- range(df_cv_results$Osservato, na.rm=TRUE)

par(mfrow=c(1,2)) # Due grafici affiancati
plot(df_cv_results$Osservato, df_cv_results$Pred_Iso,
     main = paste("CV Isotropo\nRMSE:", round(rmse_iso,3)),
     xlab = "Osservato (dS2S)", ylab = "Predetto (LOO)",
     pch = 19, col = rgb(0,0,1,0.3), xlim=lims, ylim=lims, asp=1)
abline(0,1, col="red", lwd=2)

plot(df_cv_results$Osservato, df_cv_results$Pred_Ani,
     main = paste("CV Anisotropo\nRMSE:", round(rmse_ani,3)),
     xlab = "Osservato (dS2S)", ylab = "Predetto (LOO)",
     pch = 19, col = rgb(0,0.5,0,0.3), xlim=lims, ylim=lims, asp=1)
abline(0,1, col="red", lwd=2)
par(mfrow=c(1,1)) # Reset

# ---------------------------------------------------------------------
# 22. 10-FOLD CROSS-VALIDATION
# ---------------------------------------------------------------------

cat("\n--- INIZIO 10-FOLD CROSS-VALIDATION ---\n")

run_kfold_cv <- function(raw_data, reference_sites, model_signal, coords_cols, k = 10) {
  
  library(lme4)
  library(dplyr)
  
  raw_data <- as.data.frame(raw_data)
  reference_sites <- as.data.frame(reference_sites)
  
  stazioni_target <- reference_sites$NET_STA
  n_stazioni <- length(stazioni_target)
  preds_cv <- rep(NA, n_stazioni)
  
  set.seed(123) 
  folds <- sample(rep(1:k, length.out = n_stazioni))
  sta_fold_map <- data.frame(NET_STA = stazioni_target, fold = folds)
  
  for(f in 1:k) {
    sta_train <- sta_fold_map$NET_STA[sta_fold_map$fold != f]
    sta_test  <- sta_fold_map$NET_STA[sta_fold_map$fold == f]
    
    df_train <- raw_data[raw_data$NET_STA %in% sta_train, ]
    
    mod_cv <- try(lmer(res ~ 1 + (1 | EVENT), data = df_train, REML = TRUE), silent = TRUE)
    if(inherits(mod_cv, "try-error")) next
    
    offset_train <- fixef(mod_cv)[1]
    
    dBe_train <- ranef(mod_cv)$EVENT
    dBe_train$EVENT <- rownames(dBe_train)
    colnames(dBe_train)[1] <- "dBe"
    
    mod_sep_cv <- try(lmer(res ~ 1 + (1 | EVENT) + (1 | NET_STA), data = df_train, REML = TRUE), silent=TRUE)
    var_intra_train <- if(!inherits(mod_sep_cv, "try-error")) sigma(mod_sep_cv)^2 else sigma(mod_cv)^2
    
    site_train <- df_train %>%
      left_join(dBe_train, by = "EVENT") %>%
      mutate(dWs = res - offset_train - dBe) %>%
      group_by(NET_STA) %>%
      summarize(
        dS2S = mean(dWs, na.rm = TRUE),
        meas_error_var = var_intra_train / n()
      ) %>%
      ungroup() %>%
      left_join(reference_sites %>% select(NET_STA, all_of(coords_cols)), by = "NET_STA")
    
    site_test <- reference_sites %>% filter(NET_STA %in% sta_test)
    if(nrow(site_test) == 0) next
    
    coords_train <- as.matrix(site_train[, coords_cols])
    dist_train <- as.matrix(dist(coords_train))
    
    cov_train_list <- variogramLine(model_signal, dist_vector = as.vector(dist_train), covariance = TRUE)
    Sigma_train <- matrix(cov_train_list$gamma, nrow = nrow(coords_train))
    
    R_err <- diag(site_train$meas_error_var, nrow = nrow(site_train))
    Sigma_tot <- Sigma_train + R_err
    
    N_train <- nrow(site_train)
    A_mat <- rbind(cbind(Sigma_tot, 1), 1)
    A_mat[N_train+1, N_train+1] <- 0
    A_inv <- try(solve(A_mat), silent = TRUE)
    
    if(inherits(A_inv, "try-error")) next
    
    for(i_t in 1:nrow(site_test)) {
      target_sta <- site_test$NET_STA[i_t]
      
      dist_pt <- sqrt((site_train[[coords_cols[1]]] - site_test[[coords_cols[1]]][i_t])^2 + 
                        (site_train[[coords_cols[2]]] - site_test[[coords_cols[2]]][i_t])^2)
      
      cov_target <- variogramLine(model_signal, dist_vector = dist_pt, covariance = TRUE)$gamma
      b_vec <- c(cov_target, 1)
      weights <- (A_inv %*% b_vec)[1:N_train]
      
      idx_save <- which(stazioni_target == target_sta)
      preds_cv[idx_save] <- sum(weights * site_train$dS2S)
    }
  }
  cat("\nCalcolo K-Fold completato.\n")
  return(preds_cv)
}

# --- ESECUZIONE 10-FOLD CV ---

# 1. Modello ISOTROPO (Coordinate UTM)
cat("Esecuzione 10-Fold CV Isotropa...")
cv_kfold_iso <- run_kfold_cv(
  raw_data = itacentrale_df,            # <--- NUOVO: Passiamo i dati grezzi
  reference_sites = site_terms_all,     # <--- NUOVO: Passiamo le stazioni per l'allineamento
  model_signal = fit_site, 
  coords_cols = c("x_utm", "y_utm"),
  k = 10
)

# 2. Modello ANISOTROPO (Coordinate Trasformate)
cat("Esecuzione 10-Fold CV Anisotropa...")
cv_kfold_ani <- run_kfold_cv(
  raw_data = itacentrale_df,            # <--- NUOVO
  reference_sites = site_terms_all,     # <--- NUOVO
  model_signal = fit_site_ani_equiv, 
  coords_cols = c("x_ani", "y_ani"),
  k = 10
)

# --- CONFRONTO RISULTATI K-FOLD ---

df_kfold_results <- data.frame(
  Osservato = site_terms_all$dS2S,
  Pred_Iso  = cv_kfold_iso,
  Pred_Ani  = cv_kfold_ani
)

# RMSE
rmse_kfold_iso <- sqrt(mean((df_kfold_results$Osservato - df_kfold_results$Pred_Iso)^2, na.rm=TRUE))
rmse_kfold_ani <- sqrt(mean((df_kfold_results$Osservato - df_kfold_results$Pred_Ani)^2, na.rm=TRUE))

cat("\n--- RISULTATI 10-FOLD CV ---\n")
cat(sprintf("RMSE Isotropo:   %.4f\n", rmse_kfold_iso))
cat(sprintf("RMSE Anisotropo: %.4f\n", rmse_kfold_ani))
cat(sprintf("Miglioramento:   %.2f%%\n", (rmse_kfold_iso - rmse_kfold_ani)/rmse_kfold_iso*100))

# Plot rapido
par(mfrow=c(1,2))
lims_k <- range(df_kfold_results$Osservato, na.rm=TRUE)

plot(df_kfold_results$Osservato, df_kfold_results$Pred_Iso, main="10-Fold Isotropo",
     xlab="Osservato", ylab="Predetto", pch=19, col=rgb(0,0,1,0.3), xlim=lims_k, ylim=lims_k)
abline(0,1, col="red")

plot(df_kfold_results$Osservato, df_kfold_results$Pred_Ani, main="10-Fold Anisotropo",
     xlab="Osservato", ylab="Predetto", pch=19, col=rgb(0,0.5,0,0.3), xlim=lims_k, ylim=lims_k)
abline(0,1, col="red")
par(mfrow=c(1,1))

# ---------------------------------------------------------------------
# 23. TEST DI STAZIONARIETÀ (MOVING WINDOW STATISTICS)
# ---------------------------------------------------------------------
# Obiettivo: Verificare se Media e Varianza dipendono dalla posizione.
# Se vediamo trend chiari, l'ipotesi di stazionarietà (necessaria per Kriging Ordinario) cade.

library(sf)
library(ggplot2)
library(dplyr)

print("Inizio analisi di stazionarietà locale...")

# 1. DEFINIZIONE GRIGLIA DI ANALISI
#    Creiamo dei "quadratoni" (tiles) di dimensione fissa (es. 40x40 km)
#    Per avere statistica, ogni quadrato deve contenere almeno N stazioni.

GRID_SIZE_METRI <- 40000 # 40 km
MIN_OBS_PER_TILE <- 5    # Minimo stazioni per calcolare varianza

# Estrai i dati unici per stazione (usiamo le medie di sito o i residui grezzi)
# Usiamo i dWs grezzi per avere più punti.
data_for_test <- itacentrale_df %>%
  select(x_utm, y_utm, dWs)

# Crea oggetto SF
data_sf <- st_as_sf(data_for_test, coords = c("x_utm", "y_utm"), crs = 32633)

# Crea la griglia (Tessellazione)
grid_test <- st_make_grid(data_sf, cellsize = GRID_SIZE_METRI) %>% st_as_sf()
grid_test$ID <- 1:nrow(grid_test)

# 2. CALCOLO STATISTICHE PER CELLA (SPATIAL JOIN)
#    Associa ogni residuo al suo quadrato della griglia

joined <- st_join(data_sf, grid_test)

stats_grid <- joined %>%
  st_drop_geometry() %>%
  group_by(ID) %>%
  summarize(
    Media_Locale = mean(dWs, na.rm = TRUE),
    Varianza_Locale = var(dWs, na.rm = TRUE),
    N_Obs = n()
  ) %>%
  filter(N_Obs >= MIN_OBS_PER_TILE) # Rimuovi celle vuote o povere

# Riunisci i dati alla geometria della griglia
grid_results <- left_join(grid_test, stats_grid, by = "ID") %>%
  filter(!is.na(Media_Locale))

# 3. PLOT DIAGNOSTICI

# A. Mappa della Media Locale (Test Stazionarietà 1° Ordine)
#    Se è stazionario, deve essere tutto vicino a 0 (colore neutro).
#    Se vedi Rosso a Nord e Blu a Sud, c'è un Trend non rimosso.

limit_mean <- max(abs(grid_results$Media_Locale), na.rm=TRUE)

p1 <- ggplot() +
  geom_sf(data = grid_results, aes(fill = Media_Locale), color = "black") +
  scale_fill_distiller(palette = "RdBu", direction = -1, limits=c(-limit_mean, limit_mean)) +
  labs(title = "Stazionarietà della Media", 
       subtitle = "Cerca pattern spaziali evidenti (Trend)") +
  theme_minimal()

# B. Mappa della Varianza Locale (Test Stazionarietà 2° Ordine)
#    Se è stazionario, il colore deve essere uniforme.
#    Se vedi zone molto scure e zone molto chiare, c'è Eteroschedasticità.

p2 <- ggplot() +
  geom_sf(data = grid_results, aes(fill = Varianza_Locale), color = "black") +
  scale_fill_viridis_c(option = "magma", direction = -1) + # Magma: Giallo=Alto, Nero=Basso
  labs(title = "Stazionarietà della Varianza", 
       subtitle = "Cerca differenze di 'energia' del segnale") +
  theme_minimal()

# Visualizza

#library(gridExtra)
#grid.arrange(p1, p2, ncol = 2)
plot(p1); plot(p2)

# 4. TEST DI CORRELAZIONE LINEARE (Trend Test)
#    C'è una correlazione tra coordinata e media?

coords <- st_coordinates(st_centroid(grid_results))
df_stat <- data.frame(
  X = coords[,1],
  Y = coords[,2],
  Mean = grid_results$Media_Locale,
  Var = grid_results$Varianza_Locale
)

cor_x <- cor.test(df_stat$X, df_stat$Mean)
cor_y <- cor.test(df_stat$Y, df_stat$Mean)

cat("\n--- Risultati Test Correlazione Trend ---\n")
cat("Correlazione Longitudine (X) vs Media Locale: p-value =", round(cor_x$p.value, 4), "\n")
cat("Correlazione Latitudine (Y) vs Media Locale: p-value =", round(cor_y$p.value, 4), "\n")

if(cor_x$p.value < 0.05 || cor_y$p.value < 0.05) {
  cat("ATTENZIONE: Possibile Trend Spaziale rilevato (Non Stazionario nella Media).\n")
} else {
  cat("OK: Nessun trend lineare significativo rilevato.\n")
}

# ------------------------------------------------------------------------------
# 24. OMOSCHEDASTICITÀ
# ------------------------------------------------------------------------------

# Estraiamo i residui puri (epsilon_0) dal modello mod_sep
# Nota: resid(mod_sep) estrae i residui condizionati (cioè l'errore epsilon)
itacentrale$residuals_pure <- resid(mod_sep)

# Impostiamo una finestra grafica 2x2
par(mfrow = c(2,2))

# 1. Residui vs Magnitudo (Mw)
plot(itacentrale$Mw_RCMT, itacentrale$residuals_pure,
     xlab = "Magnitudo (Mw)", ylab = "Residui (epsilon_0)",
     main = "Residui vs Mw", pch = 20, col = rgb(0,0,0,0.3))
abline(h = 0, col = "red", lwd = 2)
# Se vedi un "imbuto" (più stretti a Mw alte), c'è eteroschedasticità.

# 2. Residui vs Distanza (RR)
plot(itacentrale$RR, itacentrale$residuals_pure,
     xlab = "Distanza (km)", ylab = "Residui (epsilon_0)",
     main = "Residui vs Distanza", pch = 20, col = rgb(0,0,0,0.3))
abline(h = 0, col = "red", lwd = 2)

# 3. Residui vs Vs30 (Proxy per il sito)
plot(itacentrale$vs30, itacentrale$residuals_pure,
     xlab = "Vs30 (m/s)", ylab = "Residui (epsilon_0)",
     main = "Residui vs Vs30", pch = 20, col = rgb(0,0,0,0.3), log="x")
abline(h = 0, col = "red", lwd = 2)
# Questo ti dice se i suoli soffici sono più "rumorosi" delle rocce.

# 4. Residui vs Valori Predetti
plot(fitted(mod_sep), itacentrale$residuals_pure,
     xlab = "Valori Predetti (Fitted)", ylab = "Residui",
     main = "Tukey-Anscombe Plot", pch = 20, col = rgb(0,0,0,0.3))
abline(h = 0, col = "red", lwd = 2)

par(mfrow = c(1,1)) # Ripristina finestra

library(dplyr)

# Calcola la deviazione standard dei residui per ogni stazione
station_noise <- itacentrale %>%
  group_by(NET_STA) %>%
  summarise(
    sd_noise = sd(residuals_pure),
    n_records = n()
  ) %>%
  filter(n_records >= 5) # Consideriamo solo stazioni con almeno 5 dati per robustezza

# Plot della variabilità per stazione
plot(station_noise$n_records, station_noise$sd_noise^2,
     xlab = "Numero di Eventi per Stazione",
     ylab = "Varianza stimata del rumore (S^2)",
     main = "Verifica Omoschedasticità tra Stazioni",
     pch = 19)
abline(h = phi_0_sq, col = "red", lwd = 2, lty = 2) # La tua assunzione globale



################################################################################
######################### COSE INUTILI O DA RIVEDERE ###########################
################################################################################


#---------------------------------------------------------------------
# 7. KRIGING (PREVISIONE SPAZIALE) SU UN SINGOLO EVENTO
#---------------------------------------------------------------------
# Ora usiamo il modello 'v.fit' (calcolato su TUTTI gli eventi)
# per fare predizioni (Kriging) su UN evento specifico.

# Scegliamo un evento per cui fare la predizione
EVENTO_TARGET <- "EMSC-20090406_0000100"
cat(paste("\nEsecuzione Kriging per l'evento:", EVENTO_TARGET), "\n")

# 1. Estrai i dati di osservazione SOLO per quell'evento
dati_osservati_evento <- filter(itacentrale_utm, EVENT == EVENTO_TARGET)

# 2. Crea l'oggetto gstat per questo evento
#    NOTA: Usiamo il modello 'v.fit' fittato globalmente!
g.evento <- gstat(
  formula = dWs ~ 1,              # ~ 1 = Ordinary Kriging (media sconosciuta)
  data = dati_osservati_evento,
  model = fit_exp)                # dai plot il fit esponenziale sembra il migliore

# 3. Crea una griglia di predizione (CON PADDING)

# 3a. Definisci il padding (es. 20km)
padding_metri <- 20000 

# 3b. Calcola il Bounding Box originale
bbox_originale <- st_bbox(dati_osservati_evento)

# 3c. Crea un Bounding Box allargato (MODO CORRETTO)

# 3c-1. Calcola il vettore numerico (invariato)
bbox_vector <- c(
  xmin = bbox_originale["xmin"] - 2*padding_metri,
  ymin = bbox_originale["ymin"] - 0.5*padding_metri,
  xmax = bbox_originale["xmax"] + 2*padding_metri,
  ymax = bbox_originale["ymax"] + 0.5*padding_metri
)

# 3c-2. Trasforma il vettore in un oggetto 'bbox' (Metodo Classico)
# Dato che st_as_bbox() non è disponibile, impostiamo la classe manualmente.
class(bbox_vector) <- "bbox"

# 3c-3. Assegna il CRS all'oggetto appena creato
st_crs(bbox_vector) <- st_crs(dati_osservati_evento) 

# 3c-4. Riassegna al nome della variabile che usiamo dopo
bbox_allargato <- bbox_vector

# 3d. Crea la griglia usando il BBOX ALLARGATO
griglia_evento <- st_make_grid(
  st_as_sfc(bbox_allargato),  # <--- Usa il bbox allargato
  n = c(30, 30))
griglia_evento_sf <- st_as_sf(griglia_evento)

# 4. Esegui la predizione (Ordinary Kriging)
ok_predizione <- predict(g.evento, newdata = griglia_evento_sf)

# 5. Plotta i risultati
global_limits <- range( c(ok_predizione$var1.pred, dati_osservati_evento$dWs), na.rm = TRUE)

print(
  ggplot() +
    geom_sf(data = ok_predizione, aes(fill = var1.pred), color = NA) +
    scale_fill_viridis_c(name = "dWs", limits = global_limits) + 
    geom_sf(data = dati_osservati_evento, aes(color = dWs), size = 4) +
    scale_color_viridis_c(name = "dWs", limits = global_limits) +
    labs(title = paste("Mappa di Kriging (dWs) per", EVENTO_TARGET)) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank() ) )

print(
  ggplot() +
    geom_sf(data = ok_predizione, aes(fill = var1.var), color = NA) +
    scale_fill_viridis_c(name = "var(dWs)") + 
    geom_sf(data = dati_osservati_evento, color = "darkorange", size = 3) +
    scale_color_viridis_c(name = "dWs") +
    labs(title = paste("Mappa di Kriging (varianza di dWs) per", EVENTO_TARGET)) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank() ) )

#---------------------------------------------------------------------
# 8. ESEMPI AGGIUNTIVI
#---------------------------------------------------------------------

# A. Predizione in un punto specifico (s0.new)
s0.new <- st_as_sf(data.frame(x = 343995.5, y = 4739530), 
                   coords = c("x", "y"), 
                   crs = st_crs(itacentrale_utm))

# Predizione del valore (var1.pred) e della sua varianza (var1.var)
pred_s0 <- predict(g.evento, s0.new)
print("Predizione in s0.new:")
print(pred_s0)

# B. Stima della media (BLUE = Best Linear Unbiased Estimator)
#    Nel Kriging Ordinario (dWs ~ 1), questo stima la media 'beta' (l'intercetta)
#    Dato che dWs sono residui, ci aspettiamo che sia vicino a 0.
blue_s0 <- predict(g.evento, s0.new, BLUE = TRUE)
print("Stima del trend (media) in s0.new:")
print(blue_s0)

# C. Predizione in un punto dove ABBIAMO dati
#    La varianza di predizione (var1.var) deve essere (quasi) zero.
stazione_1 <- dati_osservati_evento[1, ]
pred_staz1 <- predict(g.evento, stazione_1)
print("Predizione in una location osservata (varianza ~ 0):")
print(pred_staz1)

# ---------------------------------------------------------------------
# 12. CONFRONTO KRIGING: EMPIRICO (Data-Driven) vs TEORICO (Model-Based)
#                           (DA RIVEDERE)
# ---------------------------------------------------------------------
# NOTA: Non possiamo fare una mappa continua perché non conosciamo la
# covarianza empirica verso punti dove non ci sono stazioni.
# Facciamo quindi una "Leave-One-Out Cross Validation" sull'evento.

print("Inizio Kriging Manuale Empirico...")

# 1. PREPARAZIONE DATI (Come prima)
dati_evento <- itacentrale_df %>% filter(EVENT == EVENTO_TARGET)
stazioni_evento <- dati_evento$NET_STA
stazioni_valide <- intersect(stazioni_evento, colnames(emp_cov_matrix))
dati_validi <- dati_evento %>% filter(NET_STA %in% stazioni_valide)

# --- A. PREVISIONE EMPIRICA (MANUALE) ---
# (Questo è il loop che abbiamo scritto prima)
preds_empiriche <- numeric(length(stazioni_valide))

Sigma_evento <- emp_cov_matrix[stazioni_valide, stazioni_valide]
Sigma_evento[is.na(Sigma_evento)] <- 0 

for(k in 1:length(stazioni_valide)) {
  target_sta <- stazioni_valide[k]
  obs_sta    <- stazioni_valide[-k]
  
  N <- length(obs_sta)
  Sigma_obs <- Sigma_evento[obs_sta, obs_sta]
  A_mat <- rbind(cbind(Sigma_obs, 1), 1)
  A_mat[N+1, N+1] <- 0
  cov_target <- Sigma_evento[obs_sta, target_sta]
  b_vec <- c(cov_target, 1)
  
  soluzione <- try(solve(A_mat, b_vec), silent = TRUE)
  
  if(!inherits(soluzione, "try-error")) {
    weights <- soluzione[1:N]
    valori_osservati <- dati_validi$dWs[dati_validi$NET_STA %in% obs_sta]
    preds_empiriche[k] <- sum(weights * valori_osservati)
  } else {
    preds_empiriche[k] <- NA 
  }
}

# --- B. PREVISIONE TEORICA (AUTOMATICA CON GSTAT) ---
# Usiamo la funzione krige.cv che fa esattamente la stessa cosa del loop sopra,
# ma usando il modello variografico (fit_exp) e le distanze.

# 1. Converti i dati validi in sf per gstat
dati_validi_sf <- st_as_sf(dati_validi, coords = c("x_utm", "y_utm"), crs = 32633)

# 2. Esegui Cross-Validation
#    nfold = nrow(...) significa "Leave-One-Out" (togline uno alla volta)
cv_teorico <- krige.cv(
  dWs ~ 1,
  locations = dati_validi_sf,
  model = fit_exp,  # <--- QUI USA IL MODELLO FITTATO (Teorico)
  nfold = nrow(dati_validi_sf)
)

preds_teoriche <- cv_teorico$var1.pred

# --- C. CONSOLIDAMENTO DATI ---

df_confronto <- data.frame(
  NET_STA = stazioni_valide,
  Osservato = dati_validi$dWs,
  Pred_Empirica = preds_empiriche,
  Pred_Teorica = preds_teoriche
)

# Calcolo RMSE (Root Mean Square Error) per decretare il vincitore
rmse_emp <- sqrt(mean((df_confronto$Osservato - df_confronto$Pred_Empirica)^2, na.rm=TRUE))
rmse_teo <- sqrt(mean((df_confronto$Osservato - df_confronto$Pred_Teorica)^2, na.rm=TRUE))

print("--- RISULTATI SFIDA ---")
cat("RMSE Metodo Empirico (Covarianza Storica):", round(rmse_emp, 4), "\n")
cat("RMSE Metodo Teorico (Variogramma Distanza):", round(rmse_teo, 4), "\n")

# --- D. PLOT DI CONFRONTO FINALE ---

# Calcoliamo limiti per il grafico quadrato
lims <- range(c(df_confronto$Osservato, df_confronto$Pred_Empirica, df_confronto$Pred_Teorica), na.rm=TRUE)

plot(NA, xlim = lims, ylim = lims,
     main = paste("Sfida Kriging su", EVENTO_TARGET),
     xlab = "Valore Osservato (dWs)",
     ylab = "Valore Predetto",
     asp = 1) # Mantiene il rapporto 1:1 visivo
grid()

# Linea di riferimento perfetta
abline(0, 1, lwd = 2, col = "darkgrey", lty = 2)

# Punti Empirici (BLU)
points(df_confronto$Osservato, df_confronto$Pred_Empirica, 
       pch = 19, col = "blue", cex = 1.2)

# Punti Teorici (ROSSI)
points(df_confronto$Osservato, df_confronto$Pred_Teorica, 
       pch = 17, col = "red", cex = 1.2) # Triangoli rossi

# Aggiungi segmenti per evidenziare la differenza per ogni stazione? (Opzionale)
segments(df_confronto$Osservato, df_confronto$Pred_Empirica, 
         df_confronto$Osservato, df_confronto$Pred_Teorica, col = "grey80")

legend("topleft", 
       legend = c(paste("Empirico (RMSE:", round(rmse_emp, 3), ")"), 
                  paste("Teorico (RMSE:", round(rmse_teo, 3), ")"),
                  "Target Perfetto"),
       col = c("blue", "red", "darkgrey"), 
       pch = c(19, 17, NA), 
       lty = c(NA, NA, 2),
       bg = "white",
       cex = 0.8)

#---------------------------------------------------------------------
# 22. PLOT DELLA MEDIANA SENZA BINNING
#            (PESANTISSIMO)
#---------------------------------------------------------------------

library(dplyr)

# 1. Calcolo dei dati (Mediana per distanza)
df_mediana <- df_cut %>%
  mutate(h_round = round(h_m, 0)) %>% # Arrotonda al metro
  group_by(h_round) %>%
  summarize(
    dist_media = mean(h_m, na.rm = TRUE),
    gamma_mediana = median(gamma, na.rm = TRUE)
  ) %>%
  ungroup()

# 2. Plot (Base R)
plot(
  df_mediana$dist_media, 
  df_mediana$gamma_mediana,
  main = "Variogramma Station-Pair (Mediana)",
  xlab = "Distanza h (metri)",
  ylab = "Mediana di Gamma",
  pch = 20, cex = 0.4, col = rgb(0, 0, 0, 0.3), las = 1)
grid()

#---------------------------------------------------------------------
# 4.B ALTERNATIVA: VARIOGRAMMA BASATO SU COPPIE DI STAZIONI (NO BINNING)
#                         (INSOSTENIBILE)
#---------------------------------------------------------------------
print("Inizio calcolo variogramma per coppie di stazioni (Metodo A)...")

# 1. Raggruppamento per distanza esatta (Collasso delle "colonne")
#    Arrotondiamo la distanza a 1 metro per gestire minime differenze numeriche
#    e raggruppare la stessa coppia fisica.

soglia = 0.2

vario_pairs_df <- df_cut %>%
  mutate(h_round = round(h_m, 0)) %>% # Arrotonda al metro
  group_by(h_round) %>%
  summarize(
    dist = mean(h_m, na.rm = TRUE),     # Distanza esatta media
    gamma = median(gamma, na.rm = TRUE), # Mediana robusta della "colonna"
    np = as.numeric(n())                       # Numero di eventi per questa coppia
  ) %>%
  ungroup() %>%
  filter(np > 0) %>% # Rimuovi eventuali errori
  filter(gamma <= soglia)

cat("Numero di coppie uniche di stazioni trovate:", nrow(vario_pairs_df), "\n")

# 2. Formattazione per gstat
#    Creiamo un oggetto 'gstatVariogram' con migliaia di punti invece di 15 bin.

vario_pairs <- as.data.frame(vario_pairs_df)
vario_pairs$h_round <- NULL # Rimuovi colonna ausiliaria
class(vario_pairs) <- c("gstatVariogram", "data.frame")
attr(vario_pairs, "model") <- vgm(model = "Nug") # Attributo necessario per gstat

plot(
  vario_pairs$dist,
  vario_pairs$gamma,
  main = "Variogramma per Coppie di Stazioni (Station-Pair Averaging)",
  xlab = "Distanza h (metri)",
  ylab = "Semivarianza (gamma mediana)",
  pch = 20, cex = 0.4, col = rgb(0, 0, 0, 0.3), ylim = c(0, max(vario_pairs$gamma) * 1.1))
grid()

# Non si capisce niente, il resto della sezione è inutile

# 3. Fit dei Modelli su 'vario_pairs'
#    Nota: Il fitting sarà più lento perché deve ottimizzare su molti più punti.

print("Fitting dei modelli sulla nuvola delle coppie...")

# Usiamo gli stessi parametri iniziali definiti prima
nugget_start <- 0.02
sill_start <- 0.08
range_start_metri <- 200000 

# Definizioni modelli
vgm_sph <- vgm(psill = sill_start-nugget_start, model = "Sph", range = range_start_metri, nugget = nugget_start)
vgm_exp <- vgm(psill = sill_start-nugget_start, model = "Exp", range = range_start_metri/3, nugget = nugget_start)
vgm_mat <- vgm(psill = sill_start-nugget_start, model = "Mat", range = range_start_metri/4, nugget = nugget_start, kappa = 1)

# Esecuzione Fit (potrebbe impiegare qualche secondo in più)
fit_sph_pairs <- try(fit.variogram(vario_pairs, model = vgm_sph, fit.method = 7), silent = TRUE)
fit_exp_pairs <- try(fit.variogram(vario_pairs, model = vgm_exp, fit.method = 7), silent = TRUE)
fit_mat_pairs <- try(fit.variogram(vario_pairs, model = vgm_mat, fit.method = 7), silent = TRUE)

print("--- Modello Esponenziale (Pairs) ---")
print(fit_exp_pairs)

# 4. Plot Confronto
#    Attenzione: ci sono molti punti, usiamo un simbolo piccolo (pch='.')

maxdist <- max(vario_pairs$dist, na.rm = TRUE)

# Aggiungi le linee fittate
if(!inherits(fit_sph_pairs, "try-error")) {
  lines(variogramLine(fit_sph_pairs, maxdist = maxdist, n = 200), col = "red", lwd = 2)}

if(!inherits(fit_exp_pairs, "try-error")) {
  lines(variogramLine(fit_exp_pairs, maxdist = maxdist, n = 200), col = "blue", lwd = 2)}

if(!inherits(fit_mat_pairs, "try-error")) {
  lines(variogramLine(fit_mat_pairs, maxdist = maxdist, n = 200), col = "forestgreen", lwd = 2)}

legend(
  "bottomright",
  legend = c("Coppie Stazioni (Mediana)", "Sferico", "Esponenziale", "Matérn"),
  pch = c(20, NA, NA, NA),
  lty = c(NA, 1, 1, 1),
  lwd = c(NA, 2, 2, 2),
  col = c("black", "red", "blue", "forestgreen"),
  bg = "white"
)

#---------------------------------------------------------------------
# 4.C ALTERNATIVA: MAXIMUM LIKELIHOOD ESTIMATION (MLE)
#                         (INSOSTENIBILE)
#---------------------------------------------------------------------

library(geoR)
print("Inizio stima MLE (Metodo B)...")

# 1. Preparazione e Campionamento Dati
#    MLE inverte matrici enormi. Non usare con > 1500 punti.
#    Prendiamo un sottocampione casuale dal tuo dataset proiettato.

N_SAMPLES <- 1000 # Numero massimo di punti da usare
set.seed(123)     # Per riproducibilità

if(nrow(itacentrale_df) > N_SAMPLES) {
  cat("Dataset troppo grande per MLE puro. Eseguo subsampling a", N_SAMPLES, "punti.\n")
  data_mle_subset <- itacentrale_df[sample(1:nrow(itacentrale_df), N_SAMPLES), ]
} else {
  data_mle_subset <- itacentrale_df
}

# Crea oggetto 'geodata' (specifico di geoR)
geodata_obj <- as.geodata(data_mle_subset, coords.col = c("x_utm", "y_utm"), data.col = "dWs")

# 2. Stima dei Parametri (likfit)
#    Bisogna fornire dei valori iniziali (ini.cov.pars) per aiutare l'algoritmo.
#    ini.cov.pars = c(partial_sill, range)

print("Esecuzione likfit (potrebbe impiegare minuti)...")

# Stima modello Esponenziale
mle_fit_exp <- likfit(
  geodata = geodata_obj, 
  ini.cov.pars = c(0.06, 50000), # Sill parziale 0.06, Range 50km (tentativo iniziale)
  nugget = 0.02,                 # Nugget iniziale
  cov.model = "exponential",     # Modello
  fix.nugget = FALSE,            # Lascia che il nugget venga stimato
  messages = TRUE
)

print("--- Risultati MLE (Esponenziale) ---")
summary(mle_fit_exp)

# 3. Plot di confronto
#    MLE non produce punti, ma una curva. La sovrapponiamo al variogramma binned
#    calcolato nella Sezione 4 classica per vedere se ha senso.

# Ricalcoliamo velocemente il binned per averlo come sfondo
# (Assumiamo che 'exp_vario' esista dalla Sezione 4 precedente)

plot(exp_vario$dist, exp_vario$gamma, 
     main = "Confronto: Variogramma Binned vs MLE Curve",
     xlab = "Distanza (m)", ylab = "Semivarianza",
     pch = 19, ylim = c(0, max(exp_vario$gamma)*1.2))

# Estrai parametri dal fit MLE per disegnare la linea
mle_nugget <- mle_fit_exp$nugget
mle_sill   <- mle_fit_exp$cov.pars[1]
mle_range  <- mle_fit_exp$cov.pars[2]

# Funzione per disegnare la curva esponenziale MLE
curve(mle_nugget + mle_sill * (1 - exp(-x/mle_range)), 
      add = TRUE, col = "purple", lwd = 3, lty = 2)

legend("bottomright", legend = c("Dati Binned (Riferimento)", "Fit MLE (Likfit)"),
       col = c("black", "purple"), pch = c(19, NA), lty = c(NA, 2), lwd = c(NA, 3))


#---------------------------------------------------------------------
# 4.D ALTERNATIVA: CLOUD FITTING (Minimi Quadrati su tutti i punti)
#                         (INSOSTENIBILE)
#---------------------------------------------------------------------
print("Inizio Cloud Fitting (Metodo C)...")

# 1. Preparazione Dati
#    Recuperiamo 'df_for_vario' dalla Sezione 3.
#    Dobbiamo trasformarlo in un oggetto che gstat accetta come variogramma.

# Filtriamo solo per distanza ragionevole per facilitare il fit
MAX_DIST_CLOUD <- 230000 
cloud_data <- df_for_vario %>% 
  filter(h_m > 0 & h_m <= MAX_DIST_CLOUD)

# Creiamo l'oggetto gstatVariogram "finto"
# gstat si aspetta colonne: np, dist, gamma.
# Qui ogni riga è una coppia, quindi np = 1.

cloud_gstat <- data.frame(
  np = 1, 
  dist = cloud_data$h_m, 
  gamma = cloud_data$gamma,
  dir.hor = 0, dir.ver = 0, id = "var1"
)
class(cloud_gstat) <- c("gstatVariogram", "data.frame")
attr(cloud_gstat, "model") <- vgm(model = "Nug") 

cat("Numero di punti nella nuvola:", nrow(cloud_gstat), "\n")

# 2. Fitting
#    Usiamo fit.method = 7 (Minimi quadrati pesati standard).
#    Nota: Senza binning, i pesi basati su N/h^2 sono meno efficaci, 
#    il fit sarà guidato dalla massa dei punti.

print("Fitting modelli sulla nuvola (potrebbe essere lento)...")

# Modelli iniziali (stessi di prima)
vgm_start <- vgm(psill = 0.06, model = "Exp", range = 60000, nugget = 0.02)

fit_exp_cloud <- try(
  fit.variogram(cloud_gstat, model = vgm_start, fit.method = 7), 
  silent = TRUE
)

print("--- Risultati Cloud Fitting ---")
print(fit_exp_cloud)

# 3. Plot della Nuvola e del Fit
#    Usiamo un plot denso per visualizzare la nuvola.

# Per il plot, campioniamo la nuvola se è troppo grande (solo per visualizzazione)
if(nrow(cloud_data) > 10000) {
  plot_data <- cloud_data[sample(1:nrow(cloud_data), 10000), ]
} else {
  plot_data <- cloud_data
}

plot(plot_data$h_m, plot_data$gamma,
     main = "Variogram Cloud Fitting",
     xlab = "Distanza (m)", ylab = "Semivarianza",
     pch = ".", col = rgb(0, 0, 0, 0.2), # Puntini trasparenti
     ylim = c(0, quantile(plot_data$gamma, 0.99))) # Taglia gli outlier estremi Y

grid()

# Aggiungi la linea fittata
if(!inherits(fit_exp_cloud, "try-error")) {
  lines(variogramLine(fit_exp_cloud, maxdist = MAX_DIST_CLOUD, n = 200), 
        col = "orange", lwd = 3)
}

legend("topleft", legend = c("Variogram Cloud (campionato)", "Fit su TUTTA la Cloud"),
       col = c("grey", "orange"), pch = c(46, NA), lty = c(NA, 1), lwd = c(NA, 3))

################################################################################
########################## SALVATAGGIO PER PONTE ###############################
################################################################################

cat("\n>>> Preparazione esportazione per il Ponte...\n")

# Verifica che la 10-Fold CV sia stata eseguita
if (!exists("rmse_kfold_ani") | !exists("rmse_kfold_iso")) {
  stop("ERRORE: Devi prima eseguire la Sezione 22 (10-Fold CV) per calcolare rmse_kfold_iso e rmse_kfold_ani.")
}

# --- 1. AGGIUNTA CELL_ID ALLA GRIGLIA ---
if(!"cell_id" %in% names(grid_plain)) {
  grid_plain$cell_id <- grid_plain$serpentine_idx
}

# --- 2. RECUPERO INFORMAZIONI GEOGRAFICHE STAZIONI ---
# Ci servono le Lat/Lon originali da itacentrale
info_geo_staz <- itacentrale %>%
  distinct(NET_STA, .keep_all = TRUE) %>%
  select(NET_STA, st_lon, st_lat)

# Creiamo stazioni_obs unendo le info e trasformando in km
stazioni_per_ponte <- site_terms_all %>%
  mutate(X_km = x_utm / 1000, Y_km = y_utm / 1000) %>%
  left_join(info_geo_staz, by = "NET_STA")

# Funzione per trovare il cell_id più vicino alla stazione
trova_cella_vicina <- function(x_km, y_km, grid_df) {
  distanze_sq <- (grid_df$x_utm/1000 - x_km)^2 + (grid_df$y_utm/1000 - y_km)^2
  return(grid_df$cell_id[which.min(distanze_sq)])
}

stazioni_per_ponte$cell_id <- mapply(trova_cella_vicina,
                                     stazioni_per_ponte$X_km,
                                     stazioni_per_ponte$Y_km,
                                     MoreArgs = list(grid_df = grid_plain))

# --- 3. SELEZIONE MODELLO VINCITORE E COSTRUZIONE MATRICE GLOBALE ---
if (rmse_kfold_ani < rmse_kfold_iso) {
  cat("=> Il modello ANISOTROPO ha vinto (RMSE 10-Fold:", rmse_kfold_ani, "). Esporto questo.\n")
  
  # Mappa
  mappa_out <- data.frame(
    cell_id = grid_plain$cell_id,
    X_km = df_map_ani$x / 1000,
    Y_km = df_map_ani$y / 1000,
    Pred_Staz = df_map_ani$pred,
    SE_Staz = df_map_error_ani$std_dev
  )
  
  # Matrice Covarianza Stazioni
  cov_out <- Sigma_signal_ani
  
  # CALCOLO MATRICE GLOBALE (ANISOTROPA)
  cat("Calcolo Matrice Covarianza Globale Anisotropa...\n")
  # Utilizziamo la griglia trasformata (stirata e ruotata) per calcolare le distanze
  dist_grid_ani <- as.matrix(dist(grid_coords_ani))
  cov_globale_list <- variogramLine(fit_site_ani_equiv, dist_vector = as.vector(dist_grid_ani), covariance = TRUE)
  C_globale_out <- matrix(cov_globale_list$gamma, nrow = nrow(grid_coords_ani))
  
  rmse_out <- rmse_kfold_ani
  nome_vincitore <- "Stazionario (Anisotropo)"
  
} else {
  cat("=> Il modello ISOTROPO ha vinto (RMSE 10-Fold:", rmse_kfold_iso, "). Esporto questo.\n")
  
  # Mappa
  mappa_out <- data.frame(
    cell_id = grid_plain$cell_id,
    X_km = df_map_weighted$x / 1000,
    Y_km = df_map_weighted$y / 1000,
    Pred_Staz = df_map_weighted$pred,
    SE_Staz = df_map_error$std_dev
  )
  
  # Matrice Covarianza Stazioni
  cov_out <- Sigma_signal_pure
  
  # CALCOLO MATRICE GLOBALE (ISOTROPA)
  cat("Calcolo Matrice Covarianza Globale Isotropa...\n")
  # Utilizziamo la griglia standard in metri
  dist_grid_iso <- as.matrix(dist(grid_coords))
  cov_globale_list <- variogramLine(fit_site, dist_vector = as.vector(dist_grid_iso), covariance = TRUE)
  C_globale_out <- matrix(cov_globale_list$gamma, nrow = nrow(grid_coords))
  
  rmse_out <- rmse_kfold_iso
  nome_vincitore <- "Stazionario (Isotropo)"
}

# Impostiamo nomi per la matrice stazioni
rownames(cov_out) <- site_terms_all$NET_STA
colnames(cov_out) <- site_terms_all$NET_STA
stazioni_out <- site_terms_all$NET_STA

# --- 4. ESPORTAZIONE ---

# Sostituisci questo percorso con quello PGV se stai modificando lo script PGV
cartella_export <- "C:/Users/miche/OneDrive/Desktop/prova3/INGV - pgv" 

if(!dir.exists(cartella_export)) {
  dir.create(cartella_export, recursive = TRUE)
}

# Impacchettamento nel formato "Universale" per il Ponte
risultati_stazionario <- list(
  modello_scelto = nome_vincitore,
  
  # Variabili originali
  mappa = mappa_out,
  cov = cov_out,
  stazioni = stazioni_out,
  rmse = rmse_out,
  
  # Aggiunte per abilitare la Figura 2 nel Ponte!
  cov_globale = C_globale_out,
  stazioni_obs = stazioni_per_ponte[, c("NET_STA", "X_km", "Y_km", "cell_id", "st_lon", "st_lat")]
)

path_export <- file.path(cartella_export, "Risultati_Stazionario.rds")
saveRDS(risultati_stazionario, file = path_export)
print(paste("Risultati Stazionario esportati con successo in:", path_export))
