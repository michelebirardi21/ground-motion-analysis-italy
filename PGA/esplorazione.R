

setwd("C:/Users/miche/OneDrive/Desktop/prova3/INGV")

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
names(data)
head(data)

lat_min <- 41.5
lat_max <- 44
lon_min <- 10.5
lon_max <- 15

#plot(data$ev_lon,data$ev_lat, pch = 16, col = 'blue')
#abline(h = lat_min, v = lon_min, col = 'darkgreen')
#abline(h = lat_max, v = lon_max, col = 'darkgreen')
#abline(a=50,b=-0.7, col = 'darkred')
#abline(a=54,b=-0.7, col = 'darkred')

library(ggplot2)
library(maps)
library(ggforce)

italy_map <- map_data("world", region = "Italy")

# grafico eventi completo
ggplot() +
  geom_polygon(data = italy_map, aes(x = long, y = lat, group = group),
               fill = "gray95", color = "black") +
  geom_point(data = data, aes(x = ev_lon, y = ev_lat),
             color = "blue", size = 2) +
  coord_quickmap() +
  theme_minimal()

itacentraletemp <- subset(data, ev_lat >= lat_min & ev_lat <= lat_max &
                        ev_lon >= lon_min & ev_lon <= lon_max)
nrow(itacentraletemp)
head(itacentraletemp)

# grafico eventi itacentraletemp
ggplot() +
  geom_polygon(data = italy_map, aes(x = long, y = lat, group = group),
               fill = "gray95", color = "black") +
  geom_point(data = itacentraletemp, aes(x = ev_lon, y = ev_lat),
             color = "blue", size = 2) +
  coord_quickmap() +
  theme_minimal()

x1 = 14.4; y1 = 42.8; x2 = 12.6; y2 = 44
m = (y1 - y2)/(x1 - x2)
q = y1 - m*x1

# zoom italia centrale
ggplot() +
  geom_polygon(data = italy_map, aes(x = long, y = lat, group = group), 
               fill = "gray95", color = "gray50") +
  geom_point(data = itacentraletemp, aes(x = ev_lon, y = ev_lat), 
             color = "blue", size = 1.1, alpha = 0.7) +
  scale_x_continuous(breaks = seq(from = floor(lon_min), to = ceiling(lon_max), by = 0.2)) +
  scale_y_continuous(breaks = seq(from = floor(lat_min), to = ceiling(lat_max), by = 0.2)) +
  coord_quickmap(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
    panel.grid.minor = element_blank() ) +
  geom_abline(intercept = q, slope = m, color = "red", linetype = "dashed") +
  geom_circle(aes(x0 = 13.4, y0 = 43.6, r = 0.05), inherit.aes = FALSE, color = "red", linetype = "dashed")

itacentrale <- subset(itacentraletemp, (ev_lat < m*ev_lon + q) |
                                   ( (ev_lon - 13.4)^2 + (ev_lat - 43.6)^2 < 0.05^2) )
itacentrale$EVENT <- droplevels(itacentrale$EVENT)
nrow(itacentrale)

#zoom eventi puliti
p_eventi = ggplot() +
  geom_polygon(data = italy_map, aes(x = long, y = lat, group = group), 
               fill = "gray95", color = "gray50") +
  geom_point(data = itacentrale, aes(x = ev_lon, y = ev_lat), 
             color = "blue", size = 1.1, alpha = 0.7) +
  scale_x_continuous(breaks = seq(from = floor(lon_min), to = ceiling(lon_max), by = 0.5)) +
  scale_y_continuous(breaks = seq(from = floor(lat_min), to = ceiling(lat_max), by = 0.5)) +
  coord_quickmap(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max)) +
  labs(x = "Longitude (°)", y = "Latitude (°)") +
  theme_bw() +
  theme(
    # --- DIMENSIONE DEI NUMERINI (Etichette degli assi) ---
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 16), # Numerini asse X
    axis.text.y = element_text(size = 16),                          # Numerini asse Y
    
    # --- DIMENSIONE DEI TITOLI DEGLI ASSI ---
    axis.title.x = element_text(size = 18),   # Titolo "Longitude (°)"
    axis.title.y = element_text(size = 18),   # Titolo "Latitude (°)"
    
    panel.grid.minor = element_blank()
  )

print(p_eventi)

#plot stazioni
stazioni <- itacentrale |>
  dplyr::distinct(NET_STA, st_lon, st_lat)

p_stazioni = ggplot() +
  geom_polygon(data = italy_map,
               aes(x = long, y = lat, group = group),
               fill = "gray95", color = "gray50") +
  geom_point(data = stazioni,
             aes(x = st_lon, y = st_lat),
             color = "red", size = 1.5) +
  scale_x_continuous(breaks = seq(from = floor(lon_min)-4, to = ceiling(lon_max)+4, by = 2)) +
  scale_y_continuous(breaks = seq(from = floor(lat_min)-4, to = ceiling(lat_max)+4, by = 2)) +
  coord_quickmap() +
  labs(x = "Longitude (°)", y = "Latitude (°)") +
  theme_bw() +
  theme(
    # --- DIMENSIONE DEI NUMERINI ---
    axis.text.x = element_text(size = 16),          # Numerini asse X (qui non avevi l'angolo)
    axis.text.y = element_text(size = 16),          # Numerini asse Y
    
    # --- DIMENSIONE DEI TITOLI DEGLI ASSI ---
    axis.title.x = element_text(size = 18),         # Titolo "Longitude (°)"
    axis.title.y = element_text(size = 18)          # Titolo "Latitude (°)"
  )

print(p_stazioni)

library(ggpubr)

ggarrange(p_stazioni, p_eventi, ncol = 2, nrow = 1)

options(max.print = 2000)
unique(itacentrale$EVENT)
unique(itacentrale$NET_STA)

IM = 0
M = itacentrale$Mw_RCMT
R = itacentrale$RR
vs30 =itacentrale$vs30
sf = itacentrale$sf

ita18 = ITA18_RJB(IM, M, R, vs30, sf)

log_oss = log10(itacentrale$pga)
log_pred = log10(ita18$pred)
l2b <- norm(matrix(log_oss - log_pred, ncol = 1), type = "2")

res = log_oss - log_pred
itacentrale$res = res

library(lme4)
library(lattice)

# res = offset + dBe|evento + dWs

mod <- lmer(res ~ 1 + (1 | EVENT), data = itacentrale, REML = TRUE)
dotplot(ranef(mod))
dBe <- ranef(mod)$EVENT[, "(Intercept)"]
fitted_vals <- fitted(mod)
offset <- fixef(mod)[1]
dWs <- res - fitted_vals
itacentrale$dWs = dWs

################################################################################
################################ VERSIONE GSTAT ################################
################################################################################

# -----------------------------------------------------------------
# 1. LIBRERIE NECESSARIE
# -----------------------------------------------------------------
library(dplyr)
library(sp)
library(gstat) # Ci serve per fit.variogram, vgm, e plot

# -----------------------------------------------------------------
# 2. PREPARAZIONE DATI E TIPI (Come prima)
# -----------------------------------------------------------------

itacentrale <- itacentrale %>%
  mutate(
    NET_STA = as.character(NET_STA),
    EVENT   = as.character(EVENT)
  ) %>%
  filter(!is.na(st_lat) & !is.na(st_lon) & !is.na(dWs))

# -----------------------------------------------------------------
# 3. CREAZIONE COPPIE E CALCOLI
# -----------------------------------------------------------------

print("Inizio creazione coppie...")
df_pairs <- inner_join(
  itacentrale, 
  itacentrale, 
  by = "EVENT",
  suffix = c("_i", "_j"),
  relationship = "many-to-many"
) %>%
  filter(NET_STA_i < NET_STA_j)

print("Inizio calcolo h e gamma...")
coords_i <- as.matrix(df_pairs[, c("st_lon_i", "st_lat_i")])
coords_j <- as.matrix(df_pairs[, c("st_lon_j", "st_lat_j")])

df_for_vario <- data.frame(
  h = spDists(coords_i, coords_j, longlat = TRUE, diagonal = TRUE),
  gamma = 0.5 * (df_pairs$dWs_i - df_pairs$dWs_j)^2
) %>%
  na.omit()

print("Coppie e calcoli completati.")

# -----------------------------------------------------------------
# 4. BINNING VELOCE
# -----------------------------------------------------------------
print("Inizio binning (metodo veloce con bin larghi)...")

# --- MODIFICA QUI ---
# Usiamo bin da 25km (o 20) invece di 10km per smussare i dati
summary(df_for_vario$h)
quantile(df_for_vario$h, probs = c(0.90, 0.95, 0.99))
df_cut <- df_for_vario %>% filter(h <= 230)
boundaries <- quantile(df_cut$h, probs = seq(0, 1, length.out = 11))

exp_vario_df <- df_cut %>%
  mutate(dist_bin = cut(h, breaks = boundaries, include.lowest = TRUE)) %>%
  group_by(dist_bin) %>%
  summarize(
    dist = mean(h, na.rm = TRUE),
    gamma = median(gamma, na.rm = TRUE),
    np = as.numeric(n())
  ) %>%
  ungroup()

# -----------------------------------------------------------------
# 5. IL TRUCCO: TRASFORMARE IL DPLYR in GSTAT (Come prima)
# -----------------------------------------------------------------
print("Applico il 'trucco' per gstat...")

exp_vario <- as.data.frame(exp_vario_df)
exp_vario$dist_bin <- NULL
class(exp_vario) <- c("gstatVariogram", "data.frame")
attr(exp_vario, "model") <- vgm(model = "Nug")

print("Variogramma sperimentale (binned):")
print(exp_vario) # Ora dovresti vedere meno righe (una per ogni bin da 25km)

# -----------------------------------------------------------------
# 6. FIT E PLOT (CON ESTRAZIONE PARAMETRI SUPER-ROBUSTA)
# -----------------------------------------------------------------
print("Inizio fit (metodo gstat)...")

# (La parte di FIT non cambia)

#nugget_start <- min(exp_vario$gamma)
#sill_start <- max(exp_vario$gamma)
#range_start <- median(df_for_vario$h)

nugget_start <- 0.02
sill_start <- 0.08
range_start <- 180

vgm_sph <- vgm(psill = sill_start - nugget_start, model = "Sph", range = range_start, nugget = nugget_start)
vgm_exp <- vgm(psill = sill_start - nugget_start, model = "Exp", range = range_start / 3, nugget = nugget_start)

fit_sph <- fit.variogram(exp_vario, model = vgm_sph, fit.method = 7)
fit_exp <- fit.variogram(exp_vario, model = vgm_exp, fit.method = 7)

print("--- Modello Sferico Fittato ---")
print(fit_sph)
print("--- Modello Esponenziale Fittato ---")
print(fit_exp)

# ----- SOSTITUTO: PLOT CON variogramLine() (usa i modelli fittati) -----

# Assicurati che exp_vario abbia tipi numeric
exp_vario$dist  <- as.numeric(as.character(exp_vario$dist))
exp_vario$gamma <- as.numeric(as.character(exp_vario$gamma))
exp_vario$np    <- as.numeric(as.character(exp_vario$np))

# distanza massima per tracciare le curve (usa max osservato o cutoff)
maxdist <- max(exp_vario$dist, na.rm = TRUE)

# Plot punti sperimentali
plot(
  exp_vario$dist,
  exp_vario$gamma,
  main = "Confronto Modelli Variogramma (Plot Base)",
  xlab = "Distanza h (km)",
  ylab = "Semivarianza (gamma)",
  pch = 19,
  col = "black"
)

# Genera le curve fittate con variogramLine:
# preferiamo la curva derivata da fit.variogram (se disponibile),
# altrimenti usiamo il modello iniziale vgm_sph/vgm_exp

vline_sph <- try(variogramLine(
  if(!inherits(fit_sph, "try-error") && !is.null(fit_sph)) fit_sph else vgm_sph,
  maxdist = maxdist,
  n = 200
), silent = TRUE)

vline_exp <- try(variogramLine(
  if(!inherits(fit_exp, "try-error") && !is.null(fit_exp)) fit_exp else vgm_exp,
  maxdist = maxdist,
  n = 200
), silent = TRUE)

# Traccia le curve se la generazione è riuscita
if(!inherits(vline_sph, "try-error") && !is.null(vline_sph)) {
  lines(vline_sph$dist, vline_sph$gamma, col = "red", lwd = 2)
} else {
  message("Attenzione: non è stato possibile costruire la curva Sferica con variogramLine().")
}

if(!inherits(vline_exp, "try-error") && !is.null(vline_exp)) {
  lines(vline_exp$dist, vline_exp$gamma, col = "blue", lwd = 2)
} else {
  message("Attenzione: non è stato possibile costruire la curva Esponenziale con variogramLine().")
}

# legenda
legend(
  "topleft",
  legend = c("Sperimentale", "Sferico Fittato (gstat)", "Esponenziale Fittato (gstat)"),
  pch = c(19, NA, NA),
  lty = c(NA, 1, 1),
  lwd = c(NA, 2, 2),
  col = c("black", "red", "blue")
)

################################################################################
############################## VERSIONE NLS ####################################
################################################################################

# -----------------------------------------------------------------
# 1. LIBRERIE NECESSARIE
# -----------------------------------------------------------------
library(dplyr)
library(sp)

# -----------------------------------------------------------------
# 2. PREPARAZIONE DATI E TIPI
# -----------------------------------------------------------------

itacentrale <- itacentrale %>%
  mutate(
    NET_STA = as.character(NET_STA), # <-- Correzione per errore '<' su factor
    EVENT   = as.character(EVENT) )
# -----------------------------------------------------------------
# 3. CREAZIONE DELLE COPPIE (WITHIN-EVENT)
# -----------------------------------------------------------------

df_pairs <- inner_join(
  itacentrale, 
  itacentrale, 
  by = "EVENT",
  suffix = c("_i", "_j"),
  relationship = "many-to-many" # Silenzia l'avviso (è corretto)
) %>%
  # Mantiene solo coppie uniche (es. A-B) e scarta (B-A) e (A-A)
  filter(NET_STA_i < NET_STA_j)

print(paste("Numero di coppie uniche create:", nrow(df_pairs)))

# -----------------------------------------------------------------
# 4. CALCOLO 'h' (DISTANZA) E 'gamma' (SEMIVARIANZA)
# -----------------------------------------------------------------

# Calcolo distanza 'h'
coords_i <- as.matrix(df_pairs[, c("st_lon_i", "st_lat_i")])
coords_j <- as.matrix(df_pairs[, c("st_lon_j", "st_lat_j")])
df_pairs$h <- spDists(coords_i, coords_j, longlat = TRUE, diagonal = TRUE)

# Calcolo 'gamma' (la nostra misura di dissimilarità)
df_pairs$gamma <- 0.5 * (df_pairs$dWs_i - df_pairs$dWs_j)^2

# Creiamo il dataframe finale pulito, rimuovendo coppie con calcoli falliti
df_for_vario <- df_pairs[, c("h", "gamma")] %>%
  na.omit()

# -----------------------------------------------------------------
# 5. BINNING MANUALE (VELOCE) CON DPLYR
# -----------------------------------------------------------------

# Definiamo i 'boundaries' (i bin di distanza in km)
# Regola 'seq(0, 100, by = 10)' in base al tuo range di distanze

df_cut <- df_for_vario %>% filter(h <= 230)
boundaries <- quantile(df_cut$h, probs = seq(0, 1, length.out = 11))

exp_vario_df <- df_cut %>%
  mutate(dist_bin = cut(h, breaks = boundaries, include.lowest = TRUE)) %>%
  group_by(dist_bin) %>%
  summarize(
    dist = mean(h, na.rm = TRUE),
    gamma = median(gamma, na.rm = TRUE),
    np = as.numeric(n())
  ) %>%
  ungroup()

print("Variogramma sperimentale (binned):")
print(exp_vario_df)

# -----------------------------------------------------------------
# 6. Sostituto NLS: FIT via ottimizzazione (nlminb)
# -----------------------------------------------------------------

# Assume exp_vario_df già costruito (come nella sezione precedente)
# Usa np come pesi (w = np)

# Sicurezza: conversione tipi numerici
exp_vario_df <- exp_vario_df %>%
  mutate(
    dist  = as.numeric(dist),
    gamma = as.numeric(gamma),
    np    = as.numeric(np)
  )

# -----------------------------------------------------------
# a. Definizione dei modelli locali (non globali)
# -----------------------------------------------------------

sferico_local <- function(h, nugget, psill, range) {
  h <- as.numeric(h)
  if (range <= 0) return(rep(Inf, length(h)))
  a <- h / range
  gamma <- ifelse(
    h <= range,
    nugget + psill * (1.5 * a - 0.5 * a^3),
    nugget + psill
  )
  return(gamma)
}

esponenziale_local <- function(h, nugget, psill, range) {
  h <- as.numeric(h)
  if (range <= 0) return(rep(Inf, length(h)))
  phi <- range / 3.0
  gamma <- nugget + psill * (1 - exp(-h / phi))
  return(gamma)
}

# -----------------------------------------------------------
# b. Funzioni obiettivo (Weighted SSE)
# -----------------------------------------------------------

w <- exp_vario_df$np

obj_wSSE_sph <- function(par, dist, obs, w) {
  if (any(par < 0) || par[3] <= 0) return(1e20)
  pred <- sferico_local(dist, par[1], par[2], par[3])
  sum(w * (obs - pred)^2, na.rm = TRUE)
}

obj_wSSE_exp <- function(par, dist, obs, w) {
  if (any(par < 0) || par[3] <= 0) return(1e20)
  pred <- esponenziale_local(dist, par[1], par[2], par[3])
  sum(w * (obs - pred)^2, na.rm = TRUE)
}

# -----------------------------------------------------------
# c. Stime iniziali
# -----------------------------------------------------------

nug_start   <- max(0, exp_vario_df$gamma[1] * 0.9)
sill_emp    <- max(exp_vario_df$gamma, na.rm = TRUE)
psill_start <- max(1e-8, sill_emp - nug_start)
range_start <- median(df_for_vario$h, na.rm = TRUE)

start_par <- c(
  nugget = nug_start,
  psill  = psill_start,
  range  = range_start
)

cat("Start params (nugget, psill, range):", round(start_par, 6), "\n")

# -----------------------------------------------------------
# d. Ottimizzazione (nlminb)
# -----------------------------------------------------------

# --- Sferico ---
opt_sph <- nlminb(
  start     = start_par,
  objective = obj_wSSE_sph,
  dist      = exp_vario_df$dist,
  obs       = exp_vario_df$gamma,
  w         = w,
  lower     = c(0, 1e-12, 1e-6),
  control   = list(iter.max = 1e4, eval.max = 1e4)
)

params_sph_opt <- opt_sph$par
names(params_sph_opt) <- c("nugget", "psill", "range")

cat("OPT Sferico - params:\n"); print(params_sph_opt)
cat("OPT Sferico - wSSE:", opt_sph$objective, "\n")

# --- Esponenziale ---
opt_exp <- nlminb(
  start     = start_par,
  objective = obj_wSSE_exp,
  dist      = exp_vario_df$dist,
  obs       = exp_vario_df$gamma,
  w         = w,
  lower     = c(0, 1e-12, 1e-6),
  control   = list(iter.max = 1e4, eval.max = 1e4)
)

params_exp_opt <- opt_exp$par
names(params_exp_opt) <- c("nugget", "psill", "range")

cat("OPT Esponenziale - params:\n"); print(params_exp_opt)
cat("OPT Esponenziale - wSSE:", opt_exp$objective, "\n")

# -----------------------------------------------------------
# e. Plot comparativo
# -----------------------------------------------------------

plot(
  exp_vario_df$dist,
  exp_vario_df$gamma,
  xlab = "Distanza media h (km)",
  ylab = "Semivarianza (gamma)",
  main = "Variogramma: experimental + optim(Sph/Exp) + gstat(if any)",
  pch  = 19,
  col  = "black"
)

xvec <- seq(0, max(exp_vario_df$dist, na.rm = TRUE), length.out = 300)

# Linee dai modelli ottimizzati
lines(
  xvec,
  sferico_local(xvec, params_sph_opt["nugget"], params_sph_opt["psill"], params_sph_opt["range"]),
  col = "purple",
  lwd = 2
)

lines(
  xvec,
  esponenziale_local(xvec, params_exp_opt["nugget"], params_exp_opt["psill"], params_exp_opt["range"]),
  col = "darkgreen",
  lwd = 2
)

# Se esistono i fit gstat, li plottiamo (linee tratteggiate)
if (exists("fit_sph") && !inherits(fit_sph, "try-error")) {
  v_sph <- try(variogramLine(fit_sph, maxdist = max(exp_vario_df$dist), n = 200), silent = TRUE)
  if (!inherits(v_sph, "try-error"))
    lines(v_sph$dist, v_sph$gamma, col = "red", lwd = 2, lty = 2)
}

if (exists("fit_exp") && !inherits(fit_exp, "try-error")) {
  v_exp <- try(variogramLine(fit_exp, maxdist = max(exp_vario_df$dist), n = 200), silent = TRUE)
  if (!inherits(v_exp, "try-error"))
    lines(v_exp$dist, v_exp$gamma, col = "blue", lwd = 2, lty = 2)
}

legend(
  "topleft",
  legend = c(
    "experimental",
    "optim Sferico",
    "optim Esponenziale",
    "gstat Sph (dashed)",
    "gstat Exp (dashed)"
  ),
  col  = c("black", "purple", "darkgreen", "red", "blue"),
  lwd  = c(NA, 2, 2, 2, 2),
  pch  = c(19, NA, NA, NA, NA),
  lty  = c(NA, 1, 1, 2, 2),
  bg   = "white",
  cex = 0.8
)

# -----------------------------------------------------------
# f. Tabella riassuntiva risultati
# -----------------------------------------------------------

# Crea un data frame riassuntivo con i parametri stimati e la wSSE
# per ciascun modello (ottimizzato e gstat, se esistono)

# res_tab <- data.frame(
#   method = c("optim_sph", "optim_exp", "gstat_sph", "gstat_exp"),
#   nugget = c(
#     params_sph_opt["nugget"],
#     params_exp_opt["nugget"],
#     if (exists("fit_sph") && !inherits(fit_sph, "try-error"))
#       extract_safe(fit_sph$nugget[fit_sph$model == "Nug"]) else NA,
#     if (exists("fit_exp") && !inherits(fit_exp, "try-error"))
#       extract_safe(fit_exp$nugget[fit_exp$model == "Nug"]) else NA
#   ),
#   psill = c(
#     params_sph_opt["psill"],
#     params_exp_opt["psill"],
#     if (exists("fit_sph") && !inherits(fit_sph, "try-error"))
#       extract_safe(fit_sph$psill[fit_sph$model != "Nug"]) else NA,
#     if (exists("fit_exp") && !inherits(fit_exp, "try-error"))
#       extract_safe(fit_exp$psill[fit_exp$model != "Nug"]) else NA
#   ),
#   range = c(
#     params_sph_opt["range"],
#     params_exp_opt["range"],
#     if (exists("fit_sph") && !inherits(fit_sph, "try-error"))
#       extract_safe(fit_sph$range[fit_sph$model != "Nug"]) else NA,
#     if (exists("fit_exp") && !inherits(fit_exp, "try-error"))
#       extract_safe(fit_exp$range[fit_exp$model != "Nug"]) else NA
#   ),
#   wSSE = c(
#     opt_sph$objective,
#     opt_exp$objective,
#     if (exists("fit_sph") && !inherits(fit_sph, "try-error"))
#       sum(
#         exp_vario_df$np *
#           (exp_vario_df$gamma - sferico_local(
#             exp_vario_df$dist,
#             extract_safe(fit_sph$nugget[fit_sph$model == "Nug"]),
#             extract_safe(fit_sph$psill[fit_sph$model != "Nug"]),
#             extract_safe(fit_sph$range[fit_sph$model != "Nug"])
#           ))^2,
#         na.rm = TRUE
#       ) else NA,
#     if (exists("fit_exp") && !inherits(fit_exp, "try-error"))
#       sum(
#         exp_vario_df$np *
#           (exp_vario_df$gamma - esponenziale_local(
#             exp_vario_df$dist,
#             extract_safe(fit_exp$nugget[fit_exp$model == "Nug"]),
#             extract_safe(fit_exp$psill[fit_exp$model != "Nug"]),
#             extract_safe(fit_exp$range[fit_exp$model != "Nug"])
#           ))^2,
#         na.rm = TRUE
#       ) else NA
#   )
# )

# Stampa a video la tabella riepilogativa
# print(res_tab)

################################################################################
######################## VARIOGRAMMI DIREZIONALI ###############################
################################################################################

# pacchetti utili (caricane già alcuni in testa del tuo script)
library(dplyr)
library(sp)      # hai già usato spDists
library(ggplot2) # per il plot

# ---- Assumiamo che df_pairs esista e contenga:
# st_lon_i, st_lat_i, st_lon_j, st_lat_j, dWs_i, dWs_j
# se non è così, ricrea df_pairs come nel tuo script prima di eseguire questo pezzo

# 1) calcola distanza (km) con spDists (come già facevi) e bearing con formula (no geosphere)
coords_i <- as.matrix(df_pairs[, c("st_lon_i", "st_lat_i")])
coords_j <- as.matrix(df_pairs[, c("st_lon_j", "st_lat_j")])

# distanza in km (vector)
h_km <- spDists(coords_i, coords_j, longlat = TRUE, diagonal = TRUE)

# bearing (formula geodetica, risultato in gradi, riferimento: 0 = North, increases clockwise)
# implementazione vettoriale:
deg2rad <- pi/180
rad2deg <- 180/pi
phi1 <- df_pairs$st_lat_i * deg2rad
phi2 <- df_pairs$st_lat_j * deg2rad
dlon <- (df_pairs$st_lon_j - df_pairs$st_lon_i) * deg2rad

# evitare NaN: se dlon/phi NA -> set 0 temporaneamente
phi1[is.na(phi1)] <- 0
phi2[is.na(phi2)] <- 0
dlon[is.na(dlon)] <- 0

# formula dell'initial bearing (forward azimuth)
yb <- sin(dlon) * cos(phi2)
xb <- cos(phi1) * sin(phi2) - sin(phi1) * cos(phi2) * cos(dlon)
bearing_rad <- atan2(yb, xb)  # rad
bearing_deg <- (bearing_rad * rad2deg + 360) %% 360  # 0..360, 0 = North, clockwise

# converti alla convenzione usata prima (gstat-style): 0 = East, 90 = North
gstat_angle <- (90 - bearing_deg) %% 360

# 2) assegna la direzione target (0,45,90,135) con tolleranza
targets <- c(0, 45, 90, 135)
tol <- 22.5 # tolleranza angolare (±22.5°)

# funzione differenza angolare simmetrica
angdiff <- function(a, b) {
  d <- abs(((a - b + 180) %% 360) - 180)
  return(d)
}

assign_dir <- function(angle, targets, tol) {
  d <- angdiff(angle, targets)
  minidx <- which.min(d)
  if (d[minidx] <= tol) return(as.character(targets[minidx]))
  else return(NA_character_)
}

dirs <- vapply(gstat_angle, assign_dir, FUN.VALUE = character(1), targets = targets, tol = tol)

# 3) crea dataframe con h, gamma e dir
df_dir <- data.frame(
  h = as.numeric(h_km),
  gamma = 0.5 * (df_pairs$dWs_i - df_pairs$dWs_j)^2,
  dir = dirs,
  stringsAsFactors = FALSE
)

# filtra coppie con direzione assegnata e con h non NA
cutoff <- 230  # km come prima
df_dir <- df_dir %>% filter(!is.na(dir) & is.finite(h) & h <= cutoff)

# 4) binning (stesso approccio quantile: 10 bin fino al cutoff)
boundaries <- quantile(df_dir$h, probs = seq(0, 1, length.out = 11), na.rm = TRUE)

df_dir <- df_dir %>%
  mutate(dist_bin = cut(h, breaks = boundaries, include.lowest = TRUE)) %>%
  group_by(dir, dist_bin) %>%
  summarize(
    dist = mean(h, na.rm = TRUE),
    gamma_med = median(gamma, na.rm = TRUE),
    np = n(),
    .groups = "drop"
  )

# 5) mostra sintesi del numero di coppie per direzione
cat("Numero coppie per direzione (totali):\n")
print(df_dir %>% group_by(dir) %>% summarize(total_pairs = sum(np)) )

# 6) plot: una linea per direzione, punti dimensionati con np
ggplot(df_dir, aes(x = dist, y = gamma_med, color = dir, group = dir)) +
  geom_point(aes(size = np), alpha = 0.6) +
  geom_line() +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Distanza h (km)", y = "Semivarianza (gamma, mediana)", color = "Direzione (°)", size = "n coppie") +
  theme_minimal() +
  ggtitle("Variogrammi direzionali (0,45,90,135°) — coppie same-event")


################################################################################
############################### MAPPA DWS ######################################
################################################################################

library(dplyr)
df_plot <- itacentrale %>%
  filter(!is.na(st_lon) & !is.na(st_lat) & !is.na(dWs)) %>%
  group_by(NET_STA, st_lon, st_lat) %>%
  summarize(dWs = median(dWs, na.rm = TRUE), .groups = "drop") %>%
  rename(lon = st_lon, lat = st_lat)

ggplot() +
  # base map
  geom_polygon(
    data = italy_map,
    aes(x = long, y = lat, group = group),
    fill = "gray95",
    color = "gray60",
    linewidth = 0.3
  ) +
  # bolle dei dWs: colore variabile, dimensione fissa
  geom_point(
    data = df_plot,
    aes(x = lon, y = lat, color = dWs),
    size = 3,
    alpha = 0.7
  ) +
  # scala colore monocromatica attorno a steelblue
  scale_color_gradient(
    name = "dWs",
    low  = "red",
    high = "blue"
  ) +
  scale_size_continuous(name = "|dWs|", range = c(1, 8)) +
  coord_quickmap(xlim = c(9, 17), ylim = c(39, 46)) +
  labs(
    title = "Distribuzione spaziale dei dWs",
    x = "Longitudine",
    y = "Latitudine"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold"))

################################################################################
################################ MAPPA PGA #####################################
################################################################################

library(dplyr)

events_df <- itacentrale %>%
  filter(!is.na(ev_lon) & !is.na(ev_lat) & !is.na(pga)) %>%
  group_by(EVENT, ev_lon, ev_lat) %>%
  summarize(
    pga_event    = median(pga, na.rm = TRUE),   # mediana PGA per evento (robusta)
    n_stations   = n(),                         # numero di stazioni che hanno registrato l'evento
    .groups = "drop")

events_df <- events_df %>% arrange(pga_event)

ggplot() +
  # base map
  geom_polygon(
    data = italy_map,
    aes(x = long, y = lat, group = group),
    fill = "gray95",
    color = "gray60",
    linewidth = 0.3
  ) +
  # bolle degli eventi: colore variabile, dimensione fissa
  geom_point(
    data = events_df,
    aes(x = ev_lon, y = ev_lat, color = pga_event),
    size = 3,
    alpha = 0.75
  ) +
  # scala colore monocromatica attorno a steelblue
  scale_color_gradient(
    name = "PGA",
    low  = "blue",
    high = "darkred"
  ) +
  scale_size_continuous(
    name = "|PGA|",
    range = c(2, 12)
  ) +
  coord_quickmap(xlim = c(10, 16), ylim = c(40, 45)) +
  labs(
    title = "Bubble plot: PGA per evento (coordinate epicentrali)",
    subtitle = "Colore ~ PGA",
    x = "Longitudine",
    y = "Latitudine"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold")
  )

################################################################################
############################ MORAN'S I TEST ####################################
################################################################################

library(spdep)
library(dplyr)
library(ggplot2)
library(sf)

# --- AGGIUNTA FONDAMENTALE PER LE COORDINATE ---

# 1. Creiamo un oggetto sf temporaneo usando le coordinate geografiche della stazione
# Assumiamo che nel CSV le colonne si chiamino "st_lon" e "st_lat". 
# Se hanno nomi diversi (es. "rec_lon"), cambiali qui sotto.
ita_sf_temp <- st_as_sf(itacentrale, coords = c("st_lon", "st_lat"), crs = 4326, remove = FALSE)

# 2. Trasformiamo in UTM Zona 33N (EPSG:32633) per avere i metri
# (32633 è lo standard per l'Italia centrale/orientale, coerente con il tuo codice successivo)
ita_utm_temp <- st_transform(ita_sf_temp, 32633)

# 3. Estraiamo le coordinate metriche e convertiamo in KM
coords_utm <- st_coordinates(ita_utm_temp)
itacentrale$X_km <- coords_utm[,1] / 1000  # Converte Est in km
itacentrale$Y_km <- coords_utm[,2] / 1000  # Converte Nord in km

# Ora 'itacentrale' ha le colonne X_km e Y_km pronte per l'aggregazione.

# --- 1. PREPARAZIONE DATI (Aggregazione per Stazione) ---

# Raggruppiamo per stazione e calcoliamo la media di dWs
# (Ci serve un solo punto per ogni coordinata fisica)
station_data <- itacentrale %>%
  group_by(NET_STA) %>%
  summarise(
    dWs_mean = mean(dWs, na.rm = TRUE),
    X_km = first(X_km),  # Le coordinate sono fisse per la stazione
    Y_km = first(Y_km),
    st_lon = first(st_lon),
    st_lat = first(st_lat)
  ) %>%
  filter(!is.na(dWs_mean)) # Rimuoviamo eventuali NA

# Convertiamo in oggetto sf (per i plot)
station_sf <- st_as_sf(station_data, coords = c("st_lon", "st_lat"), crs = 4326)
# O se preferisci usare UTM per il calcolo delle distanze (meglio):
station_sf_utm <- st_as_sf(station_data, coords = c("X_km", "Y_km"), crs = 32633) # Assumo X_km siano già proiettati ma qui servono i numeri puri

print(paste("Numero di stazioni uniche analizzate:", nrow(station_data)))

# --- 2. DEFINIZIONE DEI VICINI (Matrice dei Pesi W) ---

# Estraiamo la matrice delle coordinate direttamente dall'oggetto sf
# (Non serve convertire in "Spatial" o caricare la libreria 'sp')
coords <- st_coordinates(station_sf_utm)

# Metodo: K-Nearest Neighbors (k=5)
# "Chi sono i miei vicini?" -> Le 5 stazioni più vicine
k <- 5
nb <- knn2nb(knearneigh(coords, k = k))

# Creiamo la lista dei pesi (W)
# style = "W" significa "Row-standardized" (la somma dei pesi per ogni riga fa 1)
listw <- nb2listw(nb, style = "W")

# --- 3. MORAN'S I GLOBALE ---

# Test statistico

global_moran <- moran.test(station_data$dWs_mean, listw)
print("RISULTATI MORAN'S I GLOBALE")
print(global_moran)

# INTERPRETAZIONE RAPIDA:
# p-value < 0.05 -> L'autocorrelazione è significativa (non è casuale).
# Moran I statistic > 0 -> Clustering (Zone rosse vicine a rosse, blu a blu).
# Moran I statistic < 0 -> Dispersione (Scacchiera).

# --- 4. MORAN'S I LOCALE (LISA) ---

# Calcolo locale per ogni stazione
local_m <- localmoran(station_data$dWs_mean, listw)

# Uniamo i risultati al dataset originale
station_data$Ii <- local_m[,1] # Indice locale
station_data$E_Ii <- local_m[,2] # Valore atteso
station_data$Var_Ii <- local_m[,3] # Varianza
station_data$Z_Ii <- local_m[,4] # Z-score
station_data$Pr <- local_m[,5] # p-value

# --- CLASSIFICAZIONE QUADRANTI (HH, LL, HL, LH) ---
# Standardizziamo la variabile (dWs) e il Lag spaziale (media dei vicini)
# serve per capire in che quadrante siamo
scale_x <- scale(station_data$dWs_mean) # Z-score del valore alla stazione
lag_x <- lag.listw(listw, station_data$dWs_mean) # Valore medio dei vicini
scale_lag <- scale(lag_x) # Z-score dei vicini

# Creiamo la colonna Cluster
station_data$Cluster <- "Not Sig"

# Soglia di significatività (es. 95%)
alpha <- 0.05

# Logica di assegnazione (solo se p-value < alpha)
# High-High (Rosso): Stazione alta, Vicini alti (Hotspot)
station_data$Cluster[which(scale_x > 0 & scale_lag > 0 & station_data$Pr < alpha)] <- "High-High"
# Low-Low (Blu): Stazione bassa, Vicini bassi (Coldspot)
station_data$Cluster[which(scale_x < 0 & scale_lag < 0 & station_data$Pr < alpha)] <- "Low-Low"
# High-Low (Viola): Stazione alta, Vicini bassi (Outlier spaziale)
station_data$Cluster[which(scale_x > 0 & scale_lag < 0 & station_data$Pr < alpha)] <- "High-Low"
# Low-High (Azzurro): Stazione bassa, Vicini alti (Outlier spaziale)
station_data$Cluster[which(scale_x < 0 & scale_lag > 0 & station_data$Pr < alpha)] <- "Low-High"

# Ordiniamo i fattori per il grafico
station_data$Cluster <- factor(station_data$Cluster, levels = c("High-High", "Low-Low", "High-Low", "Low-High", "Not Sig"))


# --- 5. PLOT DELLA MAPPA LISA (Hotspots) ---

# Definizione colori standard per LISA maps
lisa_colors <- c("High-High" = "red", 
                 "Low-Low" = "blue", 
                 "High-Low" = "pink", 
                 "Low-High" = "lightblue", 
                 "Not Sig" = "grey80")

# Setup per i bordi (uso la tua variabile ita_border_utm se esiste, altrimenti scarica)
# Assicurati di avere ita_border_utm caricato come nel codice precedente
# Se non ce l'hai, scommenta queste righe:
library(rnaturalearth)
ita_border <- ne_countries(scale = "medium", country = "Italy", returnclass = "sf")
ita_border_utm <- st_transform(ita_border, 32633)

# Calcolo Limiti Zoom (come nel codice RDD)
limit_x_m = range(station_data$X_km * 1000)
limit_y_m = range(station_data$Y_km * 1000)
buffer_x = diff(limit_x_m) * 0.1
buffer_y = diff(limit_y_m) * 0.1
final_xlim = c(limit_x_m[1] - buffer_x, limit_x_m[2] + buffer_x)
final_ylim = c(limit_y_m[1] - buffer_y, limit_y_m[2] + buffer_y)


# GRAFICO
p_lisa <- ggplot() +
  # Bordi Italia
  geom_sf(data = ita_border_utm, fill = "white", color = "black", size = 0.3) +
  
  # Punti Stazioni Colorati per Cluster
  geom_point(data = station_data, 
             aes(x = X_km * 1000, y = Y_km * 1000, color = Cluster, size = (Cluster != "Not Sig")), 
             alpha = 0.8) +
  
  scale_color_manual(values = lisa_colors) +
  scale_size_manual(values = c("TRUE" = 2, "FALSE" = 1), guide = "none") + # Punti significativi più grandi
  
  coord_sf(xlim = final_xlim, ylim = final_ylim, expand = FALSE) +
  
  labs(title = "LISA Map (Local Moran's I) - Residui dWs",
       subtitle = "Cluster spaziali significativi (p < 0.05)",
       x = "Est (m)", y = "Nord (m)") +
  theme_minimal()

print(p_lisa)


################################################################################
############################# ANALISI DUPLICATI ################################
################################################################################

library(dplyr)

# 1. Creiamo un dataframe univoco delle stazioni
# Prendiamo solo Nome, Latitudine e Longitudine originali
stazioni_check <- itacentrale %>%
  select(NET_STA, st_lat, st_lon) %>%
  distinct() # Rimuove le righe duplicate (stessa stazione, stesso evento)

# 2. Arrotondamento
# A volte le coordinate differiscono per millimetri (es. 42.1000001 vs 42.1000000).
# Arrotondiamo alla 4a cifra decimale (~11 metri) o 5a (~1 metro) per raggrupparle.
stazioni_check <- stazioni_check %>%
  mutate(
    lat_round = round(st_lat, 5),
    lon_round = round(st_lon, 5)
  )

# 3. Troviamo i duplicati spaziali
duplicati_spaziali <- stazioni_check %>%
  group_by(lat_round, lon_round) %>%
  summarise(
    n_nomi = n_distinct(NET_STA),       # Quanti nomi diversi ci sono qui?
    elenco_nomi = paste(unique(NET_STA), collapse = " | ") # Elencali
  ) %>%
  filter(n_nomi > 1) %>%
  ungroup()

# 4. Stampiamo il risultato
if(nrow(duplicati_spaziali) > 0) {
  print("ATTENZIONE: Trovate stazioni con nomi diversi nelle stesse coordinate:")
  print(duplicati_spaziali)
} else {
  print("Tutto ok: Nessuna stazione sovrapposta con nomi diversi.")
}

# 1. Recupera le stazioni "sospette" (quelle con stesse coordinate arrotondate)
# (Assumendo che tu abbia già calcolato 'duplicati_spaziali' dal messaggio precedente)
coordinate_duplicate <- duplicati_spaziali %>% 
  select(lat_round, lon_round)

# 2. Estrai i dati di queste stazioni
dati_sospetti <- itacentrale %>%
  mutate(
    lat_round = round(st_lat, 9),
    lon_round = round(st_lon, 9)
  ) %>%
  inner_join(coordinate_duplicate, by = c("lat_round", "lon_round"))

# 3. Controllo incrociato: Condividono eventi?
analisi_overlap <- dati_sospetti %>%
  group_by(lat_round, lon_round) %>%
  summarise(
    nomi_stazioni = paste(unique(NET_STA), collapse = " vs "),
    # Contiamo quanti eventi unici ci sono in totale in quel punto
    eventi_totali = n_distinct(EVENT),
    # Contiamo quante righe ci sono (se > eventi_totali, c'è sovrapposizione)
    righe_totali = n(),
    # Calcoliamo la sovrapposizione esatta
    eventi_condivisi = sum(duplicated(EVENT))
  ) %>%
  mutate(
    DIAGNOSI = case_when(
      eventi_condivisi == 0 ~ "STESSA STAZIONE (Cambio nome temporale)",
      eventi_condivisi > 0  ~ "DUE STRUMENTI DIVERSI (Colocated)"
    )
  )

print(analisi_overlap)

################################################################################
#################### DISTRIBUZIONE STAZIONI PER EVENTO #########################
################################################################################

library(dplyr)
library(ggplot2)

# 1. Calcola quante stazioni hanno registrato ogni singolo evento
stats_eventi <- itacentrale %>%
  group_by(EVENT) %>%
  summarise(n_stazioni = n()) # Usa n() perché ogni riga è una registrazione singola

# 2. Statistiche rapide numeriche
print(summary(stats_eventi$n_stazioni))
cat(paste("Numero totale di eventi:", nrow(stats_eventi), "\n"))
cat(paste("Eventi con meno di 5 stazioni:", sum(stats_eventi$n_stazioni < 5), "\n"))

# 3. Genera l'istogramma
p_hist <- ggplot(stats_eventi, aes(x = n_stazioni)) +
  geom_histogram(binwidth = 5, fill = "firebrick", color = "white", alpha = 0.8) +
  geom_vline(aes(xintercept = mean(n_stazioni)), color="blue", linetype="dashed", size=1) +
  labs(title = "Distribuzione del Numero di Stazioni per Evento",
       subtitle = paste("Media:", round(mean(stats_eventi$n_stazioni), 1), 
                        "| Mediana:", median(stats_eventi$n_stazioni)),
       x = "Numero di Stazioni che hanno registrato l'evento",
       y = "Conteggio Eventi (Frequenza)") +
  theme_minimal() +
  theme(plot.title = element_text(face="bold"))

print(p_hist)

# Distribuzione eventi per stazione

p_hist <- ggplot(stats_stazioni, aes(x = n_eventi)) +
  geom_histogram(binwidth = 5, fill = "forestgreen", color = "black", alpha = 0.8, linewidth = 0.3) +
  # Zoom sull'asse X per tagliare la coda vuota
  coord_cartesian(xlim = c(0, 150)) +
  # Imposta tacche personalizzate per maggiore chiarezza
  scale_x_continuous(breaks = seq(0, 150, by = 25)) +
  labs(x = expression(paste("Number of recorded events per station (", N[i], ")")),
       y = "Frequency (Number of stations)") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    panel.grid.minor = element_blank()
  )

print(p_hist)


# 1. Calcola quanti eventi ha registrato ogni singola stazione (N_i)
stats_stazioni <- itacentrale |>
  group_by(NET_STA) |>
  summarise(n_eventi = n()) 

# 2. Genera l'istogramma con la tua estetica
p_hist <- ggplot(stats_stazioni, aes(x = n_eventi)) +
  geom_histogram(binwidth = 5, fill = "forestgreen", color = "white", alpha = 0.8) +
  geom_vline(aes(xintercept = median(n_eventi)), color = "firebrick", linetype = "dashed", linewidth = 0.5) +
  coord_cartesian(xlim = c(0, 150)) +
  scale_x_continuous(breaks = seq(0, 150, by = 25)) +
  labs(x = expression(paste("Number of recorded events per station (", N[i], ")")),
       y = "Frequency (number of stations)") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  )

print(p_hist)
