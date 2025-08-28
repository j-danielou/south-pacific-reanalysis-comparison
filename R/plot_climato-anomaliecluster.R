library(dplyr)


# Charger les climatologies et anomalies
climatologies = list(
  cmems = readRDS("C:/Users/jdanielou/Desktop/rds/cmems_climatologies_sss.rds"),
  bran = readRDS("C:/Users/jdanielou/Desktop/rds/bran_climatologies_sss.rds"),
  hycom = readRDS("C:/Users/jdanielou/Desktop/rds/hycom_climatologies_sss.rds"),
  glorys = readRDS("C:/Users/jdanielou/Desktop/rds/glorys_climatologies_sss.rds")
)

anomalies = list(
  cmems = readRDS("C:/Users/jdanielou/Desktop/rds/cmems_anomalies_sss.rds"),
  bran = readRDS("C:/Users/jdanielou/Desktop/rds/bran_anomalies_sss.rds"),
  hycom = readRDS("C:/Users/jdanielou/Desktop/rds/hycom_anomalies_sss.rds"),
  glorys = readRDS("C:/Users/jdanielou/Desktop/rds/glorys_anomalies_sss.rds")
)

clusters = paste0("cluster", 1:6)
colors = c("black", "blue", "darkgreen", "red")
model_names = names(climatologies)

# Ouvrir la fenêtre graphique
x11(width = 18, height = 16)
par(mfrow = c(6, 2), mar = c(1.1, 2.2, 1.1, 2.2), oma = c(1.3, 5.5, 1, 5), 
    cex.axis = 1.8, cex.lab = 2, las=1)

for (i in 1:6) {
  cluster = clusters[i]

  ## --- PLOT CLIMATOLOGIE ---
  # Calcul des limites
  clim_vals = unlist(lapply(climatologies, function(m) m[[cluster]]$climatology))
  ylim_clim = range(clim_vals, na.rm = TRUE) + c(-0.2, 0.2)
  months_numeric = as.numeric(climatologies$cmems[[cluster]]$month)
  month_labels = levels(climatologies$cmems[[cluster]]$month)

  plot(x = months_numeric,
       y = climatologies$cmems[[cluster]]$climatology,
       type = "l", col = colors[1], lwd = 1.5,
       ylim = ylim_clim,
       ylab = "PSU",
       xaxt = ifelse(i < 7, "n", "s"),
       yaxt = "n",
       xlab = "")
  
  axis(2, las = 1, at = pretty(ylim_clim, n = 3), labels = pretty(ylim_clim, n = 3))
  
  
  for (j in 2:length(model_names)) {
    lines(climatologies[[j]][[cluster]]$climatology, col = colors[j], lwd = 1.5)
  }
  if (i == 1) {
    usr <- par("usr")  # coordonnées du plot : c(x1, x2, y1, y2)
    
    # Légende gauche
    legend(x = usr[1] + 0.52 * diff(usr[1:2]),
           y = usr[4] - 0.02 * diff(usr[3:4]),
           legend = c("CMEMS", "BRAN"),
           col = colors[1:2], lwd = 2, bty = "n", cex = 1.5, x.intersp = 0.5)
    
    # Légende droite
    legend(x = usr[2] - 0.28 * diff(usr[1:2]),
           y = usr[4] - 0.02 * diff(usr[3:4]),
           legend = c("HYCOM", "GLORYS"),
           col = colors[3:4], lwd = 2, bty = "n", cex = 1.5, x.intersp = 0.5)
  }
  if (i == 6) {
    axis(1, at = 1:12, labels = month_labels, cex.axis = 1.8)
  }
  
  ## --- PLOT ANOMALIES ---
  # Plage temporelle
  date_min = as.Date("1994-01-01")
  date_max = as.Date("2015-12-31")
  
  # Référence OISST filtrée
  cmems_df = anomalies$cmems[[cluster]] %>%
    filter(time >= date_min & time <= date_max)
  
  # Définir les limites y
  ylim_anom = range(unlist(lapply(anomalies, function(m) {
    m[[cluster]] %>%
      filter(time >= date_min & time <= date_max) %>%
      pull(anomaly)
  })), na.rm = TRUE) + c(-0.2, 0.2)
  
  # Plot OISST
  plot(x = cmems_df$time,
       y = cmems_df$anomaly,
       type = "l", col = colors[1], lwd = 1.5,
       ylim = ylim_anom,
       xaxt = ifelse(i < 7, "n", "s"),
       yaxt = "n",
       xlab = "", ylab = "")
  
  axis(2, las = 1, at = pretty(ylim_anom, n = 3), labels = pretty(ylim_anom, n = 3))
  
  # Tracer les autres modèles et stocker r/RMSE
  r_text = "r-values : "
  rmse_text = "RMSE : "

  for (j in 2:length(model_names)) {
    mod_df = anomalies[[j]][[cluster]] %>%
      filter(time >= date_min & time <= date_max)

    lines(mod_df$time, mod_df$anomaly, col = colors[j], lwd = 1.5)

    common = inner_join(cmems_df, mod_df, by = "time", suffix = c("_ref", "_mod"))
    r = round(cor(common$anomaly_ref, common$anomaly_mod, use = "complete.obs"), 2)
    rmse = round(sqrt(mean((common$anomaly_ref - common$anomaly_mod)^2, na.rm = TRUE)), 2)

    abbrev = substr(toupper(model_names[j]), 1, 1)  # G, B, H
    r_text = paste0(r_text, abbrev, "=", r, ifelse(j < length(model_names), "; ", ""))
    rmse_text = paste0(rmse_text, abbrev, "=", rmse, ifelse(j < length(model_names), "; ", ""))
  }

  #Affichage des valeurs r/RMSE en haut du plot
  #mtext(paste(r_text, "  ", rmse_text), side = 4, line = 0.2, cex = 1.1)
  
  # Axe x seulement pour dernière ligne
  if (i == 6) {
    tick_dates = seq(date_min, date_max, by = "2 years")
    axis(1, at = tick_dates, labels = format(tick_dates, "%Y"), las = 1)
  }
  
}  
dev.copy(png, file = "C:/Users/jdanielou/Desktop/climato-anomalies-cluster-sss.png", width = 18, height = 10, units = "in", res = 150)
dev.off() 
