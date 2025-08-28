# -------------------------
# PACKAGES
# -------------------------
library(ncdf4)
library(gts)
library(zoo)

# -------------------------
# FONCTION : Calcul de l’ONI + Classification
# -------------------------
calculate_oni = function(model_file, ersst_file, date_start_model, ersst_range = NULL) {
  model_data = read.table(model_file)$V1
  ersst = read.table(ersst_file)$V1
  dates_model = seq(from = as.Date(date_start_model), by = "month", length.out = length(model_data))
  dates_ersst = seq(from = as.Date("1854-01-01"), by = "month", length.out = length(ersst))
  
  if (is.null(ersst_range)) {
    ersst_sub = ersst[which(dates_ersst >= min(dates_model) & dates_ersst <= max(dates_model))]
  } else {
    ersst_sub = ersst[ersst_range[1]:ersst_range[2]]
  }
  
  model_3mois = rollmean(model_data, k = 3, align = "center")
  ersst_3mois = rollmean(ersst_sub, k = 3, align = "center")
  dates_3mois = dates_model[2:(length(dates_model) - 1)]
  
  get_climatology = function(date) {
    year = as.numeric(format(date, "%Y"))
    year_group_start = floor((year - 1) / 5) * 5 + 1
    clim_start = year_group_start - 30
    clim_end   = year_group_start - 1
    date_start = as.Date(paste0(clim_start, "-01-01"))
    date_end   = as.Date(paste0(clim_end, "-12-31"))
    idx = which(dates_ersst >= date_start & dates_ersst <= date_end)
    
    if (length(idx) < 12 * 20) {
      warning(paste("Climatologie insuffisante pour", year, ":", clim_start, "-", clim_end))
      return(rep(NA, 12))
    }
    
    clim_data = ersst[idx]
    mois = rep(1:12, length.out = length(clim_data))
    climatologie = tapply(clim_data, mois, mean)
    return(climatologie)
  }
  
  oni_model = numeric(length(model_3mois))
  oni_ersst = numeric(length(ersst_3mois))
  
  for (i in seq_along(dates_3mois)) {
    date_i = dates_3mois[i]
    mois_i = as.numeric(format(date_i, "%m"))
    clim = get_climatology(date_i)
    oni_model[i] = model_3mois[i] - mean(clim[c(ifelse(mois_i-1==0,12,mois_i-1), mois_i, ifelse(mois_i+1==13,1,mois_i+1))])
    oni_ersst[i] = ersst_3mois[i] - mean(clim[c(ifelse(mois_i-1==0,12,mois_i-1), mois_i, ifelse(mois_i+1==13,1,mois_i+1))])
  }
  
  # classify_oni = function(x) {
  #   list(
  #     VS_S = which(x >= 1.5),
  #     M_W = which(x >= 0.5 & x < 1.5),
  #     neutral = which(x > -0.5 & x < 0.5),
  #     W_M = which(x >= -1.5 & x < -0.5)
  #   )
  # }
  
  classify_oni = function(x) {
    list(
      VS = which(x >= 2),
      S = which(x >= 1.5 & x < 2),
      M = which(x >= 1.0 & x < 1.5),
      W = which(x >= 0.5 & x < 1.0),
      neutral = which(x > -0.5 & x < 0.5),
      W_M = which(x >= -1.5 & x < -0.5)
    )
  }
  
  return(list(
    oni_model = oni_model,
    oni_ersst = oni_ersst,
    dates_3mois = dates_3mois,
    categories_model = classify_oni(oni_model),
    categories_ersst = classify_oni(oni_ersst)
  ))
}

# -------------------------
# CHARGEMENT DES DONNÉES
# -------------------------

ersst_path = "C:/Users/jdanielou/Desktop/plots_internship/csv/ts_csv/ersst/ersst_zone_3.4.csv"

oni_bran = calculate_oni(
  model_file = "C:/Users/jdanielou/Desktop/plot_internship/csv/ts_csv/bran/bran_zone_3.4.csv",
  ersst_file = ersst_path,
  date_start_model = "1993-01-01",
  ersst_range = c(1669, 2016)
)

oni_glorys = calculate_oni(
  model_file = "C:/Users/jdanielou/Desktop/plot_internship/csv/ts_csv/glorys/glorys_zone_3.4.csv",
  ersst_file = ersst_path,
  date_start_model = "1993-01-01",
  ersst_range = c(1669, 2016)
)

oni_hycom = calculate_oni(
  model_file = "C:/Users/jdanielou/Desktop/plot_internship/csv/ts_csv/hycom/hycom_zone_3.4.csv",
  ersst_file = ersst_path,
  date_start_model = "1994-01-01",
  ersst_range = c(1681, 1944)
)

# -------------------------
# FONCTION PLOT DENSITÉ PAR CATÉGORIE
# -------------------------
plot_density_per_category = function(oni_bran, oni_glorys, oni_hycom) {
  categories = c("VS", "S", "M", "W", "neutral", "W_M")
  title = c("Very Strong", "Strong", "Moderate", "Weak", "Neutral", "Weak - Moderate")
  i = 0
  
  colors = c(BRAN = "gold", GLORYS = "forestgreen", HYCOM = "dodgerblue", ERSST = "purple")
  
  x11(width = 18, height = 10)
  par(mfrow = c(2, 3), oma = c(0, 0, 0, 0))
  
  for (cat in categories) {
    i=i+1
    val_bran   = oni_bran$oni_model[oni_bran$categories_model[[cat]]]
    val_glorys = oni_glorys$oni_model[oni_glorys$categories_model[[cat]]]
    val_hycom  = oni_hycom$oni_model[oni_hycom$categories_model[[cat]]]
    
    val_ersst = oni_bran$oni_ersst[oni_bran$categories_ersst[[cat]]]
    
    # Choix conditionnels des labels
    xlab_val <- if (cat %in% c("neutral")) "ONI (°C)" else ""
    ylab_val <- if (cat %in% c("W", "VS")) "Density" else ""
    
    # Étendre légèrement les bornes X
    all_vals <- c(val_bran, val_glorys, val_hycom, val_ersst)
    range_val <- range(all_vals, na.rm = TRUE)
    delta <- diff(range_val) * 0.5
    range_elargie <- c(range_val[1] - delta, range_val[2] + delta)
    
    plot(density(val_ersst, na.rm = TRUE),
         col = colors["ERSST"], lwd = 2, lty = 2,
         xlab = xlab_val, ylab = ylab_val,
         ylim = c(0, 4), xlim = range_elargie, main = "",
         las = 1, cex.lab = 1.5, cex.axis = 1.3)
    
    title(main = title[i], line = 0.5, cex.main = 1.6)
    
    lines(density(val_bran, na.rm = TRUE), col = colors["BRAN"], lwd = 2)
    lines(density(val_glorys, na.rm = TRUE), col = colors["GLORYS"], lwd = 2)
    lines(density(val_hycom, na.rm = TRUE), col = colors["HYCOM"], lwd = 2)
    
    # Légende uniquement pour le panneau "M_W"
    if (cat == "M") {
      legend("topright", legend = c("BRAN", "GLORYS", "HYCOM", "ERSST"),
             col = colors, lty = c(1, 1, 1, 2), lwd = 1.8, bty = "n", cex = 1.2, inset = c(0.05, 0))
    }
  }
}

# -------------------------
# APPEL FINAL
# -------------------------
plot_density_per_category(oni_bran, oni_glorys, oni_hycom)

dev.copy(png, file = "C:/Users/jdanielou/Desktop/oni_comparison_density.png", width = 14, height = 10, units = "in", res = 150)
dev.off()
