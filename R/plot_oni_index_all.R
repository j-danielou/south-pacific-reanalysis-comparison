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
  
  classify_oni = function(x) {
    list(
      VS_S = which(x >= 1.5),
      M_W = which(x >= 0.5 & x < 1.5),
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
# CHEMINS & CALCULS
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
# PLOT COMPARATIF MULTIPANEL
# -------------------------

plot_oni_comparison = function(oni_model, oni_ersst, dates_3mois, model_name, model_color) {
  # el_nino_colors = c("#FFCCCC", "#FF9999", "#FF6666", "#FF3333", "#FF0000")
  # la_nina_colors = c("#CCCCFF", "#9999FF", "#6666FF", "#3333FF", "#0000FF")
  el_nino_thresholds = c(0.5, 1.0, 1.5, 2.0, 3.0)
  la_nina_thresholds = c(-0.5, -1.0, -1.5, -2.0, -3.0)
  
  xlim_val = if (model_name != "HYCOM") c(8790, 18620) else c(9060, 16450)
  
  plot(dates_3mois, oni_model, type = "n", ylim = c(-2.79, 2.79), xlim = xlim_val,
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
  
  # title(main = paste(model_name, "(RMSE =", round(sqrt(mean((oni_model - oni_ersst)^2, na.rm = TRUE)), 3), ")"),
  #       cex.main = 1.8, line = 0.35) 
  
   # for (i in seq_along(el_nino_thresholds)) {
   #   ybottom = if (i == 1) 0 else el_nino_thresholds[i - 1]
   #   ytop = el_nino_thresholds[i]
   #   rect(min(dates_3mois), ybottom, max(dates_3mois), ytop, col = el_nino_colors[i], border = NA)
   # }
   # for (i in seq_along(la_nina_thresholds)) {
   #   ytop = if (i == 1) 0 else la_nina_thresholds[i - 1]
   #   ybottom = la_nina_thresholds[i]
   #   rect(min(dates_3mois), ybottom, max(dates_3mois), ytop, col = la_nina_colors[i], border = NA)
   # }

  abline(h = seq(-2.5, 2.5, 0.5), col = "gray", lty = 2)
  abline(h = 0, col = "black", lwd = 2)
  
  lines(dates_3mois, oni_model, col = model_color, lwd = 3)
  lines(dates_3mois, oni_ersst, col = "purple", lwd = 3, lty = 2)
  
  axis(2, at = seq(-3, 3, 1), las = 1, cex.axis = 1.6)
  if (model_name == "GLORYS") mtext("ONI Index (°C)", side = 2, line = 3, cex = 1.2)

  years = seq(from = as.Date(format(min(dates_3mois), "%Y-01-01")),
              to = as.Date(format(max(dates_3mois) , "%Y-01-01")),
              by = "year")
  years_subset = c(years[seq(1, length(years), by = 3)], 18993)
  axis(1, at = years_subset, labels = FALSE)
  
  text(x = years_subset, y = par("usr")[3] - 0.3,
       labels = format(years_subset, "%Y"), srt = 0, adj = 0.5,
       xpd = TRUE, cex = 1.6)
  
  # axis(1, at = years, labels = FALSE)
  # text(x = years, y = par("usr")[3] - 0.2,
  #      labels = format(years, "%Y"), srt = 0, adj = 1,
  #      xpd = TRUE, cex = 1.6)
  # 
  
  if (model_name == "BRAN"){
    legend("topright", legend = c(paste("ONI Model"), "ONI ERSST"),
         text.col = "gray99", col = c(model_color, "purple"),
         lty = c(1, 2), lwd = 2.5, bty = "0", bg = "gray", cex = 1.4)
  }
}
# -------------------------
# VISUALISATION
# -------------------------

x11(width = 15, height = 12)
par(mfrow = c(3, 1), oma = c(0, 0.5, 0, 0))

plot_oni_comparison(oni_model = oni_bran$oni_model, oni_ersst = oni_bran$oni_ersst, dates_3mois = oni_bran$dates_3mois, model_name = "BRAN", model_color = "gold")
plot_oni_comparison(oni_glorys$oni_model, oni_glorys$oni_ersst, oni_glorys$dates_3mois, "GLORYS", "gold")
plot_oni_comparison(oni_hycom$oni_model, oni_hycom$oni_ersst, oni_hycom$dates_3mois, "HYCOM", "gold")


dev.copy(png, file = "C:/Users/jdanielou/Desktop/oni_comparison_all.png", width = 14, height = 13, units = "in", res = 150)
dev.off()

