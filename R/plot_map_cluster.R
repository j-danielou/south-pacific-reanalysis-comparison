#---------------------
#Packages & Sources
library(dplyr)
library(ggplot2)
source("R./function_taylor.R")

#---------------------
#Function
reorder_and_plot_clusters = function(df_all, clust1_path, clust2_path, prefix, colors = NULL) {
  library(dplyr)
  library(ggplot2)
  source("R./function_taylor.R")
  # Charger les clusters
  clust1 = readRDS(clust1_path)
  clust2 = readRDS(clust2_path)
  
  # Appliquer le 1er clustering
  df_all$cluster = clust1$clusters
  df_all$cluster_final = clust2$clusters[df_all$cluster]
  
  if (prefix=="full"){
    # Charger les centroïdes
    centroids_df = as.data.frame(clust2$centroids)[, 1:9]
    colnames(centroids_df) = c("glorys_crmsd", "glorys_R", "glorys_sd", "bran_crmsd", "bran_R", "bran_sd", "hycom_crmsd", "hycom_R", "hycom_sd")
    centroids_df$cluster = 1:nrow(centroids_df)
    
    # Distance à l'idéal
    target = c(crmsd = 0, R = 1, sd = 1)
    centroids_df$distance_to_target = sqrt(
      (centroids_df[[1]] - target["crmsd"])^2 +
        (centroids_df[[2]] - target["R"])^2 +
        (centroids_df[[3]] - target["sd"])^2 +
        (centroids_df[[4]] - target["crmsd"])^2 +
        (centroids_df[[5]] - target["R"])^2 +
        (centroids_df[[6]] - target["sd"])^2 +
        (centroids_df[[7]] - target["crmsd"])^2 +
        (centroids_df[[8]] - target["R"])^2 +
        (centroids_df[[9]] - target["sd"])^2
    )
  }else{
    # Charger les centroïdes
    centroids_df = as.data.frame(clust2$centroids)[, 1:3]
    colnames(centroids_df) = paste0(prefix, c("_crmsd", "_R", "_sd"))
    centroids_df$cluster = 1:nrow(centroids_df)
    
    # Distance à l'idéal
    target = c(crmsd = 0, R = 1, sd = 1)
    centroids_df$distance_to_target = sqrt(
      (centroids_df[[1]] - target["crmsd"])^2 +
        (centroids_df[[2]] - target["R"])^2 +
        (centroids_df[[3]] - target["sd"])^2
    )
  }
  
  # Réordonner
  centroids_df = centroids_df |>
    arrange(distance_to_target) |>
    mutate(new_cluster = row_number())
  
  # Recode
  recode_vector = centroids_df$new_cluster
  names(recode_vector) = centroids_df$cluster
  df_all$cluster_final_reordered = recode_vector[as.character(df_all$cluster_final)]
  
  cluster_labels = c("1" = "SPG", "2" = "TSW", "3" = "POI",
                     "4" = "TSE", "5" = "TSP", "6" = "PNG")
  
  # Couleurs
  n_clust = nrow(centroids_df)
  if (is.null(colors)) {
    colors = RColorBrewer::brewer.pal(n = n_clust, name = "Set3")
    length(colors) = n_clust
  }
  names(colors) = as.character(1:n_clust)
  
  # if (prefix != "full"){
  #   col = colorful::divergencePalette(n=6, col=c("dodgerblue3", "firebrick4"), intensity=1)
  #   #col[7] = "grey10"
  #   length(col) = 6
  #   names(col) = as.character(1:n_clust)
  #   
  #   # Taylor Diagram
  #   x11(width = 10, height = 10)
  #   taylor.diagram(ref = c(1, 1.1), model = c(1, 0.9), label = "Ref", add = FALSE)
  #   pch_set = 16
  #   for (j in 1:n_clust) {
  #     crmsd = centroids_df[[1]][j]
  #     R = centroids_df[[2]][j]
  #     sd = centroids_df[[3]][j]
  #     x = sd * R
  #     y = sd * sin(acos(R))
  #     points(x, y, col = col[as.character(j)], pch = pch_set, cex = 1.3)
  #   }
  #   legend("topright", legend = prefix,
  #          col = "black", pch = pch_set, pt.cex = 2, inset = c(-0.05, -0.05))
  #   legend("topright", legend = paste("Regions", 1:n_clust),
  #          col = col, pch = 15, pt.cex = 2, inset = c(-0.05, 0.08))
  #   dev.copy(png, file = paste0("C:/Users/jdanielou/Desktop/taylor-diagram-cluster-supp",prefix,".png"), width = 10, height = 10, units = "in", res = 150)
  #   dev.off()
  # }else{
  #   col = colorful::divergencePalette(n=8, col=c("dodgerblue3", "firebrick4"), intensity=1)
  #   col[7] = "grey10"
  #   length(col) = 7
  #   names(col) = as.character(1:n_clust)
  #   
  #   # Taylor Diagram
  #   x11(width = 10, height = 10)
  #   taylor.diagram(ref = c(1, 1.1), model = c(1, 0.9), label = "Ref", add = FALSE)
  #   pch_set = c(16,17,18) 
  #   # Boucle sur les 3 produits
  #   for (i in 1:3) {
  #     
  #     # Sélectionner les colonnes appropriées
  #     R_col = centroids_df[[ (i-1)*3 + 2 ]]
  #     sd_col = centroids_df[[ (i-1)*3 + 3 ]]
  #     
  #     # Ajouter les points pour chaque cluster
  #     for (j in 1:n_clust) {
  #       
  #       # Ajouter point
  #       points(sd_col[j] * R_col[j], sd_col[j] * sin(acos(R_col[j])), col = col[as.character(j)], pch = pch_set[i], cex = 1.2)
  #       
  #     }
  #   }
  #   legend("topright", legend = c("GLORYS", "BRAN", "HYCOM"),
  #          col = "black", pch = pch_set, pt.cex = 2, inset = c(-0.05, -0.05))
  #   legend("topright", legend = paste("Cluster", 1:n_clust),
  #          col = col, pch = 15, pt.cex = 2, inset = c(-0.05, 0.08))
  #   dev.copy(png, file = paste0("C:/Users/jdanielou/Desktop/taylor-diagram-cluster-supp-",prefix,".png"), width = 10, height = 10, units = "in", res = 150)
  #   dev.off()
  # }
  
  # Carte
  c2t = gts:::coord2text
  cl = gts:::checkLongitude
  
  x11(width = 16, height = 8)
  p = ggplot() +
    geom_polygon(data = map_data("world2"), aes(x = long, y = lat, group = group),
                 fill = "grey80", color = "black", linewidth = 0.1, inherit.aes = FALSE) +
    geom_point(data = df_all, aes(x = lon, y = lat, color = as.factor(cluster_final_reordered)),
               size = 1) +
    coord_fixed(xlim = range(df_all$lon), ylim = range(df_all$lat), expand = FALSE) +
    scale_color_manual(values = colors, labels = cluster_labels) +
    scale_x_continuous(breaks = pretty(df_all$lon, n=6), labels = c2t(cl(pretty(df_all$lon, n=6)), "lon"), expand = c(0, 0))+
    scale_y_continuous(breaks = pretty(df_all$lat, n=4), labels = c2t(pretty(df_all$lat, n=4), "lat"), expand = c(0,0))+
    labs(color = "Region", title = "") +
    guides(color = guide_legend(override.aes = list(size = 10), label.theme = element_text(size = 15),
                                title.theme = element_text(size = 15))) + 
    theme_minimal() +
    theme(plot.title = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15))
  
  print(p)
  dev.copy(png, file = paste0("C:/Users/jdanielou/Desktop/map-cluster-",prefix,".png"), width = 16, height = 7, units = "in", res = 150)
  dev.off() 
  #Sys.sleep(1)
  
  # Retourne les résultats utiles si besoin
  return(list(df = df_all, centroids = centroids_df, colors = colors, clust1 = clust1))
}

glorys_df = read.csv("C:/Users/jdanielou/Desktop/plot_internship/csv/metric_csv/glorys/taylor_metrics_pixel_glorys.csv")
bran_df = read.csv("C:/Users/jdanielou/Desktop/plot_internship/csv/metric_csv/bran/taylor_metrics_pixel_bran.csv")
hycom_df = read.csv("C:/Users/jdanielou/Desktop/plot_internship/csv/metric_csv/hycom/taylor_metrics_pixel_hycom.csv")

# Renommage des colonnes
glorys_df = glorys_df |>
  rename(glorys_crmsd = crmsd, glorys_R = R, glorys_sd = sd, glorys_rmse = rmse, glorys_biais = biais)

bran_df = bran_df |>
  rename(bran_crmsd = crmsd, bran_R = R, bran_sd = sd, bran_rmse = rmse, bran_biais = biais)

hycom_df = hycom_df |>
  rename(hycom_crmsd = crmsd, hycom_R = R, hycom_sd = sd, hycom_rmse = rmse, hycom_biais = biais)

# Merge des données
df_all = glorys_df %>%
  left_join(bran_df, by = c("lon", "lat")) %>%
  left_join(hycom_df, by = c("lon", "lat")) %>%
  filter(complete.cases(.))

col = colorful::divergencePalette(n=8, col=c("dodgerblue3", "firebrick4"), intensity=0.3)
length(col) = 7


result = reorder_and_plot_clusters(
  df_all = df_all,
  clust1_path = "C:/Users/jdanielou/Desktop/clust_centroid_full.rds",
  clust2_path = "C:/Users/jdanielou/Desktop/clust_centroid_2_full.rds",
  prefix = "full",
  colors =  col #c("seagreen", "lightseagreen", "lightblue1", "khaki", "orange", "brown4", "grey10")
)

saveRDS(result$df, file ="C:/Users/jdanielou/Desktop/df_all.rds")


##################################
#######   NEW VERSION   ##########
##################################
#---------------------
#Packages & Sources
library(dplyr)
library(ggplot2)
#---------------------
df_all = readRDS("C:/Users/jdanielou/Desktop/df_all.rds")

# Ton ordre défini
ordre_défini = c(1, 4, 2, 3, 5, 6)

# Création d’un vecteur de mapping : nom = ancien cluster, valeur = nouveau cluster
mapping = setNames(seq_along(ordre_défini), ordre_défini)

# Ajout de la colonne cluster_ordered
glorys_mat$cluster_ordered = mapping[as.character(glorys_mat$cluster)]

# Réordonner
glorys_mat = glorys_mat |>
  arrange(cluster_ordered)


cluster_labels = c("1" = "SPG", "2" = "TSW", "3" = "POI",
                   "4" = "TSE", "5" = "TSP", "6" = "PNG")
colors = colorful::divergencePalette(n=8, col=c("dodgerblue3", "firebrick4"), intensity=0.5)
length(colors) = 6

# Carte améliorée
c2t = gts:::coord2text
cl = gts:::checkLongitude

x11(width = 16, height = 8)
p = ggplot() +
  # Fond de carte plus clair
  geom_polygon(data = map_data("world2"), aes(x = long, y = lat, group = group),
               fill = "grey95", color = "grey70", linewidth = 0.2) +
  
  # Points colorés par cluster
  geom_point(data = glorys_mat, aes(x = lon, y = lat, color = as.factor(cluster_ordered)),
             size = 0.85) +
  
  # Axes fixes
  coord_fixed(xlim = range(glorys_mat$lon), ylim = range(glorys_mat$lat), expand = FALSE) +
  
  # Couleurs et labels
  scale_color_manual(values = colors, labels = cluster_labels, name = "Region") +
  
  # Axes longitude/latitude
  scale_x_continuous(
    breaks = pretty(glorys_mat$lon, n = 6),
    labels = c2t(cl(pretty(glorys_mat$lon, n = 6)), "lon"),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = pretty(glorys_mat$lat, n = 4),
    labels = c2t(pretty(glorys_mat$lat, n = 4), "lat"),
    expand = c(0, 0)
  ) +
  
  # Légende & axes
  guides(color = guide_legend(
    override.aes = list(size = 6),
    label.theme = element_text(size = 14),
    title.theme = element_text(size = 15)
  )) +
  
  # Thème graphique
  theme_minimal(base_size = 15) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.text = element_text(size = 13),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.key = element_rect(fill = NA),
    panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA),
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 10)
  )

print(p)

# Export PNG haute résolution
dev.copy(png, file = "C:/Users/jdanielou/Desktop/map-cluster.png", width = 16, height = 7, units = "in", res = 300)
dev.off()



###########################
#######   TEST   ##########
###########################

prefix="glorys"
centroids_df = as.data.frame(result$clust1$centroids)[, 1:3]
colnames(centroids_df) = paste0(prefix, c("_crmsd", "_R", "_sd"))
centroids_df$cluster = 1:nrow(centroids_df)

# Distance à l'idéal
target = c(crmsd = 0, R = 1, sd = 1)
centroids_df$distance_to_target = sqrt(
  (centroids_df[[1]] - target["crmsd"])^2 +
    (centroids_df[[2]] - target["R"])^2 +
    (centroids_df[[3]] - target["sd"])^2
)

# Réordonner
centroids_df = centroids_df |>
  arrange(distance_to_target) |>
  mutate(new_cluster = row_number())

# Recode
recode_vector = centroids_df$new_cluster
names(recode_vector) = centroids_df$cluster


col = colorful::divergencePalette(n=1000, col=c("dodgerblue3", "firebrick4"), intensity=0.2)
result_df_all = result$df
result_df_all$cluster_final_reordered = recode_vector[as.character(result_df_all$cluster)]

c2t = gts:::coord2text
cl = gts:::checkLongitude

x11(width = 16, height = 12)
ggplot() +
  geom_polygon(data = map_data("world2"), aes(x = long, y = lat, group = group),
               fill = "grey80", color = "black", linewidth = 0.1, inherit.aes = FALSE) +
  geom_point(data = result_df_all, aes(x = lon, y = lat, color = as.factor(cluster_final_reordered)),
             size = 1) +
  coord_fixed(xlim = range(result_df_all$lon), ylim = range(result_df_all$lat), expand = FALSE) +
  scale_color_manual(values = col) +
  scale_x_continuous(breaks = pretty(result_df_all$lon, n=6), labels = c2t(cl(pretty(result_df_all$lon, n=6)), "lon"), expand = c(0, 0))+
  scale_y_continuous(breaks = pretty(result_df_all$lat, n=4), labels = c2t(pretty(result_df_all$lat, n=4), "lat"), expand = c(0,0))+
  labs(color = "Cluster", title = "") +
  guides(color = guide_legend(override.aes = list(size = 10), label.theme = element_text(size = 15),
                              title.theme = element_text(size = 15))) + 
  theme_minimal() +
  theme(plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.position = "none")
        
dev.copy(png, file = paste0("C:/Users/jdanielou/Desktop/map-cluster-supp",prefix,".png"), width = 16, height = 7, units = "in", res = 200)
dev.off() 
