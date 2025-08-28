library(ncdf4)
library(fields)
library(ncdf4)
library(rgl)

nc = nc_open("C:/Users/jdanielou/Desktop/glorys-v1-southpacific-to-monthly-199301-199303.nc")

relief_3d = function(nc, varid, longitude = "longitude", latitude = "latitude", depths = "depth", scale_factor = 0.3, load = TRUE){
  
  # Ouvrir le fichier NetCDF
  if (load == TRUE){
    var = ncvar_get(nc, varid)  # [lon, lat, depth, time]
    lon = ncvar_get(nc, longitude)
    lat = ncvar_get(nc, latitude)
    depth = ncvar_get(nc, depths)
  }else{
    var = data$var  # [lon, lat, depth, time]
    lon = data$lon
    lat = data$lat
    depth = data$depth
  }
  
  nc_close(nc)
  
  # Sélection d’un mois
  temp_month = var[,,,1]   # [lon, lat, depth, time]
  
  # Dimensions
  nx = length(lon)
  ny = length(lat)
  
  # Matrice pour stocker la profondeur maximale non NA
  last_valid_depth = matrix(NA, nrow = nx, ncol = ny)
  
  for (i in 1:nx) {
    for (j in 1:ny) {
      profile = temp_month[i, j, ]
      if (all(is.na(profile))) next
      last_index = max(which(!is.na(profile)))
      last_valid_depth[i, j] = depth[last_index]
    }
  }
  
  # Préparer les matrices pour plot3d (rgl attend des grilles)
  lon_mat = matrix(rep(lon, each = ny), nrow = nx)
  lat_mat = matrix(rep(lat, times = nx), nrow = nx)
  depth_mat = last_valid_depth
  
  # Nettoyage : retirer les NA
  valid_mask = !is.na(depth_mat)
  lon_mat[!valid_mask] = NA
  lat_mat[!valid_mask] = NA
  depth_mat[!valid_mask] = NA
  
  # Palette inversée : rouge pour les faibles profondeurs, bleu pour les grandes
  col_scale = viridis(100, option = "C", direction = -1)
  col_idx = round(scales::rescale(depth_mat, to = c(1, 100)))
  col_mat = col_scale[col_idx]
  
  # Plot 3D relief
  open3d()
  surface3d(
    x = lon,
    y = lat,
    z = -depth_mat,   # Inverser l'axe des profondeurs
    color = t(col_mat),
    back = "lines"
  )
  aspect3d(1, 1, scale_factor)
  axes3d()
  title3d("Relief 3D - derniere profondeur avec temperature", 
          xlab = "Longitude", 
          ylab = "Latitude", 
          zlab = "Profondeur")
  
 return(list(var = var, lon = lon, lat = lat, depth = depth))
}

data = relief_3d(nc, varid = "to", load = TRUE)
