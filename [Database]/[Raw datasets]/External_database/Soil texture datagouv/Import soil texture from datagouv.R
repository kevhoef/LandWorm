library(raster)

# Convert the gps_x and gps_y columns to numeric type (if not already converted)
mydata$gps_x <- as.numeric(mydata$gps_x)
mydata$gps_y <- as.numeric(mydata$gps_y)


# Specify the paths to your TIFF files
argile_0_5_path <- "[Database]/[Raw datasets]/External_database/Soil texture datagouv/argile.0_5.tif"
argile_5_15_path <- "[Database]/[Raw datasets]/External_database/Soil texture datagouv/argile.5_15.tif"
limon_0_5_path <- "[Database]/[Raw datasets]/External_database/Soil texture datagouv/limon.0_5.tif"
limon_5_15_path <- "[Database]/[Raw datasets]/External_database/Soil texture datagouv/limon.5_15.tif"
sable_0_5_path <- "[Database]/[Raw datasets]/External_database/Soil texture datagouv/sable.0_5.tif"
sable_5_15_path <- "[Database]/[Raw datasets]/External_database/Soil texture datagouv/sable.5_15.tif"

# List of full file paths to TIFF files
tif_file_paths <- c(argile_0_5_path, argile_5_15_path, limon_0_5_path, limon_5_15_path, sable_0_5_path, sable_5_15_path)

# Create an empty dataframe to store soil values
soil_values <- data.frame(gps_x = mydata$gps_x, gps_y = mydata$gps_y)

# Loop to process each TIFF file and add corresponding columns to the dataframe
for (tif_file_path in tif_file_paths) {
  # Read the TIFF file as a raster
  raster_data <- raster(tif_file_path)
  
  # Create a SpatialPoints object with GPS coordinates in geographic coordinate system (WGS84)
  gps_coords <- data.frame(x = mydata$gps_x, y = mydata$gps_y)
  gps_coords_sp <- SpatialPoints(gps_coords, proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  # Project the GPS coordinates to the raster's coordinate system
  gps_coords_proj <- spTransform(gps_coords_sp, crs(raster_data))
  
  # Extract the soil values for these projected GPS coordinates
  soil_values[[basename(tif_file_path)]] <- extract(raster_data, gps_coords_proj)
}

# Add the extracted soil values as new columns to mydata
mydata <- cbind(mydata, soil_values)

