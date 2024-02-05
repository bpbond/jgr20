
# /----------------------------------------------------------------------------#
#/   FUNCTION REPROJECTING 

### FIX THE CORNER DUPLICATE RASTERS IN ROBINSON PROJ

# Function converting format & proj 
# Question: is the the one that prevents exceeded area in robinson proj?
WGSraster2dfROBIN <- function(r){
    library(terra)
    crs(r) <- crs('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
    # r <- terra::rast(r)
    r_robin <- terra::project(r, '+proj=robin', method='near', mask=T)
    r_robin_df <- as.data.frame(r_robin, xy=TRUE, na.rm=TRUE) 
    return(r_robin_df) }

