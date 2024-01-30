# /-----------------------------------------------------------------------------
#/    get polygon of continent outline                                     -----
# FIXED - SO THERE'S NO LINE STRETCHING BETWEEN ALASKA-KAMATCHAKA when plotting map.
# natearth_dir <- "../../chap5_global_inland_fish_catch/data/gis/nat_earth"
natearth_dir <- "../data/nat_earth"



sf_use_s2(FALSE) # https://stackoverflow.com/questions/69227019/odd-sfst-crop-behaviour

# crs(com_ext) <- 
robin_crs <- st_crs("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m")
library(rgdal)
library(raster)

com_ext <- sf::st_bbox(c(xmin = -180, xmax = 180, ymin = -60, ymax = 90), crs = 4326) %>% 
    sf::st_as_sfc() %>%
    sf::st_transform(crs = 4326)



# read and reproject outside box
bbox <- st_read(natearth_dir, "ne_110m_wgs84_bounding_box") %>% sf::st_transform(crs = 4326)
bbox <- st_intersection(bbox, com_ext, dimension = "polygon")  # Set smaller extent, excl. Antarctica
bbox_robin <- st_transform(bbox, robin_crs)  # reproject bounding box



#/    get polygon of continent outline                                     -----
library(rworldmap)
library(sp)
data(coastsCoarse)

coastsCoarse <- st_as_sf(coastsCoarse) 
coastsCoarse_robin <- st_intersection(coastsCoarse, bbox) %>% 
    st_transform(robin_crs)


# read and reproject countries  -  and ticks to  Robinson 
countries <- st_read(natearth_dir, "ne_110m_admin_0_countries")
countries_robin <- st_intersection(countries, bbox) %>% st_transform(robin_crs)

