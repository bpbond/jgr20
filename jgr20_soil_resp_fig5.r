# Script plotting Figure 5 for the paper for JGR 20-year anniversary:
# Bond-Lamberty et al. - Twenty years of progress, challenges, and opportunities in measuring and understanding soil respiration


# /----------------------------------------------------------------------------#
#/   Import libraries
library(colorspace)
library(cowplot)
library(dplyr)
library(here)
library(ggplot2)
library(lubridate)
library(sf)
library(StatCompLab)
library(terra)
library(tidyverse)


# /----------------------------------------------------------------------------#
#/  Get SRDB data, and filter records
srdb <- read.csv('../data/srdb/srdb-data.csv') %>% 
        filter(Study_midyear >= 1991) %>% 
        filter(Duplicate_record == '') %>% 
        filter(Latitude >= -60) %>% 
        dplyr::select(Record_number, Longitude, Latitude) %>% 
        filter(complete.cases(.))

# Convert SRDB to spatial object using lat/lon coordinates
srdb <- srdb %>%  
        st_as_sf(coords = c("Longitude","Latitude")) %>% 
        st_set_crs(st_crs(4326))


# /----------------------------------------------------------------------------#
#/  Get GSWP3 temperature data
# https://data.isimip.org/datasets/ae8ab033-cea9-4fca-a275-5c46533ec970/
# The GSWP3-W5E5 dataset is based on GSWP3 v1.09 (Kim 2017) and W5E5 v2.0 (Cucchi et al. 2020, Lange et al. 2021). 
# The GSWP3 dataset is a dynamically downscaled and bias-adjusted version of the 
# Twentieth Century Reanalysis version 2 (20CRv2; Compo et al. 2011). 

# Read GSWP global temperature grids; one file per decade; 
gswp3_1991_2000 <- rast('../data/gswp3/gswp3-w5e5_obsclim_tas_global_daily_1991_2000.nc')
gswp3_2001_2010 <- rast('../data/gswp3/gswp3-w5e5_obsclim_tas_global_daily_2001_2010.nc')
gswp3_2011_2019 <- rast('../data/gswp3/gswp3-w5e5_obsclim_tas_global_daily_2011_2019.nc')

# Stack decadal raster files of daily temperature
gswp3_daily_1991_2019 <- c(gswp3_1991_2000, gswp3_2001_2010, gswp3_2011_2019)


# Create a sequence of years between the start and end dates
yr_seq <- year(seq(ymd("1991-01-01"), ymd("2019-12-31"), by = "days"))

# Calculate mean temperature per year 
gswp3_1991_2019_annual <- tapp(gswp3_daily_1991_2019, yr_seq, fun=mean, na.rm=T, cores=12)

# Clean env
rm(gswp3_1991_2000, gswp3_2001_2010, gswp3_2011_2019, gswp3_daily_1991_2019)


# /----------------------------------------------------------------------------#
#/  Apply landmask to GSWP3 temperature grids
landseamask <- rast('../data/gswp3/landseamask.nc')
landseamask[landseamask == 0] <- NA

# Apply landmask
gswp3_1991_2019_annual <- gswp3_1991_2019_annual * landseamask

# Crop out Antarctica
gswp3_1991_2019_annual <- crop(gswp3_1991_2019_annual, ext(-180, 180, -60, 90))


# /----------------------------------------------------------------------------#
#/  Get pixel-wise trend slope of GSWP3 temperature 
years <- seq(1991, 2019)
bfun <- function(x) { if (is.na(x[1])){ c(NA, NA) } else {lm(x ~ years)$coefficients}} 
lm_coef <- app(gswp3_1991_2019_annual, bfun)
names(lm_coef) <- c('intercept', 'slope')

# Convert to slope decG/yr to degC/decade
global_tas_slope <- lm_coef[[2]] *10

# Compute gridcell area
global_tas_slope <- c(global_tas_slope, cellSize(global_tas_slope, unit='km'))

global_tas_slope_df <- data.frame(global_tas_slope) %>% 
                       mutate(set = 'Global land area') %>% 
                       filter(!is.na(slope))


# /----------------------------------------------------------------------------#
#/  Extract slope at SRDB points locations

srdb_tas_slope <- extract(global_tas_slope, srdb) %>% 
                  mutate(set = 'SRDB locations')

# /----------------------------------------------------------------------------#
#/  Run KS test to evaluate if the distribution of temperature trend at SRDB locations 
#   is part of the global background slope distribution.

ks.test(global_tas_slope_df$slope, srdb_tas_slope$slope, alternative = c("two.sided"))


# /----------------------------------------------------------------------------#
#/  Figure 5A - cumulative distribution of SRDB and global temperature trend

# Bind rows from SRDB and global trends to the same dataframe for plotting
tas_slope <- bind_rows(srdb_tas_slope, global_tas_slope_df)

# Import custom figure theme
source('themes/line_plot_theme.r')


cdf_srdb_tas_slope <-
    ggplot()+
    geom_vline(xintercept=0, color='grey85', linewidth=0.4) +
    geom_density(data=tas_slope, aes(x=slope, group=set, color=set), 
                 size=0.5, outline.type = 'upper', n=200, adjust=1.75) +
    xlab(expression(paste('Temperature change over 1991-2019 (°C ', decade^-1, ')'))) +
    ylab('Density') + 
    scale_color_manual(values = c('black','red')) +
    scale_y_continuous(expand=c(0,0))+
    line_plot_theme() + 
    theme(legend.position = c(0.9, 0.15),
          plot.margin=unit(c(4, 30, 2, 30), 'mm'))


# /----------------------------------------------------------------------------#
#/  Figure 5B-  Map of SRDB location over global background warming trend

# Import background polygon data for global map
source('fcns/get_country_bbox_shp_for_ggplot_map.r')

### Reproject GSWP3 trend to Robinson 
# Import function reprojecting raster data to Robinson projection and converting to dataframe
source('fcns/wgsRaster2dfRobin.r')
# Apply function
global_tas_slope_robin_df <- WGSraster2dfROBIN(global_tas_slope)

### Reproject SRDB point locations to Robinson
srdb_robin <- st_transform(srdb, st_crs("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m"))
srdb_robin$set <- 'SRDB locations'

# Import custom theme
source('themes/raster_map_theme.r')



tas_slope_map <- 
    ggplot()+
    # Add GSWP3 trend layer
    geom_raster(data=global_tas_slope_robin_df, aes(x, y, fill=slope)) +
    # Add coastline line
    geom_sf(data=coastsCoarse_robin, fill=NA, color='grey40', linewidth=0.15)+
    # Add SRDB points
    geom_sf(data=srdb_robin, size=0.25, aes(shape=set)) +
    # Add bounding box 
    geom_sf(data=bbox_robin, fill=NA, color='black', linewidth=0.15)+

    scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0, rev=T) + 
    raster_map_theme() +
    theme(legend.position = 'bottom',
          plot.margin=unit(c(3, -11, 2, -11), 'mm')) +

    guides(shape=guide_legend(order = 1, title="", label.position ='top'),
           fill = guide_colorbar(
                        nbin=10, raster=F, barheight = 0.4, barwidth=6,
                        frame.colour=c('black'), frame.linewidth=0.25,
                        ticks.colour='black',  direction='horizontal',
                        title.position='top',
                        title = expression(paste('Temperature change over 1991-2019 (°C ', decade^-1, ')'))))


# /----------------------------------------------------------------------------#
#/  Figure 5C - Latitudinal band distribution of SRDB location

# Preprocess data for plotting: Get number of SRDB counts per latitude band
srdb_df <- srdb %>% 
           mutate(lat_band=round(Latitude, -1)) %>% 
           group_by(lat_band) %>% 
           summarize(frac=n() / nrow(.)) %>% 
           ungroup() %>%
           add_row(lat_band = -60, frac = 0) %>% 
           arrange(lat_band) %>% 
           mutate(lat_band_fac=as.factor(lat_band)) 

# Preprocess data for plotting: Group GSWP3 pixels per  
global_tas_slope_df <- as.data.frame(global_tas_slope, xy=TRUE, na.rm=TRUE) %>% 
                        mutate(lat_band=round(y, -1)) %>% 
                        add_row(lat_band = 90, slope = 0) %>% 
                        arrange(lat_band) %>% 
                        mutate(lat_band_fac=as.factor(lat_band))


#  Plot latitudinal distribution of SRDB and global warming trend
lat_band_plot <- 
    ggplot()+
    # Add 0 line
    geom_hline(yintercept=0, color='grey85', linewidth=0.4) +
    # Global warming as boxplot
    geom_boxplot(data= global_tas_slope_df, aes(x=lat_band_fac, y=slope), 
                 color='red', width=0.6, size=0.3, outlier.shape=NA) +
    # SRDB as points
    geom_point(data=srdb_df, aes(x=lat_band_fac, y=frac), color='black', size=1.2)+
    line_plot_theme() +  
    scale_y_continuous(
        name = expression(paste('Temperature change over 1991-2019 (°C ', decade^-1, ')')),
        sec.axis = sec_axis( trans=~.*100, name="Percentage of SRDB samples (%)")) +
    xlab('Latitude band') +
    theme( axis.line.x.bottom = element_line(color = "red"), 
           axis.ticks.x.bottom = element_line(color = "red"),
           axis.text.x.bottom = element_text(color = "red"),
           axis.title.x.bottom = element_text(color = "red")) +
    coord_flip() +
    theme(legend.position = c(0.9, 0.15),
          plot.margin=unit(c(4, 30, 2, 30), 'mm'))




# /----------------------------------------------------------------------------#
#/  Combine panels into a single multi-panel figure 5

#
srdb_tas_slope_plot <- cowplot::plot_grid( cdf_srdb_tas_slope, lat_band_plot, tas_slope_map, 
                                           ncol=1,
                                           rel_heights = c(0.7, 0.9, 1),
                                           labels=c('a','b', 'c'))

# Fill white background
srdb_tas_slope_plot<- 
    cowplot::ggdraw(srdb_tas_slope_plot) + 
    theme(plot.background = element_rect(fill="white", color = NA))


# Save figure as PNG
ggsave("../output/figures/srdb_tas_slope_plot_v4.png",
       srdb_tas_slope_plot, width=180, height=210, dpi=500, units="mm")

# Save figure as PDF
ggsave("../output/figures/srdb_tas_slope_plot_v4.pdf",
       srdb_tas_slope_plot, width=180, height=210,  units="mm") 
