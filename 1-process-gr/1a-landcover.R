######################
# 1a - Process MODIS landcover (MCD12Q1) - create forest mask for greenup
#
# Runtime: < 15 min total
######################


# Raw data DL instructions ------------------------------------------------

#go to: https://lpdaacsvc.cr.usgs.gov/appeears/task/area
#upload cropped shp file
#select MCD12Q1 -> all layers
#select GeoTIFF and geographic projection
#DL MCD12Q1
#select MCD12Q2 -> all layers
#select GeoTIFF and geographic projection
#DL MCD12Q2


# specify dir -------------------------------------------------------------

dir <- '~/Google_Drive/R/pheno_trends/'
run_date <- '2020-08-06'


# load packages -----------------------------------------------------------

library(raster)
library(sp)
library(ggplot2)
library(rgdal)


# MCD12Q1 - Landcover -----------------------------------------------------

#LCCS1 Land Cover Layer
#from user guide:
#1 - barren
#2 - perm snow and ice
#3 - water bodies
#11 - evergreen needleleaf forests
#12 - evergreen broadleaf forests
#13 - deciduous needleleaf forests
#14 - deciduous broadleaf forests
#15 - mixed broadleaf/needleleaf forests
#16 - mixed broadleaf evergreen/deciduous forests
#21 - open forests
#22 - sparse forests
#31 - dense herbaceous
#32 - sparse herbaceous
#41 - dense shrublands
#42 - shrubland/grassland mosaics
#43 - sparse shrublands
#255 - unclassified


# #read in  tif - 2018 LCCS1
# setwd(paste0(dir, 'pheno_trends/Data/environment/RAW/MCD12Q1/'))
# lc1 <- raster::raster('MCD12Q1.006_LC_Prop1_doy2018001_aid0001.tif')

#read in  tif - 2017 LCCS1
setwd(paste0(dir, '/Data/environment/RAW/MCD12Q1/'))
lc1 <- raster::raster('MCD12Q1.006_LC_Prop1_doy2017001_aid0001.tif')


# filter by cover type ----------------------------------------------------

# ever_needle <- raster::match(lc1, 11, nomatch = NA)
# ever_broad <- raster::match(lc1, 12, nomatch = NA)
# decid_needle <- raster::match(lc1, 13, nomatch = NA)
# decid_broad <- raster::match(lc1, 14, nomatch = NA)
# mix_bn <- raster::match(lc1, 15, nomatch = NA)
# mix_broad <- raster::match(lc1, 16, nomatch = NA)
# open_for <- raster::match(lc1, 21, nomatch = NA)
# sparse_for <- raster::match(lc1, 22, nomatch = NA)
# blank <- raster::match(lc1, c(1, 2, 11, 12, 13, 
#                               14, 15, 16, 21, 22, 
#                               31, 32, 41, 42, 43, 
#                               255), nomatch = NA)


# write to KMLs -----------------------------------------------------------

# setwd(paste0(dir, 'pheno_trends/Data/environment/landcover_kml'))
# raster::KML(ever_needle, file = 'ever_needle', col = 'darkolivegreen4', 
#     maxpixels = ncell(lc1)/100)
# raster::KML(ever_broad, file = 'ever_broad', col = 'darkolivegreen2', 
#     maxpixels = ncell(lc1)/100)
# raster::KML(decid_needle, file = 'decid_needle', col = 'darkolivegreen3', 
#     maxpixels = ncell(lc1)/100)
# raster::KML(decid_broad, file = 'decid_broad', col = 'darkolivegreen1', 
#     maxpixels = ncell(lc1)/100)
# raster::KML(mix_bn, file = 'mix_bn', col = 'gold4', 
#     maxpixels = ncell(lc1)/100)
# raster::KML(mix_broad, file = 'mix_broad', col = 'gold3', 
#     maxpixels = ncell(lc1)/100)
# raster::KML(open_for, file = 'open_for', col = 'gold2', 
#     maxpixels = ncell(lc1)/100)
# raster::KML(sparse_for, file = 'sparse_for', col = 'gold1', 
#     maxpixels = ncell(lc1)/100)
# raster::KML(blank, file = 'blank', col = 'white', 
#     maxpixels = ncell(lc1)/100)


# create mask of just forest land cover type ------------------------------

lc_for_mask <- lc1 %in% c(11, 12, 13, 14, 15, 16, 21, 22)


# save to RDS -------------------------------------------------------------

#create output dir if it doesn't exist
ifelse(!dir.exists(paste0(dir, 'Data/environment/processed/', run_date)),
       dir.create(paste0(dir, 'Data/environment/processed/', run_date)),
       FALSE)

setwd(paste0(dir, 'Data/environment/processed/', run_date))
saveRDS(lc_for_mask, paste0('lc_for_mask-', run_date,'.rds'))

#copy script to dir
system(paste0('cp ', dir, 'Scripts/1-process-gr/1a-landcover.R ', dir, 'Data/environment/processed/', run_date, '/1a-landcover-', run_date, '.R'))
