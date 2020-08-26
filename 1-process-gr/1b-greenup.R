######################
# 1b - Process MODIS greenup (MCD12Q2)
#
# Runtime: ~ 1/2 day
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

dir <- '~/Google_Drive/R/pheno_trends'
run_date <- '2020-08-06'


# load packages -----------------------------------------------------------

library(raster)
library(sp)
library(rgdal)
library(dggridR)


# veg cycles --------------------------------------------------------------

#check number of veg cycles - vast majority of cells have 1 - only going to use phenometrics from first veg cycle (product comes with greenup from two veg cycles, if applicable - see user manual)
# setwd(paste0(dir, 'pheno_trends/Data/environment/RAW/MCD12Q2/'))
# nc <- raster::raster(paste0('MCD12Q2.006_NumCycles_doy2001001_aid0001.tif'))
# nc1 <- nc == 1
# plot(nc1)
# nc2 <- nc >= 2
# plot(nc2)


# Greenup -----------------------------------------------------------------

tt <- proc.time()

#hexgrid
setwd(paste0(dir, 'Data/hex_grid_crop/'))
hexgrid <- rgdal::readOGR('hex_grid_crop.shp', verbose = FALSE)
hexgrid_cells <- as.numeric(as.character(hexgrid@data$cell))

#read in forest mask
setwd(paste0(dir, 'Data/environment/processed/', run_date))
lc_for_mask <- readRDS(paste0('lc_for_mask-', run_date ,'.rds'))

#extract greenup for each year, filter by QA score and forest land cover type
#'Greenup' = onset (15%)
#'MidGreenup' = mid (50%)
#forest = TRUE -> filter by forest land cover type
gr_pro_fun <- function(gr_type = 'Greenup', years = 2001:2017, forest = TRUE)
{
  setwd(paste0(dir, 'Data/environment/RAW/MCD12Q2/'))
  # define functions to convert binary QA score to interpretable scale - modified from MCD12Q2 User Manual
  # output score: 0 - best, 1 - good, 2 - fair, 3 - poor
  # output in order: Greenup, MidGreenup, Maturity, Peak, Senescence, MidGreendown, Dormancy
  # onset: ST = 1, FN = 2
  # mid: ST = 3, FN = 4
  # if (gr_type == 'Greenup')
  # {
  #   ST <- 1
  #   FN <- 2
  # }
  # if (gr_type == 'MidGreenup')
  # {
  #   ST <- 3
  #   FN <- 4
  # }
  # 
  # UnpackDetailedQA <- function(x, ST, FN)
  # {
  #   bits <- as.integer(intToBits(x))
  #   #only first 4 bits (Greenup and MidGreenup)
  #   quals <- sapply(seq(ST, FN, by = 2), function(i) sum(bits[i:(i+1)] * 2^c(0, 1)))
  #   return(quals)
  # }
  lna <- length(years) * length(hexgrid_cells)
  out_df <- data.frame(year = rep(NA, lna),
                       cell = rep(NA, lna),
                       gr_mn = rep(NA, lna),
                       gr_sd = rep(NA, lna),
                       gr_ncell = rep(NA, lna),
                       gr_pcell = rep(NA, lna),
                       gr_type = gr_type)
  counter <- 1
  #process rasters
  for (i in 1:length(years))
  {
    #i <- 1
    print(paste0('Processing: ', years[i]))
    
    #QA overall
    qao <- raster::raster(paste0('MCD12Q2.006_QA_Overall_0_doy', years[i], '001_aid0001.tif'))
    
    #tQA <- raster::calc(qa, fun = function(x) UnpackDetailedQA(x, ST = ST, FN = FN))
    
    #create mask based on QA scores
    #only 'best' and 'good'
    QA_mask <- qao %in% c(0,1)
    rm(qao)
    
    #greenup (onset or mid)
    gr <- raster::raster(paste0('MCD12Q2.006_', gr_type, '_0_doy', years[i], '001_aid0001.tif'))
    
    #filter by QA score
    gr_f <- raster::mask(gr, QA_mask, maskvalue = 1,
                         inverse = TRUE)
    rm(QA_mask)
    
    if (forest == TRUE)
    {
      #filter by forest mask
      gr_f2 <- raster::mask(gr_f, lc_for_mask, maskvalue = 1,
                            inverse = TRUE)
      rm(gr_f)
    } else {
      gr_f2 <- gr_f
      rm(gr_f)
    }
    
    #convert filtered greenup raster to spdf
    green_spdf <- raster::rasterToPoints(gr_f2, spatial = TRUE)
    #unfiltered greenup - to cal num of pixels in each cell
    ogr_spdf <- raster::rasterToPoints(gr, spatial = TRUE)
    #get mean greenup within each hex cell
    gr_mn <- sp::over(hexgrid, green_spdf, na.rm = TRUE, fn = mean)
    #get sd greenup within each hex cell
    gr_sd <- sp::over(hexgrid, green_spdf, na.rm = TRUE, fn = sd)
    #convet to julian day
    jday <- as.numeric(julian(as.Date(gr_mn[,1], 
                                      origin = as.Date('1970-01-01')), 
                              origin = as.Date(paste0(years[i], '-01-01'))))
    #number of cells with greenup within each hex cell
    gr_nc <- sp::over(hexgrid, green_spdf, 
                       fn = function(x) sum(!is.na(x)))
    #total number of cells within each hex cell
    gr_tc <- sp::over(hexgrid, ogr_spdf, 
                      fn = function(x) sum(!is.na(x)))
    
    nc <- length(hexgrid_cells)
    out_df$year[counter:(counter + nc - 1)] <- rep(years[i], (nc))
    out_df$cell[counter:(counter + nc - 1)] <- hexgrid_cells
    out_df$gr_mn[counter:(counter + nc - 1)] <- round(jday, 2)
    out_df$gr_sd[counter:(counter + nc - 1)] <- round(gr_sd[,1], 2)
    out_df$gr_ncell[counter:(counter + nc - 1)] <- gr_nc[,1]
    out_df$gr_pcell[counter:(counter + nc - 1)] <- round(gr_nc[,1] / gr_tc[,1], 2)
    
    counter <- (counter + nc)
    
    rm(gr)
    rm(gr_f2)
    rm(gr_mn)
    rm(gr_sd)
    rm(jday)
    rm(gr_nc)
    rm(gr_tc)
    rm(green_spdf)
    rm(ogr_spdf)
    gc()
  }
  return(out_df)
}


# run function ------------------------------------------------------------

#filter by forest cover
# gr_out <- gr_pro_fun(gr_type = 'Greenup', forest = TRUE)
mid_gr_out <- gr_pro_fun(gr_type = 'MidGreenup', forest = TRUE)

#no land cover type filtering
# gr_out_nf <- gr_pro_fun(gr_type = 'Greenup', forest = FALSE)
mid_gr_out_nf <- gr_pro_fun(gr_type = 'MidGreenup', forest = FALSE)



# add lat/lon -------------------------------------------------------------

hexgrid6 <- dggridR::dgconstruct(res = 6)
cellcenters <- dggridR::dgSEQNUM_to_GEO(hexgrid6, mid_gr_out$cell)

# gr_out$cell_lat <- cellcenters$lat_deg
# gr_out$cell_lng <- cellcenters$lon_deg
mid_gr_out$cell_lat <- cellcenters$lat_deg
mid_gr_out$cell_lng <- cellcenters$lon_deg

# gr_out_nf$cell_lat <- cellcenters$lat_deg
# gr_out_nf$cell_lng <- cellcenters$lon_deg
mid_gr_out_nf$cell_lat <- cellcenters$lat_deg
mid_gr_out_nf$cell_lng <- cellcenters$lon_deg

proc.time() - tt


# save to RDS -------------------------------------------------------------

setwd(paste0(dir, 'Data/environment/processed/', run_date))
# saveRDS(gr_out, paste0('Greenup-', run_date, '-forest.rds'))
saveRDS(mid_gr_out, paste0('MidGreenup-', run_date, '-forest.rds'))
# saveRDS(gr_out_nf, paste0('Greenup-', run_date, '-all.rds'))
saveRDS(mid_gr_out_nf, paste0('MidGreenup-', run_date, '-all.rds'))

#copy script to dir
system(paste0('cp ', dir, 'Scripts/1-process-gr/1b-greenup.R ', dir, 'Data/environment/processed/', run_date, '/1b-greenup-', run_date, '.R'))


