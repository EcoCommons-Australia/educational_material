# Environmental data

library(raster)
require(rgdal)

# this is a soil wetness index from 2013 downloaded from the EcoCommons site 
# https://api.data-ingester.app.ecocommons.org.au/api/data/69b045b3-29d1-5068-a5e6-f0783c68c1b8/download/awap_run26j_2013_ann-TempMin.tif
# Simply enter the above url, download the file, and load it from the directory path you saved it
#- ideally you would match the time the layer was made to the time 
# your occurrence data was collected, and you might average over the last ten years, but this is just an example of the kinds of data you can use
# the EcoCommons point-and-click environment does this for you when in the SDM workflow
#you can select to make all the variables the finest or coarsest resolution

direct<- "/Users/s2992269/Documents/Use_cases"
folder <- "/SDM_in_R"

#this sets your working director for all subsequent chunks of code in your R Markdown script
setwd(paste0(direct,folder))

AWAP1<-raster("raw_data/awap_WRel1End.tif")
plot(AWAP1)
AWAP1

base <- raster("data/base_LiPe.asc")
crs(base) <- "+proj=longlat +datum=WGS84 +no_defs"

# reduce the resolution of this layer on soil wetness, and match the extent of base
# note there are a variety of ways to this. In some cases these decisions on how 
# to reduce or increase resolution matter to the result.  Often the defaults are fine.
# know your data, and dig deeper into the methods underlying these functions
#try typing 
# ??raster::projectRaster
AWAP2 <-projectRaster(AWAP1,base)
AWAP3<- mask(AWAP2,base)
plot(AWAP3)

AWAP <- AWAP3

#now we want to be sure the name of the variable AWAP is the same in memory, and in a folder
# full of all predictors so we can predict our fitted model into the entire area
# possibly using other point occurrence data (independent data)
writeRaster(AWAP,"predictors/AWAP.asc")
writeRaster(AWAP,"predictors_future/AWAP.asc")

# Lets now upload the Vegetation classification layer NVIS from EcoCommons
# https://api.data-ingester.app.ecocommons.org.au/api/data/ffcf7885-95e0-5873-8178-4a5ecc14123b/download/nvis-2020-90m_aus6_0e_mvg_amvg.tif

NVIS1<-raster("raw_data/nvis-2020-90m_aus6_0e_mvg_amvg.tif")
# each grid cell value should correspond to one of these categegories
# https://www.awe.gov.au/sites/default/files/env/pages/ba1d4b30-d46f-42f7-bec2-fac391f26072/files/mvsg60-sort-order.pdf
freq(NVIS1)
freq(NVIS1,value=44)
# habitat of frog https://environment.des.qld.gov.au/wildlife/animals/a-z/striped-marsh-frog
# it is hard to see one of these classes being super helpful, class 44 freshwater did not return anything
# if we did use it we need to be careful changing the resolution because these are categorical variables 
# (use nearest neighbor when changing resolution)

# Randome useful tidbit
# A handy function to clip points to a polygon
# clipped points <- spatial_point_file[polygon_file, ]

# get a wetland file
# these steps take too long to run here, so we supply the wetland layer, steps are below
# download an Australia wetland shapefile from here
# https://ecat.ga.gov.au/geonetwork/srv/eng/catalog.search#/metadata/83135
# select a 0.1 degree grid that covers Australia
# read in Wetland polygones
# wetland <- st_read("SurfaceHydrologyPolygonsNational.gdb")
# then turn polygons into raster
# the argument getCover in the rasterize function calculates the area of each cell covered by wetland
# wetland <- rasterize(wetland,Point1_degree_grid,getCover = TRUE)
# for this variable it would have made more sense to source a polygon layer of freshwater only wetland

setwd("~/Documents/Use_cases/SDM_in_R/Raw_data")
wetland1 <- raster('wet_cov.asc')
crs(wetland1) <- "+proj=longlat +datum=WGS84 +no_defs"
wetland2 <-projectRaster(wetland1,base)
wetland3<- mask(wetland2,base)
plot(wetland3)

# Here we use the focal function (a neighborhood function) to transform each grid cell
# value to the sum of the central cell and all cells 5 cells away surrounding that central cell
# this is a surrogate of connectivity, isolated wetlands will have lower values

wetland_connectivity <- focal(wetland3, w=matrix(1, nrow = 11, ncol = 11),fun=sum,na.rm=TRUE)
plot(wetland_connectivity)
setwd("~/Documents/Use_cases/SDM_in_R")

writeRaster(wetland_connectivity,"predictors/wetland_connectivity.asc")
writeRaster(wetland_connectivity,"predictors_future/wetland_connectivity.asc")

# data comes in many formats, and uploading often requires different mehtods
# NetSDF files are increasingly common, below is a sript to bring in that file type

#require(sf)
#require(ncdf4)
#require(rasterVis)
#require(raster)

# MXtemp <- brick("Terraclim_EY_NSW.nc", varname="tmax") # this NetCDF file includes many 
#days of maximum temperature data with that subset of data having the varname = "tmax"
# MXtemp_mean <- mean(MXtemp) # if we just want one layer which is the mean of those daily totals

# we will now download a vegetation greenness index from EcoCommons
#https://api.data-ingester.app.ecocommons.org.au/api/data/34ab5ea6-650f-503f-9446-88d0ae9effe1/download/ndlc-2004-250m_trend-evi-mean.tif

EVI1<- raster("https://api.data-ingester.app.ecocommons.org.au/api/data/34ab5ea6-650f-503f-9446-88d0ae9effe1/download/ndlc-2004-250m_trend-evi-mean.tif")

# process it as previous layers

EVI2 <-projectRaster(EVI1,base)
EVI3<- mask(EVI2,base)
plot(EVI3)
EVI <- EVI3
writeRaster(EVI,"predictors/EVI.asc")
writeRaster(EVI,"predictors_future/EVI.asc")

# Here we just show how to bring in current and future data from EcoCommons

# after looking at correlations between bioclim variables this one was dropped
#Bioclim01_1 <- raster("https://api.data-ingester.app.ecocommons.org.au/api/data/90317596-ddef-5666-91c5-9cbc25c24fbc/download/current_1976-2005_bioclim-01.tif")
#Bioclim01_2 <-projectRaster(Bioclim01_1,base)
#Bioclim01_3 <- mask(Bioclim01_2 ,base)
#plot(Bioclim01_3)
#Bioclim01 <- Bioclim01_3
#writeRaster(Bioclim01,"predictors/Bioclim01.asc")

# we will not write the bioclim data to the predictions_future folder because
# we have separate future climate predictions we can use.  We do not have future
# predictions for EVI or Wetlands so we are assuming those things will stay the same
# in the future, the names of the future variables need to be the same as the
# current time variable names in order for the predict function to work, so just
# be sure to keep the different variables with the same name in different folders.

Bioclim05_1 <- raster("https://api.data-ingester.app.ecocommons.org.au/api/data/cc081daa-f524-58c2-939e-166a2b2e79eb/download/current_1976-2005_bioclim-05.tif")
Bioclim05_2 <-projectRaster(Bioclim05_1,base)
Bioclim05_3 <- mask(Bioclim05_2 ,base)
plot(Bioclim05_3)
Bioclim05 <- Bioclim05_3
writeRaster(Bioclim05,"predictors/Bioclim05.asc")

Bioclim06_1 <- raster("https://api.data-ingester.app.ecocommons.org.au/api/data/476e4343-99f2-578e-b44a-951a55c6b7c2/download/current_1976-2005_bioclim-06.tif")
Bioclim06_2 <-projectRaster(Bioclim06_1,base)
Bioclim06_3 <- mask(Bioclim06_2 ,base)
plot(Bioclim06_3)
Bioclim06 <- Bioclim06_3
writeRaster(Bioclim06,"predictors/Bioclim06.asc")

Bioclim12_1 <- raster("https://api.data-ingester.app.ecocommons.org.au/api/data/b2d70413-6b74-5366-b5d5-82ed919ded93/download/current_1976-2005_bioclim-12.tif")
Bioclim12_2 <-projectRaster(Bioclim12_1,base)
Bioclim12_3 <- mask(Bioclim12_2 ,base)
plot(Bioclim12_3)
Bioclim12 <- Bioclim12_3
writeRaster(Bioclim12,"predictors/Bioclim12.asc")

rm(Bioclim05_1,Bioclim05_2,Bioclim05_3,Bioclim06_1,Bioclim06_2,Bioclim06_3,Bioclim12_1,Bioclim12_2,Bioclim12_3)

# after looking at correlations between bioclim variables this one was dropped
#Bioclim14_1 <- raster("https://api.data-ingester.app.ecocommons.org.au/api/data/e64dfa1d-ece8-50e7-9ed1-d603bfb7a101/download/current_1976-2005_bioclim-14.tif")
#Bioclim14_2 <-projectRaster(Bioclim14_1,base)
#Bioclim14_3 <- mask(Bioclim14_2 ,base)
#plot(Bioclim14_3)
#Bioclim14 <- Bioclim14_3
#writeRaster(Bioclim14,"predictors/Bioclim14.asc")

#create a raster stack of current data, be careful not to include the future variables below in this stack (they have the same names)

preds_current <- stack(AWAP,wetland_connectivity,EVI,Bioclim05,Bioclim06,Bioclim12)
plot(preds_current)
names(preds_current)
preds_current2 <- setNames(preds_current,c("AWAP","wetland_connectivity","EVI","Bioclim05","Bioclim06","Bioclim12"))
names(preds_current2)

rm(preds_current)

#Repeat with future climate data, from EcoCommons
# Australia, Climate Projection, SRESA1B based on INM-CM30, 30 arcsec (~1km) - 2085

# after looking at correlations between bioclim variables this one was dropped
#Bioclim01_1 <- raster("https://api.data-ingester.app.ecocommons.org.au/api/data/4e5534eb-c743-52af-a9c3-d137df8ba8c1/download/SRESA1B_inm-cm30_2085_bioclim-01.tif")
#Bioclim01_2 <-projectRaster(Bioclim01_1,base)
#Bioclim01_3 <- mask(Bioclim01_2 ,base)
#plot(Bioclim01_3)
#Bioclim01 <- Bioclim01_3
#writeRaster(Bioclim01,"predictors_future/Bioclim01.asc")

Bioclim05_1 <- raster("https://api.data-ingester.app.ecocommons.org.au/api/data/34fe63e3-4152-5406-96f7-dd67f54be5f5/download/SRESA1B_inm-cm30_2085_bioclim-05.tif")
Bioclim05_2 <-projectRaster(Bioclim05_1,base)
Bioclim05_3 <- mask(Bioclim05_2,base)
plot(Bioclim05_3)
Bioclim05 <- Bioclim05_3
writeRaster(Bioclim05,"predictors_future/Bioclim05.asc")

Bioclim06_1 <- raster("https://api.data-ingester.app.ecocommons.org.au/api/data/36cf3700-d238-511c-bb8b-bc3c3677309d/download/SRESA1B_inm-cm30_2085_bioclim-06.tif")
Bioclim06_2 <-projectRaster(Bioclim06_1,base)
Bioclim06_3 <- mask(Bioclim06_2,base)
plot(Bioclim06_3)
Bioclim06 <- Bioclim06_3
writeRaster(Bioclim06,"predictors_future/Bioclim06.asc")

Bioclim12_1 <- raster("https://api.data-ingester.app.ecocommons.org.au/api/data/449cee85-5f6a-5b25-aca8-0fd2838be09f/download/SRESA1B_inm-cm30_2085_bioclim-12.tif")
Bioclim12_2 <-projectRaster(Bioclim12_1,base)
Bioclim12_3 <- mask(Bioclim12_2,base)
plot(Bioclim12_3)
Bioclim12 <- Bioclim12_3
writeRaster(Bioclim12,"predictors_future/Bioclim12.asc")

rm(Bioclim05_1,Bioclim05_2,Bioclim05_3,Bioclim06_1,Bioclim06_2,Bioclim06_3,Bioclim12_1,Bioclim12_2,Bioclim12_3)

# after looking at correlations between bioclim variables this one was dropped
#Bioclim14_1 <- raster("https://api.data-ingester.app.ecocommons.org.au/api/data/178ee1bc-0aa1-58a0-abaf-ff5e19fcf363/download/SRESA1B_inm-cm30_2085_bioclim-14.tif")
#Bioclim14_2 <-projectRaster(Bioclim14_1,base)
#Bioclim14_3 <- mask(Bioclim14_2 ,base)
#plot(Bioclim14_3)
#Bioclim14 <- Bioclim14_3
#writeRaster(Bioclim14,"predictors_future/Bioclim14.asc")

#again just be sure to run this in order so you are not getting current bioclim data mixed with the future data which has the same name

preds_future <- stack(AWAP,wetland_connectivity,EVI,Bioclim05,Bioclim06,Bioclim12)
names(preds_future)
preds_future2 <- setNames(preds_current,c("AWAP","wetland_connectivity","EVI","Bioclim05","Bioclim06","Bioclim12"))
names(preds_future2)
preds_current2

rm(preds_future)

# here we will show you how to take a mean of many months of data if you have many months of data in one folder (you may wa#
# perhaps you are exploring what a suitable niche looks like during a drought
# You may then want climate like averages like precipitation, or average temperature during the drought months

# one way to do this if you have separate rasters fir each month in one folder
#raster_list <- list.files(pattern='.tif$', all.files=TRUE)  # if there are 30 files of tifs from January, this creates a list of those files
#all_chl <- stack(raster_list) # this reads each of those files in, and makes them into a raster stack
#mean_chl <- mean(all_chl) # and this produces one raster where each cell value is the mean value of all the rasters that were in the that folder

# NETCDF files also are a common way for large volumes of data to be
#this is an example loop using Terraclimate data
# read NetCDF data file obtained from "TerraClimate"

# import netCDF file - data from TerraClimate python download

# require(sf)
# require(ncdf4)
# require(rasterVis)
# require(raster)

# download the file here: https://drive.google.com/file/d/1d4aCrdwjWRgENFWqkOucwEZhyDsXpnSV/view?usp=sharing 
# setwd("~/Documents/Use_cases/EY_frogs/data")  # this is the directory where this NETCDF file is stored
# MXtemp <- brick("Terraclim_EY_NSW.nc", varname="tmax")
# Mxtemp_mean <- mean(MXtemp)
# plot(Mxtemp_mean)
# writeRaster(Mxtemp_mean,"MXtemp_TERRA_Sydney_region.asc", overwrite=TRUE)

# Rain <- brick("Terraclim_EY_NSW.nc", varname="ppt")
# Rain_mean <- mean(Rain)
# plot(Rain_mean)
# writeRaster(Rain_mean,"Rain_TERRA_Sydney_region.asc", overwrite=TRUE)

# MNtemp <- brick("Terraclim_EY_NSW.nc", varname="tmin")
# MNtemp_mean <- mean(MNtemp)
# plot(MNtemp_mean)
# writeRaster(MNtemp_mean,"MNtemp_TERRA_Sydney_region.asc", overwrite=TRUE)

# Soil <- brick("Terraclim_EY_NSW.nc", varname="soil")
# Soil_mean <- mean(Soil)
# plot(Soil_mean)
# writeRaster(Soil_mean,"Soil_wetness_TERRA_Sydney_region.asc", overwrite=TRUE)
# crs(Soil_mean)
