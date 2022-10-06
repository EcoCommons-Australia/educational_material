

##########################################################################
###      EcoCommons running species distribution models in R           ###
###   Step 1, set up your workspace, download and clean occurrence data ###
##########################################################################
##
# Author details: EcoCommons Platform 
# Contact details: comms@ecocommons.org.au
# This script is the product of EcoCommons platform.
# 
# please cite:
#   
#   EcoCommons Australia 2022. EcoCommons Australia – the platform of choice for ecological and environmental modelling, EcoCommons, Griffith University, Nathan, 
# Queensland. Available: https://data–explorer.app.ecocommons.org.au/ (Accessed: Date[ e.g., January 19, 2022]).  https://doi.org/10.47486/PL108
# 
# Date: March 2022
# 
# Script and data info:
#   The script and data in this file draws on several sources including 
# 
# 1. The ATLAS of Living Australia (ALA), script modified from script written by Jenna Wraith, ALA Training & Outreach Coordinator
# 
# Stevenson M, Westgate M, Newman P (2022). _galah: Atlas of Living Australia (ALA) Data and
# Resources in R_. R package version 1.3.1, <URL: https://CRAN.R-project.org/package=galah>.
# 
# See for vingnettes on how to use the package: https://atlasoflivingaustralia.github.io/galah/
#   
#   2. Global Biodiversity Information Facility (GBIF)
# 
# source of base code which has been modified here: https://www.r-bloggers.com/2021/03/downloading-and-cleaning-gbif-data-with-r/
#   
#   Chamberlain S, Barve V, Mcglinn D, Oldoni D, Desmet P, Geffert L, Ram K (2022). _rgbif:
#   Interface to the Global Biodiversity Information Facility API_. R package version 3.6.0, 
# <URL: https://CRAN.R-project.org/package=rgbif>.

# Setting up your workspace

# First, if you are new to using R, we strongly suggest you visit this website and fo through that material 
# before trying this: https://datacarpentry.org/R-ecology-lesson/
#   
#   You should be able to run all this code in EcoCommons' coding cloud using either a jupyter notebook (R Kernal) or 
# by importing this notebook into R studio as an R Markdown, i.e. https://rmarkdown.rstudio.com/docs/reference/convert_ipynb.html
# 
# Below is code to install and load the needed packages, and to set up your working directories and subfolders

# Install packages if needed, check first to see if they are installed, if they are not run the required line of 
# code for the needed package

# install.packages("galah")
# install.packages("rgbif")
# install.packages("maps")
# install.packages("tidyr")
# install.packages("jpeg")
# install.packages("raster")
# install.packages("rgeos")
# install.packages("sp")


# load all required packages
# Use the install packages line above for any packages not already installed. 
# install.packages("tidyr") # installs the tidyr package
# library(tidyr) #loads that package into the current session

library(galah)
library(rgbif)
library(maps)
library(tidyverse)
library(jpeg)
library(raster)
library(rgeos)
library(sp)
library(knitr)
library(rmarkdown)


#confirm which packages are loaded
(.packages())



# identify your working directory
getwd()

direct<- "/Users/s2992269/Documents/Use_cases"
folder <- "/SDM_in_R"

# Set the directory
setwd(paste0(direct,folder))

# double check your working directory
getwd()

# create a variety of folders in your working directory, if statement ignores this if the folder exists already

raw_data_folder <- paste0(getwd(),"/raw_data")
  if (!file.exists(raw_data_folder)) {
    dir.create(raw_data_folder)}

data_folder     <- paste0(getwd(),"/data")
  if (!file.exists(data_folder)) {
    dir.create(data_folder)}

results_folder  <- paste0(getwd(),"/results")
  if (!file.exists(results_folder)) {
    dir.create(results_folder)}
      
results_folder  <- paste0(getwd(),"/results1")
  if (!file.exists(results_folder)) {
    dir.create(results_folder)}
    
results_folder  <- paste0(getwd(),"/results2")
  if (!file.exists(results_folder)) {
    dir.create(results_folder)}
    
results_folder  <- paste0(getwd(),"/results3")
  if (!file.exists(results_folder)) {
    dir.create(results_folder)}
      
results_folder  <- paste0(getwd(),"/results4")
  if (!file.exists(results_folder)) {
    dir.create(results_folder)}
      
results_folder  <- paste0(getwd(),"/results_brt")
  if (!file.exists(results_folder)) {
    dir.create(results_folder)}
      
results_folder  <- paste0(getwd(),"/results_glm")
  if (!file.exists(results_folder)) {
    dir.create(results_folder)}
    
results_folder  <- paste0(getwd(),"/results_gam")
  if (!file.exists(results_folder)) {
    dir.create(results_folder)}
      
scripts_folder  <- paste0(getwd(),"/scripts")
    if (!file.exists(scripts_folder)) {
        dir.create(scripts_folder)}

raw_data_folder <- paste0(getwd(),"/predictors")
  if (!file.exists(raw_data_folder)) {
    dir.create(raw_data_folder)}

raw_data_folder <- paste0(getwd(),"/predictors_future")
  if (!file.exists(raw_data_folder)) {
    dir.create(raw_data_folder)}

# Set galah_config by adding your email
# galah_config(email = "your-email@email.com") # This email needs to be registered with ALA
# you can register here: https://auth.ala.org.au/userdetails/registration/createAccount 

galah_config(email = "r.clemens@griffith.edu.au") # again put the email you registered to ALA with

# select an ATLAS
## These configuration options will only be saved for this session
## Set preserve = TRUE to preserve them for future sessions

show_all_atlases()

# we will be using the Australian ATLAS data

galah_config(atlas="Australia")

# show all the fields or columns of data within in the ALA

show_all_fields()

# find the kinds of data in each field replace the field name of interest in quotes

search_fields("australian states and territories")

search_fields("coordinateUncertaintyinMeters")

search_fields("occurrenceID")

# look at the field values 

search_field_values("cl22")

# in this example we are looking at frogs
search_taxa("Amphibia")

search_taxa(genus = "Limnodynastes")

search_taxa(species = "Limnodynastes peroni")

# download some occurence data

# the galah_call function starts a query, then galah_identify selects the taxa you are 
# interested in, galah_filter selects records for specified columns and atlas_occurrences 
# retrieves the specified occurrence records

# look for records submitted by FrogID

search_field_values("datasetName")  # note you can see more rows, click lower right below

# first check the number of records your query will return, there are 67,000+ records in ALA

galah_call() %>%
  galah_identify("Limnodynastes peroni")%>%
  atlas_counts()

# Then filter those records so only those records from "FrogID" are returned for this species
# & remove records with a coordinate uncertainty greater than 100m, there are 36,000+ records from 
# the FrogID dataset, and with only coordinates < 100m precision
galah_call() %>%
  galah_identify("Limnodynastes peroni")%>%
  galah_filter(datasetName == "FrogID")%>%
  galah_filter(stateProvince == "Queensland")%>%
  galah_filter(coordinateUncertaintyInMeters < 100)%>%
  atlas_counts()

# select the occurrence records, galah_select returns wanted columns
# We could also filter by year, but FrogID data is all pretty recent
# Often you will want to ensure you are not including really old records in your data
# i.e. galah_filter(year >= 2020)

LiPe <- galah_call() %>%
  galah_identify("Limnodynastes peroni")%>%
  galah_filter(datasetName == "FrogID")%>%
  galah_filter(coordinateUncertaintyInMeters < 100)%>%
  galah_filter(stateProvince == "Queensland")%>%
  atlas_occurrences()
  
# get familiar with data - this will return the column names or fields

head(LiPe)

# generate a summary of those data

names(LiPe)

write.csv(LiPe, "raw_data/Limnodynastes peroni.csv")

summary(LiPe)

# drop columns not needed in analyses
# 
# Note you will often want to filter on additional fields not shown here

LiPe<-LiPe[,c("decimalLatitude","decimalLongitude","eventDate","scientificName")]  # This returns columns 2,3,4 & 5, if we had 2:5 before the comma, it would return rows 2,3,4, & 5
# if we had [,c(2,3,5)] we would return rows 2, 3 & 5

# look at the top rows of the new dataset
head(LiPe)

# this is the earliest date, 2017 is fairly recent, and we want to make sure our predictor variables are available for these years

min(LiPe$eventDate)

#This is the most recent date, 2020
max(LiPe$eventDate)

# save your filtered data to another folder, again you would often filter on other fields as well

write.csv(LiPe, "data/Limnodynastes peroni.csv")

# double check the directory within the chunks is correct
# note the directory outside the code chunks might be different

getwd()

# read in data and overwrite LiPe
LiPe <- read.csv(paste0(getwd(),"/data/Limnodynastes peroni.csv"))

# FrogID records are pretty clean, but often when you map data you will see odd locations
map("world", xlim = range(LiPe$decimalLongitude),
    ylim = range(LiPe$decimalLatitude))  
points(LiPe[ , c("decimalLongitude", "decimalLatitude")], pch = ".")
# 
# When you look at these records mapped, you can see that most records are concentrated where people are near the bigger cities of eastern Australiathis kind of sampling bias is common in Australian occurrence data.  It takes so much more effort to sample in remote locations
# There are a few things you can do to reduce the problem of bias, you can break up your study area into areas near cities, and more remote areas 
# (we will not do that here)
# You can reduce the number of records close to one another.  This reduces spatial autocorrelation, and is a good step for these kinds of records
# before we show how to do spatial - thinning (reducing how close records are together) we need a base raster layer
# 
# In all the modelling we will show here, each environmental predictor needs to cover the same area (have the same extent), have grid cells that are the
# same size (same resolution), and use the same coordinate system to define a method for turning a spherical planet earth into a flat map (same coordinate system & map projection = Coordinate Reference System - CRS)
# 
# A base layer is one with a CRS, extent and resolution that all other layers will be turned into.  Generally it is a good idea to use the largest possible extent, i.e. you often get better predictions when you include the entire range of your species, and when you use the finest (smallest grid cells) possible.  Smaller grid cells only help prediction when the values vary from one grid cell to the next grid cell. I recommed choosing the variable that is most closely related to the ecology of your species and has the smallest resolution (grid cell size).  Which variable do you think will predict best where your species are found? 

# this is Enhanced Vegetation Index (EVI) which captures greenness index while correcting ndvi issues, from 2004 downloaded from the EcoCommons site - 
# https://api.data-ingester.app.ecocommons.org.au/api/data/34ab5ea6-650f-503f-9446-88d0ae9effe1/download/ndlc-2004-250m_trend-evi-mean.tif

EVI<-raster("raw_data/ndlc-2004-250m_trend-evi-mean.tif")

plot(EVI)

EVI

LiPe_pts <- SpatialPoints(coords = cbind(LiPe$decimalLongitude, LiPe$decimalLatitude),CRS(as.character("+proj=longlat +datum=WGS84 +no_defs")) )

require(rgeos)
LiPe_convH<- rgeos::gConvexHull(LiPe_pts)

# check if your convex hull looks correct by plotting
# useful information on working with spatial data in R https://cengel.github.io/R-spatial/intro.html 

plot(EVI)
lines(LiPe_convH)

#increase the size of your extent (study area to beyone the extent of your point data)
# our CRS is in decimal degrees so 1 degree ~ 111km larger

LiPe_extent <- rgeos::gBuffer(LiPe_convH, width = 1)

# check if your convex hull or extent is now larger

plot(EVI)
lines(LiPe_extent)

EVI_LiPe<- raster::crop(EVI,LiPe_extent)
plot(EVI_LiPe)
EVI_LiPe

EVI_LiPe_mask <- raster::mask(EVI_LiPe,LiPe_extent)
plot(EVI_LiPe_mask)

# create your reference raster, all other rasters will be set to this resolution, extent and CRS, dividing raster by itself sets all values to 1

base <- EVI_LiPe_mask/EVI_LiPe_mask

crs(base)<-crs(EVI)

plot(base)

writeRaster(base,"data/base_LiPe.asc", overwrite=TRUE)

getwd()

# run these lines if you are coming back to the script to load your base layer
require(raster)
base <- raster("data/base_LiPe.asc")
#check that it looks correct
plot(base)

### spatially thin occurrence records so that only one record is selected randomly from within a Grid cell that is 4 times larger than "base"

#aggregate reduce resolution (make grid cells larger) (factor = 4) notice grid cells are now 4 times larger
large_base <- aggregate(base, fact=4)
res(base)
res(large_base)

LiPe <- read.csv("data/Limnodynastes peroni.csv")
LiPe_pts <- sp::SpatialPoints(coords = cbind(LiPe$decimalLongitude, LiPe$decimalLatitude),CRS(as.character("+proj=longlat +datum=WGS84 +no_defs")) )

cell_no<- raster::extract(large_base,LiPe_pts,cellnumbers=TRUE)
head(cell_no)

LiPe_cells<- cbind(LiPe,cell_no)

head(LiPe_cells)

require(dplyr)
LiPe_thinned <- LiPe_cells %>% 
  group_by(cells) %>% 
  slice_sample(n = 1)

length(LiPe_thinned$cells)

# plot to look at results & write results to file

LiPe_pts2 <- SpatialPoints(coords = cbind(LiPe_thinned$decimalLongitude, LiPe_thinned$decimalLatitude),CRS(as.character("+proj=longlat +datum=WGS84 +no_defs")) )
plot(base)
points(LiPe_pts2,pch=20,cex=0.2)

write.csv(LiPe_thinned,"data/LiPe_thinned.csv")

# download all FrogID survey records to capture survey effort, spatial sampling bias, and to zero-fill

galah_call() %>%
  galah_identify("Amphibia")%>%
  galah_filter(datasetName == "FrogID")%>%
  galah_filter(coordinateUncertaintyInMeters < 100)%>%
  galah_filter(stateProvince == "Queensland")%>%
  galah_select("datasetName")%>%
  atlas_counts()

frogs <- galah_call() %>%
  galah_identify("Amphibia")%>%
  galah_filter(datasetName == "FrogID")%>%
  galah_filter(coordinateUncertaintyInMeters < 100)%>%
  galah_filter(stateProvince == "Queensland")%>%
  galah_select(datasetName,occurrenceID)%>%
  atlas_occurrences()

head(frogs)

length(frogs$occurrenceID)

frogs$unique_visit<- paste0(frogs$decimalLatitude,frogs$decimalLongitude,frogs$eventDate)

length(unique(frogs$unique_visit))

length(unique(frogs$eventDate))

frogs$visitID <- as.numeric(as.factor(frogs$unique_visit))

length(unique(frogs$visitID))

head(frogs$visitID, 100)

head(frogs)

names(frogs)

frogs2<-frogs[,c(1,3:6,11)]

head(frogs2)

getwd()

write.csv(frogs2,"raw_data/FrogID_all_ALA_Mar2022.csv")

require(dplyr)
names(frogs2)

#Richness on each visit
frogs3 <- frogs2 %>%
  group_by(decimalLatitude, decimalLongitude, visitID) %>%
  summarise(no_spp = length(unique(scientificName)))

head(frogs3)

length(frogs3$decimalLatitude)

summary(frogs3$no_spp)

# calculate the number of visits to the same lat long location
frogs4 <- frogs3 %>%
  group_by(decimalLatitude, decimalLongitude) %>%
  summarise(no_visits = length(visitID))

length(frogs4$decimalLatitude)

summary(frogs4$no_visits)
# note there were 448 visits to one lat long location, but most locations only had one visit

require(sp)
visits_pts<-SpatialPoints(coords = cbind(frogs4$decimalLongitude, frogs4$decimalLatitude),CRS(as.character("+proj=longlat +datum=WGS84 +no_defs")))

#woops we need to include only those points within our study area
plot(base)
points(visits_pts,pch=20,cex=0.2)

# This creates a raster that totals the total number of visits at each 
# Lat / long location
b1 <- rasterize(visits_pts, base, frogs4$no_visits, fun=sum, background=0)
b1

# This is one way to generate a bias file, simply make a layer where the value in each cell = the number of surveys 
# done in that cell. Maxent will not work with values of zero in the bias file, so below we add 1 to each cell value 
# in our study area (all values in base = 1). Here you can see that there are too few points to show up on a map.  
#If unable to add the bias layer directly into Maxent (args = c("biaslayer = bias")) - this method does not work for 
# me and unfortunately means that when using R it is not using the FACTORBIASOUT argument, but a work around is to account 
# for bias within Maxent by using your bias grid to determine the probability of sampling that grid cell to generate your 
# background.  Here most of the values in the grid are the same value of 1 which in our case means no sampling was done there.  
# So in this example we have a number of sites with an increased probability of being selected as a background point in Maxent, 
# but most of the study area has the same low probability of being selected as background.  (Some people might choose to use this)

numb_visits_bias <- b1+base

# Below we highlight an alternative method to generate a bias grid. If we assume areas close to areas where surveys have been conducted 
# are more likely to have been surveyed, so there is a general spatial sampling bias, we may want to identify those spatial locations 
# where data was more likely to have come from. you can do this kind of thing with a focal function (see step 2 script),
#  but furthe below we generate a smoothed spatial surface where the probability of an area being sampled is a function of distance and 
# sampling intensity. Using a Kernel density function is one way to do this. 

bb <- bbox(b1)

visit_locations<- raster::extract(b1, visits_pts, cellnumbers=TRUE)

# these are the points that only occur in the study area
visit_locations2 <- as.data.frame(na.omit(visit_locations))

cellID <- unique(visit_locations2$cells)

xy_visits <- xyFromCell(b1,cell = cellID)

lg<-c(bb[1,1],bb[1,2])
lt<-c(bb[2,1],bb[2,2])
coordss<-cbind(lg,lt)
v_xy2<-rbind(xy_visits,coordss)

# For some of these functions kde2d, and xyFromCell I had to increase my local R-environ memory size
# Step 1: Open terminal,
# 
# Step 2:
# 
# cd ~
# touch .Renviron
# open .Renviron
# Step 3: Save the following as the first line of .Renviron:
# 
# R_MAX_VSIZE=100Gb 

require(MASS)

dens <- kde2d(v_xy2[,1], v_xy2[,2], n = c(nrow(b1), ncol(b1)))
b5 <- raster(dens)
plot(b5)

b6<-projectRaster(b5,base)

b7<-mask(b6,base)
plot(b7)

writeRaster(b7,"data/Bias_LiPe_kd.asc",overwrite=TRUE)

# Now lets clean up our global environment, and get rid some of these big files

rm(b5)
rm(b6)
rm(EVI)
rm(EVI_LiPe)
rm(EVI_LiPe_mask)
rm(frogs2)
rm(frogs3)
rm(visit_locations)

# Zero filling
# In order to identify locations that were surveyed but where LiPe was not detected.
# 
# The FrogID project attempts to identify all frog species that were heard during a period of sound 
# recording. While there is a possibility that frogs were present but not detected it is correct to 
# assume that if the target species was not recorded on a certain date and time at a specific lat/long 
# the there was a zero detection.  You can then make some assumptions around the number of surveys needed 
# to be sure there were no frogs detected, while understanding that the spatial resolution of the modelling 
# you are doing may include many frog habitats some of which may not have been surveyed.  Still this method 
# of generating pseudo-absences is at least putting zeros in areas where surveys were done, but the species 
# was not detected. If you relax the assumptions further, then species for which counts are usually inclusive 
# of all the species detected (like bird lists), those cells where 10 (the number of surveys will vary by 
# species dataset etc) surveys have been done but which did not detect the target species can be assumed 
# to have zero of that species.
# 
# The raster layer we created and named b1 has the number of surveys at all the locations where frogs were 
# surveyed.  First, we are going to check how many cells have a value of 3 or higher.  So is there enough 
# data to produce pseudo-absence zeros if we assume that when three visits were done in a grid cell we should 
# have identified our target species (here striped marsh frog).  In a perfect world we may have needed 10 
# surveys at each of the habitats within each of the grid cells in order to be certain that the frog is truly 
# absent from this location. Here we have some evidence that the frog was not present, and we have the B1 raster 
# which captures survey effort which is a useful covariate in occupancy modelling (something we won't go into here).

# reclassify raster with number of surveys to 1 if more than 2 surveys were done, and zero otherwise

m <- c(0, 2.9, 0,  2.9, 3247, 1)
reclass <- matrix(m, ncol= 3, byrow= TRUE)
rc <- reclassify(b1, reclass)

visits_3or_more <- mask(rc,base)
# since all values are 1 or zero and ones are where there were three or more visits, cellStats gives us the number of grid cells where FrogID surveys were conducted (9196) in our case
freq(visits_3or_more)

#create a vector of 1's of the same length as the LiPe_pres
LiPe_pres<- rep(1,length(LiPe$decimalLatitude))

# Then create a raster with a value of 1 for each gridcell where a LiPe was recorded
LiPe_pres_raster<-rasterize(LiPe_pts,base,LiPe_pres, fun = min, background=0)

LiPe_pres_raster2<-mask(LiPe_pres_raster,base)

#How many grid cells of the pre-thinned data had LIPE recorded in them
freq(LiPe_pres_raster2)

# confirm that your two raters are the same extent, CRS, resolution etc
compareRaster(LiPe_pres_raster2,visits_3or_more,extent=TRUE, rowcol=TRUE, crs=TRUE, res=TRUE, orig=TRUE,rotation=TRUE, values=FALSE)

Zero_LiPe <- visits_3or_more - LiPe_pres_raster2
# note there are 3179 locations with a value of -1 indicating the number of locations where LiPe was the 
# only species observed in that grid cell, while 5617 cells have a value of 1, these are cells where frogs 
# other than LiPe were recorded and we will consider these zeros from here on.  This is not enough zeros to 
# use as a targeted background in Maxent (usually best to go with 10,000 background locations), but it is 
# enough for most other methods that require zeros BRT, GLM etc.

freq(Zero_LiPe)

#reclassify so a zero raster has the value of 1 for all LiPe zeros, and NA or 0 for all other values
m <- c(-2, 0.1, 0,  0.1, 2, 1)
reclass2 <- matrix(m, ncol= 3, byrow= TRUE)
rc2 <- reclassify(Zero_LiPe, reclass2)
Zero_LiPe2<-mask(rc2,base)
freq(Zero_LiPe2)

# extract the cell numbers from the 0 grid where the value ==1
cell_vals_0<-Which(Zero_LiPe2 ==1,cells=TRUE)
# these are the lat / longs for locations where at least three surveys were done, but zero LiPe were detected - these are our pseudo absences
xy_zero_LiPe <- xyFromCell(Zero_LiPe2,cell = cell_vals_0)
#plot the absence locations
plot(base)
points(xy_zero_LiPe, pch=20,cex=0.2)

# plot the presence locations
plot(base)
points(LiPe[ , c("decimalLongitude", "decimalLatitude")], pch = ".")

# we will use these pseudo absences for BRT & GLM later

LiPe_zeros<-as.data.frame(xy_zero_LiPe)
colnames(LiPe_zeros)<-c("Long","Lat")
LiPe_zeros$spp<-"Limnodynastes_peroni"
LiPe_zeros$pres<- 0
write.csv(LiPe_zeros,"data/LiPe_zero_locations_3vis.csv")

writeRaster(Zero_LiPe2,"data/LiPe_absences1_3vis.asc",overwrite=TRUE)

zeros_3visit <- read.csv("data/LiPe_zero_locations_3vis.csv")

# reclassify raster with number of surveys to 1 if more than 1 surveys was done, and zero otherwise

m <- c(0, 0.9, 0,  0.9, 3247, 1)
reclass3 <- matrix(m, ncol= 3, byrow= TRUE)
rc3 <- reclassify(b1, reclass3)
freq(rc3)

# It is a good idea to check your work as you go.  Plotting data, looking at freq are some ways to do this.  
# Here we double check that the number of cells with values greater than 1 in our b1 raster match the number 
# of 1's in the reclassification

t1<-as.data.frame(freq(b1))
sum(t1$count)-t1[1,2]

# sure enough freq of rc3 returns 29161 cells with a value of 1, and the sum of cell counts with values 
# greater than 1 in the original b1 raster are also 29161

# Why am I showing this?  Because my first attempt at generating the rc3 raster did not have equal numbers 
# using a different method, and a typo in my first reclassification call kept these numbers separate.  
# It is just a good idea to always verify that what you think you did in your code do not have any errors.

visits_1or_more <- mask(rc3,base)
freq(visits_1or_more)

# confirm that your two raters are the same extent, CRS, resolution etc
compareRaster(LiPe_pres_raster2,visits_1or_more,extent=TRUE, rowcol=TRUE, crs=TRUE, res=TRUE, orig=TRUE,rotation=TRUE, values=FALSE)

Zero_LiPe1visit <- visits_1or_more - LiPe_pres_raster2
freq(Zero_LiPe1visit)
# here we have 22134 grid cells where another frog(s) was recorded but not our target spp.  We only need 
# 10,000 cells (usually) for a targeted background in Maxent, so we can afford to do some spatial thinning. 
# (If we don't have close to 10,000 grid cells after spatial thinning, we will just use the kde bias layer)

cell_vals_0_1 <- Which(Zero_LiPe1visit==1,cells=TRUE)

xy_zero_LiPe_1vis <- as.data.frame(xyFromCell(Zero_LiPe2,cell = cell_vals_0_1))

colnames(xy_zero_LiPe_1vis)<-c("Long","Lat")
xy_zero_LiPe_1vis$spp<-"Limnodynastes_peroni"
xy_zero_LiPe_1vis$pres<- 0
write.csv(xy_zero_LiPe_1vis,"data/LiPe_zero_locations_1vis.csv")

LiPe_pts2 <- SpatialPoints(coords = cbind(xy_zero_LiPe_1vis$Long, xy_zero_LiPe_1vis$Lat),CRS(as.character("+proj=longlat +datum=WGS84 +no_defs")) )
cell_no2<- as.data.frame(raster::extract(large_base,LiPe_pts2,cellnumbers=TRUE))

head(cell_no2)

# notice here I did not need to use a subsequent colnames function

LiPe_cells2<- as.data.frame(cbind(Long = xy_zero_LiPe_1vis$Long, Lat = xy_zero_LiPe_1vis$Lat,cell_no = cell_no2$cells))

head(LiPe_cells2)

require(dplyr)
LiPe_zeros_thinned_1vis <- LiPe_cells2 %>% 
  group_by(cell_no) %>% 
  slice_sample(n = 1)

length(LiPe_zeros_thinned_1vis$cell_no)

# note this in this object there is still more, than 10,000 areas where frog surveys were done, but no LiPe 
# were recorded. We might see how many cells we would have if there were two or more visits, and then thin 
# that result, but we will try to use these and the points from the bias layer to fit our models.  There are 
# trade-offs with every modelling decision, some are important for your result, some do not really impact the 
# result.  Your decisions on what is best for your model depends on your question, the available data, the 
# ecology of your species, and an understanding of what has worked or is important according to the literature 
# for your species and your kind of question.

LiPe_zeros_thinned_1vis$spp<-"Limnodynastes_peroni"
LiPe_zeros_thinned_1vis$pres<- 0
head(LiPe_zeros_thinned_1vis)
write.csv(LiPe_zeros_thinned_1vis,"data/LiPe_zero_locations_thinned_1vis.csv")

# note if you just use head(LiPe_zeros_thinned_1vis), you won't see all the decimal points in Long and Lat, 
# print.data.frame shows all the decimal

print.data.frame(head(LiPe_zeros_thinned_1vis))

# It is often a good idea to break your data into randomly selected training and testing data.  
# Often you would randomly remove 20% or more of your data and withold that from the model building process.  
# Once your model was finsihed using your training data, you would then test your model with this training data. 
# If you have very little data, bootstrapping can be used to see if removal of a small percentage of data 
# repeatedly changes results.  This gives you a good understanding of the confidence intervals around your 
# results and can reduce the impact of outliers on your final result.  Cross-validation requires more data.  
# Ideally for cross-validation you break your data into 10 folds (subsets) and compare results between folds, 
# or summarise variability between folds. Some people select folds spatially, with each fold an independent 
# spatial block of data, but random folds seem to work just as well or better. If you can afford to set aside 
# 20% + of your data, withholding a test data set is good practice.  Ideally, you will test your model with 
# completely independent data.
