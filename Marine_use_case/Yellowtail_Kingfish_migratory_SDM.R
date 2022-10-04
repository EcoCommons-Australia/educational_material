#  Prepare Marine data  Yellowtail Kingfish example

require(galah)
require(terra)
require(tidyverse)
require(googledrive)
require(purrr)

# we recommend going to this google link and downloading the whole
# of the "Marine folder" and put it in a local directory
#  https://drive.google.com/drive/folders/19YemtBtbZbdtwDYcorKOyzlBCNHvRv9h?usp=sharing

setwd("~/Downloads/Marine") # set this to your local directory

# You can download records directly from the Atlas of Living Australia and this
# includes IMOS records. The exact data we used to generate the use case is available through the google line:41 below.
# Keep in mind with stochastic learning algorithms like Maxent, results may vary somewhat from run to run, even
# if using the same data.
# Again, we recommend downloading directly from a large database each time you start a new model, but here
# we got some additional records from AODN which are not in the ALA data, so use the google link (line 41) to replicate our 
# results. Note our ALA download includes WA near-coastal records, and does not include SA near-coastal data which
# were in the cleaned IMOS AODN data.

# ALA download example
# There are a few steps you need to take to download records from ALA
# First, you need to register with ALA
# you can register here: https://auth.ala.org.au/userdetails/registration/createAccount 
# then run the config function in R
# galah_config(email = "your-email@email.com") # This email needs to be registered with ALA

galah_config(email = "r.clemens@griffith.edu.au") # replace my email with the email you registered to ALA with

YT <- galah_call() %>%
  galah_identify("Seriola lalandi")%>%
  galah_filter(basisOfRecord == "MACHINE_OBSERVATION" | basisOfRecord == "HUMAN_OBSERVATION")%>%
  galah_filter(coordinateUncertaintyInMeters < 501)%>%
  galah_filter(year > 2011)%>%
  galah_select("month")%>%
  atlas_occurrences()

# use atlas_counts instead of atlas_occurrences() if you want just a count of # of records in your query

## to use the same occurrence data as we did in the marine use case, download the data from Google, link below.
## use direct links to get a google drive file url which will download:  https://sites.google.com/site/gdocs2direct/
YT <- read.csv("https://drive.google.com/uc?export=download&id=1CRv5ZWIS-jWwnz9VH4kng1A_Ts053wA_", stringsAsFactors = FALSE)
# when checking the data, there is a negative latitude which plots to an odd place, so we will filter out that data.
YT_2 <- YT[YT$lon >0,]

YT_pts <- terra::vect(YT_2, geom=c("lon", "lat")) # this turns the dataframe into a spatial points layer with attributes

# look at your data
YT_pts
plot(YT_pts)

# Now we are going to upload some bathemetry data
## bath
bath <- rast("https://drive.google.com/uc?export=download&id=1uSOShdrsuz4OBmzGu0yvyP7XjRlrQnfS")
crs(bath)
plot(bath)
points(YT_pts)

# here we use a focal function on the edge of the raster to estimate bathemetry values at near-coastal areas
# where bathymetry data is missing (na.policy="only",na.rm=T), by averaging values within 7 grid cells
wt <- matrix (data = 1, nrow = 7, ncol = 7)
bath2 <- terra::focal(bath, w=wt, fun="mean", na.policy="only", na.rm=T) 
plot(bath2)
# verify that we now have bathemetry data at all the sensor locations, including those in shallow water
YT_bath <- terra::extract(bath2,YT[,c(4,3)])
summary(YT_bath)

dir.create("data1")
writeRaster(bath2,paste0(getwd(),"/",("data1/bathemetry.tif")))


#####
# note we dropped ucur variable after initial analysis (it added little to the result)
# we also dropped the distance from the coast & bathymetry (both resulted in underprediction far from the coast)
# so we show how to do it here, but we did not generate monthly 

####################################
###################################
### example of how to geneerate a BIAS layer
## AODN sensors are deployed at shallow areas around the Australian continent
## The sensors detect any tagged fish that swims by a sensor

stations <- read.csv("IMOS_array_summarised.csv")

st_xy <- cbind(stations$lon,stations$lat)
colnames(st_xy)<- c("Long","Lat")

days <- stations$days_deployed

b1 <- terra::rasterize(st_xy,bath2,days,fun=sum,background=1)

plot(b1)
b1
summary(b1)

b2 <- terra::project(b1,bath2)
plot(b2)
b3<-mask(b2,bath2)
plot(b3)  # if there had been enough sampling in enouch cells this could be the bias layer
terra::ext(b3)
lg<-c(109.2333,163.2)
lt<-c(-47.2,-8.875)
coordss<-cbind(lg,lt)
st_xy2<-rbind(st_xy,coordss) # this step ensures we have coordinates
#  needed to calculate a kernal density for the entire extent

require(MASS)

dens <- kde2d(st_xy2[,1], st_xy2[,2], n = c(nrow(b2), ncol(b2)))
b4<-raster::raster(dens)
plot(b4)
b4a<-terra::rast(b4)
b5 <- terra::project(b4a,bath2)
plot(b5)

b6 <- mask(b5,bath2)
plot(b6)

dir.create("bias")
writeRaster(b6,paste0(getwd(),"/",("bias/Bias_AODN.tif")))

#######################
#######################
#######################
#######################
#######################
#######################
#######################
#######################
## Generate monthly data  

#chl

dir.create("monthly_predictors")

setwd("~/Downloads/Marine/raw_data/chl")

t1<-list.files(recursive=TRUE)
#get only unique charater strings near end of string in new list
t2<-unique(str_sub(t1,-6,-5))

for (i in 1:length(unique(t2))){
  txt_files=list.files(pattern=paste("*\\-",unique(t2)[i],'.*\\.tif',sep=""),recursive=TRUE)
  stack<-terra::rast(txt_files)
  filename<-"chl"
  filen<-mean(stack)
  mean_chl_1 <-project(filen,bath2)
  wt2 <- matrix (data = 1, nrow = 33, ncol = 33)
  mean_chl_2 <- focal(mean_chl_1, wt2, fun = mean, na.rm = TRUE, NAonly = TRUE)
  mean_chl_3 <- raster::mask(mean_chl_2, bath2)
  writeRaster(mean_chl_3,paste0("~/Downloads/Marine/monthly_predictors/",unique(t2)[i],"/chl.tif"),overwrite=TRUE)
}

## Generate monthly data  gsla

setwd("~/Downloads/Marine/raw_data/gsla")

t1<-list.files(recursive=TRUE)
#get only unique charater strings near end of string in new list
t2<-unique(str_sub(t1,-6,-5))

for (i in 1:length(unique(t2))){
  txt_files=list.files(pattern=paste("*\\-",unique(t2)[i],'.*\\.tif',sep=""),recursive=TRUE)
  stack<-terra::rast(txt_files)
  filename<-"gsla"
  filen<-mean(stack)
  mean_gsla_1 <-project(filen,bath2)
  wt2 <- matrix (data = 1, nrow = 15, ncol = 15)
  mean_gsla_2 <- focal(mean_gsla_1, wt2, fun = mean, na.rm = TRUE, NAonly = TRUE)
  mean_gsla_3 <- raster::mask(mean_gsla_2, bath2)
  writeRaster(mean_gsla_3,paste0("~/Downloads/Marine/monthly_predictors/",unique(t2)[i],"/gsla.tif"),overwrite=TRUE)
}



## Generate monthly data  sst

setwd("~/Downloads/Marine/raw_data/sst")

t1<-list.files(recursive=TRUE)
#get only unique charater strings near end of string in new list
t2<-unique(str_sub(t1,-6,-5))

for (i in 1:length(unique(t2))){
  txt_files=list.files(pattern=paste("*\\-",unique(t2)[i],'.*\\.tif',sep=""),recursive=TRUE)
  stack<-rast(txt_files)
  all_sst_2 <-project(stack,bath2)
  all_sst_3 <- raster::mask(all_sst_2, bath2)
  filen<- mean(all_sst_3)
  wt2 <- matrix (data = 1, nrow = 35, ncol = 35)
  mean_sst_2 <- focal(filen, wt2, fun = mean, na.rm = TRUE, NAonly = TRUE)
  mean_sst_3 <- raster::mask(mean_sst_2, bath2)
  writeRaster(mean_sst_3,paste0("~/Downloads/Marine/monthly_predictors/",unique(t2)[i],"/sst.tif"),overwrite=TRUE)
}


## bath & dist

# bath<-raster("~/Documents/Use_cases/Marine/output/bath/bath.tif")
# plot(bath)
# dist<-raster("~/Documents/Use_cases/Marine/output/dist_to_land/dist.tif")
# plot(dist)

# wt <- matrix (data = 1, nrow = 7, ncol = 7)
# bath2 <- focal(bath, wt, fun = mean, na.rm = TRUE, NAonly = TRUE)

# dist2 <- projectRaster(dist,bath2)
# wt <- matrix (data = 1, nrow = 7, ncol = 7)
# dist3 <- focal(dist2, wt, fun = mean, na.rm = TRUE, NAonly = TRUE)
# dist4 <- raster::mask(dist3, bath2)

# note we dropped distance to coast and bathymetery after initial analysis
# so we show how to do it here, but don't run the four lines
# below if you want to replicate our results

#for (i in 1:length(unique(t2))){
#  writeRaster(dist4,paste0("~/Documents/Use_cases/Marine/Yellowtail_R/data/monthly_predictors/",unique(t2)[i],"/dist.asc"))
#  writeRaster(bath2,paste0("~/Documents/Use_cases/Marine/Yellowtail_R/data/monthly_predictors/",unique(t2)[i],"/bath.asc"))
#}


## vcur

setwd("~/Downloads/Marine/raw_data/vcur")

t1<-list.files(recursive=TRUE)
#get only unique charater strings near end of string in new list
t2<-unique(str_sub(t1,-6,-5))

for (i in 1:length(unique(t2))){
  txt_files=list.files(pattern=paste("*\\-",unique(t2)[i],'.*\\.tif',sep=""),recursive=TRUE)
  stack<-rast(txt_files)
  filename<-"vcur"
  filen<-mean(stack)
  mean_vcur_1 <-project(filen,bath2)
  wt2 <- matrix (data = 1, nrow = 15, ncol = 15)
  mean_vcur_2 <- focal(mean_vcur_1, wt2, fun = mean, na.rm = TRUE, NAonly = TRUE)
  mean_vcur_3 <- raster::mask(mean_vcur_2, bath2)
  writeRaster(mean_vcur_3,paste0("~/Downloads/Marine/monthly_predictors/",unique(t2)[i],"/vcur.tif"))
}


# here we generate a unique raster for each month of the year
# with the values corresponding to the numeric month value

# Create a base layer where each value is equal to one
base <- bath2/bath2

for (i in 1:length(unique(t2))){
  month<- as.numeric(unique(t2)[i])
  mo<- base*month
  writeRaster(mo,paste0("~/Downloads/Marine/monthly_predictors/",unique(t2)[i],"/mo.tif"))
}

########################################################################
########################################################################
########################################################################
########################################################################

# Yellowtail Kingfish monthly predictions

library(dismo)
library(raster)
library(sf)
library(rJava)
library(jpeg)
require(stringr)

# set your working directory
#setwd("~/Documents/Use_cases/Marine/Yellowtail_R")

## read occurrence data from Yellowtailed Kingfish
#YT <- read.csv("~/Documents/Use_cases/Marine/Yellowtail_R/data/ala_aodn_occurrences_filtered.csv", stringsAsFactors = FALSE)
#head(YT)
#length(YT$Long)
# the sprintf function returns a character vector containing a formatted combination of text and variable values
YT_2$month <- sprintf("%02d", as.numeric(YT$month))
YT_2$mo <- as.factor(as.numeric(YT$month))
head(YT_2)

# check if Maxent is in the folder
jar <- paste(system.file(package = "dismo"), "/java/maxent.jar", sep = '') 
if (file.exists(jar)) {
  cat("can continue, maxent is available")
} else {
  cat('cannot run this because maxent is not available')
}

setwd("~/Downloads/Marine/monthly_predictors/01")
rast_lst <- list.files(pattern='.tif$', all.files=TRUE)
rast_lst<-rast_lst[-3]
YT_env_vars <- raster::stack(rast_lst)
crs(YT_env_vars)<-crs(bath)
YT_env_vars


###########  Check for correlation in grids

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use="pairwise.complete.obs"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
rpoints<-randomPoints(YT_env_vars,1000)
samp<-raster::extract(YT_env_vars,rpoints)
pairs(samp,lower.panel=panel.smooth,upper.panel=panel.cor)



x_YTKF<- YT_2$lon
y_YTKF<-YT_2$lat
xy_YTKF<-cbind(x_YTKF,y_YTKF)
xy.sp_YTKF<-SpatialPoints(xy_YTKF)
crs(xy.sp_YTKF) <- crs(YT_env_vars)
crs(xy.sp_YTKF)

plot(YT_env_vars[[1]])
points(xy.sp_YTKF,pch=20,cex=0.2)

mask_r <- YT_env_vars[[2]]
ext <- raster::extent(YT_env_vars)



############# extract monthly data

bias <- raster("~/Downloads/Marine/bias/Bias_AODN.tif")

dir.create("data")
dir.create("~Downloads/Marine/data/env_dframe")



t1<-unique(YT_2$month)

t1<-sort(t1)

#defaultW <- getOption("warn") 

#options(warn = -1) 

for (i in 1:length(unique(t1))){
  subsetYT<-YT[which(YT$month==unique(t1)[i]), ]
  longitude<-subsetYT_2$lon
  latitude<-subsetYT_2$lat
  xy<-cbind(longitude,latitude)
  xy.sp<-SpatialPoints(xy)
  raster::crs(xy.sp) <- raster::crs(YT_env_vars)
  
  chl<-raster(paste0("~/Downloads/Marine/monthly_predictors/",unique(t2)[i],"/chl.tif"))
  
  gsla<-raster(paste0("~/Downloads/Marine/monthly_predictors/",unique(t2)[i],"/gsla.tif"))
  
  mo<-raster(paste0("~/Downloads/Marine/monthly_predictors/",unique(t2)[i],"/mo.tif"))
  
  sst<-raster(paste0("~/Downloads/Marine/monthly_predictors/",unique(t2)[i],"/sst.tif"))
  
  vcur<-raster(paste0("~/Downloads/Marine/monthly_predictors/",unique(t2)[i],"/vcur.tif"))
  
  EVstack<-raster::stack(chl,gsla,mo,sst,vcur)
  raster::crs(EVstack) <- raster::crs(YT_env_vars)
  
  presEV <- as.data.frame(raster::extract(EVstack,xy.sp))
  pres_p <- cbind(longitude,latitude,presEV)
  pres_p$pres <- 1
  
  
  rpoints<-as.data.frame(xyFromCell(bias, sample(which(!is.na(values(bias))), 2000, prob=values(bias)[!is.na(values(bias))])))
  colnames(rpoints)<-c("longitude","latitude")
  xy.bk<-SpatialPoints(rpoints)
  
  bkgEV <- as.data.frame(raster::extract(EVstack,xy.bk))
  bkg_p <- cbind(rpoints,bkgEV)
  bkg_p$pres <- 0
  
  dat1 <- rbind(pres_p,bkg_p)
  write.csv(dat1,paste0("~/Downloads/Marine/data/env_dframe/env_",unique(t2)[i],".csv"))
}

files <- list.files("~/Downloads/Marine/data/env_dframe/")
setwd("~/Downloads/Marine/data/env_dframe/")

data1 <- do.call("rbind",lapply(files,read.csv,header=TRUE))
data1$mo <- as.factor(data1$mo)
pa<- data1$pres

data1$pres<-NULL
data1$X<-NULL
data1$longitude<-NULL
data1$latitude<-NULL

write.csv(data1,"~/Downloads/Marine/data/maxent_env_dat.csv")
write.csv(pa,"~/Downloads/Marine/data/pa_vector.csv")

## This is the overall Maxent model using data from all months and locations
maxent_args <- c('removeduplicates=TRUE','jackknife=TRUE','responsecurves=TRUE','plots=TRUE','betamultiplier=3')
#setwd("~/Documents/Use_cases/Marine/Yellowtail_R/")
Maxent_model1 <- dismo::maxent(x=data1,p=pa,path = "results_April22",args= maxent_args)

#plot the response curves
response(Maxent_model1)

# look in results_April22 folder for more results

# plot monthly predictions

for (i in 1:length(unique(t1))){
  subsetYT<-YT_2[which(YT_2$month==unique(t1)[i]), ]
  longitude<-subsetYT$lon
  latitude<-subsetYT$lat
  xy<-cbind(longitude,latitude)
  xy.sp<-SpatialPoints(xy)
  raster::crs(xy.sp) <- raster::crs(YT_env_vars)
  
  chl<-raster(paste0("~/Downloads/Marine/monthly_predictors/",unique(t2)[i],"/chl.tif"))
  
  gsla<-raster(paste0("~/Downloads/Marine/monthly_predictors/",unique(t2)[i],"/gsla.tif"))
  
  mo<-raster(paste0("~/Downloads/Marine/monthly_predictors/",unique(t2)[i],"/mo.tif"))
  
  sst<-raster(paste0("~/Downloads/Marine/monthly_predictors/",unique(t2)[i],"/sst.tif"))
  
  vcur<-raster(paste0("~/Downloads/Marine/monthly_predictors/",unique(t2)[i],"/vcur.tif"))
  
  Vstack<-stack(chl,gsla,mo,sst,vcur)
  raster::crs(EVstack) <- raster::crs(YT_env_vars)
  
  filenameM <- paste0("/results_April22/map_pred_",unique(t1)[i],".tif")
  map_predictions <- predict(Maxent_model1, Vstack,filename=paste0(getwd(),filenameM), args= c('outputformat=cloglog','outputgrids=TRUE','applythresholdrule=Maximum training sensitivity plus specificity'))
  
  jpeg(paste0(getwd(),"/results_April22/YTKF_predicted",unique(t1)[i],".jpeg"))
  plot(map_predictions, main=paste0("Yellowtail Kingfish distribution month =",unique(t1)[i]))
  points(xy.sp,pch=20,cex=0.2)
  dev.off()
  
  results<-read.csv(paste0(getwd(),"/results_April22/maxentResults.csv"))
  
  thresh<-results$Maximum.training.sensitivity.plus.specificity.area
  m <- c(0, thresh, 0,  thresh, 1, 1)
  reclass <- matrix(m, ncol= 3, byrow= TRUE)
  rc <- reclassify(map_predictions, reclass)
  
  jpeg(paste0(getwd(),"/results_April22/YTKF_thresholded",unique(t1)[i],".jpeg"))
  plot(rc, main=paste0("Yellowtail Kingfish distribution month =",unique(t1)[i]))
  points(xy.sp,pch=20,cex=0.2)
  dev.off()
}
