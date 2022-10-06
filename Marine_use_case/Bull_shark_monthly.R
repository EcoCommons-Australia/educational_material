# Bull Shark monthly predictions

library(dismo)
library(raster)
library(sf)
library(rJava)
library(jpeg)
require(stringr)
library(galah)

galah_config(email = "r.clemens@griffith.edu.au")

# Bull shark

CL <- galah_call() %>%
  galah_identify("Carcharhinus leucas")%>%
  galah_filter(coordinateUncertaintyInMeters < 1000)%>%
  galah_select("month")%>%
  atlas_occurrences()

CL2<-CL[,c(1:3)]

CL3 <- na.omit(CL2)

colnames(CL3)<-c("month","Lat","Long")

CL<-CL3

head(CL)
length(CL$Long)
# the sprintf function returns a character vector containing a formatted combination of text and variable values
CL$month <- sprintf("%02d", as.numeric(CL$month))
CL$mo <- as.numeric(CL$month)
head(CL)

# check if Maxent is in the folder
jar <- paste(system.file(package = "dismo"), "/java/maxent.jar", sep = '') 
if (file.exists(jar)) {
  cat("can continue, maxent is available")
} else {
  cat('cannot run this because maxent is not available')
}


x_CL<- CL$Long
y_CL<-CL$Lat
xy_CL<-cbind(x_CL,y_CL)
xy.sp_CL<-SpatialPoints(xy_CL)
crs(xy.sp_CL) <- crs(YT_env_vars)
crs(xy.sp_CL)

plot(EVstack[[1]])
points(xy.sp_CL,pch=20,cex=0.2)

mask_r <- CL_env_vars[[2]]
ext <- raster::extent(CL_env_vars)


############# extract monthly data


t1<-unique(CL$month)

t1<-sort(t1)
t1

#defaultW <- getOption("warn") 

#options(warn = -1) 

for (i in 1:length(unique(t1))){
  subsetCL<-CL[which(CL$month==unique(t1)[i]), ]
  longitude<-subsetCL$Long
  latitude<-subsetCL$Lat
  xy<-cbind(longitude,latitude)
  xy.sp<-SpatialPoints(xy)
  crs(xy.sp) <- crs(bath2)
  
  chl<-raster(paste0("~/Documents/Use_cases/Marine/Yellowtail_R/data/monthly_predictors/",unique(t2)[i],"/chl.asc"))
  
  gsla<-raster(paste0("~/Documents/Use_cases/Marine/Yellowtail_R/data/monthly_predictors/",unique(t2)[i],"/gsla.asc"))
  
  mo<-raster(paste0("~/Documents/Use_cases/Marine/Yellowtail_R/data/monthly_predictors/",unique(t2)[i],"/mo.asc"))
  
  sst<-raster(paste0("~/Documents/Use_cases/Marine/Yellowtail_R/data/monthly_predictors/",unique(t2)[i],"/sst.asc"))
  
  vcur<-raster(paste0("~/Documents/Use_cases/Marine/Yellowtail_R/data/monthly_predictors/",unique(t2)[i],"/vcur.asc"))
  
  EVstack<-stack(chl,gsla,mo,sst,vcur)
  crs(EVstack) <- crs(bath2)
  
  presEV <- as.data.frame(raster::extract(EVstack,xy.sp))
  pres_p <- cbind(longitude,latitude,presEV)
  pres_p$pres <- 1
  
  # I did not use a bias layer for background points for Bull shark, just random points
  rpoints<-as.data.frame(randomPoints(EVstack, 2000))
  colnames(rpoints)<-c("longitude","latitude")
  xy.bk<-SpatialPoints(rpoints)
  
  bkgEV <- as.data.frame(raster::extract(EVstack,xy.bk))
  bkg_p <- cbind(rpoints,bkgEV)
  bkg_p$pres <- 0
  
  dat1 <- rbind(pres_p,bkg_p)
  write.csv(dat1,paste0("~/Documents/Use_cases/Marine/Bull_shark/data/env_dframe/env_",unique(t2)[i],".csv"))
}

files <- list.files("~/Documents/Use_cases/Marine/Bull_shark/data/env_dframe/")
setwd("~/Documents/Use_cases/Marine/Bull_shark/data/env_dframe/")

data1 <- do.call("rbind",lapply(files,read.csv,header=TRUE))
data1$mo <- as.factor(data1$mo)
pa<- data1$pres

data1$pres<-NULL
data1$X<-NULL
data1$longitude<-NULL
data1$latitude<-NULL

write.csv(data1,"~/Documents/Use_cases/Marine/Bull_shark/data/maxent_env_dat.csv")
write.csv(pa,"~/Documents/Use_cases/Marine/Bull_shark/data/pa_vector.csv")

maxent_args <- c('removeduplicates=TRUE','jackknife=TRUE','responsecurves=TRUE','plots=TRUE','betamultiplier=3')
setwd("~/Documents/Use_cases/Marine/Bull_shark/")
Maxent_model1 <- dismo::maxent(x=data1,p=pa,path = "results_April22",args= maxent_args)

#plot the response curves
response(Maxent_model1)

# look in results_April22 folder for more results

# plot monthly predictions

for (i in 1:length(unique(t1))){
  subsetCL<-CL[which(CL$month==unique(t1)[i]), ]
  longitude<-subsetCL$Long
  latitude<-subsetCL$Lat
  xy<-cbind(longitude,latitude)
  xy.sp<-SpatialPoints(xy)
  crs(xy.sp) <- crs(bath2)
  
  chl<-raster(paste0("~/Documents/Use_cases/Marine/Yellowtail_R/data/monthly_predictors/",unique(t2)[i],"/chl.asc"))
  
  gsla<-raster(paste0("~/Documents/Use_cases/Marine/Yellowtail_R/data/monthly_predictors/",unique(t2)[i],"/gsla.asc"))
  
  mo<-raster(paste0("~/Documents/Use_cases/Marine/Yellowtail_R/data/monthly_predictors/",unique(t2)[i],"/mo.asc"))
  
  sst<-raster(paste0("~/Documents/Use_cases/Marine/Yellowtail_R/data/monthly_predictors/",unique(t2)[i],"/sst.asc"))
  
  vcur<-raster(paste0("~/Documents/Use_cases/Marine/Yellowtail_R/data/monthly_predictors/",unique(t2)[i],"/vcur.asc"))
  
  Vstack<-stack(chl,gsla,mo,sst,vcur)
  crs(Vstack) <- crs(bath2)
  
  filenameM <- paste0("/results_April22/map_pred_",unique(t1)[i],".asc")
  map_predictions <- predict(Maxent_model1, Vstack,filename=paste0(getwd(),filenameM), args= c('outputformat=cloglog','outputgrids=TRUE','applythresholdrule=Maximum training sensitivity plus specificity'))
  
  jpeg(paste0(getwd(),"/results_April22/Bshark_predicted",unique(t1)[i],".jpeg"))
  plot(map_predictions, main=paste0("Bull Shark distribution month =",unique(t1)[i]))
  points(xy.sp,pch=20,cex=0.2)
  dev.off()
  
  results<-read.csv(paste0(getwd(),"/results_April22/maxentResults.csv"))
  
  thresh<-results$Maximum.training.sensitivity.plus.specificity.area
  m <- c(0, thresh, 0,  thresh, 1, 1)
  reclass <- matrix(m, ncol= 3, byrow= TRUE)
  rc <- reclassify(map_predictions, reclass)
  
  jpeg(paste0(getwd(),"/results_April22/Bshark_thresholded",unique(t1)[i],".jpeg"))
  plot(rc, main=paste0("Bull Shark distribution month =",unique(t1)[i]))
  points(xy.sp,pch=20,cex=0.2)
  dev.off()
}
