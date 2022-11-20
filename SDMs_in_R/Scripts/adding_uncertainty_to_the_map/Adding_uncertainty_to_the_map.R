##  Mapping uncertainty example script

require(dismo)
require(sp)
require(raster)
library(rJava)
library(biomod2)

# check if maxent is available

jar <- paste(system.file(package = "dismo"), "/java/maxent.jar", sep = '') 
if (file.exists(jar)) {
  cat("can continue, maxent is available")
} else {
  cat('cannot run this because maxent is not available')
}

# set your working directory to where you downloaded the data
# setwd("~/Your_path")
Bias <- raster("Bias.tif")
Forest <- raster("Forest_cover.tif")
Height <- raster("Forest_Height.tif")
NDVI <- raster("NDVI.tif")
Remnant <- raster("Remnant_cover.tif")
Stream <- raster("Stream_density.tif")

preds <- stack(Forest,Height,NDVI,Remnant,Stream)

# An example of a way to download data directly from EcoCommons
DEM <- raster("https://api.data-ingester.app.ecocommons.org.au/api/data/fde33c41-b4ec-53a0-8fcb-44caeb07eece/download/digital_elevation_model_9s_2008_9s.tif")

DEM2<-projectRaster(DEM,preds) # woops, made 5 copies
DEM3 <- mask(DEM2,preds)
DEM4 <- DEM3[[1]]

# Grab all species X data for the region (includes non breeding sightings)
# set WD to where you downloaded the species records
all <- read.csv("2020_speciesX_records.csv")
longlat4 <- SpatialPoints(all)
crs(longlat4) <- crs(preds)

longlat5 <- crop(longlat4,preds[[3]])

plot(preds[[3]])
points(longlat5)

extr <- extract(preds[[3]],longlat5)
ex2<-as.data.frame(extr)
coord <- as.data.frame(geom(longlat5))
t2 <- cbind(ex2$extr,coord$x,coord$y)
t3 <- as.data.frame(t2)
t4 <- t3[complete.cases(t3),]

longitude <- as.numeric(t4$V2)
latitude <- as.numeric(t4$V3)
lonlat3 <- cbind(longitude,latitude)
longlat3 <- SpatialPoints(lonlat3)
crs(longlat3) <- crs(preds)

df2<-data.frame(data=rep(1,6208))
spdf2<-SpatialPointsDataFrame(longlat3, df2)

preds2 <- stack(preds,DEM4)

# the biomod2 package is generating some odd errors, try running the bm_SRE functions below
# if that doesn't work, upload the SRE results

sre.100 <- bm_SRE(resp.var = spdf2, 
    expl.var = preds2, 
    new.env = preds2, 
    quant = 0)

plot(sre.100)

sre.90 <- 
  bm_SRE(
    resp.var = spdf2, 
    expl.var = preds2, 
    new.env = preds2, 
    quant = 0.05
  )

plot(sre.90)

sre.80 <- 
  bm_SRE(
    resp.var = spdf2, 
    expl.var = preds2, 
    new.env = preds2, 
    quant = 0.1
  )

plot(sre.80)

sre.70 <- 
  bm_SRE(
    resp.var = spdf2, 
    expl.var = preds2, 
    new.env = preds2, 
    quant = 0.15
  )

plot(sre.70)

sre.60 <- 
  bm_SRE(
    resp.var = spdf2, 
    expl.var = preds2, 
    new.env = preds2, 
    quant = 0.2
  )

plot(sre.60)

## If the above functions were not working upload the resulting files and stack them

sre.100 <- raster("sre.100.tif")
sre.90 <- raster("sre.90.tif")
sre.80 <- raster("sre.80.tif")
sre.70 <- raster("sre.70.tif")
sre.60 <- raster("sre.60.tif")


SRE <- stack(sre.100,sre.90,sre.80,sre.70,sre.60)

# when we add these together, the areas that most closely match the available environmental space have high values
SRE_v <- sum(SRE)

plot(SRE_v)
# by subtracting each value from the maximum value the places with the least certainty of containing environmental space get the largest values
inv_SRE <- 6 - SRE_v
inv_SRE2 <- scale(inv_SRE,center = FALSE, scale = TRUE)
SRE_base <- inv_SRE2/inv_SRE2
SRE_max <- SRE_base*1.326596
SRE_min <- SRE_base*0.221099
SRE_scale <- (inv_SRE2-SRE_min)/(SRE_max-SRE_min)

##### Here we take the maximum value in the bias layer and subtract bias values
#  essentially this sets the locations with the fewest bird surveys with the highest values

inv_Bias <- 280 - Bias
# we then scale the results
inv_Bias2 <- scale(inv_Bias,center = FALSE, scale = TRUE)
Bias_base <- inv_Bias2/inv_Bias2
Bias_max <- Bias_base*1.270755
Bias_min <- Bias_base*0.006793233
Bias_scale <- (inv_Bias2-Bias_min)/(Bias_max-Bias_min)

# this takes the statistical union of these two rasters
uncertainty <- (1-((1-Bias_scale)*(1-SRE_scale))) 
plot(uncertainty)
# this step just gives more separation of the colors for the map
uncertainty2 <- ((uncertainty)^5)
plot(uncertainty2)

bias_geo_env <- 5 - uncertainty2


##### select background points using the values in the bias layer to define the probability of selecting any cell


bg <- xyFromCell(bias_geo_env, sample(which(!is.na(values(bias_geo_env))), 10000, prob=values(bias_geo_env)[!is.na(values(bias_geo_env))]))

#### thinned occurrence data

base <- DEM4/DEM4

large_base <- aggregate(base, fact=8,expand=TRUE)
cell_no<- raster::extract(large_base,longlat3,cellnumbers=TRUE)

cells_df<- as.data.frame(cbind(geom(longlat3),cell_no))
head(cells_df)
# by selecting only one value in these much larger grid cells we are thinning out areas where
# there are high densities of records

require(dplyr)
thinned <- cells_df %>% 
  group_by(cells) %>% 
  slice_sample(n = 1)

LatLon_thin <- SpatialPoints(coords = cbind(thinned$x, thinned$y),CRS(as.character("+proj=longlat +datum=WGS84 +no_defs")) )

Mod_thin <- dismo::maxent(x=preds, p=LatLon_thin,a=bg)
pred_maxent_thin <- predict(Mod_thin,preds)
plot(pred_maxent_thin)


# Attempt to create bivariate raster

library(classInt)
library(raster)
library(rgdal)
library(dismo)
library(XML)
library(maps)
library(sp)

# these functions are great at generating a biavariate map of two rasters

# https://rfunctions.blogspot.com/2015/03/bivariate-maps-bivariatemap-function.html 

# Below we use estimated suitability on the y axis, and uncertainty on the x axis

# This first function creates that legend grid
colmat<-function(nquantiles=10, upperleft=rgb(0,150,235, maxColorValue=255), upperright=rgb(130,0,80, maxColorValue=255), bottomleft="grey", bottomright=rgb(255,230,15, maxColorValue=255), xlab="x label", ylab="y label"){
  my.data<-seq(0,1,.01)
  my.class<-classIntervals(my.data,n=nquantiles,style="quantile")
  my.pal.1<-findColours(my.class,c(upperleft,bottomleft))
  my.pal.2<-findColours(my.class,c(upperright, bottomright))
  col.matrix<-matrix(nrow = 101, ncol = 101, NA)
  for(i in 1:101){
    my.col<-c(paste(my.pal.1[i]),paste(my.pal.2[i]))
    col.matrix[102-i,]<-findColours(my.class,my.col)}
  plot(c(1,1),pch=19,col=my.pal.1, cex=0.5,xlim=c(0,1),ylim=c(0,1),frame.plot=F, xlab=xlab, ylab=ylab,cex.lab=1.3)
  for(i in 1:101){
    col.temp<-col.matrix[i-1,]
    points(my.data,rep((i-1)/100,101),pch=15,col=col.temp, cex=1)}
  seqs<-seq(0,100,(100/nquantiles))
  seqs[1]<-1
  col.matrix<-col.matrix[c(seqs), c(seqs)]}
#This line runs the above function, have a play with different colors
col.matrix<-colmat(nquantiles=10, upperleft="blue", upperright="grey", bottomleft="yellow", bottomright="black", xlab="Uncertainty", ylab="Suitability")

# This function generates the new colors for your mapped area based on two rasters with the same extent, resolution and CRS
bivariate.map<-function(rasterx, rastery, colormatrix=col.matrix, nquantiles=10){
  quanmean<-getValues(rasterx)
  temp<-data.frame(quanmean, quantile=rep(NA, length(quanmean)))
  brks<-with(temp, quantile(temp,na.rm=TRUE, probs = c(seq(0,1,1/nquantiles))))
  r1<-within(temp, quantile <- cut(quanmean, breaks = brks, labels = 2:length(brks),include.lowest = TRUE))
  quantr<-data.frame(r1[,2]) 
  quanvar<-getValues(rastery)
  temp<-data.frame(quanvar, quantile=rep(NA, length(quanvar)))
  brks<-with(temp, quantile(temp,na.rm=TRUE, probs = c(seq(0,1,1/nquantiles))))
  r2<-within(temp, quantile <- cut(quanvar, breaks = brks, labels = 2:length(brks),include.lowest = TRUE))
  quantr2<-data.frame(r2[,2])
  as.numeric.factor<-function(x) {as.numeric(levels(x))[x]}
  col.matrix2<-colormatrix
  cn<-unique(colormatrix)
  for(i in 1:length(col.matrix2)){
    ifelse(is.na(col.matrix2[i]),col.matrix2[i]<-1,col.matrix2[i]<-which(col.matrix2[i]==cn)[1])}
  cols<-numeric(length(quantr[,1]))
  for(i in 1:length(quantr[,1])){
    a<-as.numeric.factor(quantr[i,1])
    b<-as.numeric.factor(quantr2[i,1])
    cols[i]<-as.numeric(col.matrix2[b,a])}
  r<-rasterx
  r[1:length(r)]<-cols
  return(r)}

# this generates the data
bivmap<-bivariate.map(uncertainty2,pred_maxent_thin,colormatrix=col.matrix, nquantiles=10)
# this plots the bivariate map
plot(bivmap,frame.plot=F,axes=F,box=F,add=F,legend=F,col=as.vector(col.matrix))
# map(interior=T,add=T) # this adds the country boundary to the plot

## this bivariate mapping package works well with polygon data
library(remotes)
remotes::install_github(repo = "lydialucchesi/Vizumap", build_vignettes = TRUE, force = TRUE)
library(Vizumap)
vignette("Vizumap")

cmBivPal <- build_palette(name = "CyanMagenta")
view(cmBivPal)


## Targetted background

# GBIF.org (26 October 2022) GBIF Occurrence Download https://doi.org/10.15468/dl.zk8ekz # Southern Boobook
# GBIF.org (26 October 2022) GBIF Occurrence Download https://doi.org/10.15468/dl.z57gut # Masked Owl
# GBIF.org (26 October 2022) GBIF Occurrence Download https://doi.org/10.15468/dl.thremx #Greater Sooty Owl
