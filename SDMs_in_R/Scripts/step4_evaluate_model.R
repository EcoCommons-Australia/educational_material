# step 4 evaluate the model, be sure to check all these packages are installed

library(caret)
library(dismo)
library(gbm)
library(raster)
library(sp)
library(jpeg)
library(galah)
library(tidyr)

# There are many ways to evaluate your model.  With statistical models you can evaluate residual plots to check that your model meets assumptions.
#
## It is often a good idea to break your data into randomly selected training and testing data.  Often you would randomly remove 20% or more of 
# your data and withold that from the model building process.  Once your model was finsihed using your training data, you would then test your 
# model with this training data. If you have very little data, bootstrapping can be used to see if removal of a small percentage of data repeatedly 
# changes results.  This gives you a good understanding of the confidence intervals around your results and can reduce the impact of outliers on your 
# inal result.  Cross-validation requires more data.  In a perfect world for cross-validation you break your data into 10 folds (subsets) and compare results between 
# folds, or summarise variability between folds.  If you can afford to set aside 20% + of your data, withholding a test data set 
# is good practice.  Ideally, you will test your model with completely independent data.

# here we evaluate the BRT model with the training data used to construct the model

# First we run the BRT model again without additional variables identified in our simplify function

LiPe_mod_new <- gbm.step(data = LiPe_preds3, gbm.x = brt_drops$pred.list[[2]], gbm.y = 3, family = "bernoulli", tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.75)

LiPe_1s <- LiPe_preds3[LiPe_preds3$pres==1,]
x_1s<- LiPe_1s$Long
y_1s<-LiPe_1s$Lat
xy_1s<-cbind(x_1s,y_1s)
xy.sp_1s<-SpatialPoints(xy_1s)
crs(xy.sp_1s) <- crs(LiPe_predictors)
crs(xy.sp_1s)

LiPe_0s <- LiPe_preds3[LiPe_preds3$pres==0,]
x_0s<- LiPe_0s$Long
y_0s<-LiPe_0s$Lat
xy_0s<-cbind(x_0s,y_0s)
xy.sp_0s<-SpatialPoints(xy_0s)
crs(xy.sp_0s) <- crs(LiPe_predictors)
crs(xy.sp_0s)

# look at which predictors are in our predictors stack
names(LiPe_predictors)

# Here we drop the variables that were dropped from the LiPe_mod_new BRT model
# This is so we can predict all locations using the subset of variables used in the model
# In any model, the variables you use, need to match the variables you use to predict
# column names, or raster stack names need to match exactly

LiPe_preds_brt_new <- LiPe_predictors[[which(c(TRUE,TRUE,TRUE,TRUE,FALSE,FALSE))]]

brt_pred <- predict(LiPe_preds_brt_new,LiPe_mod_new,type="link") # note we usually use type link, but the evaluate function is based on the "link" scale

brt_eval_training <- dismo::evaluate(xy.sp_1s, xy.sp_0s, LiPe_mod_new, LiPe_preds_brt_new)

#There are a variety of ways to calculate a threshold, maximising "kappa" is one well supported method, but there are many others  
# in Maxent we saw a table with a variety of ways to calculate thresholds

# 
thresh_brt <- threshold(brt_eval_training,stat = "kappa")

# another approach to thresholding - there are many choose the one that best fits in your model
thresh_brt2 <- threshold(brt_eval_training,stat = "spec_sens")

m <- c(-5, thresh_brt, 0,  thresh_brt, 2.1, 1)
reclass <- matrix(m, ncol= 3, byrow= TRUE)
rc_brt <- reclassify(brt_pred, reclass)

plot(rc_brt)
points(xy.sp_LiPe,pch=20,cex=0.2)


# look at training evaluation
brt_eval_training

x_s<- LiPe_preds3$Long
y_s<-LiPe_preds3$Lat
xy_s<-cbind(x_s,y_s)
xy.sp_s<-SpatialPoints(xy_s)
crs(xy.sp_s) <- crs(LiPe_predictors)

predicted <- raster::extract(rc_brt,xy.sp_s)

test1 <- as.data.frame(cbind(predicted = predicted, reference = LiPe_preds3$pres))

test1$reference <- factor(test1$reference, levels = c("1","0"))
test1$predicted <- factor(test1$predicted, levels = c("1","0"))

conf_matrix <- confusionMatrix(test1$predicted,test1$reference)

# this gives us the confusion matrix
conf_matrix$table

TP <- conf_matrix$table[1,1]
FP <- conf_matrix$table[1,2]
FN <- conf_matrix$table[2,1]
TN <- conf_matrix$table[2,2]

FPR <- FP/(FP + TN)
FNR <- FN/(FN + TP)
TPR <- TP/(TP + FN)
TNR <- TN/(TN + FP)

TSS <- TPR + TNR - 1  # if we were predicting perfectly our TSS would = 1
TSS # we are far from perfect in this model, and a much lower number that AUC
# a perfect AUC value would also be 1
# TSS often gives a better indication of prediction because it takes into
# account anbalanced errors between positives and negatives

Precision <- TP / (TP + FP)
Recall <- TP / (TP + FN)
F1 <- (2*Precision*Recall)/(Precision+Recall)#  F1 is another statistic, now thought to be
# better than TSS


F1 <- TP /(TP + 0.5*(FP + FN)) # this is another way to calculate the same value
# but the literature refers to Precision and Recall above (same thing though)
# again perfect score would be 1

F1 # as you can see in our case F1 falls between AUC and TSS scores

# now we are going to test our results from ALA records that were not from the FrogID project
# for the most part we are just repeating code from step1

galah_config(email = "xxxxx@griffith.edu.au") #insert your own email address

LiPe_ALA<-galah_call() %>%
  galah_identify("Limnodynastes peroni")%>%
  galah_filter(datasetName != "FrogID")%>%
  galah_filter(coordinateUncertaintyInMeters < 200)%>%
  galah_filter(year>1999)%>%
  galah_select("datasetName","year")%>%
  atlas_occurrences()

LiPe_ALA<-na.omit(LiPe_ALA)

LiPe_pts_ALA <- SpatialPoints(coords = cbind(LiPe_ALA$decimalLongitude, LiPe_ALA$decimalLatitude),CRS(as.character("+proj=longlat +datum=WGS84 +no_defs")) )

LiPe_pres_ALA<- rep(1,length(LiPe$decimalLatitude))

# Then create a raster with a value of 1 for each gridcell where a LiPe was recorded
LiPe_pres_raster_ALA<-rasterize(LiPe_pts_ALA,base,LiPe_pres_ALA, fun = min, background=0)

LiPe_pres_raster2_ALA<-mask(LiPe_pres_raster,base)

cell_no_ALA<- raster::extract(large_base,LiPe_pts_ALA,cellnumbers=TRUE)

LiPe_cells_ALA<- cbind(LiPe_ALA,cell_no_ALA)

require(dplyr)
LiPe_thinned_ALA <- LiPe_cells_ALA %>% 
  group_by(cells) %>% 
  slice_sample(n = 1)


frogs_ALA<-galah_call() %>%
  galah_identify("Amphibia")%>%
  galah_filter(datasetName != "FrogID")%>%
  galah_filter(coordinateUncertaintyInMeters < 100)%>%
  galah_filter(year>1999)%>%
  galah_select("datasetName","year")%>%
  atlas_occurrences()

frogs_ALA$unique_visit<- paste0(frogs_ALA$decimalLatitude,frogs_ALA$decimalLongitude,frogs_ALA$eventDate)
frogs_ALA$visitID <- as.numeric(as.factor(frogs_ALA$unique_visit))

frogs_ALA2 <- frogs_ALA %>%
  group_by(decimalLatitude, decimalLongitude, visitID) %>%
  summarise(no_spp = length(unique(scientificName)))

frogs_ALA3 <- frogs_ALA2 %>%
  group_by(decimalLatitude, decimalLongitude) %>%
  summarise(no_visits = length(visitID))

frogs_ALA3 <- na.omit(frogs_ALA3)

visits_pts_ALA<-SpatialPoints(coords = cbind(frogs_ALA3$decimalLongitude, frogs_ALA3$decimalLatitude),CRS(as.character("+proj=longlat +datum=WGS84 +no_defs")))

b1_ALA <- rasterize(visits_pts_ALA, base, frogs_ALA3$no_visits, fun=sum, background=0)

bb_ALA <- bbox(b1_ALA)

visit_locations_ALA<- raster::extract(b1_ALA, visits_pts_ALA, cellnumbers=TRUE)

visit_locations2_ALA <- as.data.frame(na.omit(visit_locations_ALA))

cellID_ALA <- unique(visit_locations2_ALA$cells)

xy_visits_ALA <- raster::xyFromCell(b1_ALA,cell = cellID_ALA)

cellStats(b1_ALA ,"max")
m <- c(0, 2.9, 0,  2.9, 880, 1)
reclass <- matrix(m, ncol= 3, byrow= TRUE)
rc_ALA <- reclassify(b1_ALA, reclass)

visits_3or_more_ALA <- mask(rc_ALA,base)

Zero_LiPe_ALA <- visits_3or_more_ALA - LiPe_pres_raster2_ALA

freq(Zero_LiPe_ALA)

m <- c(-2, 0.1, 0,  0.1, 2, 1)
reclass2 <- matrix(m, ncol= 3, byrow= TRUE)
rc2_ALA <- reclassify(Zero_LiPe_ALA, reclass2)
Zero_LiPe2_ALA<-mask(rc2_ALA,base)
freq(Zero_LiPe2_ALA)

# extract the cell numbers from the 0 grid where the value ==1
cell_vals_0_ALA<-Which(Zero_LiPe2_ALA ==1,cells=TRUE)
# these are the lat / longs for locations where at least three surveys were done, but zero LiPe were detected - these are our pseudo absences
xy_zero_LiPe_ALA <- xyFromCell(Zero_LiPe2_ALA,cell = cell_vals_0_ALA)
#plot the absence locations
plot(base)
points(xy_zero_LiPe_ALA, pch=20,cex=0.2)

xy_zero_LiPe_ALA <- as.data.frame(xy_zero_LiPe_ALA)

## Now lets turn our presence and absence points into one dataset

head(LiPe_thinned_ALA)

ALA_preds1 <- as.data.frame(cbind(Long = LiPe_thinned_ALA$decimalLongitude, Lat = LiPe_thinned_ALA$decimalLatitude,pres = 1))
ALA_preds2 <- as.data.frame(cbind(Long = xy_zero_LiPe_ALA$x, Lat = xy_zero_LiPe_ALA$y,pres = 0))

ALA_preds <- rbind(ALA_preds1,ALA_preds2)

x_a<- ALA_preds$Long
y_a<-ALA_preds$Lat
xy_a<-cbind(x_a,y_a)
xy.sp_a<-SpatialPoints(xy_a)
crs(xy.sp_a) <- crs(LiPe_predictors)

## now we just extract the 1's & 0's using our lat longs from the new ALA data 
# extracted from the same thresholded brt predictions from above

predicted_a <- raster::extract(rc_brt,xy.sp_a)

test2 <- as.data.frame(cbind(predicted = predicted_a, reference = ALA_preds$pres))

test2$reference <- factor(test2$reference, levels = c("1","0"))
test2$predicted <- factor(test2$predicted, levels = c("1","0"))

conf_matrix2 <- confusionMatrix(test2$predicted,test2$reference)
conf_matrix2

TP2 <- conf_matrix2$table[1,1]
FP2 <- conf_matrix2$table[1,2]
FN2 <- conf_matrix2$table[2,1]
TN2 <- conf_matrix2$table[2,2]

FPR2 <- FP2/(FP2 + TN2)
FNR2 <- FN2/(FN2 + TP2)
TPR2 <- TP2/(TP2 + FN2)
TNR2 <- TN2/(TN2 + FP2)

TSS2 <- TPR2 + TNR2 - 1  
TSS2 # notice this is a big drop from our original TSS statistic
TSS

Precision2 <- TP2 / (TP2 + FP2)
Recall2 <- TP2 / (TP2 + FN2)
F1_b <- 2*((Precision2*Recall2)/(Precision2+Recall2))
F1_b # notice this is a very slight drop from our original F1 statsitic
F1

# When we used data to test the model that was not associated with
# the FrogID project, TSS, accuracy, & precision scores all fell markedly
# Encouragingly the F1 score did not fall that much, and as a user
# you will need to decidd on what levels of accuracy in your
# confusion matrix are good enough for the intended use of the model

# If I was looking to predict well in order to select what reserves
# to make for this frog, I would want better performance.

# The first step I would take to improve this model, would be to select
# pseudo_absence locations at only those locations where similarly
# distributed frogs were seen but our target species was not seen
# on three visits. Our psuedo absence points are still including
# places in the arid interior where this species would never occur.

# Second, as mentioned before, I would also look for better freshwater
# wetland layers, perhaps look at how often an area is wet over the last
# decoade.




# OTHER places to learn how to generate SDMs in R
# http://www.earthskysea.org/best-practices-in-species-distribution-modeling-a-workshop-in-r/
# a series of lectures on SDMs https://www.youtube.com/watch?v=obuMW5NAtJE 

