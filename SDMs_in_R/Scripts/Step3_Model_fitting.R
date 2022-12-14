# title: "Step3_fit_model"

direct<- "/Users/s2992269/Documents/Use_cases"
folder <- "/SDM_in_R"
setwd(paste0(direct,folder))

getwd()

# fit a Maxent model and then fit a BRT model
# java will need to be installed in your local machine
# and the Maxent files will need to be in your local directory
# dowload Maxent from here: https://biodiversityinformatics.amnh.org/open_source/maxent/
# ensure these files are in your local directory: maxent.jar  maxent.bat  maxent.sh

library(dismo)
library(gbm)
library(raster)
library(sf)
library(rJava)
library(boot)
library(jpeg)
library(MASSExtra)
library(mgcv)

## read in LiPe_predictors stack

#LiPe_predictors

## read occurrence data which was thinned for LiPe

LiPe_thinned <- read.csv("data/LiPe_thinned.csv", stringsAsFactors = FALSE)
head(LiPe_thinned)

# check if Maxent is in the folder
jar <- paste(system.file(package = "dismo"), "/java/maxent.jar", sep = '') 
if (file.exists(jar)) {
  cat("can continue, maxent is available")
} else {
  cat('cannot run this because maxent is not available')
}

setwd("~/Documents/Use_cases/SDM_in_R/predictors")
rast_lst <- list.files(pattern='.asc$', all.files=TRUE)
rast_lst
LiPe_predictors <- stack(rast_lst)
crs(LiPe_predictors) <- "+proj=longlat +datum=WGS84 +no_defs"
LiPe_predictors 

# earlier runs of this function indicated BioClim01 and BioClim14 were highly correlated with other variables, so they were dropped
# With this set of variables we still see high correlations between BioClim 12, AWAP and BioClim 05
# however, all correlation are below 0.80, and for the purposes of prediction these correlations are not too high
# for looking at comparisons of the importance of variables these correlations are still probably
# a bit high 

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
rpoints<-randomPoints(LiPe_predictors,1000)
samp<-extract(LiPe_predictors,rpoints)
pairs(samp,lower.panel=panel.smooth,upper.panel=panel.cor)

x_LiPe<- LiPe_thinned$decimalLongitude
y_LiPe<-LiPe_thinned$decimalLatitude
xy_LiPe<-cbind(x_LiPe,y_LiPe)
xy.sp_LiPe<-SpatialPoints(xy_LiPe)
crs(xy.sp_LiPe) <- crs(LiPe_predictors)
crs(xy.sp_LiPe)

plot(LiPe_predictors[[1]])
points(xy.sp_LiPe,pch=20,cex=0.2)

mask_r <- LiPe_predictors[[2]]
ext <- raster::extent(LiPe_predictors)

#might not need this step
#setwd("~/Documents/Use_cases/Marine/Yellowtail_R")


#output_file <- paste0(getwd(), "/output")
#if (!file.exists(output_file)) {
#  dir.create(output_file)
#}

#maxent_args <- c('removeduplicates=TRUE','jackknife=TRUE','replicates=3','replicatetype=crossvalidate','replicatetype=crossvalidate','betamultiplier=1','responsecurves=TRUE', 'pictures=TRUE','plots=TRUE','defaultprevalence=0.5','outputformat=raw')

#maxent_args <- c('removeduplicates=TRUE','jackknife=TRUE','responsecurves=TRUE', 'outputformat=cloglog','plots=TRUE','outputgrids=TRUE','applythresholdrule=Maximum training sensitivity plus specificity')

# we then choose which arguments we want to use, there are many scientific publications that look at which arguments to use when
# a good place to start to understand what Maxent is doing can be found here: https://biodiversityinformatics.amnh.org/open_source/maxent/Maxent_tutorial2017.pdf
# and this is an excellent overview of the process and the things to consider https://doi.org/10.1111/ecog.04960 
# how you fit your model depends on your question & your data
#*remove duplicates is good to make sure is true, just in case, but our thinning step should have removed all duplicates
# *jackknife gives an indication of how much each variable contributes to the results
# cloglog format is usually what gives optimal predicted values *outputformat = cloglog https://www.researchgate.net/figure/Comparison-of-logistic-and-cloglog-transforms-of-Maxent-output-For-simplicity-no_fig2_315501215
# * plots=TRUE gives a variety of plots in the output (results1) folder
# a complete list of args that can be used can be found below
# This is an argument that is explored alot in the literature *'betamultiplier=1', #multiply all automatic regularization parameters by this number. A higher number gives a more spread-out distribution.

maxent_args <- c('removeduplicates=TRUE','jackknife=TRUE','responsecurves=TRUE','plots=TRUE')

# this is the model fitting function, and path can be tricky, but when set correctly you get many results written to this folder

Mod1_LiPe <- dismo::maxent(LiPe_predictors, xy.sp_LiPe, path = "results", args= maxent_args)

# if you open up the results 1 folder and click on "maxent.html"  A webpage should open up showing some results
# The ROC curve and AUC score < 0.9 suggests we are not predicting our training super accurately - 
# if the goal of this exercise is prediction I would probably stop here and make sure My occurrence data is accurate, and think about what
# other variables might help us predict the presence of this species. Also, I would run the model multiple times
# Machine / stochatic learning algorithms often result in different results each time
# It is a good idea to run the model multiple times to get a sense of the variation
# We might get better predictions with a wetland layer that only included freshwater wetland, 
# perhaps averaging AWAP soil wetness layers over a decade would help? Perhaps there is a national layer that captures small wetlands better?
# Jackknife results, and variable importance do not indicate in my results that EVI is a good predictor, so we might look to drop that one from future models
# note those Jackknife results and variable response plots are also in the plots folder in the results folder
# Try deleting the results in results1 folder, and running the model again, Are results exactly the same?  WHY NOT?

# now lets predict a gridcell value for each cell in our extent

map_predictions <- predict(Mod1_LiPe, LiPe_predictors,args=maxent_args)
plot(map_predictions)
points(xy.sp_LiPe,pch=20,cex=0.2)

# if you want to save a jpg of the mapped predictions use the code below

require("jpeg")
setwd("~/Documents/Use_cases/SDM_in_R") # make sure your current directory is your working directory, use getwd() to ensure you have asigned the directory correctly
jpeg(paste0(getwd(),"/results/LiPe_predicted.jpeg"))
plot(map_predictions)
points(xy.sp_LiPe,pch=20,cex=0.2)
dev.off()

# This threshold maximises the cases where the model incorrectly assigns unsuitable habitat (true negative) and misses suitable habitat (false positive)
# keep in mind how you set your threshold will depend on your use of the resulting map, or your question

results1<-read.csv("results/maxentResults.csv")
thresh<-results1$Maximum.training.sensitivity.plus.specificity.area
thresh

m <- c(0, thresh, 0,  thresh, 1, 1)
reclass <- matrix(m, ncol= 3, byrow= TRUE)
rc <- reclassify(map_predictions, reclass)

plot(rc)
points(xy.sp_LiPe,pch=20,cex=0.2)

#  What if someone gave us $1million to buy some LiPe habitat?
# We might want to be really sure the species is present, if that is the case we might want to find a model with AUC > 0.9
# but just to illustrate the result lets try a much higher threshold
# first we want to rescale the data so values run between 0 and 1, then select a threshold of 0.8
# There are many reasons and ways to set thresholds, this is just one example
# of selecting a more conservative threshold that identifies locations where we are
# relatively more sure are suitable

min.R <- cellStats(map_predictions,"min")
max.R <- cellStats(map_predictions,"max")

map_predictions2 <- ((map_predictions - min.R) / (max.R - min.R) - 0 ) * 1

thresh2 <- 0.8

m2 <- c(0, thresh2, 0,  thresh2, 1, 1)
reclass2 <- matrix(m2, ncol= 3, byrow= TRUE)
rc2 <- reclassify(map_predictions2, reclass2)

plot(rc2)

# there are far fewer locations where we are pretty sure someone will find LiPe, and what if our original model is not that good, how can we assess this?

#######################
#######################
#######################
#######################

# now we are going to explore how consistent results are with a three fold-cross validation
# 3-folds or replicates = 3, when using cross validation as the replicate type

#*replicatetype=crossvalidate, *replicates=3

maxent_args2 <- c('removeduplicates=TRUE','jackknife=TRUE','responsecurves=TRUE', 'outputformat=cloglog','plots=TRUE','replicatetype=crossvalidate','replicates=3')

# this is the model fitting function, and path can be tricky, but when set correctly you get many results written to this folder

Mod2_LiPe <- dismo::maxent(LiPe_predictors, xy.sp_LiPe, path = "results2", args= maxent_args2)

# notice that in the folder results2 the differences between the three models (maxent_0.html, maxent_2.html, maxent_2.html) are very similar (likely due to large numer of occurrence records)
# what other way could you compare the performance of these three models, and how?

#######################
#######################
#######################
#######################
#######################
#######################
#######################
#######################
#######################
#######################

# Now let briefly explore how to account for sampling bias

bias_raster <- raster("data/Bias_LiPe_kd.asc")
crs(bias_raster) <- "+proj=longlat +datum=WGS84 +no_defs"  #always double check you have defined the CRS correctly
plot(bias_raster)  # and a good idea to double check that it looks like it should

bg <- xyFromCell(bias_raster, sample(which(!is.na(values(bias_raster))), 10000, prob=values(bias_raster)[!is.na(values(bias_raster))]))

# args = c('biasfile=bias_raster')  # some documentation indicates that you could use this argument to allow Maxent to correct for bias (looks like a different method though) I could not get this argument to work though

maxent_args <- c('removeduplicates=TRUE','jackknife=TRUE','responsecurves=TRUE', 'outputformat=cloglog','plots=TRUE')
Mod3_LiPe <- dismo::maxent(x=LiPe_predictors, p=xy.sp_LiPe,a=bg, path = "results3", args= maxent_args)

map_predictions3 <- predict(Mod3_LiPe, LiPe_predictors)

plot(map_predictions3)
points(xy.sp_LiPe,pch=20,cex=0.2)

# compare these results, what else might we do to compare results between models?

# there are some pretty big changes, and the model fitting does not give as high an AUC, Why? What if we had not used thinned data? Would the AUC have been higher? Why?

#######################
#######################
#######################
#######################
#######################
#######################
#######################
#######################
#######################
#######################

# Now lets try using targeted background points

#LiPe_zeros_thinned_1vis
#LiPe_thinned

zeros_1visit <- read.csv("data/LiPe_zero_locations_thinned_1vis.csv") # these are the locations where we assume there were no LiPe in that cell because at least one

x_bk2<- zeros_1visit$Long
y_bk2<-zeros_1visit$Lat
xy_bk2<-cbind(x_bk2,y_bk2)
xy.sp_bk2<-SpatialPoints(xy_bk2)
crs(xy.sp_bk2) <- crs(LiPe_predictors)
crs(xy.sp_bk2)

plot(LiPe_predictors[[2]])
points(xy.sp_bk2,pch=20,cex=0.2) #notice there has been alot of inland locations where frogs were surveyed for, but which did not record LiPe

maxent_args <- c('removeduplicates=TRUE','jackknife=TRUE','responsecurves=TRUE', 'outputformat=cloglog','plots=TRUE')
Mod4_LiPe <- dismo::maxent(x=LiPe_predictors, p=xy.sp_LiPe, a=xy.sp_bk2, path = "results4", args= maxent_args)

map_predictions4 <- predict(Mod4_LiPe, LiPe_predictors)

plot(map_predictions4)
points(xy.sp_LiPe,pch=20,cex=0.2)

# Notice our AUC values have fallen again and while the response curves look similar, variable importance has shifted.  There is a clear difference
# between the overall background, and a targeted background that only uses environmental variables where someone has done a frog survey
# should we consider excluding desert areas, where LiPe is not observed from our study area, or extent? Will we get more realistic estimates of
# model performance, and of the subtle differences between the areas LiPe occurs and where is does not if we focus our extent to only occur areas 
# where they might occur.  The map does not look radically different to the map predictions from model 1, but the low AUC is highligting how the
# it is hard to distinguish between presence locations and nearby locations where LiPe was not observed.  This is due to the high spatial autocorrelation
# in our predictor values and in our sampling effort.  This is a better AUC, and indicates there is much work to do to predict well.

# what if we predict our current model into the future?

setwd("~/Documents/Use_cases/SDM_in_R/predictors_future")
rast_lst2 <- list.files(pattern='.asc$', all.files=TRUE)
rast_lst2
LiPe_predictors_future <- stack(rast_lst2)
crs(LiPe_predictors_future) <- "+proj=longlat +datum=WGS84 +no_defs"
LiPe_predictors_future

map_predictions5 <- predict(Mod4_LiPe, LiPe_predictors_future)

plot(map_predictions5)
points(xy.sp_LiPe,pch=20,cex=0.2)

# any differences expected in the future map (slight difference with retraction away from inland and northern areas, possible suitability increase at higher elevations  
# How confident are we?  Was it a problem that we used one source of climate model to represent current climate
# and another source to represent future climate?

####################
####################
####################
####################
####################
####################
####################
####################
## Now we will fit a Boosted Regression Tree (BRT) model, which is essentially a generalized boosting model (gbm)

### First, lets add our data where no LiPe were found on three visits

LiPe_zero_3vis<-read.csv("data/LiPe_zero_locations_3vis.csv")
head(LiPe_zero_3vis)
LiPe_zero_3vis$X <- NULL
LiPe_zero_3vis$spp <- NULL
head(LiPe_zero_3vis)
LiPe_thin <- as.data.frame(cbind(Long = LiPe_thinned$decimalLongitude, Lat = LiPe_thinned$decimalLatitude, pres = LiPe_thinned$layer))
LiPe_thin2 <- LiPe_thin[complete.cases(LiPe_thin), ]
head(LiPe_thin2)
LiPe_p_a <- rbind(LiPe_zero_3vis,LiPe_thin2)

x_pa<- LiPe_p_a$Long
y_pa<-LiPe_p_a$Lat
xy_pa<-cbind(x_pa,y_pa)
xy.sp_pa<-SpatialPoints(xy_pa)
crs(xy.sp_pa) <- crs(LiPe_predictors)
crs(xy.sp_pa)

LiPe_preds <- extract(LiPe_predictors,xy.sp_pa)
LiPe_preds2 <- cbind(LiPe_p_a,LiPe_preds)
LiPe_preds3 <- as.data.frame(LiPe_preds2[complete.cases(LiPe_preds2), ])
summary(LiPe_preds3)
names(LiPe_preds3) # note in gbm.step you specify the column index location, so pres is indexed as 3 in names

# these models use stochastic learning, so you can expect a different answer each time unless you set the random seed, or make them
# deterministic by setting bag fraction to 1

LiPe_mod5_brt <- gbm.step(data = LiPe_preds3, gbm.x = c(4:9), gbm.y = 3, family = "bernoulli", tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.75)

# notice the plot that is produced automatically that shows the drop in holdout deviance as the number of trees increases, we are looking to explain more deviance with each tree
# ideally this curve would not decrease quickly at the start, there would be gentler curve.  Steep drops early in the model are often less stable models
# reducing the learning rate from 0.01 to 0.003 might lessen the drop, but this one is not too bad, and it would take much longer to run if we slow the learning rate down

# Elith et al. provides general guidelines on where to start, and suggests a variety of ways to select optimal settings.

# keep in mind a tree complexity of 3 suggests you suspect that three way interactions make sense ecologically. 
# perhaps the suitability of frog ponds relates to how hot it has been, what the soil moisture was like and how much rain occurred?

# also, simply running this model again, and comparing results will give you some indication of how stable your model is

# there are alot of interesting things stored in the model object, and you can browse through them

names(LiPe_mod5_brt)
LiPe_mod5_brt$gbm.call # this shows all the other default settings as well as those specified in the call

# lets predict and plot results, this step takes a couple hours
# this predict function will use the n.trees = to the best trees number found in the original model check it using:

pred_map_brt <- predict(LiPe_predictors,LiPe_mod5_brt,type = "response")
LiPe_mod5_brt$gbm.call$best.trees

# 3000 trees will be used, you can choose another value if you have a good reason too.

require("jpeg")
setwd("~/Documents/Use_cases/SDM_in_R") # make sure your current directory is your working directory, use getwd() to ensure you have asigned the directory correctly
jpeg(paste0(getwd(),"/results_brt/LiPe_brt_predicted.jpeg"))
plot(pred_map_brt)
points(xy.sp_LiPe,pch=20,cex=0.2)
dev.off()

# write raster
# writeRaster(pred_map_brt, "results_brt/LiPe_brt_predicted.asc")

# there are alot more options to explore with BRT: https://rspatial.org/raster/sdm/9_sdm_brt.html 
# curious what a function does

??gbm.plot

# gbm.plot generates response curves

jpeg(paste0(getwd(),"/results_brt/LiPe_brt_var_response_curves.jpeg"))
gbm.plot(LiPe_mod5_brt,rug=TRUE)
dev.off()

gbm.plot.fits(LiPe_mod5_brt)

gbm.interactions(LiPe_mod5_brt)

BRT_var_import<-as.data.frame(LiPe_mod5_brt$contributions)

# write.csv(BRT_var_import,"results_brt/LiPe_brt_var_importance.csv")

# var importance suggests two variables don't add much, lets look at possibly dropping those two

brt_drops <- gbm.simplify(LiPe_mod5_brt, n.drops = 4)

# if we look at the graph produced by gbm.simplify we can see a red verticle line at 2 variables removed 
# gbm uses boosting to focus on improving fit (deviance) on each iteration, but if it is focussing on a
# a variable that does not help, it actually increases deviance by having that extra predictor in the model
# A better predictive model might result from dropping 1 or 2 variables

names(LiPe_preds3)[brt_drops$pred.list[[1]]]

#compare that to the full list used in the model

names(LiPe_preds3)[c(4:9)]

# Not surprisingly this indicates that the variable EVI which contributes least can be dropped, Ecologically it might be a stretch
# to think that wetlands for frogs will be related to how green vegetation is in an area, so it might make sense to drop it
# Note few people worry about this step, but ecologists like simpler models, and keeping EVI actually increases predictive deviance, so dropping makes sense

names(LiPe_preds3)[brt_drops$pred.list[[2]]]

# if we were to drop two variables it would be EVI and wetland connectivity.  We know wetland connectivity includes salty wetlands which are no good for frogs
# and the layer does not include most of the small wetlands used by frogs, so ecologically might make sense to drop it, and a slight increase in predictive deviance
# when included.  Given my objective is to maximise predictions, I would drop both.

# a short cut for running another model with those two variables dropped is below, notice the diffent way to call gbm.x = 

# LiPe_mod6_brt <- gbm.step(data = LiPe_preds3, gbm.x = brt_drops$pred.list[[2]], gbm.y = 3, family = "bernoulli", tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.75)

# other plotting methods, evaluation, bootstrapping or ways to optimise the selection of the number of trees can be found here:
# https://rspatial.org/raster/sdm/9_sdm_brt.html
# Elith, J., Leathwick, J.R., and Hastie, T. (2008). Boosted regression trees - a new technique for modelling ecological data. Journal of Animal Ecology

####################
####################
####################
####################
####################
####################
####################
####################
## Now we will very quickly fit a GLM
# from the Mass package require(MASSExtra)

glm1<- glm(pres ~ AWAP + Bioclim05 + Bioclim06 + Bioclim12 + EVI+ wetland_connectivity, family = binomial, data = LiPe_preds3)

glm1
summary(glm1)

# this is a great new function which simplifies variable selection (remember the days of forward, backward) stepAIC is often good, stepBIC uses a different metric

glm1a <- step_BIC(glm1)

# this suggests our best model includes AWAP, Bioclim 05, and Bioclim06

# check residual plots, interpretation of residual plots for binomial family is different to other residual plot diagnostics, hear there appears to be little
# violation of linearity, and there are some influencial data outliers
#  from library(boot)

model_diags <- glm.diag(glm1a)    #residual diagnostics
glm.diag.plots(glm1a,model_diags)

# lets see how the predictions look

names(LiPe_predictors)
LiPe_preds_subset <- LiPe_predictors[[which(c(TRUE,TRUE,TRUE,FALSE,FALSE,FALSE))]]

glm_pred <- predict(LiPe_preds_subset,glm1a,type="response")

setwd("~/Documents/Use_cases/SDM_in_R") # make sure your current directory is your working directory, use getwd() to ensure you have asigned the directory correctly
jpeg(paste0(getwd(),"/results_glm/LiPe_glm_predicted.jpeg"))
plot(glm_pred)
points(xy.sp_LiPe,pch=20,cex=0.2)
dev.off()

####################
####################
####################
####################
####################
####################
####################
####################
## finally we will git a GAM

# from the mgcv package require(mgcv)
# a cr spline is defined by having knots spread throughout the co-variate values and which the amount of smooting is optimised, so it can become less wiggly

gam1<- gam(pres ~ s(AWAP, bs="cr") + s(Bioclim05, bs="cr") + s(Bioclim06,bs="cr"), family = binomial, data = LiPe_preds3)
gam.check(gam1)

# if we compare AIC values it looks like the GAM picks up some important non-linear relationships and is better than a GLM for this variable set anyway.
# BIC also suggests gam1 is better than either glm model

AIC(glm1a)
AIC(gam1)

plot(gam1)

gam_pred <- predict(LiPe_preds_subset,gam1,type="response")

# generate a jpeg map of the GAM predictions

setwd("~/Documents/Use_cases/SDM_in_R") # make sure your current directory is your working directory, use getwd() to ensure you have asigned the directory correctly
jpeg(paste0(getwd(),"/results_gam/LiPe_glm_predicted.jpeg"))
plot(gam_pred)
points(xy.sp_LiPe,pch=20,cex=0.2)
dev.off()

####################
####################
####################
####################
####################
####################
####################
####################
## list of Maxent args

#maxent_args <- c(

<!--   #duplicate records -->
  <!--   'removeduplicates=TRUE', #remove duplicate presence records. If environmental data are in grids, duplicates are records in the same grid cell, otherwise, duplicates are records with identical coordinates. -->

<!--   #background records -->
  <!--   'maximumbackground=10000', #if the number of background points/grid cells is larger than this number, then this number of cells is chosen randomly for background points. -->
<!--   'addsamplestobackground=TRUE', #add to the background any sample for which has a combination of environmental values that isn't already present in the background -->
<!--   'addallsamplestobackground=FALSE', #add all samples to the background, even if they have combinations of environmental values that are already present in the background -->

<!--   #missing data -->
  <!--   'allowpartialdata=FALSE', #during model training, allow use of samples that have nodata values for one or more environmental variables -->

<!--   #variable importance -->
  <!--   'jackknife=TRUE', #NB: default=FALSE; measure importance of each environmental variable by training with each environmental variable first omitted, then used in isolation. -->

<!--   #random seed -->
  <!--   'randomseed=FALSE', #if selected, a different random seed will be used for each run, so a different random test/train partition will be made and a different random subset of the background will be used, if applicable. -->

<!--   #prevalence -->
  <!--   'defaultprevalence=0.5', #default prevalence of the species: probability of presence at ordinary occurrence points. See Elith et al. Diversity and Distributions, 2011 for details -->

<!--   #train/test settings -->
  <!--   'randomtestpoints=0', #percentage of presence localities to be randomly set aside as test points, used to compute AUC, omission, etc. -->
<!--   'replicates=1', #number of replicate runs to do when cross-validating, bootstrapping or doing sampling with replacement runs. -->
<!--   'replicatetype=crossvalidate', #if replicates > 1, do multiple runs of this type:  -->
<!--   #crossvalidate: samples divided into replicates fods; each fold in turn used for test data -->
  <!--   #bootstrap: replicate sample sets chosen by sampling with replacement -->
  <!--   #subsample: replicate sample sets chosen by removing random test percentage without replacement to be used for evaluation -->
  <!--   'maximumiterations=500', #stop training after this many iterations of the optimization algorithm -->
<!--   'convergencethreshold=0.00001', #stop training when the drop in log loss per iteration drops below this number  -->

<!--   #feature selection -->
  <!--   'autofeature=TRUE', #automatically select which feature classes to use, based on number of training samples -->
<!--   'linear=TRUE', #allow linear features to be used -->
<!--   'quadratic=TRUE', #allow quadratic features to be used -->
<!--   'product=TRUE', #allow product features to be used -->
<!--   'threshold=TRUE', #allow threshold features to be used -->
<!--   'hinge=TRUE', #allow hinge features to be used -->

<!--   #feature settings -->
  <!--   'lq2lqptthreshold=80', #number of samples at which product and threshold features start being used -->
<!--   'l2lqthreshold=10', #number of samples at which quadratic features start being used -->
<!--   'hingethreshold=15', #number of samples at which hinge features start being used -->

<!--   #regularization settings -->
  <!--   'betamultiplier=1', #multiply all automatic regularization parameters by this number. A higher number gives a more spread-out distribution. -->
<!--   'beta_threshold=-1', #regularization parameter to be applied to all threshold features; negative value enables automatic setting -->
<!--   'beta_categorical=-1', #regularization parameter to be applied to all categorical features; negative value enables automatic setting -->
<!--   'beta_lqp=-1', #regularization parameter to be applied to all linear, quadratic and product features; negative value enables automatic setting -->
<!--   'beta_hinge=-1', #regularization parameter to be applied to all hinge features; negative value enables automatic setting -->

<!--   #outputs - NB. These are not shown in the UI, so unable to be changed by user -->
  
  <!--   'responsecurves=TRUE', #NB. default=FALSE; create graphs showing how predicted relative probability of occurrence depends on the value of each environmental variable -->
<!--   'responsecurvesexponent=FALSE', #instead of showing the logistic value for the y axis in response curves, show the exponent (a linear combination of features). -->
<!--   'pictures=TRUE', #create a .png image for each output grid -->
<!--   'outputformat=raw', #representation of probabilities used in writing output grids, see Help for details -->
<!--   'writeclampgrid=TRUE', #write a grid that shows the spatial distribution of clamping. At each point, the value is the absolute difference between prediction values with and without clamping. -->
<!--   'writemess=TRUE', #a multidimensional environmental similarity surface (MESS) shows where novel climate conditions exist in the projection layers. The analysis shows both the degree of novelness and the variable that is most out of range at each point. -->
<!--   'writeplotdata=FALSE', #write output files containing the data used to make response curves, for import into external plotting software. -->
<!--   'outputgrids=TRUE', #write output grids. Turning this off when doing replicate runs causes only the summary grids (average, std, deviation, etc) to be written, not those for the individual runs. -->
<!--   'plots=TRUE', #write various plots for inclusion in .html output -->
<!--   'logfile=maxent.log', #file name to be used for writing debugging information about a run in output directory -->
<!--   #'applythresholdrule=Fixed cumulative value 1', #apply a threshold rule, generating a binary outputgrid in addition to the regular prediction grid. Use the full name of the threshold rule in Maxent's html output as the argument. For example 'applyThresholdRule=Fixed cumulative value 1'. -->
  <!--   'logscale=TRUE', #if selected, all pictures of models will use a logarithmic scale for color-coding -->
<!--   'writebackgroundpredictions=FALSE', #write .csv file with predictions at background points -->
<!--   'fadebyclamping=FALSE', #reduce prediction at each point in projections by the difference between clamped and non-clamped output at that point -->

<!--   #projection settings NB. These are not shown in the UI, so unable to be changed by user -->
  
  <!--   'extrapolate=TRUE', #predict to regions of environmental space outside the limits encountered during training -->
<!--   'doclamp=TRUE' #apply clamping when projecting -->

<!-- #) -->