### Author: Maria Petrova
# Date: May, 2014
# Purpose: Boundary-line analysis of the data with a moving 'window' 
# Bootstrap on the boundary line to access robustness of the line


###################################################################################################
###################################################################################################
#STEP #1: ADMINISTRATIVE TASKS

#Libraries
library(ggplot2)
library(xlsx)
library(reshape)
library(plyr)
library(boot)
library(nlme)

### START data import and preparation ###
#setwd('C:\\Users\\Mirra\\Dropbox\\Flatwoods Community Ecology\\')
#Not needed for read.table()

###################################################################################################
###################################################################################################
#STEP #2: OPEN AND ADJUST BIOMASS DATA


#########################################################
#2A The steps below create a more appropriate dataset by combining some categories.


### Data - Biomass
plant <- read.table(
  "C:/usfs_sef_data_output/sef_Ecology_BiomassPlotMatrix_Ep1_OriginalOulierX_2014-05-05_15.00.18.csv", 
  header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE,
  stringsAsFactors = F)
str(plant)

#Convert file name
plotBiomass <- plant

#########################################################
#2B The steps below create a more appropriate dataset by combining some categories.

#Combine like categories
#Dead forb and live forb > these are being combined because some sites were sampled in 
#winter and some in summer so live or dead forb loadings are an indication of season
#rather then site differences.
forb <- plotBiomass$live.forb + plotBiomass$dead.forb

#Dead shrub and dead tree
#Dead standing woody should be one category
dead.woody <- plotBiomass$dead.shrub + plotBiomass$dead.tree

#Dead palmetto and live palmetto
#Same reasoning for combining live and dead forbs
palmetto <- plotBiomass$live.palmetto + plotBiomass$dead.palmetto

#Bluestem and grasses
#The grass will be difficult to interpret in the analysis because it contains a mixture of 
#primarily rhizomatous grasses but also bunch grasses that were not wiregrass or bluestem.

#Remove original categories that were combined into new categories in the steps above.
plotBiomass2 <- plotBiomass[,!colnames(plotBiomass) %in% c("live.forb", "dead.forb", 
                                                           "dead.shrub", "dead.tree", 
                                                           "live.palmetto", 
                                                           "dead.palmetto")]

#Add in new category columns
plotBiomass3 <- data.frame(plotBiomass2, forb = forb, dead.woody = dead.woody, palmetto = palmetto)

###################################################################################################
###################################################################################################
#STEP #3: OPEN AND ADJUST COVER DATA


### Data - Cover
plotCover <- read.table(
  "C:/usfs_sef_data_output/sef_Ecology_CoverPlotMatrix_Ep1_OriginalOulierX_2014-08-21_17.24.49.csv", 
  header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE,
  stringsAsFactors = F)
str(plotCover)

###################################################################################################
###################################################################################################
#STEP #4: OPEN AND ADJUST ENVIRONMENTAL DATA

#########################################################
#4a: Open environmental matrix

### Environmental matrix
env <- read.table(
  "C:/usfs_sef_data_output/2014.03.13_EnvironmentalMatrix.csv", 
  header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE,
  stringsAsFactors = F)
str(env)

#########################################################
#4b:subset to only have the variables needed for modeling - FireRotation
env <- env[c("Site", "FireRotation")]

###################################################################################################
###################################################################################################
#STEP #5: MERGE DATA

#Do you want to look at biomass or cover data
plant <- plotBiomass3#biomass
plant <- plotCover#

### merge the Environmental and species (cover or biomass) datasets 
data <- merge(env, plant, by.y = 'siteName', by.x = 'Site'); dim(data); str(data)

###################################################################################################
###################################################################################################
#STEP #6: BOUNDARY LINE REGRESSION


#########################################################
#6a: Create a boundary line regression function
boundaryLine <- function(data, species)
{
  
  #subset the data to contain only FireRotation and the specified species
  data <- data[ c('FireRotation', species) ]
  # for each unique predictor value, point falling within a certain range were subset
  # and the max identified
  data$window <- data$FireRotation + 0.5 
  data$max <- 0
  for (i in 1:nrow(data))
  data$max[i] <- max( data[which(data$FireRotation <= data$window[i]), 'FireRotation'] )
  
  

  ## fit the boundary line(log model) to the max points
  # 1. for the plot rename columns of data with x and y
  # 2. same for the maxpairs dataset consisting of max for each window
  tempdf <- data.frame(x = data['FireRotation'], y = data[species])
  names(tempdf) <- c('x', 'y') ## names in the datasets have to be the same or errors
  maxPairs <- data.frame(x = data['max'], y = data[species])
  names(maxPairs) <- c('x', 'y') ## errors if names are not x, y in the formula and data
  theme_set(theme_bw())
  p <- ggplot(tempdf, aes(x, y)) + geom_point() + xlab("Fire Rotation (years)") + ylab(species)
  p <- p + stat_smooth(method="lm", data = maxPairs, formula = y ~ log(x), se = T, 
                       colour="black", fullrange = T, fill="gray91")
  
  print(p)
  #return (lm(maxPairs[,'y'] ~ log(maxPairs[,'x']) )$coef) # linear -log model
  return(maxPairs)
  
}

#boundaryLine(data, 'forb')
#model <- lm(tempdf[,'y'] ~ log(tempdf[,'x']) ) # linear -log model
#summary(model)

#J Cronan >> code below is testing linear mixed effects model which is needed
#to deal with pseudoreplication (plots) and probably also regional differences.
#model <- lme(y ~ log(x), random = ~1|Site, data = templme)
#templme <- data.frame(x = data['x'], y = data['y'], site = data['Site'])
#head(templme)
#summary(model)

#########################################################
#6b: Create a function for bootstrapping data.

### Bootstrapping - resample the original data with replacement to create new datasets with 
### same number of observations
### same as the boundary line function but without plotting 
### fit a boundary line to each set 
boundaryLineBoot <- function(data, species, indices)
{
  
  # for loop that creats a group label for each level of FireRotation 
  data <- data[indices, ] # indices is the random indexes for the bootstrap sample
  data <- data[ c('FireRotation', species) ]
  data$group <- NA
  for (i in 3:9)
    data$group[which(data$FireRotation > i - 1 & data$FireRotation < i)] <- paste('group',
                                                                                  i-2, sep = '')
  
  
  sortnames <- c('group', species) 
  dataSorted <- data[do.call('order', data[sortnames]), ] ## do.call to sort over two columns
  dataSorted <- dataSorted[, c('group', 'FireRotation',species)]
  
  # cumsum is the cumulative sum over the count of elements in each group
  # returns the index of the max values pair since data is sorted based on both columns
  maxPairs <- dataSorted[cumsum(count(dataSorted, 'group')$freq), c('FireRotation', species)]
  
  ## fit the boundary line(log model) to the max points
  #x = maxdata.bootairs[ ,'FireRotation']; y = maxPairs[ , species]; d <- data.frame(x,y)
  #return(lm(maxPairs[,species] ~ log(maxPairs[,'FireRotation']), data = maxPairs)$coef[2])  
  # return the beta coefficient
  loglin <- lm(maxPairs[,species] ~ log(maxPairs[,'FireRotation']) )
  return(coef(loglin)) # linear-log model
}

  


### Fake data to test CI
data <- merge(env, plant, by.y = 'siteName', by.x = 'Site'); dim(data); str(data)
test <- rbind(data, data, data, data) # to make a bigger dataset
test[, 'forb'] <- jitter(test[,'forb']) # Add a small amount of noise to species data
test$forb[which(test$forb < 0)] <- 0 
hist(data$forb);hist(test$forb)


# bootstrap
# data: data frame to which bootstrap resampling is to be applied.
# statistic: A function that returns the (possibly vector-valued) statistic to be bootstrapped. 
# R is how many bootstrap samples
data.boot = boot(data=test, statistic=boundaryLineBoot, species = 'forb', R=3000) 

data.boot
### average prediction and 95% confidence interval were estimated from the results and plottted
CI.inter <- quantile(data.boot$t[,1], c(.025, .975) )#0.05, 0.95 for 90% CI
CI.slope <- quantile(data.boot$t[,2], c(.025, .975) )#0.025, 0.975 for 90% CI
mean(data.boot$t[,1]); mean(data.boot$t[,2])
boundaryLine(data, 'forb')


hist(data.boot$t[,1], breaks = 100, main = 'Histogram of bootstrap resluts\n(2500 iterations)'
     , xlab = 'Intercept values', col = 'gray91')
abline(v=CI.inter[[1]][1], lwd = 2 , lty = 3) # lower intercept
abline(v=CI.inter[[2]][1], lwd = 2 , lty = 3) # upper intercept

hist(data.boot$t[,2], breaks = 100, main = 'Histogram of bootstrap results\n (2500 iterations)'
     , xlab = 'Slope values', col = 'gray91')
abline(v=CI.slope[[1]][1], lwd = 2 , lty = 3) # lower slope value
abline(v=CI.slope[[2]][1], lwd = 2 , lty = 3) # upper slope value




## plot the 95%confidence interval
maxPairs <- boundaryLine(data, 'forb')
maxPairsLow <- boundaryLine(data, 'forb')## data to plot the lower CI limits
maxPairsUp <- boundaryLine(data, 'forb')## data to plot the upper CI limits
## recalculate y based in the CI instercept and slope;
## formula is y = intersect + slope*log(x)
maxPairsLow$y = CI.inter[[1]][1] + CI.slope[[1]][1]*log(maxPairs$x)
maxPairsUp$y= CI.inter[[2]][1] + CI.slope[[2]][1]*log(maxPairs$x)
names(data)[which(names(data) == 'FireRotation' | names(data) == 'forb')] <- c('x','y')

p <- ggplot(data, aes(x, y)) + geom_point() + xlab('Fire Rotation (years)') + ylab("Forb Loading")
#p + stat_smooth(method="lm", data = maxPairs, formula = y ~ log(x), geom = "point") + 
#    stat_smooth(method="lm", data = maxPairs, formula = y ~ log(x), geom = "errorbar")

p <- p + stat_smooth(method="lm", data = maxPairs, formula = y ~ log(x), se = T, 
                     colour="black", fullrange = T, fill="gray91")

p <- p + stat_smooth(method="lm", data = maxPairsLow, formula = y ~ log(x), se = F, 
                colour="black", fullrange = T, linetype="dashed")

p + stat_smooth(method="lm", data = maxPairsUp, formula = y ~ log(x), se = F, 
                colour="black", fullrange = T, linetype="dashed")



head(exp, data$x,)




### plot log curve with base graphing package
x = maxPairs[ ,'FireRotation']; y = maxPairs[ , species]; d <- data.frame(x,y)
logmodel <- lm(y ~ log(x),data=d)
plot(dataSorted[ ,'FireRotation'], dataSorted[ , species]
     , xlab = "FireRotation", ylab = species, pch = 16)

xvec <- seq(2,8, length=nrow(data)*3) ### fake vector ro smooth out the curve
logpred <- predict(logmodel, newdata=data.frame(x=xvec))
lines(xvec,logpred)

