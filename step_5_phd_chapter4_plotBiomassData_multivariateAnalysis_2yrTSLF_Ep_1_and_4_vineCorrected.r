###################################################################################################
############################-----start-----########################################################
#The purpose of this script is to explore the environmental data for the community ecology
#multivariate analysis.

#NOTE
#This script was modified in October, 2022 and name was changed.
#Script was maintained on GitHub (LongleafCommunityEcology repository) and previous 
#modification occurred on 1-Jan-2020. Script name was changed from: 
#sef_2014.09.25_Analysis_OriginalData_MeadowOutliersX_mod
#to:
#

#Which computer are you using?
#Forest Service>>>>>>> FS
#Personal >>>>>>>>>>>> JC
computer <- "JC"

###################################################################################################
###################################################################################################
#STEP #1: ADMINISTRATIVE TASKS

#Reset functions
#rm(list=ls())#do not run this function b/c it will erase biostats script
dev.off()

#Libraries
library(Hmisc)
library(PerformanceAnalytics)
library(vegan)
library(pastecs)
#library(simba) -- package no longer maintained on CRAN (Cronan - 2022-10-11)
library(fields)#for set.panel()
library(labdsv)#Brooke's recommendation


#Set working directory
if(computer == "FS")
{
  setwd("C:/Users/jcronan/Box/01. james.cronan Workspace/Research/UW_PHD/Dissertation/4_Chapter_4/Data")
} else
{
  setwd("C:/Users/james/Box/01. james.cronan Workspace/Research/UW_PHD/Dissertation/4_Chapter_4/Data")
}

###################################################################################################
###################################################################################################
#STEP #2: OPEN AND ADJUST ENVIRONMENTAL DATA

#########################################################
#2a: Open environmental matrix

#Open data
siteEnv <- read.table(paste("Fire_History/", "phd_chapter4_environmentalData.csv", sep = ""), 
                      header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE, stringsAsFactors = F, fill = T)

#########################################################
#2b: Re-assign column 1 (site names) to  a row name.
#Also remove column 2 (site numbers)
siteEnv2 <- siteEnv[,3:(length(siteEnv[1,]))]
rownames(siteEnv2) <- siteEnv[,1]

#########################################################
#2c: Environmental data, remove data you will not analyze
siteEnv3 <- siteEnv2[,-c(11,12)]
#Column 11 is soil drainage. This data is from USDA Soil Survey and not accurate
#Column 12 is region. Only two regions, at this point I don't believe this adds 
#substance to the analysis because there are no physical properties for each region
#that can confound how fire characteristics drive understory plant composition
#siteEnv3 <- siteEnv2[,-c(11)]
#maintain region column for restricted permutation tests

###################################################################################################
#STEP 3: CREATE A BOUNDARY LAYER REGRESSION FUNCTION for Fire Rotation

###################################################################################################
############################################FUNCTION###############################################

#Create a function that will conduct a boundary layer regression with y variable as exponent,
#print summary of regression, and plot data
blr <- function(a,x,y,z)
{
  cf <- 0.8
  mv <- vector(mode = "numeric")
  fv <- vector(mode = "numeric")
  le <- vector(mode = "numeric")
  
  fri <- 2:6
  for(i in fri)
  {
    mv[i-1] <- max(y[x > i & x < i + 1])
    fv[i-1] <- min(x[x > i & x < i + 1 & y == mv[i-1]])
    le[i-1] <- length(y[x > i & x < i + 1])
    
    mv[i-1] <- mv[i-1] + 0.0001
  }
  
  x2 <- fv
  y2 <- mv  
  
  d <- data.frame(x2,y2)
  logmodel <- lm(y2~a(x2),data=d)
  
  
  ### fake vector
  xvec <- seq(0,8, length=101)
  logpred <- predict(logmodel, newdata=data.frame(x2=xvec))
  
  #Plot with boundary layer regression
  plot(x, y, xlab = "Fire Rotation (years)", ylab = "biomass (Mg/ha)", main = z)
  points(fv, mv, pch = 22)
  lines(xvec,logpred)
  text(4.1, (max(y) - max(y)/16), paste("Intercept", round(logmodel$coefficients[1],4)), cex = cf)
  text(5.9, (max(y) - max(y)/16), paste("Slope", round(logmodel$coefficients[2],4)), cex = cf)
  ms <- summary(logmodel)
  text(4.1, (max(y) - max(y)/8), paste("r-squared", round(ms$r.squared,2)), cex = cf)
  text(5.9, (max(y) - max(y)/8), paste("P-value", round(ms$coefficients[8],4)), cex = cf)
  print(summary(logmodel))
}

###################################################################################################
############################################FUNCTION###############################################

###################################################################################################
###################################################################################################
#STEP #2: OPEN AND ADJUST BIOMASS DATA

#########################################################
#2A: Open biomass data

#Plot-level biomass data. This data has been modified so that vine species are consistent
#across all sites.
#Outlier sites have not been removed.
plotBiomass <- read.table(paste("Understory_Vegetation_FlatFiles/stage_4_aggregate/outputs/", 
                                "phd_chapter4_biomassPlotMatrix_2yr_Ep_1_4_Species_2022-10-17_15.34.40.csv",
                                sep = ""), header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE,
                          stringsAsFactors = F)

#########################################################
#2C: calculate site means for untransformed plot-level data.
#USED FOR DATA SCREENING ONLY - DATA NOT USED IN ANALYSIS.

#Calculate untransformed site means from plot-level data 
siteBiomass <- summarize(X = plotBiomass[,4:length(plotBiomass[1,])], by = plotBiomass$siteName, 
                         colMeans, stat.name = colnames(plotBiomass)[4])
colnames(siteBiomass)[1] <- "siteName"

#Convert object back into array format
siteBiomass2 <- as.data.frame.array(siteBiomass)

#Remove site names from matrix columns and assign as row names.
#All biomass >> UNTRANSFORMED
siteBiomass3 <- siteBiomass2[,2:(length(siteBiomass2[1,]))]
rownames(siteBiomass3) <- siteBiomass2[,1]

#########################################################
#2d: Conduct log transformation on plot-level data and then calculate site means.

#Conduct log-transformation on plot-level data
#>>>>Justification for log-transformation?
plotBiomassLogTrans <- data.trans(plotBiomass[,4:length(plotBiomass[1,])], method = 'log', 
                                  plot = F)

#Calculate site means from log-transformed plot-level data (USED IN ANALYSIS)
siteBiomassLogTrans <- summarize(X = plotBiomassLogTrans, by = plotBiomass$siteName, colMeans, 
                                 stat.name = colnames(plotBiomass)[4])
colnames(siteBiomassLogTrans)[1] <- "siteName"


#Convert object back into array format
siteBiomassLogTrans2 <- as.data.frame.array(siteBiomassLogTrans)

#Remove site names from matrix columns and assign as row names.
#All biomass >> LOG-TRANSFORMED
siteBiomassLogTrans3 <- siteBiomassLogTrans2[,2:(length(siteBiomassLogTrans2[1,]))]
rownames(siteBiomassLogTrans3) <- siteBiomassLogTrans2[,1]

###################################################################################################
###################################################################################################
#STEP #5: DATASET ADJUSTMENTS

###################################################################################################
#5a Standardize row order of data sets (i.e. order site names (row names) of biological data 
#according to order in environmental data set) and remove columns with no values > 0. 

#Reorder biological data (biomass) so sites are in same order as environmental data
#Step not necessary for cover data... sites are already in same order as environmental data.
ro <- data.frame(siteNo = 1:length(rownames(siteEnv3)), siteName = I(rownames(siteEnv3)))
ro2 <- ro[order(ro[,2]),]

#Untranformed biomass data
siteBiomass4 <- data.frame(so = ro2[,1], siteBiomass3)
siteBiomass5 <- siteBiomass4[order(siteBiomass4$so),]
siteBiomass6 <- siteBiomass5[,!colnames(siteBiomass5) %in% c("so")]

#Log-tranformed biomass data
siteBiomassLogTrans4 <- data.frame(so = ro2[,1], siteBiomassLogTrans3)
siteBiomassLogTrans5 <- siteBiomassLogTrans4[order(siteBiomassLogTrans4$so),]
siteBiomassLogTrans6 <- siteBiomassLogTrans5[,!colnames(siteBiomassLogTrans5) %in% c("so")]

#Remove taxa with no occurences, this is necessary before you can use foa.plots below
biomassOrig <- drop.var(siteBiomass5, min.fo = 1)
biomassLog <- drop.var(siteBiomassLogTrans5, min.fo = 1)
#Show taxa that were dropped:
cbind(colnames(siteBiomass5), match(colnames(siteBiomass5), colnames(biomassOrig)))
#dropped seven genus. From 40 to 33.
#dropped eight 'genus'species. From 50 to 42.

###################################################################################################
###################################################################################################
#STEP #6: DATA SCREENING

#########################################################
#6a: Summary stats for biomass data

#Summary stats
round(stat.desc(biomassOrig),2)

#Look at how data is distributed
foa.plots(biomassOrig)

#Cum. Number of Species (aka Genus) vs Frequency of Occurrence
#>> There are 33 genus of plants. Half (17) of those occur on 14 of the sites
#   The other half are uncommon and occur at less than 7 (total n = 21) sites.

#Cum. Number of Species vs Percent Occurrence
#Not sure what "percent occurrence" (y-axis) means. How is this different than the prior chart?
#>> Similar trend as above
#   21 genus occur less than 50%
#   12 genus have occurence greater than 70, with all but one higher than 90%

#Histogram of Species Occurrence
#>> Bimodal with largest peak in the 0-5 sites range and a smaller, trough in the 10-15 
#   sites range and secondary peak in the 20-25 sites range.

#Histogram of Log(Species Occurrence)
#>> Not sure what this histogram is displaying but distrubution is more or less flat.

#Cumulative Distribution of Species Mean Abundance
#>> Exponentially increasing
# This shows that about 5 genus are dominant (genus on steep curve)
# and the remaining 28 genus have low biomass.

#Species Occurrence vs Mean Abundance
#Genus with highest abundance have higher occurrence, but there are also many genus with high occurrence
#and low abundance.

#Species Occurrence vs Log(Mean Abunandance)
#>> There is a weak (probably significant) upward trend of increasing abundance with frequency
#   of occurrence.

#Cumulative Distribution of Plot (Site) Richness.
#>> Site richness doubles from 12 to 21 and increases on a steady slope across all sites. I.e., there are no
#distinct groups of low and high richness plots, but rather a gradual transition from low to high.

#Cumulative Distribution of Plot Total Abundance
#>> There is a linear upward trend. Not sure what abundance is, loading?
#   Looks like it is showing total loading by site arranged from lowest to highest
#   but given the chart title (cumulative) I believe I am misreading both this, and the
#   above chart.

#Plot Richness vs Total Abundance
#>> There is a weak and possibly insignificant upward trend.
#   Plots along the line from low abundance/low richness to high:
#   A32, E508A, E505, E807D, E807B, E100B-E
#   Outliers include:
#   Low plot richness and high abundance: S209
#   Medium plot richness and low abundance: E100B-W

#Drop rare species < 50% of sites.
spo <- drop.var(biomassOrig, min.po=5)
spl <- drop.var(biomassLog, min.po=5)

#Check to see how many species were dropped for untransformed biomass data
length(biomassOrig[1,])#
length(spo[1,])#
length(biomassLog[1,])#
length(spl[1,])
#42 to 13 variables (29 drops) for species-level data

str(spo)#22 obs and 13 variables
str(spl)#22 obs and 13 variables

#Convert data.frame to a matrix (needed for uv.plots)
so2 <- data.matrix(frame = spo, rownames.force = NA)
sl2 <- data.matrix(frame = spl, rownames.force = NA)

#uv.plots displays histogram, box and whisker, cumulative distribution and normal q-q plots
#in one pane for each variable.
uv.plots(so2)#Many species are right skewed with a long right tail. To reduce the effect of these outliers use log-transformed data.
uv.plots(sl2)#Distributions are not skewed when log transformation is applied. Use this data for PCA.

#What percent of values are zero, if it is over 50% you should consider changing the data
#to presence/absence
length(sl2[sl2 == 0])/(length(sl2[1,])*length(sl2[,1]))
#6% zero values >>> no need to convert data to presence/absence

#Look at correlation between species variables.
#generate table (you need to run correlation_matrix() function for this to work:
#https://www.r-bloggers.com/2020/07/create-a-publication-ready-correlation-matrix-with-significance-levels-in-r/
scm <- correlation_matrix(as.data.frame(so2), type = "pearson", show_significance = T, 
                         digits = 2, use = "lower", replace_diagonal = T)
chart.Correlation(so2, method = "pearson")#performanceAlanlytics
#there is only one correlation significant at the P < 0.05 level
#St. Andrews Cross and Darrow's blueberry

#########################################################
#6c: Summary stats for environmental data

#Show how environmental variables are correlated.
chart.Correlation(siteEnv3)
ecm <- correlation_matrix(as.data.frame(siteEnv3), type = "pearson", show_significance = T, 
                          digits = 2, use = "lower", replace_diagonal = T)

###################################################################################################
###################################################################################################
#STEP #7: PRELIMINARY ANALYSIS

#########################################################
#7a: Determine species with the highest and second highest biomass at each site
#and then display occurrence of dominant and co-dominant species with a histogram.

#Create a histogram of primary and secondary dominance across sites by species.

#Remove non-living categories (dead woody) from the biomass data
last_col <- length(biomassOrig[1,])
nextLast_col <- length(biomassOrig[1,]) - 1

#List of dominant species
prim <- mapply(function(x) {colnames(biomassOrig)[order(biomassOrig[x,])][last_col]}, 
               1:length(biomassOrig[,1]))
seco <- mapply(function(x) {colnames(biomassOrig)[order(biomassOrig[x,])][nextLast_col]}, 
               1:length(biomassOrig[,1]))

dTable <- data.frame(Sites = I(rownames(biomassOrig)), Primary = I(prim), Secondary = I(seco))

udom <- unique(c(unique(dTable[,2]),unique(dTable[,3])))

hprim <- mapply(function(x) {length(dTable[,2][dTable[,2] == udom[x]])}, 1:length(udom))
hseco <- mapply(function(x) {length(dTable[,3][dTable[,3] == udom[x]])}, 1:length(udom))

hdom <- data.frame(Species = I(udom), Primary = hprim, Secondary = hseco)
h2dom <- hdom[order(hdom[,2], decreasing = T),]

par(mai = c(2.5,1.2,1,1))
barplot(t(cbind(h2dom$Primary,h2dom$Secondary)), main = "", 
        xaxt = "n", ylab = "Number of sites", xlab = "", axes = F, beside = T, 
        col = c("white", "dark grey"))
axis(2)
text(matrix(c(seq(2,nrow(h2dom)*3,3), rep(-0.25,nrow(h2dom))), nrow = nrow(h2dom), 
            ncol = 2, byrow = F), srt = 60, adj = 1, xpd = T, labels = paste(h2dom$Species), 
     cex = 0.95)

legend(15,6, c("Highest biomass", "Second highest biomass"), fill = c("white", "dark grey"), bty = "n")

#########################################################
#7b: Arrange species on bar chart with highest to lowest average biomass

#Make a table showing mean/sd biomass
bm <- apply(biomassOrig,2,mean)
bs <- apply(biomassOrig,2,sd)

biomass <- data.frame(Mean = round(bm,2), SD = round(bs,2))

###################################################################################################
###################################################################################################
#STEP #8: Gradient length diagnostics with DCA for untransformed biomass data, log-transformed
#biomass data, and untransformed cover data.

#Create a bindary dataset for species presence absence
#Biomass - untransformed
speocc <- data.trans(so2, method = "power", exp = 0, plot = F)

decorana(speocc, ira = 0)

#Biomass - log-transformed
speocc <- data.trans(sl2, method = "power", exp = 0, plot = F)

decorana(speocc, ira = 0)

#DCA 1 axis length is 0.59714
#This is way below the 3-4 cut-off proposed by ter Braak and Smilauer (2002) used as a threshold 
#for selecting a multivariate analysis technique (< 3 PCA/RDA; data are suited for linear method, and
#> 4 CA/CCA; data are suited for unimodal methods).

#This value indicates a short gradient length (i.e., use linear analysis). This is consistent 
#with sites selection which focused on a narrow set of environmental conditions. 
#I.e. mesic flatwoods with a history of management with frequent rx fire and no other evidence of 
#major disturbance.

###################################################################################################
###################################################################################################
#STEP #9: PCA of the covariance matrix
pca <- rda(decostand(sl2, method = "hellinger", scale = T))#scale = T should be PCA of the correlation matrix
pca
#Plot pca
biplot(pca, scaling = 1)
#Use brokan stick method to determine how many PCA axes to retain
screeplot(pca, bstick = T, type = "l", main = NULL)
#Extract eigenvalues
eigenvals(pca)
summary(eigenvals(pca))

#Fitting environmental data
set.seed(42)#use this to make this permutation analysis reproducibe, otherwise it will be different
#every time, especially when number of permutations is low.
ev <- envfit(pca ~ ., data = siteEnv3, choices = 1:2, scaling = 'symmetric', permutations = 1000)
#options for scaling: "species", "sites", "symmetric" (scales both for species and sites), "none" (raw scores).
#cca() hill = T, rda(): correlation = T (standardizes scores among species).
ev

plot(pca, display = "sites", type = "n", scaling = "symmetric")
points(pca, display = "sites", scaling = "symmetric")
plot(ev, add = T)

surf <- ordisurf(pca ~ FineWD, data = siteEnv3, knots = 10, isotropic = T, main = NULL)
surf <- ordisurf(pca ~ CoarseWD, data = siteEnv3, knots = 10, isotropic = T, main = NULL)
surf <- ordisurf(pca ~ Litter, data = siteEnv3, knots = 10, isotropic = T, main = NULL)
surf <- ordisurf(pca ~ mfri_20yr, data = siteEnv3, knots = 10, isotropic = T, main = NULL)
summary(surf)



###################################################################################################
###################################################################################################
#STEP #2: OPEN AND ADJUST BIOMASS DATA

#########################################################
#2A: Open biomass data

#Plot-level biomass data. This data has been modified so that vine species are consistent
#across all sites.
#Outlier sites have not been removed.
plotBiomass <- read.table(paste("Understory_Vegetation_FlatFiles/stage_4_aggregate/outputs/", 
                                "phd_chapter4_biomassPlotMatrix_2yr_Ep_1_4_FuncGroup_2022-10-19_12.47.54.csv",
                                sep = ""), header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE,
                          stringsAsFactors = F)

#########################################################
#2C: calculate site means for untransformed plot-level data.
#USED FOR DATA SCREENING ONLY - DATA NOT USED IN ANALYSIS.

#Calculate untransformed site means from plot-level data 
siteBiomass <- summarize(X = plotBiomass[,4:length(plotBiomass[1,])], by = plotBiomass$siteName, 
                         colMeans, stat.name = colnames(plotBiomass)[4])
colnames(siteBiomass)[1] <- "siteName"

#Convert object back into array format
siteBiomass2 <- as.data.frame.array(siteBiomass)

#Remove site names from matrix columns and assign as row names.
#All biomass >> UNTRANSFORMED
siteBiomass3 <- siteBiomass2[,2:(length(siteBiomass2[1,]))]
rownames(siteBiomass3) <- siteBiomass2[,1]

#########################################################
#2d: Conduct log transformation on plot-level data and then calculate site means.

#Conduct log-transformation on plot-level data
#>>>>Justification for log-transformation?
plotBiomassLogTrans <- data.trans(plotBiomass[,4:length(plotBiomass[1,])], method = 'log', 
                                  plot = F)

#Calculate site means from log-transformed plot-level data (USED IN ANALYSIS)
siteBiomassLogTrans <- summarize(X = plotBiomassLogTrans, by = plotBiomass$siteName, colMeans, 
                                 stat.name = colnames(plotBiomass)[4])
colnames(siteBiomassLogTrans)[1] <- "siteName"


#Convert object back into array format
siteBiomassLogTrans2 <- as.data.frame.array(siteBiomassLogTrans)

#Remove site names from matrix columns and assign as row names.
#All biomass >> LOG-TRANSFORMED
siteBiomassLogTrans3 <- siteBiomassLogTrans2[,2:(length(siteBiomassLogTrans2[1,]))]
rownames(siteBiomassLogTrans3) <- siteBiomassLogTrans2[,1]

###################################################################################################
###################################################################################################
#STEP #5: DATASET ADJUSTMENTS

###################################################################################################
#5a Standardize row order of data sets (i.e. order site names (row names) of biological data 
#according to order in environmental data set) and remove columns with no values > 0. 

#Reorder biological data (biomass) so sites are in same order as environmental data
#Step not necessary for cover data... sites are already in same order as environmental data.
ro <- data.frame(siteNo = 1:length(rownames(siteEnv3)), siteName = I(rownames(siteEnv3)))
ro2 <- ro[order(ro[,2]),]

#Untranformed biomass data
siteBiomass4 <- data.frame(so = ro2[,1], siteBiomass3)
siteBiomass5 <- siteBiomass4[order(siteBiomass4$so),]
siteBiomass6 <- siteBiomass5[,!colnames(siteBiomass5) %in% c("so")]

#Log-tranformed biomass data
siteBiomassLogTrans4 <- data.frame(so = ro2[,1], siteBiomassLogTrans3)
siteBiomassLogTrans5 <- siteBiomassLogTrans4[order(siteBiomassLogTrans4$so),]
siteBiomassLogTrans6 <- siteBiomassLogTrans5[,!colnames(siteBiomassLogTrans5) %in% c("so")]

#Remove taxa with no occurences, this is necessary before you can use foa.plots below
biomassOrig <- drop.var(siteBiomass6, min.fo = 1)
biomassLog <- drop.var(siteBiomassLogTrans6, min.fo = 1)
#Show taxa that were dropped:
cbind(colnames(siteBiomass6), match(colnames(siteBiomass6), colnames(biomassOrig)))
#does not occur with functional group data. No drops

###################################################################################################
###################################################################################################
#STEP #6: DATA SCREENING

#########################################################
#6a: Summary stats for biomass data

#Summary stats
round(stat.desc(biomassOrig),2)

#Look at how data is distributed
foa.plots(biomassOrig)
dev.off()

#Drop rare species < 10% of sites.
spo <- drop.var(biomassOrig, min.po=10)
spl <- drop.var(biomassLog, min.po=10)

#Check to see how many species were dropped for untransformed biomass data
length(biomassOrig[1,])#
length(spo[1,])#
length(biomassLog[1,])#
length(spl[1,])
#no drops

str(spo)#21 obs and 13 variables
str(spl)#21 obs and 13 variables

#Convert data.frame to a matrix (needed for uv.plots)
so2 <- data.matrix(frame = spo, rownames.force = NA)
sl2 <- data.matrix(frame = spl, rownames.force = NA)

#uv.plots displays histogram, box and whisker, cumulative distribution and normal q-q plots
#in one pane for each variable.
uv.plots(so2)#Many species are right skewed with a long right tail. To reduce the effect of these outliers use log-transformed data.
uv.plots(sl2)#Distributions are not skewed when log transformation is applied. Use this data for PCA.

#What percent of values are zero, if it is over 50% you should consider changing the data
#to presence/absence
length(sl2[sl2 == 0])/(length(sl2[1,])*length(sl2[,1]))
#8% zero values >>> no need to convert data to presence/absence

#Look at correlation between species variables.
#generate table (you need to run correlation_matrix() function for this to work:
#https://www.r-bloggers.com/2020/07/create-a-publication-ready-correlation-matrix-with-significance-levels-in-r/
scm <- correlation_matrix(as.data.frame(so2), type = "pearson", show_significance = T, 
                          digits = 2, use = "lower", replace_diagonal = T)
chart.Correlation(so2, method = "pearson")#performanceAlanlytics
#there is only one correlation significant at the P < 0.05 level
#understory trees and vines

#########################################################
#6c: Summary stats for environmental data

#Show how environmental variables are correlated.
chart.Correlation(siteEnv3)
ecm <- correlation_matrix(as.data.frame(siteEnv3), type = "pearson", show_significance = T, 
                          digits = 2, use = "lower", replace_diagonal = T)

###################################################################################################
###################################################################################################
#STEP #7: PRELIMINARY ANALYSIS

#########################################################
#7a: Determine species with the highest and second highest biomass at each site
#and then display occurrence of dominant and co-dominant species with a histogram.

#Create a histogram of primary and secondary dominance across sites by species.

#Remove non-living categories (dead woody) from the biomass data
last_col <- length(biomassOrig[1,])
nextLast_col <- length(biomassOrig[1,]) - 1

#List of dominant species
prim <- mapply(function(x) {colnames(biomassOrig)[order(biomassOrig[x,])][last_col]}, 
               1:length(biomassOrig[,1]))
seco <- mapply(function(x) {colnames(biomassOrig)[order(biomassOrig[x,])][nextLast_col]}, 
               1:length(biomassOrig[,1]))

dTable <- data.frame(Sites = I(rownames(biomassOrig)), Primary = I(prim), Secondary = I(seco))

udom <- unique(c(unique(dTable[,2]),unique(dTable[,3])))

hprim <- mapply(function(x) {length(dTable[,2][dTable[,2] == udom[x]])}, 1:length(udom))
hseco <- mapply(function(x) {length(dTable[,3][dTable[,3] == udom[x]])}, 1:length(udom))

hdom <- data.frame(Species = I(udom), Primary = hprim, Secondary = hseco)
h2dom <- hdom[order(hdom[,2], decreasing = T),]

par(mai = c(2.5,1.2,1,1))
barplot(t(cbind(h2dom$Primary,h2dom$Secondary)), main = "", 
        xaxt = "n", ylab = "Number of sites", xlab = "", axes = F, beside = T, 
        col = c("white", "dark grey"))
axis(2)
text(matrix(c(seq(2,nrow(h2dom)*3,3), rep(-0.25,nrow(h2dom))), nrow = nrow(h2dom), 
            ncol = 2, byrow = F), srt = 60, adj = 1, xpd = T, labels = paste(h2dom$Species), 
     cex = 0.95)

legend(10,6, c("Highest biomass", "Second highest biomass"), fill = c("white", "dark grey"), bty = "n")

#########################################################
#7b: Arrange species on bar chart with highest to lowest average biomass

#Make a table showing mean/sd biomass
bm <- apply(biomassOrig,2,mean)
bs <- apply(biomassOrig,2,sd)

biomass <- data.frame(Mean = round(bm,2), SD = round(bs,2))

###################################################################################################
###################################################################################################
#STEP #8: Gradient length diagnostics with DCA for untransformed biomass data, log-transformed
#biomass data, and untransformed cover data.

#Create a bindary dataset for species presence absence
#Biomass - untransformed
speocc <- data.trans(so2, method = "power", exp = 0, plot = F)

decorana(speocc, ira = 0)

#Biomass - log-transformed
speocc <- data.trans(sl2, method = "power", exp = 0, plot = F)

decorana(speocc, ira = 0)

#DCA 1 axis length is 0.59714
#This is way below the 3-4 cut-off proposed by ter Braak and Smilauer (2002) used as a threshold 
#for selecting a multivariate analysis technique (< 3 PCA/RDA; data are suited for linear method, and
#> 4 CA/CCA; data are suited for unimodal methods).

#This value indicates a short gradient length (i.e., use linear analysis). This is consistent 
#with sites selection which focused on a narrow set of environmental conditions. 
#I.e. mesic flatwoods with a history of management with frequent rx fire and no other evidence of 
#major disturbance.

###################################################################################################
###################################################################################################
#STEP #9: PCA of the covariance matrix
pca <- rda(decostand(sl2, method = "hellinger", scale = T))#scale = T should be PCA of the correlation matrix
pca
#Plot pca
biplot(pca, scaling = 1)
#Use brokan stick method to determine how many PCA axes to retain
screeplot(pca, bstick = T, type = "l", main = NULL)
#Extract eigenvalues
eigenvals(pca)
summary(eigenvals(pca))

#Fitting environmental data
set.seed(42)#use this to make this permutation analysis reproducibe, otherwise it will be different
#every time, especially when number of permutations is low.
ev <- envfit(pca ~ ., data = siteEnv3, choices = 1:2, scaling = 'symmetric', permutations = 1000)
#options for scaling: "species", "sites", "symmetric" (scales both for species and sites), "none" (raw scores).
#cca() hill = T, rda(): correlation = T (standardizes scores among species).
ev

plot(pca, display = "sites", type = "n", scaling = "symmetric")
points(pca, display = "sites", scaling = "symmetric")
plot(ev, add = T)

surf <- ordisurf(pca ~ FineWD, data = siteEnv3, knots = 10, isotropic = T, main = NULL)
surf <- ordisurf(pca ~ CoarseWD, data = siteEnv3, knots = 10, isotropic = T, main = NULL)
surf <- ordisurf(pca ~ Litter, data = siteEnv3, knots = 10, isotropic = T, main = NULL)
surf <- ordisurf(pca ~ mfri_20yr, data = siteEnv3, knots = 10, isotropic = T, main = NULL)
summary(surf)

###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
#STEP 13: BOUNDARY LAYER REGRESSION FOR FIRE ROTATION

#Old file reassignments
emat <- siteEnv3
ab2mat <- biomassOrig

ev <- emat$mfri_20yr
etext <- "mFRI (20 years)"

###################################################################################################
###################################################################################################
#PANEL PLOT BY ORIGINAL DATA __1
dev.off()
set.panel(2,4)
cf <- 0.5

#1
#Forb
#Plot with sites
genus <- ab2mat$forb
Genus <- "Forbs"

#Boundary layer regression
blr(log, ev, genus, Genus)

#2
#Bunchgrass
#Plot with sites
genus <- ab2mat$bunchgrass
Genus <- "Bunchgrass"

#Boundary layer regression
blr(log, ev, genus, Genus)

#3
#Grass
#Plot with sites
genus <- ab2mat$grass
Genus <- "Grass"

#Boundary layer regression
blr(log, ev, genus, Genus)


#4
#Palmetto
#Plot with sites
genus <- ab2mat$palmetto
Genus <- "Palmetto"

#Boundary layer regression
blr(log, ev, genus, Genus)

#5
#understory
#Plot with sites
genus <- ab2mat$understory
Genus <- "Understory Trees"

#Boundary layer regression
blr(log, ev, genus, Genus)

#6
#Vines
#Plot with sites
genus <- ab2mat$vine
Genus <- "Vines"

#Boundary layer regression
blr(log, ev, genus, Genus)

#7
#Sub-shrubs
#Plot with sites
genus <- ab2mat$subshrub
Genus <- "Sub-shrubs"

#Boundary layer regression
blr(log, ev, genus, Genus)

#8
#Shrubs
#Plot with sites
genus <- ab2mat$shrub
Genus <- "Shrubs"

#Boundary layer regression
blr(log, ev, genus, Genus)


###################################################################################################

ev <- emat$RxFireMgmt
etext <- "Duration of Management (years)"

ev <- emat$BurnFreeInterval
etext = "Time-Since-Last-Fire (years)"

ev <- emat$Season
etext = "Percent Growing Season Burns"

ev <- emat$FireRotation
etext <- "mFRI (years)"


#1-2
#Forb
#Plot with sites
genus <- ab2mat$forb
Genus <- "Forb"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)


#3-11
#Grass
#Plot with sites
genus <- ab2mat$grass + ab2mat$bluestem
Genus <- "Graminoid (ex. Aristida)"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)


#10
#Wiregrass
#Plot with sites
genus <- ab2mat$wiregrass
Genus <- "Wiregrass"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)


#12
#Yucca
#Plot with sites
genus <- ab2mat$yucca
Genus <- "Yucca"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

###################################################################################################
###################################################################################################
#PANEL PLOT BY ORIGINAL DATA __2
dev.off()
set.panel(2,4)
cf <- 0.5

#4-5
#All Palmetto
#Plot with sites
genus <- ab2mat$palmetto
Genus <- "Palmetto"

#Boundary layer regression
blr(log, ev, genus, Genus)

#23
#Carpet Oak
#Plot with sites
genus <- ab2mat$oak.groundcover
Genus <- "Grouncover Oak"

#Boundary layer regression
blr(log, ev, genus, Genus)

#25
#Gallberry
#Plot with sites
genus <- ab2mat$gallberry
Genus <- "Gallberry"

#Linear regression
logmodel <- lm(genus ~ ev)
y <- genus
plot(ev, genus, 
     xlab = etext, ylab = "Biomass (Mg/ha)", main = "Gallberry")
text(4.1, (max(y) - max(y)/16), paste("Intercept", round(logmodel$coefficients[1],4)), cex = cf)
text(5.9, (max(y) - max(y)/16), paste("Slope", round(logmodel$coefficients[2],4)), cex = cf)
ms <- summary(logmodel)
text(4.1, (max(y) - max(y)/8), paste("r-squared", round(ms$r.squared,2)), cex = cf)
text(5.9, (max(y) - max(y)/8), paste("P-value", round(ms$coefficients[8],4)), cex = cf)
abline(logmodel)

#29
#Mound Oak
#Plot with sites
genus <- ab2mat$oak.shrub
Genus <- "Shrub Oak"

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#4-5
#All Palmetto
#Plot with sites
genus <- ab2mat$palmetto
Genus <- "Palmetto"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#23
#Carpet Oak
#Plot with sites
genus <- ab2mat$oak.groundcover
Genus <- "Grouncover Oak"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#25
#Gallberry
#Plot with sites
genus <- ab2mat$gallberry
Genus <- "Gallberry"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#29
#Mound Oak
#Plot with sites
genus <- ab2mat$oak.shrub
Genus <- "Shrub Oak"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

###################################################################################################
###################################################################################################
#PANEL PLOT BY ORIGINAL DATA __3
dev.off()
set.panel(2,4)
cf <- 0.5

#30
#Rough Fetterbush
#Plot with sites
genus <- ab2mat$rough.fetterbush
Genus <- "Rough Fetterbush"

#Boundary layer regression
blr(log, ev, genus, Genus)

#19
#Darrow's blueberry
#Plot with sites
genus <- ab2mat$Darrows.blueberry
Genus <- "Darrow's blueberry"

#Boundary layer regression
blr(log, ev, genus, Genus)

#35
#Highbush blueberry
#Plot with sites
genus <- ab2mat$highbush.blueberry
Genus <- "Highbush Blueberry"

#Boundary layer regression
blr(exp, ev, genus, Genus)

#43
#Sweetbay
#Plot with sites
genus <- ab2mat$sweetbay
Genus <- "Sweetbay"

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#30
#Rough Fetterbush
#Plot with sites
genus <- ab2mat$rough.fetterbush
Genus <- "Rough Fetterbush"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#19
#Darrow's blueberry
#Plot with sites
genus <- ab2mat$Darrows.blueberry
Genus <- "Darrow's blueberry"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#35
#Highbush blueberry
#Plot with sites
genus <- ab2mat$highbush.blueberry
Genus <- "Highbush Blueberry"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#43
#Sweetbay
#Plot with sites
genus <- ab2mat$sweetbay
Genus <- "Sweetbay"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

###################################################################################################
###################################################################################################
#PANEL PLOT BY ORIGINAL DATA __4
dev.off()
set.panel(2,4)
cf <- 0.5

#18
#Huckleberry
#Plot with sites
genus <- ab2mat$huckleberry
Genus <- "Huckleberry"

#Boundary layer regression
blr(exp, ev, genus, Genus)

#38
#Titi
#Plot with sites
genus <- ab2mat$titi
Genus <- "Titi"

#Boundary layer regression
blr(log, ev, genus, Genus)

#41
#Oak Tree
#Plot with sites
genus <- ab2mat$oak.tree
Genus <- "Oak Tree"

#Boundary layer regression
blr(log, ev, genus, Genus)

#15
#Greenbriar
#Plot with sites
genus <- ab2mat$greenbriar
Genus <- "Greenbriar"

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#18
#Huckleberry
#Plot with sites
genus <- ab2mat$huckleberry
Genus <- "Huckleberry"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#38
#Titi
#Plot with sites
genus <- ab2mat$titi
Genus <- "Titi"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#41
#Oak Tree
#Plot with sites
genus <- ab2mat$oak.tree
Genus <- "Oak Tree"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#15
#Greenbriar
#Plot with sites
genus <- ab2mat$greenbriar
Genus <- "Greenbriar"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

###################################################################################################
###################################################################################################
#PANEL PLOT BY ORIGINAL DATA __5
dev.off()
set.panel(2,4)
cf <- 0.5

###################################################################################################
#14
#Yellow Jessamine
#Plot with sites
genus <- ab2mat$yellow.jessamine
Genus <- "Yellow Jessamine"

#Boundary layer regression
blr(log, ev, genus, Genus)

#37
#Sparkleberry
#Plot with sites
genus <- ab2mat$sparkleberry
Genus <- "Sparkleberry"

#Boundary layer regression
blr(log, ev, genus, Genus)

#32
#Smooth Fetterbush
#Plot with sites
genus <- ab2mat$smooth.fetterbush
Genus <- "Smooth Fetterbush"

#Boundary layer regression
blr(exp, ev, genus, Genus)

#33
#Large Gallberry
#Plot with sites
genus <- ab2mat$large.gallberry
Genus <- "Large Gallberry"

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#14
#Yellow Jessamine
#Plot with sites
genus <- ab2mat$yellow.jessamine
Genus <- "Yellow Jessamine"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#37
#Sparkleberry
#Plot with sites
genus <- ab2mat$sparkleberry
Genus <- "Sparkleberry"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#32
#Smooth Fetterbush
#Plot with sites
genus <- ab2mat$smooth.fetterbush
Genus <- "Smooth Fetterbush"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#33
#Large Gallberry
#Plot with sites
genus <- ab2mat$large.gallberry
Genus <- "Large Gallberry"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

###################################################################################################
###################################################################################################
#PANEL PLOT BY ORIGINAL DATA __6
dev.off()
set.panel(2,4)
cf <- 0.5

#13
#Muscadine Grape
#Plot with sites
genus <- ab2mat$muscadine.grape
Genus <- "Muscadine Grape"

#Boundary layer regression
blr(log, ev, genus, Genus)

#36
#Yaupon
#Plot with sites
genus <- ab2mat$yaupon
Genus <- "Yaupon"

#Boundary layer regression
blr(log, ev, genus, Genus)

#42
#Sweetleaf
#Plot with sites
genus <- ab2mat$sweetleaf
Genus <- "Sweetleaf"

#Boundary layer regression
blr(log, ev, genus, Genus)

#21
#Wax Myrtle
#Plot with sites
genus <- ab2mat$wax.myrtle
Genus <- "Wax Myrtle"

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#13
#Muscadine Grape
#Plot with sites
genus <- ab2mat$muscadine.grape
Genus <- "Muscadine Grape"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#36
#Yaupon
#Plot with sites
genus <- ab2mat$yaupon
Genus <- "Yaupon"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#42
#Sweetleaf
#Plot with sites
genus <- ab2mat$sweetleaf
Genus <- "Sweetleaf"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#21
#Wax Myrtle
#Plot with sites
genus <- ab2mat$wax.myrtle
Genus <- "Wax Myrtle"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

###################################################################################################
###################################################################################################
#PANEL PLOT BY ORIGINAL DATA __7
dev.off()
set.panel(2,4)
cf <- 0.5

#17
#Wicky
#Plot with sites
genus <- ab2mat$wicky
Genus <- "Wicky"

#Boundary layer regression
blr(exp, ev, genus, Genus)

#20
#St. Andrew's Cross
#Plot with sites
genus <- ab2mat$St..Andrews.cross
Genus <- "St. Andrew's Cross"

#Boundary layer regression
blr(log, ev, genus, Genus)

#16
#Gopher Apple
#Plot with sites
genus <- ab2mat$gopher.apple
Genus <- "Gopher Apple"

#Boundary layer regression
blr(exp, ev, genus, Genus)

#31
#Sweet Pepperbush
#Plot with sites
genus <- ab2mat$sweet.pepperbush
Genus <- "Sweet Pepperbush"

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#17
#Wicky
#Plot with sites
genus <- ab2mat$wicky
Genus <- "Wicky"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#20
#St. Andrew's Cross
#Plot with sites
genus <- ab2mat$St..Andrews.cross
Genus <- "St. Andrew's Cross"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#16
#Gopher Apple
#Plot with sites
genus <- ab2mat$gopher.apple
Genus <- "Gopher Apple"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#31
#Sweet Pepperbush
#Plot with sites
genus <- ab2mat$sweet.pepperbush
Genus <- "Sweet Pepperbush"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

###################################################################################################
###################################################################################################
#PANEL PLOT BY ORIGINAL DATA __8
dev.off()
set.panel(2,4)
cf <- 0.5

#39
#Fringe Tree
#Plot with sites
genus <- ab2mat$fringe.tree
Genus <- "Fringe Tree"

#Boundary layer regression
blr(log, ev, genus, Genus)

#24
#Chokeberry
#Plot with sites
genus <- ab2mat$chokeberry
Genus <- "Chokeberry"

#Boundary layer regression
blr(log, ev, genus, Genus)

#22
#Conradina
#Plot with sites
genus <- ab2mat$conradina
Genus <- "Conradina"

#Boundary layer regression
blr(log, ev, genus, Genus)

#40
#American Holly
#Plot with sites
genus <- ab2mat$American.holly
Genus <- "American Holly"

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#39
#Fringe Tree
#Plot with sites
genus <- ab2mat$fringe.tree
Genus <- "Fringe Tree"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#24
#Chokeberry
#Plot with sites
genus <- ab2mat$chokeberry
Genus <- "Chokeberry"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#22
#Conradina
#Plot with sites
genus <- ab2mat$conradina
Genus <- "Conradina"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#40
#American Holly
#Plot with sites
genus <- ab2mat$American.holly
Genus <- "American Holly"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

###################################################################################################
###################################################################################################
#PANEL PLOT BY ORIGINAL DATA __9
dev.off()
set.panel(2,4)
cf <- 0.5

#26
#Dewberry
#Plot with sites
genus <- ab2mat$dewberry
Genus <- "Dewberry"

#Boundary layer regression
blr(log, ev, genus, Genus)

#27
#St. Johnswort
#Plot with sites
genus <- ab2mat$St..Johnswort
Genus <- "St. Johnswort"

#Boundary layer regression
blr(log, ev, genus, Genus)

#28
#Narrowleaf Paw Paw
#Plot with sites
genus <- ab2mat$narrowleaf.pawpaw
Genus <- "Narrowleaf Paw Paw"

#Boundary layer regression
blr(log, ev, genus, Genus)

#34
#Winged Sumac
#Plot with sites
genus <- ab2mat$winged.sumac
Genus <- "Winged Sumac"

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#26
#Dewberry
#Plot with sites
genus <- ab2mat$dewberry
Genus <- "Dewberry"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#27
#St. Johnswort
#Plot with sites
genus <- ab2mat$St..Johnswort
Genus <- "St. Johnswort"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#28
#Narrowleaf Paw Paw
#Plot with sites
genus <- ab2mat$narrowleaf.pawpaw
Genus <- "Narrowleaf Paw Paw"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#34
#Winged Sumac
#Plot with sites
genus <- ab2mat$winged.sumac
Genus <- "Winged Sumac"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

###################################################################################################
###################################################################################################
#PANEL PLOT BY ORIGINAL DATA __10
dev.off()
set.panel(2,4)
cf <- 0.5

#6
#Live Shrub
#Plot with sites
genus <- ab2mat$live.shrub
Genus <- "Live shrub"

#Boundary layer regression
blr(log, ev, genus, Genus)

#7
#Dead Shrub
#Plot with sites
genus <- ab2mat$dead.shrub
Genus <- "Dead shrub"

#Boundary layer regression
blr(log, ev, genus, Genus)

#8
#Live Tree
#Plot with sites
genus <- ab2mat$live.tree
Genus <- "Live tree"

#Boundary layer regression
blr(log, ev, genus, Genus)

#9
#Dead Tree
#Plot with sites
genus <- ab2mat$dead.tree
Genus <- "Dead tree"

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#6
#Live Shrub
#Plot with sites
genus <- ab2mat$live.shrub
Genus <- "Live shrub"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#7
#Dead Shrub
#Plot with sites
genus <- ab2mat$dead.shrub
Genus <- "Dead shrub"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#8
#Live Tree
#Plot with sites
genus <- ab2mat$live.tree
Genus <- "Live tree"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#9
#Dead Tree
#Plot with sites
genus <- ab2mat$dead.tree
Genus <- "Dead tree"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)






dev.off()



###################################################################################################
#1
#live.forb
#Plot with sites
genus <- ab2mat$live.forb
Genus <- "Live forb"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#2
#dead.forb
#Plot with sites
genus <- ab2mat$dead.forb
Genus <- "Dead forb"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#3
#grass
#Plot with sites
genus <- ab2mat$grass
Genus <- "Grass (ex. bunchgrass)"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#4
#Live Palmetto
#Plot with sites
genus <- ab2mat$live.palmetto
Genus <- "Live palmetto"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#5
#Dead Palmetto
#Plot with sites
genus <- ab2mat$dead.palmetto
Genus <- "Dead palmetto"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#5-7-9
#Dead Woody
#Plot with sites
genus <- ab2mat$dead.palmetto + ab2mat$dead.shrub + ab2mat$dead.tree
Genus <- "Dead Woody"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#11
#Bluestem
#Plot with sites
genus <- ab2mat$bluestem.grass
Genus <- "Bluestem"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#Boundary layer regression
blr(log, ev, genus, Genus)



##################################################################################################
#23-29-41
#Oak
#Plot with sites
genus <- ab2mat$oak.groundcover + ab2mat$oak.shrub + ab2mat$oak.tree
Genus <- "Oak"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#Boundary layer regression
blr(log, ev, genus, Genus)

plot(0,0, main = "NA")

###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
























































































































###################################################################################################
###################################################################################################
#STEP #6: DATA SCREENING

#########################################################
#6a: Summary stats for biomass data

#Summary stats
round(stat.desc(biomassOrig),2)

#Look at how data is distributed
foa.plots(biomassOrig)

#Cum. Number of Species (aka Genus) vs Frequency of Occurrence
#>> There are 33 genus of plants. Half (17) of those occur on 14 of the sites
#   The other half are uncommon and occur at less than 7 (total n = 21) sites.

#Cum. Number of Species vs Percent Occurrence
#Not sure what "percent occurrence" (y-axis) means. How is this different than the prior chart?
#>> Similar trend as above
#   21 genus occur less than 50%
#   12 genus have occurence greater than 70, with all but one higher than 90%

#Histogram of Species Occurrence
#>> Bimodal with largest peak in the 0-5 sites range and a smaller, trough in the 10-15 
#   sites range and secondary peak in the 20-25 sites range.

#Histogram of Log(Species Occurrence)
#>> Not sure what this histogram is displaying but distrubution is more or less flat.

#Cumulative Distribution of Species Mean Abundance
#>> Exponentially increasing
# This shows that about 5 genus are dominant (genus on steep curve)
# and the remaining 28 genus have low biomass.

#Species Occurrence vs Mean Abundance
#Genus with highest abundance have higher occurrence, but there are also many genus with high occurrence
#and low abundance.

#Species Occurrence vs Log(Mean Abunandance)
#>> There is a weak (probably significant) upward trend of increasing abundance with frequency
#   of occurrence.

#Cumulative Distribution of Plot (Site) Richness.
#>> Site richness doubles from 12 to 21 and increases on a steady slope across all sites. I.e., there are no
#distinct groups of low and high richness plots, but rather a gradual transition from low to high.

#Cumulative Distribution of Plot Total Abundance
#>> There is a linear upward trend. Not sure what abundance is, loading?
#   Looks like it is showing total loading by site arranged from lowest to highest
#   but given the chart title (cumulative) I believe I am misreading both this, and the
#   above chart.

#Plot Richness vs Total Abundance
#>> There is a weak and possibly insignificant upward trend.
#   Plots along the line from low abundance/low richness to high:
#   A32, E508A, E505, E807D, E807B, E100B-E
#   Outliers include:
#   Low plot richness and high abundance: S209
#   Medium plot richness and low abundance: E100B-W

#Review distance measures
rankindex(siteEnv3, siteBiomassLogTrans5, 
          indices = c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis", "chisq", "chord"), 
          stepacross = F, method = "spearman")

#Drop rare species < 50% of sites.
bo <- drop.var(biomassOrig, min.po=50)
bl <- drop.var(biomassLog, min.po=50)

bore1 <- drop.var(bor1, min.po=80)
bore2 <- drop.var(bor2, min.po=80)
blre1 <- drop.var(blr1, min.po=80)
blre2 <- drop.var(blr2, min.po=80)


#Check to see how many species were dropped for untransformed biomass data

#All
length(biomassOrig[1,])#
length(bo[1,])#
length(biomassLog[1,])#
length(bl[1,])
#39 to 13 variables (26 drops)

#Oct 18, 2022
#Drops from 42 to 13 variables.

#Original - Egin
length(bor1[1,])#
length(bore1[1,])
#39 to 14 variables (25 drops)

#Original - Apalachicola
length(bor2[1,])#
length(bore2[1,])
#39 to 15 variables (24 drops)

#Log - Eglin
length(blr1[1,])#
length(blre1[1,])
#39 to 14 variables (25 drops)

#Log - Apalachicola
length(blr2[1,])#
length(blre2[1,])
#39 to 15 variables (24 drops)

str(bo)#22 obs and 13 variables
str(bl)#22 obs and 13 variables

#Convert data.frame to a matrix (needed for uv.plots)

#All
bo2 <- data.matrix(frame = bo, rownames.force = NA)
bl2 <- data.matrix(frame = bl, rownames.force = NA)

#Original (E)glin and (A)palachicola
boE <- data.matrix(frame = bore1, rownames.force = NA)
boA <- data.matrix(frame = bore2, rownames.force = NA)

#Log-transformed (E)glin and (A)palachicola
blE <- data.matrix(frame = blre1, rownames.force = NA)
blA <- data.matrix(frame = blre2, rownames.force = NA)


#uv.plots displays histogram, box and whisker, cumulative distribution and normal q-q plots
#in one pane for each variable.
uv.plots(bo2)
uv.plots(bl2)

uv.plots(boE)
uv.plots(boA)


#What percent of values are zero, if it is over 50% you should consider changing the data
#to presence/absence
length(biomassOrig[biomassOrig == 0])/(length(biomassOrig[1,])*length(biomassOrig[,1]))
#54% zero values >>> consider a binary transformation

length(bo2[bo2 == 0])/(length(bo2[1,])*length(bo2[,1]))
#13% zero values >>> do not do binary transformation
length(bl2[bl2 == 0])/(length(bl2[1,])*length(bl2[,1]))
#13% zero values >>> do not do binary transformation

length(boE[boE == 0])/(length(boE[1,])*length(boE[,1]))
#12% zero values >>> do not do binary transformation
length(boA[boA == 0])/(length(boA[1,])*length(boA[,1]))
#16% zero values >>> do not do binary transformation


#Look at correlation between species variables.
chart.Correlation(biomassOrig, method = "pearson")
chart.Correlation(bo2, method = "pearson")
chart.Correlation(bl2, method = "pearson")
chart.Correlation(boE, method = "pearson")
chart.Correlation(boA, method = "pearson")

dev.off()

#########################################################
#6c: Summary stats for environmental data

#Show how environmental variables are correlated.
chart.Correlation(siteEnv3[,-1])

#Effect of column-standardization on untransformed data.
siteEnv4 <- data.stand(siteEnv3[,-1], method = 'total', margin = 'column', plot = F)

#Split data by region
seE <- siteEnv4[re == 1,]
seA <- siteEnv4[re == 2,]

###################################################################################################
###################################################################################################
#STEP #7: PRELIMINARY ANALYSIS


#########################################################
#7a: Determine species with the highest and second highest biomass at each site
#and then display occurrence of dominant and co-dominant species with a histogram.

#Create a histogram of primary and secondary dominance across sites by species.

#Remove non-living categories (dead woody) from the biomass data
abhmat <- biomassOrig#send newer file name to old file name
#abhmat <- ab2mat[,-39] -- removed 2022-10-01 -- this would have removed palmetto, no 
#non-living categories present.

#List of dominant species
prim <- mapply(function(x) {colnames(abhmat)[order(abhmat[x,])][33]}, 1:length(abhmat[,1]))
seco <- mapply(function(x) {colnames(abhmat)[order(abhmat[x,])][32]}, 1:length(abhmat[,1]))

dTable <- data.frame(Sites = I(rownames(abhmat)), Primary = I(prim), Secondary = I(seco))

udom <- unique(c(unique(dTable[,2]),unique(dTable[,3])))

hprim <- mapply(function(x) {length(dTable[,2][dTable[,2] == udom[x]])}, 1:length(udom))
hseco <- mapply(function(x) {length(dTable[,3][dTable[,3] == udom[x]])}, 1:length(udom))

hdom <- data.frame(Species = I(udom), Primary = hprim, Secondary = hseco)
h2dom <- hdom[order(hdom[,2], decreasing = T),]

par(mai = c(2.5,1.2,1,1))
barplot(t(cbind(h2dom$Primary,h2dom$Secondary)), main = "", 
        xaxt = "n", ylab = "Number of sites", xlab = "", axes = F, beside = T, 
        col = c("white", "dark grey"))
axis(2)
text(matrix(c(seq(2,nrow(h2dom)*3,3), rep(-0.25,nrow(h2dom))), nrow = nrow(h2dom), 
            ncol = 2, byrow = F), srt = 60, adj = 1, xpd = T, labels = paste(h2dom$Species), 
     cex = 0.95)

legend(20,6, c("Dominant", "Co-dominant"), fill = c("white", "dark grey"), bty = "n")

#########################################################
#7b: Arrange species on bar chart with highest to lowest average biomass, also
#display average cover

#Make a table showing mean/sd biomass
bm <- apply(abhmat,2,mean)
bs <- apply(abhmat,2,sd)

biomass <- data.frame(Mean = round(bm,2), SD = round(bs,2))

###################################################################################################
###################################################################################################
#STEP #8: DATA STANDARDIZATION

#Not necessary to examine effect of column-standardizations on data. Units are the same. 

#Don't want to do a standardization by row but here's the effects on data.

###################################################################################################
###################################################################################################
#STEP #9: OUTLIERS

#Find univariate outliers
uv.outliers(biomassOrig, id = names(biomassOrig)[1]:names(biomassOrig)[length(names(biomassOrig))], sd.limit = 1)
#did not figure out this function because this is not so important in multivariate space.
#I don't think Euclidean distance is appropriate for biological data, should I be using
#a different method? At any rate, 3 is the recommended sd.limit. None of your sites
#exceeds 3 for the Euclidean method.

mv.outliers(bo2, method = 'mahalanobis', sd.limit=1)#E103B_S3
mv.outliers(bl2, method = 'mahalanobis', sd.limit=1)#A18
mv.outliers(co2, method = 'mahalanobis', sd.limit=1)#E403B

###################################################################################################
###################################################################################################
#STEP #10: PCoA


###################################################################################################
#STEP #10a: Principal Coordinate Analysis on biomass data.

lbmO1 <- vegdist(bo2, "manhattan")#untransformed data with manhattan distance measure.
lbmO2 <- vegdist(bo2, "bray")#untransformed data with bray-curtis distance measure.
lbmL1 <- vegdist(bl2, "manhattan")#log-transformed data with manhattan distance measure.
lbmL2 <- vegdist(bl2, "bray")#log-transformed data with bray-curtis distance measure.

pcoaO1 <- cmdscale(lbmO1, k = 21, eig = T, add = F)
pcoaO2 <- cmdscale(lbmO2, k = 21, eig = T, add = F)
pcoaL1 <- cmdscale(lbmL1, k = 21, eig = T, add = F)
pcoaL2 <- cmdscale(lbmL2, k = 21, eig = T, add = F)

#Percent variation explained by each principal coordinate
round(pcoaO1$eig/sum(pcoaO1$eig)*100,2)
round(pcoaO2$eig/sum(pcoaO2$eig)*100,2)
round(pcoaL1$eig/sum(pcoaL1$eig)*100,2)
round(pcoaL2$eig/sum(pcoaL2$eig)*100,2)

#Assign pcoa object to functions between  #>>>><<<<<# lines
x <- pcoaL2
y <- bl2
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<##
#Broken stick plot
plot(x$eig/sum(pcoaO1$eig)*100-bstick(21)*100,xlab = "PC", 
     ylab="Actual-random % variation explained")
abline(h=0)


#Ordination plot
ordiplot(x, choices = c(1,2), type = "text", display = "sites", xlab = "PC 1", 
         ylab = "PC 2")

#Conduct a linear correlation analysis to determine influence of plant types
vec.pcoa <- envfit(x$points, y, perm = 1000)

#Another way of writing above function
#vec.pcoa <- envfit(scores(x), y, perm = 1000)

#View significance pf variables
vec.pcoa

#Sub-objects
names(vec.pcoa$vectors)
vec.pcoa$vectors$arrows

#Ordination plot with variables
dev.off()
ordiplot(x, choices = c(1,2), type = "text", display = "sites", xlab = "PC 1", 
         ylab = "PC 2")
plot(vec.pcoa, p.max=0.1, col="blue")
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<##


###################################################################################################
#STEP #10b: Principal Coordinate Analysis on original biomass data for log-transformed data.

lcmO1 <- vegdist(co2, "manhattan")#untransformed data with manhattan distance measure.
lcmO2 <- vegdist(co2, "bray")#untransformed data with bray-curtis distance measure.

pcoaO1 <- cmdscale(lcmO1, k = 21, eig = T, add = F)
pcoaO2 <- cmdscale(lcmO2, k = 21, eig = T, add = F)

#Percent variation explained by each principal coordinate
round(pcoaO1$eig/sum(pcoaO1$eig)*100,2)
round(pcoaO2$eig/sum(pcoaO2$eig)*100,2)

#Assign pcoa object to functions between  #>>>><<<<<# lines
x <- pcoaO2
y <- co2
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<##
#Broken stick plot
plot(x$eig/sum(pcoa1$eig)*100-bstick(21)*100,xlab = "PC", 
     ylab="Actual-random % variation explained")
abline(h=0)


#Ordination plot
ordiplot(x, choices = c(1,2), type = "text", display = "sites", xlab = "PC 1", 
         ylab = "PC 2")

#Conduct a linear correlation analysis to determine influence of plant types
vec.pcoa <- envfit(x$points, y, perm = 1000)

#Another way of writing above function
#vec.pcoa <- envfit(scores(x), y, perm = 1000)

#View significance pf variables
vec.pcoa

#Sub-objects
names(vec.pcoa$vectors)
vec.pcoa$vectors$arrows

#Ordination plot with variables
dev.off()
ordiplot(x, choices = c(1,2), type = "text", display = "sites", xlab = "PC 1", 
         ylab = "PC 2")
plot(vec.pcoa, p.max=0.1, col="blue")
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<##

###################################################################################################
###################################################################################################
#STEP #11: Gradient length diagnostics with DCA for untransformed biomass data, log-transformed
#biomass data, and untransformed cover data.

#Create a bindary dataset for species presence absence
#Biomass - untransformed
speocc <- data.trans(bo2, method = "power", exp = 0, plot = F)

decorana(speocc, ira = 0)
#Axis length is 1.3968

#Biomass - log-transformed
speocc <- data.trans(bl2, method = "power", exp = 0, plot = F)

decorana(speocc, ira = 0)
#Axis length is 1.3968

#Cover - intransformed
speocc <- data.trans(co2, method = "power", exp = 0, plot = F)

decorana(speocc, ira = 0)
#Axis length is 1.2931

#Length of axis 1 is 1.46 (un-tranformed and log-transformed) and and indicates the gradient length 
#is short (linear analysis). This is consistent with sites selection which focused on a narrow set 
#of environmental conditions. I.e. mesic flatwoods with a history of management with frequent rx 
#fire and no other evidence of major disturbance.

###################################################################################################
###################################################################################################
#STEP #12: Detect and remove(?) outliers

#Univariate
#Environmental data
uv.outliers(siteEnv4, id ='Canopy:Duration', var = 'Soil_Pref', sd.limit = 2)

#multivariate
#Environmental data
mv.outliers(siteEnv4, method = "euclidean", sd.limit =2)#E100B-E & E503_S4

#biomass
mv.outliers(bo2, method = "manhattan", sd.limit =2)#S209 & E100B_S4
mv.outliers(bo2, method = "bray", sd.limit =2)#E103B_S3
mv.outliers(bl2, method = "manhattan", sd.limit =2)#S209
mv.outliers(bl2, method = "bray", sd.limit =2)#none

#cover
mv.outliers(co2, method = "manhattan", sd.limit =2)#E501B
mv.outliers(co2, method = "bray", sd.limit =2)#E100B-E


#E100B-E
ematR1 <- siteEnv4[-(row.names(siteEnv4) == "E100B-E"),]
e1matR1 <- e1mat[-(row.names(e1mat) == "E100B-E"),]
gb2matR1 <- gb2mat[-(row.names(gb2mat) == "E100B-E"),]
gb3matR1 <- gb3mat[-(row.names(gb3mat) == "E100B-E"),]
lbmatR1 <- lbmat[-(row.names(lbmat) == "E100B-E"),]
lb2matR1 <- lb2mat[-(row.names(lb2mat) == "E100B-E"),]

###################################################################################################
###################################################################################################
#STEP #13: Direct gradient analysis with CAP

#October 14, 2022
#Based on tutorial here:
#https://archetypalecology.wordpress.com/2018/02/21/distance-based-redundancy-analysis-db-rda-in-r/
dbRDA_biomassOrig <- capscale(biomassOrig ~ Canopy + Litter + CoarseWD + FineWD + mfri_20yr + 
                                StDevFRI_20YR + Season_20yr + RxFireMgmt, 
                              siteEnv3, dist = "bray")
plot(dbRDA_biomassOrig)
anova(dbRDA_biomassOrig)
anova(dbRDA_biomassOrig, by = 'axis', perm.max = 500)
anova(dbRDA_biomassOrig, by = 'terms', perm.max = 200)
scores_all <- scores(dbRDA_biomassOrig)
scores_sites <- scores_all$sites
scores_species <- scores_all$species
scores_sites_environment <- cbind(scores_sites, siteEnv3[c(1:5,7:8,10)])
correlations <- cor(scores_sites_environment)
correlations2 <- correlations[3:10,1:2]
correlations3 <- correlations2[c(2:4,8),1:2]

#October 14, 2022
#Based on tutorial here:
#https://www.davidzeleny.net/anadat-r/doku.php/en:rda_cca_examples
b.hell <- decostand(biomassOrig, 'hell')
dbRDA <- rda(b.hell ~ mfri_20yr + Season_20yr + RxFireMgmt, siteEnv3)
constrained_eig <- dbRDA$CCA$eig/dbRDA$tot.chi*100
unconstrained_eig <- dbRDA$CA$eig/dbRDA$tot.chi*100
expl_var <- c(constrained_eig, unconstrained_eig)
barplot (expl_var[1:20], col = c(rep ('red', length (constrained_eig)), rep ('black', length (unconstrained_eig))),
         las = 2, ylab = '% variation')
ordiplot(dbRDA)
head(summary(dbRDA))





#hist(emat$FireRotation, xlab = "Fire rotation (years)", main = "Distribution of fire rotation among sites", 
#     ylim = c(0,8))
#text(5,7, paste("Mean:", round(mean(emat$FireRotation),1), "years"), pos = 4)
#text(5,6.5, paste("Range:", round(min(emat$FireRotation),1),"-",round(max(emat$FireRotation),1), "years"), pos = 4)

#hist(emat$FireRotation_10yr, xlab = "Fire rotation (years)", main = "Distribution of fire rotation among sites", 
#     ylim = c(0,8))
#text(5,7, paste("Mean:", round(mean(emat$FireRotation_10yr),1), "years"), pos = 4)
#text(5,6.5, paste("Range:", round(min(emat$FireRotation_10yr),1),"-",round(max(emat$FireRotation_10yr),1), "years"), pos = 4)


#hist(emat$Season, xlab = "Growing:dormant season ratio", main = "Distribution of growing:dormant season burns among sites", 
#     ylim = c(0,8))
#text(0.4,7, paste("Mean:", round(mean(emat$Season),1), "years"), pos = 4)
#text(0.4,6.5, paste("Range:", round(min(emat$Season),1),"-",round(max(emat$Season),1), "years"), pos = 4)


#hist(emat$Season_10yr, xlab = "Growing:dormant season ratio", main = "Distribution of growing:dormant season burns among sites", 
#     ylim = c(0,12))
#text(0.4,7, paste("Mean:", round(mean(emat$Season_10yr),1), "years"), pos = 4)
#text(0.4,6.5, paste("Range:", round(min(emat$Season_10yr),1),"-",round(max(emat$Season_10yr),1), "years"), pos = 4)


#hist(emat$RxFireMgmt, xlab = "Prescribed fire management (years)", main = "Distribution of management periods among sites", 
#     ylim = c(0,15))
#text(17,12, paste("Mean:", round(mean(emat$RxFireMgmt),1), "years"), pos = 4)
#text(17,11, paste("Range:", round(min(emat$RxFireMgmt),1),"-",round(max(emat$RxFireMgmt),1), "years"), pos = 4)

#hist(emat$BurnFreeInterval, xlab = "Recent fire free period (years)", main = "Distribution of most recent fire free periods among sites", 
#     ylim = c(0,8))
#text(4,6, paste("Mean:", round(mean(emat$BurnFreeInterval),1), "years"), pos = 4)
#text(4,5.5, paste("Range:", round(min(emat$BurnFreeInterval),1),"-",round(max(emat$BurnFreeInterval),1), "years"), pos = 4)



###################################################################################################
###################################################################################################
#STEP #14: #Distance-based RDA

###################################################################################################
#STEP #14a: untransformed biomass (manhattan distance)
y <- vector('numeric', length = 1000)
for (i in 1:1000)
{
gcap1 <- capscale(bo2 ~ FireRotation + Season + RxFireMgmt, 
                  data = siteEnv4, distance = "manhattan")

gcap1s <- summary(gcap1)

x <- anova(gcap1)#
y[i] <- x[1,5]
}
hist(y)
median(y)
range(y)

anova(gcap1, by = "axis")#
anova(gcap1, by = "term")#
plot(gcap1, choices = c(1,2), display = c("wa", "sp", "bp"), scaling = 3, 
     title = "Biomass by genera (model significance: p =)", xlab = "CAP1 (%)", 
     ylab = "CAP2 (%)")

###################################################################################################
#STEP #14b-a: untransformed biomass (bray-curtis distance)
y <- vector('numeric', length = 1000)
for (i in 1:1000)
{
  
gcap2a <- capscale(bo2 ~ FireRotation + Season + RxFireMgmt, 
                  data = siteEnv4, distance = "bray")

gcap2as <- summary(gcap2a)
x <- anova(gcap2a)#
y[i] <- x[1,5]
}
hist(y)
median(y)
range(y)

anova(gcap2a, by = "axis")#
anova(gcap2a, by = "term")#
plot(gcap2a, choices = c(1,2), display = c("wa", "sp", "bp"), scaling = 3, 
     title = "Biomass by genera (model significance: p = )", xlab = "CAP1 (%)", ylab = "CAP2 (%)")

###################################################################################################
#STEP #14b-b: untransformed biomass (bray-curtis distance) > test order of entry
y <- vector('numeric', length = 1000)
for (i in 1:1000)
{
  
  gcap2b <- capscale(bo2 ~ BurnFreeInterval + Season + RxFireMgmt + FireRotation, 
                    data = siteEnv4, distance = "bray")
  
  gcap2bs <- summary(gcap2b)
  x <- anova(gcap2b)#
  y[i] <- x[1,5]
}
hist(y)
median(y)
range(y)

anova(gcap2b, by = "axis")#
anova(gcap2b, by = "term")#
plot(gcap2b, choices = c(1,2), display = c("wa", "sp", "bp"), scaling = 3, 
     title = "Biomass by genera (model significance: p = )", xlab = "CAP1 (%)", ylab = "CAP2 (%)")

###################################################################################################
#STEP #14b-c: untransformed biomass (bray-curtis distance) > test order of entry
y <- vector('numeric', length = 1000)
for (i in 1:1000)
{
  
  gcap2c <- capscale(bo2 ~ Season + RxFireMgmt + FireRotation + BurnFreeInterval, 
                     data = siteEnv4, distance = "bray")
  
  gcap2cs <- summary(gcap2c)
  x <- anova(gcap2c)#
  y[i] <- x[1,5]
}
hist(y)
median(y)
range(y)

anova(gcap2c, by = "axis")#
anova(gcap2c, by = "term")#
plot(gcap2c, choices = c(1,2), display = c("wa", "sp", "bp"), scaling = 3, 
     title = "Biomass by genera (model significance: p = )", xlab = "CAP1 (%)", ylab = "CAP2 (%)")

###################################################################################################
#STEP #14b-d: untransformed biomass (bray-curtis distance) > test order of entry
y <- vector('numeric', length = 1000)
for (i in 1:1000)
{
  
  gcap2d <- capscale(bo2 ~ Season + RxFireMgmt + FireRotation + BurnFreeInterval, 
                     data = siteEnv4, distance = "bray")
  
  gcap2ds <- summary(gcap2d)
  x <- anova(gcap2d)#
  y[i] <- x[1,5]
}
hist(y)
median(y)
range(y)

anova(gcap2d, by = "axis")#
anova(gcap2d, by = "term")#
plot(gcap2d, choices = c(1,2), display = c("wa", "sp", "bp"), scaling = 3, 
     title = "Biomass by genera (model significance: p = )", xlab = "CAP1 (%)", ylab = "CAP2 (%)")

###################################################################################################
#STEP #14b-x: untransformed biomass (bray-curtis distance) > test remaining order of entry
gcap2x <- capscale(bo2 ~ RxFireMgmt + Season + FireRotation + BurnFreeInterval, 
                     data = siteEnv4, distance = "bray")
anova(gcap2x, by = "axis")#
anova(gcap2x, by = "term")#
plot(gcap2x, choices = c(1,2), display = c("wa", "sp", "bp"), scaling = 3, 
     title = "Biomass by genera (model significance: p = )", xlab = "CAP1 (%)", ylab = "CAP2 (%)")

###################################################################################################
#STEP #14b-Brooke: untransformed biomass (bray-curtis distance) > test remaining order of entry
gcap2y <- capscale(bo2 ~ RxFireMgmt + FireRotation, 
                   data = siteEnv4, distance = "bray")
anova(gcap2y, by = "axis")#
anova(gcap2y, by = "term")#
plot(gcap2y, choices = c(1,2), display = c("wa", "sp", "bp"), scaling = 3, 
     title = "Biomass by genera (model significance: p = )", xlab = "CAP1 (%)", ylab = "CAP2 (%)")


###################################################################################################
#STEP #14c: log-transformed biomass (manhattan distance)
y <- vector('numeric', length = 1000)
for (i in 1:1000)
{
  
gcap3 <- capscale(bl2 ~ FireRotation + BurnFreeInterval + Season + RxFireMgmt, 
                  data = siteEnv4, distance = "manhattan")

gcap3s <- summary(gcap3)
x <- anova(gcap3)#
y[i] <- x[1,5]
}
hist(y)
median(y)
range(y)
anova(gcap3, by = "axis")#
anova(gcap3, by = "terms")#

plot(gcap3, choices = c(1,2), type = "points", display = "wa", scaling = 2)
plot(gcap3, choices = c(1,2), type = "n", display = "lc", scaling = 2)
text(gcap3, choices = c(1,2), labels = row.names(bl2), cex = 0.8)
round(intrasetcor(gcap3), 5)
round(intersetcor(gcap3), 5)

plot(gcap3, choices = c(1,2), display = c("wa", "sp", "bp"), scaling = 3, 
     xlab = "CAP1 (%)", ylab = "CAP2 (%)")

###################################################################################################
#STEP #14d: log-transformed biomass (bray-curtis distance)
y <- vector('numeric', length = 1000)
for (i in 1:1000)
{
  
gcap4 <- capscale(bl2 ~ FireRotation + BurnFreeInterval + Season + RxFireMgmt, 
                  data = siteEnv4, distance = "bray")

gcap4s <- summary(gcap4)
x <- anova(gcap4)#
y[i] <- x[1,5]
}
hist(y)
median(y)
range(y)

anova(gcap4, by = "axis")#
anova(gcap4, by = "terms")#
plot(gcap4, choices = c(1,2), display = c("wa", "sp", "bp"), scaling = 3, 
     xlab = "CAP1 (%)", ylab = "CAP2 (%)")

###################################################################################################
#STEP #14e: untransformed cover (manhattan distance)
y <- vector('numeric', length = 1000)
for (i in 1:1000)
{
  gcap5 <- capscale(co2 ~ FireRotation + BurnFreeInterval + Season + RxFireMgmt, 
                    data = siteEnv4, distance = "manhattan")
  
  gcap5s <- summary(gcap5)
  
  x <- anova(gcap5)#
  y[i] <- x[1,5]
}
hist(y)
median(y)
range(y)

anova(gcap5, by = "axis")#
anova(gcap5, by = "term")#
plot(gcap5, choices = c(1,2), display = c("wa", "sp", "bp"), scaling = 3, 
     title = "Biomass by genera (model significance: p =)", xlab = "CAP1 (%)", 
     ylab = "CAP2 (%)")

###################################################################################################
#STEP #14f: untransformed cover (bray-curtis distance)
y <- vector('numeric', length = 1000)
for (i in 1:1000)
{
  gcap6 <- capscale(co2 ~ FireRotation + BurnFreeInterval + Season + RxFireMgmt, 
                    data = siteEnv4, distance = "bray")
  
  gcap6s <- summary(gcap6)
  
  x <- anova(gcap6)#
  y[i] <- x[1,5]
}
hist(y)
median(y)
range(y)

anova(gcap6, by = "axis")#
anova(gcap6, by = "term")#
plot(gcap6, choices = c(1,2), display = c("wa", "sp", "bp"), scaling = 3, 
     title = "Biomass by genera (model significance: p =)", xlab = "CAP1 (%)", 
     ylab = "CAP2 (%)")

###################################################################################################
#STEP #14g: orginal biomass data (Eglin) (bray-curtis distance)
y <- vector('numeric', length = 1000)
for (i in 1:1000)
{
  gcap7 <- capscale(boE ~ FireRotation + BurnFreeInterval + Season + RxFireMgmt, 
                    data = seE, distance = "bray")
  
  gcap7s <- summary(gcap7)
  
  x <- anova(gcap7)#
  y[i] <- x[1,5]
}
hist(y)
median(y)
range(y)

anova(gcap7, by = "axis")#
anova(gcap7, by = "term")#
plot(gcap7, choices = c(1,2), display = c("wa", "sp", "bp"), scaling = 3, 
     title = "Biomass by genera (model significance: p =)", xlab = "CAP1 (%)", 
     ylab = "CAP2 (%)")

###################################################################################################
#STEP #14h: orginal biomass data (Apalachicola) (bray-curtis distance)
y <- vector('numeric', length = 1000)
for (i in 1:1000)
{
  gcap8 <- capscale(boA ~ FireRotation + BurnFreeInterval + Season + RxFireMgmt, 
                    data = seA, distance = "bray")
  
  gcap8s <- summary(gcap8)
  
  x <- anova(gcap8)#
  y[i] <- x[1,5]
}
hist(y)
median(y)
range(y)

anova(gcap8, by = "axis")#
anova(gcap8, by = "term")#
plot(gcap8, choices = c(1,2), display = c("wa", "sp", "bp"), scaling = 3, 
     title = "Biomass by genera (model significance: p =)", xlab = "CAP1 (%)", 
     ylab = "CAP2 (%)")


#How does mFRI vary with CAP2 scores
plot(gcap1s$sites[,2], emat$FireRotation[order(gcap1s$sites[,2], decreasing = T)], type = "n")
text(gcap1s$sites[,2], emat$FireRotation[order(gcap1s$sites[,2], decreasing = T)], labels = rownames(emat))

model.1 <- lm(as.vector(emat$FireRotation[order(gcap1s$sites[,2], decreasing = T)] ~ gcap1s$sites[,2]))
summary(model.1)
abline(model.1)
segments(gcap1s$sites[,2], emat$FireRotation[order(gcap1s$sites[,2], decreasing = T)], 
         gcap1s$sites[,2], predict(model.1))

#How does FineWD vary with CAP2 scores
plot(gcap1s$sites[,1], emat$FineWD[order(gcap1s$sites[,1], decreasing = T)], type = "n")
text(gcap1s$sites[,1], emat$FineWD[order(gcap1s$sites[,1], decreasing = T)], labels = rownames(emat))

model.1 <- lm(as.vector(emat$FineWD[order(gcap1s$sites[,1], decreasing = T)] ~ gcap1s$sites[,1]))
summary(model.1)
abline(model.1)
segments(gcap1s$sites[,1], emat$FineWD[order(gcap1s$sites[,1], decreasing = T)], 
         gcap1s$sites[,1], predict(model.1))



###################################################################################################
#STEP 15: CREATE A BOUNDARY LAYER REGRESSION FUNCTION for Fire Rotation

###################################################################################################
############################################FUNCTION###############################################

#Create a function that will conduct a boundary layer regression with y variable as exponent,
#print summary of regression, and plot data
blr <- function(a,x,y,z)
{
cf <- 0.8
mv <- vector(mode = "numeric")
fv <- vector(mode = "numeric")
le <- vector(mode = "numeric")
  
fri <- 2:7
for(i in fri)
  {
    mv[i-1] <- max(y[x > i & x < i + 1])
    fv[i-1] <- min(x[x > i & x < i + 1 & y == mv[i-1]])
    le[i-1] <- length(y[x > i & x < i + 1])
    
    mv[i-1] <- mv[i-1] + 0.0001
  }

x2 <- fv
y2 <- mv  
  
d <- data.frame(x2,y2)
logmodel <- lm(y2~a(x2),data=d)


### fake vector
xvec <- seq(0,8, length=101)
logpred <- predict(logmodel, newdata=data.frame(x2=xvec))

#Plot with boundary layer regression
plot(x, y, xlab = "Fire Rotation (years)", ylab = "biomass (Mg/ha)", main = z)
points(fv, mv, pch = 22)
lines(xvec,logpred)
text(4.1, (max(y) - max(y)/16), paste("Intercept", round(logmodel$coefficients[1],4)), cex = cf)
text(5.9, (max(y) - max(y)/16), paste("Slope", round(logmodel$coefficients[2],4)), cex = cf)
ms <- summary(logmodel)
text(4.1, (max(y) - max(y)/8), paste("r-squared", round(ms$r.squared,2)), cex = cf)
text(5.9, (max(y) - max(y)/8), paste("P-value", round(ms$coefficients[8],4)), cex = cf)
print(summary(logmodel))
}

###################################################################################################
############################################FUNCTION###############################################

###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
#STEP 13: BOUNDARY LAYER REGRESSION FOR FIRE ROTATION

#Old file reassignments
emat <- siteEnv3
ab2mat <- biomassOrig

###################################################################################################
###################################################################################################
#PANEL PLOT BY ORIGINAL DATA __1
dev.off()
set.panel(2,4)
cf <- 0.5

#1-2
#Forb
#Plot with sites
genus <- ab2mat$forb
Genus <- "Forb"

#Boundary layer regression
blr(log, ev, genus, Genus)

#3-11
#Grass
#Plot with sites
genus <- ab2mat$grass + ab2mat$bluestem
Genus <- "Graminoid (ex. Aristida)"

#Boundary layer regression
blr(log, ev, genus, Genus)

#10
#Wiregrass
#Plot with sites
genus <- ab2mat$wiregrass
Genus <- "Wiregrass"

#Boundary layer regression
blr(log, ev, genus, Genus)


#12
#Yucca
#Plot with sites
genus <- ab2mat$yucca
Genus <- "Yucca"

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################

ev <- emat$RxFireMgmt
etext <- "Duration of Management (years)"

ev <- emat$BurnFreeInterval
etext = "Time-Since-Last-Fire (years)"

ev <- emat$Season
etext = "Percent Growing Season Burns"

ev <- emat$FireRotation
etext <- "mFRI (years)"


#1-2
#Forb
#Plot with sites
genus <- ab2mat$forb
Genus <- "Forb"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)


#3-11
#Grass
#Plot with sites
genus <- ab2mat$grass + ab2mat$bluestem
Genus <- "Graminoid (ex. Aristida)"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)


#10
#Wiregrass
#Plot with sites
genus <- ab2mat$wiregrass
Genus <- "Wiregrass"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)


#12
#Yucca
#Plot with sites
genus <- ab2mat$yucca
Genus <- "Yucca"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

###################################################################################################
###################################################################################################
#PANEL PLOT BY ORIGINAL DATA __2
dev.off()
set.panel(2,4)
cf <- 0.5

#4-5
#All Palmetto
#Plot with sites
genus <- ab2mat$palmetto
Genus <- "Palmetto"

#Boundary layer regression
blr(log, ev, genus, Genus)

#23
#Carpet Oak
#Plot with sites
genus <- ab2mat$oak.groundcover
Genus <- "Grouncover Oak"

#Boundary layer regression
blr(log, ev, genus, Genus)

#25
#Gallberry
#Plot with sites
genus <- ab2mat$gallberry
Genus <- "Gallberry"

#Linear regression
logmodel <- lm(genus ~ ev)
y <- genus
plot(ev, genus, 
     xlab = etext, ylab = "Biomass (Mg/ha)", main = "Gallberry")
text(4.1, (max(y) - max(y)/16), paste("Intercept", round(logmodel$coefficients[1],4)), cex = cf)
text(5.9, (max(y) - max(y)/16), paste("Slope", round(logmodel$coefficients[2],4)), cex = cf)
ms <- summary(logmodel)
text(4.1, (max(y) - max(y)/8), paste("r-squared", round(ms$r.squared,2)), cex = cf)
text(5.9, (max(y) - max(y)/8), paste("P-value", round(ms$coefficients[8],4)), cex = cf)
abline(logmodel)

#29
#Mound Oak
#Plot with sites
genus <- ab2mat$oak.shrub
Genus <- "Shrub Oak"

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#4-5
#All Palmetto
#Plot with sites
genus <- ab2mat$palmetto
Genus <- "Palmetto"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#23
#Carpet Oak
#Plot with sites
genus <- ab2mat$oak.groundcover
Genus <- "Grouncover Oak"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#25
#Gallberry
#Plot with sites
genus <- ab2mat$gallberry
Genus <- "Gallberry"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#29
#Mound Oak
#Plot with sites
genus <- ab2mat$oak.shrub
Genus <- "Shrub Oak"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

###################################################################################################
###################################################################################################
#PANEL PLOT BY ORIGINAL DATA __3
dev.off()
set.panel(2,4)
cf <- 0.5

#30
#Rough Fetterbush
#Plot with sites
genus <- ab2mat$rough.fetterbush
Genus <- "Rough Fetterbush"

#Boundary layer regression
blr(log, ev, genus, Genus)

#19
#Darrow's blueberry
#Plot with sites
genus <- ab2mat$Darrows.blueberry
Genus <- "Darrow's blueberry"

#Boundary layer regression
blr(log, ev, genus, Genus)

#35
#Highbush blueberry
#Plot with sites
genus <- ab2mat$highbush.blueberry
Genus <- "Highbush Blueberry"

#Boundary layer regression
blr(exp, ev, genus, Genus)

#43
#Sweetbay
#Plot with sites
genus <- ab2mat$sweetbay
Genus <- "Sweetbay"

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#30
#Rough Fetterbush
#Plot with sites
genus <- ab2mat$rough.fetterbush
Genus <- "Rough Fetterbush"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#19
#Darrow's blueberry
#Plot with sites
genus <- ab2mat$Darrows.blueberry
Genus <- "Darrow's blueberry"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#35
#Highbush blueberry
#Plot with sites
genus <- ab2mat$highbush.blueberry
Genus <- "Highbush Blueberry"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#43
#Sweetbay
#Plot with sites
genus <- ab2mat$sweetbay
Genus <- "Sweetbay"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

###################################################################################################
###################################################################################################
#PANEL PLOT BY ORIGINAL DATA __4
dev.off()
set.panel(2,4)
cf <- 0.5

#18
#Huckleberry
#Plot with sites
genus <- ab2mat$huckleberry
Genus <- "Huckleberry"

#Boundary layer regression
blr(exp, ev, genus, Genus)

#38
#Titi
#Plot with sites
genus <- ab2mat$titi
Genus <- "Titi"

#Boundary layer regression
blr(log, ev, genus, Genus)

#41
#Oak Tree
#Plot with sites
genus <- ab2mat$oak.tree
Genus <- "Oak Tree"

#Boundary layer regression
blr(log, ev, genus, Genus)

#15
#Greenbriar
#Plot with sites
genus <- ab2mat$greenbriar
Genus <- "Greenbriar"

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#18
#Huckleberry
#Plot with sites
genus <- ab2mat$huckleberry
Genus <- "Huckleberry"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#38
#Titi
#Plot with sites
genus <- ab2mat$titi
Genus <- "Titi"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#41
#Oak Tree
#Plot with sites
genus <- ab2mat$oak.tree
Genus <- "Oak Tree"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#15
#Greenbriar
#Plot with sites
genus <- ab2mat$greenbriar
Genus <- "Greenbriar"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

###################################################################################################
###################################################################################################
#PANEL PLOT BY ORIGINAL DATA __5
dev.off()
set.panel(2,4)
cf <- 0.5

###################################################################################################
#14
#Yellow Jessamine
#Plot with sites
genus <- ab2mat$yellow.jessamine
Genus <- "Yellow Jessamine"

#Boundary layer regression
blr(log, ev, genus, Genus)

#37
#Sparkleberry
#Plot with sites
genus <- ab2mat$sparkleberry
Genus <- "Sparkleberry"

#Boundary layer regression
blr(log, ev, genus, Genus)

#32
#Smooth Fetterbush
#Plot with sites
genus <- ab2mat$smooth.fetterbush
Genus <- "Smooth Fetterbush"

#Boundary layer regression
blr(exp, ev, genus, Genus)

#33
#Large Gallberry
#Plot with sites
genus <- ab2mat$large.gallberry
Genus <- "Large Gallberry"

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#14
#Yellow Jessamine
#Plot with sites
genus <- ab2mat$yellow.jessamine
Genus <- "Yellow Jessamine"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#37
#Sparkleberry
#Plot with sites
genus <- ab2mat$sparkleberry
Genus <- "Sparkleberry"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#32
#Smooth Fetterbush
#Plot with sites
genus <- ab2mat$smooth.fetterbush
Genus <- "Smooth Fetterbush"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#33
#Large Gallberry
#Plot with sites
genus <- ab2mat$large.gallberry
Genus <- "Large Gallberry"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

###################################################################################################
###################################################################################################
#PANEL PLOT BY ORIGINAL DATA __6
dev.off()
set.panel(2,4)
cf <- 0.5

#13
#Muscadine Grape
#Plot with sites
genus <- ab2mat$muscadine.grape
Genus <- "Muscadine Grape"

#Boundary layer regression
blr(log, ev, genus, Genus)

#36
#Yaupon
#Plot with sites
genus <- ab2mat$yaupon
Genus <- "Yaupon"

#Boundary layer regression
blr(log, ev, genus, Genus)

#42
#Sweetleaf
#Plot with sites
genus <- ab2mat$sweetleaf
Genus <- "Sweetleaf"

#Boundary layer regression
blr(log, ev, genus, Genus)

#21
#Wax Myrtle
#Plot with sites
genus <- ab2mat$wax.myrtle
Genus <- "Wax Myrtle"

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#13
#Muscadine Grape
#Plot with sites
genus <- ab2mat$muscadine.grape
Genus <- "Muscadine Grape"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#36
#Yaupon
#Plot with sites
genus <- ab2mat$yaupon
Genus <- "Yaupon"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#42
#Sweetleaf
#Plot with sites
genus <- ab2mat$sweetleaf
Genus <- "Sweetleaf"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#21
#Wax Myrtle
#Plot with sites
genus <- ab2mat$wax.myrtle
Genus <- "Wax Myrtle"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

###################################################################################################
###################################################################################################
#PANEL PLOT BY ORIGINAL DATA __7
dev.off()
set.panel(2,4)
cf <- 0.5

#17
#Wicky
#Plot with sites
genus <- ab2mat$wicky
Genus <- "Wicky"

#Boundary layer regression
blr(exp, ev, genus, Genus)

#20
#St. Andrew's Cross
#Plot with sites
genus <- ab2mat$St..Andrews.cross
Genus <- "St. Andrew's Cross"

#Boundary layer regression
blr(log, ev, genus, Genus)

#16
#Gopher Apple
#Plot with sites
genus <- ab2mat$gopher.apple
Genus <- "Gopher Apple"

#Boundary layer regression
blr(exp, ev, genus, Genus)

#31
#Sweet Pepperbush
#Plot with sites
genus <- ab2mat$sweet.pepperbush
Genus <- "Sweet Pepperbush"

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#17
#Wicky
#Plot with sites
genus <- ab2mat$wicky
Genus <- "Wicky"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#20
#St. Andrew's Cross
#Plot with sites
genus <- ab2mat$St..Andrews.cross
Genus <- "St. Andrew's Cross"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#16
#Gopher Apple
#Plot with sites
genus <- ab2mat$gopher.apple
Genus <- "Gopher Apple"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#31
#Sweet Pepperbush
#Plot with sites
genus <- ab2mat$sweet.pepperbush
Genus <- "Sweet Pepperbush"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

###################################################################################################
###################################################################################################
#PANEL PLOT BY ORIGINAL DATA __8
dev.off()
set.panel(2,4)
cf <- 0.5

#39
#Fringe Tree
#Plot with sites
genus <- ab2mat$fringe.tree
Genus <- "Fringe Tree"

#Boundary layer regression
blr(log, ev, genus, Genus)

#24
#Chokeberry
#Plot with sites
genus <- ab2mat$chokeberry
Genus <- "Chokeberry"

#Boundary layer regression
blr(log, ev, genus, Genus)

#22
#Conradina
#Plot with sites
genus <- ab2mat$conradina
Genus <- "Conradina"

#Boundary layer regression
blr(log, ev, genus, Genus)

#40
#American Holly
#Plot with sites
genus <- ab2mat$American.holly
Genus <- "American Holly"

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#39
#Fringe Tree
#Plot with sites
genus <- ab2mat$fringe.tree
Genus <- "Fringe Tree"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#24
#Chokeberry
#Plot with sites
genus <- ab2mat$chokeberry
Genus <- "Chokeberry"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#22
#Conradina
#Plot with sites
genus <- ab2mat$conradina
Genus <- "Conradina"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#40
#American Holly
#Plot with sites
genus <- ab2mat$American.holly
Genus <- "American Holly"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

###################################################################################################
###################################################################################################
#PANEL PLOT BY ORIGINAL DATA __9
dev.off()
set.panel(2,4)
cf <- 0.5

#26
#Dewberry
#Plot with sites
genus <- ab2mat$dewberry
Genus <- "Dewberry"

#Boundary layer regression
blr(log, ev, genus, Genus)

#27
#St. Johnswort
#Plot with sites
genus <- ab2mat$St..Johnswort
Genus <- "St. Johnswort"

#Boundary layer regression
blr(log, ev, genus, Genus)

#28
#Narrowleaf Paw Paw
#Plot with sites
genus <- ab2mat$narrowleaf.pawpaw
Genus <- "Narrowleaf Paw Paw"

#Boundary layer regression
blr(log, ev, genus, Genus)

#34
#Winged Sumac
#Plot with sites
genus <- ab2mat$winged.sumac
Genus <- "Winged Sumac"

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#26
#Dewberry
#Plot with sites
genus <- ab2mat$dewberry
Genus <- "Dewberry"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#27
#St. Johnswort
#Plot with sites
genus <- ab2mat$St..Johnswort
Genus <- "St. Johnswort"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#28
#Narrowleaf Paw Paw
#Plot with sites
genus <- ab2mat$narrowleaf.pawpaw
Genus <- "Narrowleaf Paw Paw"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#34
#Winged Sumac
#Plot with sites
genus <- ab2mat$winged.sumac
Genus <- "Winged Sumac"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

###################################################################################################
###################################################################################################
#PANEL PLOT BY ORIGINAL DATA __10
dev.off()
set.panel(2,4)
cf <- 0.5

#6
#Live Shrub
#Plot with sites
genus <- ab2mat$live.shrub
Genus <- "Live shrub"

#Boundary layer regression
blr(log, ev, genus, Genus)

#7
#Dead Shrub
#Plot with sites
genus <- ab2mat$dead.shrub
Genus <- "Dead shrub"

#Boundary layer regression
blr(log, ev, genus, Genus)

#8
#Live Tree
#Plot with sites
genus <- ab2mat$live.tree
Genus <- "Live tree"

#Boundary layer regression
blr(log, ev, genus, Genus)

#9
#Dead Tree
#Plot with sites
genus <- ab2mat$dead.tree
Genus <- "Dead tree"

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#6
#Live Shrub
#Plot with sites
genus <- ab2mat$live.shrub
Genus <- "Live shrub"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#7
#Dead Shrub
#Plot with sites
genus <- ab2mat$dead.shrub
Genus <- "Dead shrub"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#8
#Live Tree
#Plot with sites
genus <- ab2mat$live.tree
Genus <- "Live tree"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#9
#Dead Tree
#Plot with sites
genus <- ab2mat$dead.tree
Genus <- "Dead tree"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)






dev.off()



###################################################################################################
#1
#live.forb
#Plot with sites
genus <- ab2mat$live.forb
Genus <- "Live forb"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#2
#dead.forb
#Plot with sites
genus <- ab2mat$dead.forb
Genus <- "Dead forb"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#3
#grass
#Plot with sites
genus <- ab2mat$grass
Genus <- "Grass (ex. bunchgrass)"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#4
#Live Palmetto
#Plot with sites
genus <- ab2mat$live.palmetto
Genus <- "Live palmetto"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#5
#Dead Palmetto
#Plot with sites
genus <- ab2mat$dead.palmetto
Genus <- "Dead palmetto"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#5-7-9
#Dead Woody
#Plot with sites
genus <- ab2mat$dead.palmetto + ab2mat$dead.shrub + ab2mat$dead.tree
Genus <- "Dead Woody"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#Boundary layer regression
blr(log, ev, genus, Genus)

###################################################################################################
#11
#Bluestem
#Plot with sites
genus <- ab2mat$bluestem.grass
Genus <- "Bluestem"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#Boundary layer regression
blr(log, ev, genus, Genus)



##################################################################################################
#23-29-41
#Oak
#Plot with sites
genus <- ab2mat$oak.groundcover + ab2mat$oak.shrub + ab2mat$oak.tree
Genus <- "Oak"
plot(ev, genus, type = "n", 
     xlab = etext, ylab = "Biomass (Mg/ha)")
text(ev, genus, labels = row.names(ab2mat), cex = cf)

#Boundary layer regression
blr(log, ev, genus, Genus)

plot(0,0, main = "NA")

###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

###################################################################################################
#STEP 14: CREATE A BOUNDARY LAYER REGRESSION FUNCTION for Fire Rotation

###################################################################################################
############################################FUNCTION###############################################

#Create a function that will conduct a boundary layer regression with y variable as exponent,
#print summary of regression, and plot data
blr <- function(a,x,y,z)
{   
  cf <- 0.4
  mv <- vector(mode = "numeric")
  fv <- vector(mode = "numeric")
  le <- vector(mode = "numeric")
  
  fri <- seq(0, 1.5, 0.25)
  for(i in 1:length(fri))
  {
    mv[i] <- max(y[x > fri[i] & x < fri[i] + 0.25])
    fv[i] <- min(x[x > fri[i] & x < fri[i] + 0.25 & y == mv[i]])
    le[i] <- length(y[x > fri[i] & x < fri[i] + 0.25])
    
    mv[i] <- mv[i] + 0.0001
  }
  
  x2 <- fv
  y2 <- mv  
  
  d <- data.frame(x2,y2)
  logmodel <- lm(y2~a(x2),data=d)
  
  
  ### fake vector
  xvec <- seq(0,1.8, length=101)
  logpred <- predict(logmodel, newdata=data.frame(x2=xvec))
  
  #Plot with boundary layer regression
  plot(x, y, xlab = "Fine Woody Debris (Mg/ha)", ylab = "biomass (Mg/ha)", main = z)
  points(fv, mv, pch = 22)
  lines(xvec,logpred)
  text(0.6, (max(y) - max(y)/16), paste("Intercept", round(logmodel$coefficients[1],4)), cex = cf)
  text(1.2, (max(y) - max(y)/16), paste("Slope", round(logmodel$coefficients[2],4)), cex = cf)
  ms <- summary(logmodel)
  text(0.6, (max(y) - max(y)/8), paste("r-squared", round(ms$r.squared,2)), cex = cf)
  text(1.2, (max(y) - max(y)/8), paste("P-value", round(ms$coefficients[8],4)), cex = cf)
  print(summary(logmodel))
}

###################################################################################################
############################################FUNCTION###############################################


###################################################################################################
#STEP 14: BOUNDARY LAYER REGRESSION FOR FINE WOODY DEBRIS


###################################################################################################
###################################################################################################
#PANEL PLOT BY ORIGINAL DATA __1
dev.off()
set.panel(2,4)
cf <- 0.5

#Forb
#Plot with sites
genus <- ab3mat$forb
Genus <- "Forb"

#Boundary layer regression
blr(log, emat$FineWD, genus, Genus)

#Grass
#Plot with sites
genus <- ab3mat$grass
Genus <- "Graminoid (ex. Aristida)"

#Boundary layer regression
blr(exp, emat$FineWD, genus, Genus)

#Wiregrass
#Plot with sites
genus <- ab3mat$wiregrass
Genus <- "Wiregrass"

#Boundary layer regression
blr(exp, emat$FineWD, genus, Genus)

#Greenbriar
#Plot with sites
genus <- ab3mat$greenbriar
Genus <- "Greenbriar"

#Boundary layer regression
blr(exp, emat$FineWD, genus, Genus)

###################################################################################################
#Forb
#Plot with sites
genus <- ab3mat$forb
Genus <- "Forb"
plot(emat$FineWD, genus, type = "n", 
     xlab = "Fine Woody Debris (years)", ylab = "Biomass (Mg/ha)")
text(emat$FineWD, genus, labels = row.names(ab3mat), cex = cf)

#Grass
#Plot with sites
genus <- ab3mat$grass
Genus <- "Graminoid (ex. Aristida)"
plot(emat$FineWD, genus, type = "n", 
     xlab = "Fine Woody Debris (years)", ylab = "Biomass (Mg/ha)")
text(emat$FineWD, genus, labels = row.names(ab3mat), cex = cf)

#Wiregrass
#Plot with sites
genus <- ab3mat$wiregrass
Genus <- "Wiregrass"
plot(emat$FineWD, genus, type = "n", 
     xlab = "Fine Woody Debris (years)", ylab = "Biomass (Mg/ha)")
text(emat$FineWD, genus, labels = row.names(ab3mat), cex = cf)

#Greenbriar
#Plot with sites
genus <- ab3mat$greenbriar
Genus <- "Greenbriar"
plot(emat$FineWD, genus, type = "n", 
     xlab = "Fine Woody Debris (years)", ylab = "Biomass (Mg/ha)")
text(emat$FineWD, genus, labels = row.names(ab3mat), cex = cf)

###################################################################################################
###################################################################################################
#PANEL PLOT BY ORIGINAL DATA __2
dev.off()
set.panel(2,4)
cf <- 0.5

#4-5
#All Palmetto
#Plot with sites
genus <- ab3mat$palmetto
Genus <- "Palmetto"

#Boundary layer regression
blr(log, emat$FineWD, genus, Genus)

#Carpet Oak
#Plot with sites
genus <- ab3mat$oak.groundcover
Genus <- "Grouncover Oak"

#Boundary layer regression
blr(log, emat$FineWD, genus, Genus)


#Mound Oak
#Plot with sites
genus <- ab3mat$oak.shrub
Genus <- "Shrub Oak"

#Boundary layer regression
blr(exp, emat$FineWD, genus, Genus)

#Tree Oak
#Plot with sites
genus <- ab3mat$oak.tree
Genus <- "Tree Oak"

#Boundary layer regression
blr(exp, emat$FineWD, genus, Genus)

###################################################################################################
#All Palmetto
#Plot with sites
genus <- ab3mat$palmetto
Genus <- "Palmetto"
plot(emat$FineWD, genus, type = "n", 
     xlab = "Fine Woody Debris (years)", ylab = "Biomass (Mg/ha)")
text(emat$FineWD, genus, labels = row.names(ab3mat), cex = cf)

#Carpet Oak
#Plot with sites
genus <- ab3mat$oak.groundcover
Genus <- "Grouncover Oak"
plot(emat$FineWD, genus, type = "n", 
     xlab = "Fine Woody Debris (years)", ylab = "Biomass (Mg/ha)")
text(emat$FineWD, genus, labels = row.names(ab3mat), cex = cf)


#Mound Oak
#Plot with sites
genus <- ab3mat$oak.shrub
Genus <- "Shrub Oak"
plot(emat$FineWD, genus, type = "n", 
     xlab = "Fine Woody Debris (years)", ylab = "Biomass (Mg/ha)")
text(emat$FineWD, genus, labels = row.names(ab3mat), cex = cf)

#Tree Oak
#Plot with sites
genus <- ab3mat$oak.tree
Genus <- "Tree Oak"
plot(emat$FineWD, genus, type = "n", 
     xlab = "Fine Woody Debris (years)", ylab = "Biomass (Mg/ha)")
text(emat$FineWD, genus, labels = row.names(ab3mat), cex = cf)
###################################################################################################
###################################################################################################
#PANEL PLOT BY ORIGINAL DATA __3
dev.off()
set.panel(2,4)
cf <- 0.5

#Darrow's blueberry
#Plot with sites
genus <- ab3mat$Darrows.blueberry
Genus <- "Darrow's blueberry"

#Boundary layer regression
blr(log, emat$FineWD, genus, Genus)

#Highbush blueberry
#Plot with sites
genus <- ab3mat$highbush.blueberry
Genus <- "Highbush Blueberry"

#Boundary layer regression
blr(exp, emat$FineWD, genus, Genus)

#Huckleberry
#Plot with sites
genus <- ab3mat$huckleberry
Genus <- "Huckleberry"

#Boundary layer regression
blr(log, emat$FineWD, genus, Genus)

#Gallberry
#Plot with sites
genus <- ab3mat$gallberry
Genus <- "Gallberry"

#Linear regression
logmodel <- lm(genus ~ emat$FineWD)
y <- genus
plot(emat$FineWD, genus, 
     xlab = "Fine Woody Debris (years)", ylab = "Biomass (Mg/ha)", main = "Gallberry")
text(4.1, (max(y) - max(y)/16), paste("Intercept", round(logmodel$coefficients[1],4)), cex = cf)
text(5.9, (max(y) - max(y)/16), paste("Slope", round(logmodel$coefficients[2],4)), cex = cf)
ms <- summary(logmodel)
text(4.1, (max(y) - max(y)/8), paste("r-squared", round(ms$r.squared,2)), cex = cf)
text(5.9, (max(y) - max(y)/8), paste("P-value", round(ms$coefficients[8],4)), cex = cf)
abline(logmodel)


###################################################################################################
#Darrow's blueberry
#Plot with sites
genus <- ab3mat$Darrows.blueberry
Genus <- "Darrow's blueberry"
plot(emat$FineWD, genus, type = "n", 
     xlab = "Fine Woody Debris (years)", ylab = "Biomass (Mg/ha)")
text(emat$FineWD, genus, labels = row.names(ab3mat), cex = cf)

#Highbush blueberry
#Plot with sites
genus <- ab3mat$highbush.blueberry
Genus <- "Highbush Blueberry"
plot(emat$FineWD, genus, type = "n", 
     xlab = "Fine Woody Debris (years)", ylab = "Biomass (Mg/ha)")
text(emat$FineWD, genus, labels = row.names(ab3mat), cex = cf)

#Huckleberry
#Plot with sites
genus <- ab3mat$huckleberry
Genus <- "Huckleberry"
plot(emat$FineWD, genus, type = "n", 
     xlab = "Fine Woody Debris (years)", ylab = "Biomass (Mg/ha)")
text(emat$FineWD, genus, labels = row.names(ab3mat), cex = cf)

#Gallberry
#Plot with sites
genus <- ab3mat$gallberry
Genus <- "Gallberry"
plot(emat$FineWD, genus, type = "n", 
     xlab = "Fine Woody Debris (years)", ylab = "Biomass (Mg/ha)")
text(emat$FineWD, genus, labels = row.names(ab3mat), cex = cf)

###################################################################################################
###################################################################################################
#PANEL PLOT BY ORIGINAL DATA __4
dev.off()
set.panel(2,4)
cf <- 0.5

#Dead Woody
#Plot with sites
genus <- ab3mat$dead.woody
Genus <- "Dead woody"

#Boundary layer regression
blr(exp, emat$FineWD, genus, Genus)

#Dead Woody
#Plot with sites
genus <- ab3mat$dead.woody
Genus <- "Dead woody"
plot(emat$FineWD, genus, type = "n", 
     xlab = "Fine Woody Debris (years)", ylab = "Biomass (Mg/ha)")
text(emat$FineWD, genus, labels = row.names(ab3mat), cex = cf)


















###################################################################################################
#1
#live.forb
#Plot with sites
genus <- ab3mat$live.forb
Genus <- "Live forb"
plot(emat$FineWD, genus, type = "n", 
     xlab = "Fine Woody Debris (years)", ylab = "Biomass (Mg/ha)")
text(emat$FineWD, genus, labels = row.names(ab3mat), cex = cf)

#Boundary layer regression
blr(log, emat$FineWD, genus, Genus)

###################################################################################################
#2
#dead.forb
#Plot with sites
genus <- ab3mat$dead.forb
Genus <- "Dead forb"
plot(emat$FineWD, genus, type = "n", 
     xlab = "Fine Woody Debris (years)", ylab = "Biomass (Mg/ha)")
text(emat$FineWD, genus, labels = row.names(ab3mat), cex = cf)

#Boundary layer regression
blr(log, emat$FineWD, genus, Genus)

###################################################################################################
#3
#grass
#Plot with sites
genus <- ab3mat$grass
Genus <- "Grass (ex. bunchgrass)"
plot(emat$FineWD, genus, type = "n", 
     xlab = "Fine Woody Debris (years)", ylab = "Biomass (Mg/ha)")
text(emat$FineWD, genus, labels = row.names(ab3mat), cex = cf)

#Boundary layer regression
blr(log, emat$FineWD, genus, Genus)

###################################################################################################
#4
#Live Palmetto
#Plot with sites
genus <- ab3mat$live.palmetto
Genus <- "Live palmetto"
plot(emat$FineWD, genus, type = "n", 
     xlab = "Fine Woody Debris (years)", ylab = "Biomass (Mg/ha)")
text(emat$FineWD, genus, labels = row.names(ab3mat), cex = cf)

#Boundary layer regression
blr(log, emat$FineWD, genus, Genus)

###################################################################################################
#5
#Dead Palmetto
#Plot with sites
genus <- ab3mat$dead.palmetto
Genus <- "Dead palmetto"
plot(emat$FineWD, genus, type = "n", 
     xlab = "Fine Woody Debris (years)", ylab = "Biomass (Mg/ha)")
text(emat$FineWD, genus, labels = row.names(ab3mat), cex = cf)

#Boundary layer regression
blr(log, emat$FineWD, genus, Genus)

###################################################################################################
#5-7-9
#Dead Woody
#Plot with sites
genus <- ab3mat$dead.palmetto + ab3mat$dead.shrub + ab3mat$dead.tree
Genus <- "Dead Woody"
plot(emat$FineWD, genus, type = "n", 
     xlab = "Fine Woody Debris (years)", ylab = "Biomass (Mg/ha)")
text(emat$FineWD, genus, labels = row.names(ab3mat), cex = cf)

#Boundary layer regression
blr(log, emat$FineWD, genus, Genus)

###################################################################################################
#11
#Bluestem
#Plot with sites
genus <- ab3mat$bluestem.grass
Genus <- "Bluestem"
plot(emat$FineWD, genus, type = "n", 
     xlab = "Fine Woody Debris (years)", ylab = "Biomass (Mg/ha)")
text(emat$FineWD, genus, labels = row.names(ab3mat), cex = cf)

#Boundary layer regression
blr(log, emat$FineWD, genus, Genus)



##################################################################################################
#23-29-41
#Oak
#Plot with sites
genus <- ab3mat$oak.groundcover + ab3mat$oak.shrub + ab3mat$oak.tree
Genus <- "Oak"
plot(emat$FineWD, genus, type = "n", 
     xlab = "Fine Woody Debris (years)", ylab = "Biomass (Mg/ha)")
text(emat$FineWD, genus, labels = row.names(ab3mat), cex = cf)

#Boundary layer regression
blr(log, emat$FineWD, genus, Genus)

plot(0,0, main = "NA")

###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

######################################
#Looking at bunchgrasses and their relationship to herbaceous loading/
lb$SiteName
reg <- c(2,1,1,1,2,2,1,2,1,2,1,1,1,2,2,1,2,1,2,2,2,1)
boxplot(lb$bunchgrass ~ region)
boxplot(lb$herb ~ region)
hbr <- lb$bunchgrass/(lb$bunchgrass+lb$herb)
boxplot(hbr ~ region)
hist(hbr[region == 1])
hbrL <- data.frame(site = I(lb$SiteName), hbr = round(hbr,2), bunchgrass = lb$bunchgrass, region = reg)
hbrL <- hbrL[order(hbrL$bunchgrass),]

fill = vector(mode = "character")
for(i in 1:length(hbrL$region))
{
  if(hbrL$region[i] == 1)
  {
    fill[1+length(fill)] <- "dark green" 
    fill[1+length(fill)] <- "light green"  
  }
  else
  {
    fill[1+length(fill)] <- "red"
    fill[1+length(fill)] <- "orange"  
  }
}



barplot(as.matrix(t(hbrL[,2:3])), main = "ratio of bunchgrass to herbaceous material", 
        ylab = "", xlab = "", axes = F, beside = T, col = fill, labels = F)
axis(2)
text(matrix(c(seq(1.8,(nrow(hbrL))*3.0,3.0), rep(-0.02,nrow(hbrL))), nrow = nrow(hbrL), 
            ncol = 2, byrow = F), srt = 60, adj = 1, xpd = T, labels = paste(hbrL$site), cex = 0.65)


ma <- vector(mode = "numeric")
for(i in 1:length(hbrL$bunchgrass))
{
  ma[i] <- mean(hbrL$hbr[i-2:i])
}

points(ma)

lm2 <- lm(ma)


mean(c(43,45,4848,48,4848))


######################################
#Looking at fire rotation and it's relationship to FineWD loading
dev.off()
plot(emat$FireRotation, emat$FineWD, type = "n")
text(emat$FireRotation, emat$FineWD, labels = rownames(emat))
ffm <- lm(emat$FineWD ~ emat$FireRotation)
summary(ffm)
abline(ffm)
segments(emat$FireRotation, emat$FineWD, emat$FireRotation, predict(ffm))
cor(emat$FireRotation, emat$FineWD, method = "pearson")


###################################################################################################
###################################################################################################
#TUTORIAL

###################################################################################################
###################################################################################################
#STEP #XX: PCA of the covariance matrix
pca <- rda(decostand(bl, method = "hellinger", scale = T))#scale = T should be PCA of the correlation matrix
#"Inertia is correlation instead of Inertia is variance
pca
#Plot pca
biplot(pca, scaling = 1)
#Use brokan stick method to determine how many PCA axes to retain
screeplot(pca, bstick = T, type = "l", main = NULL)
#Extract eigenvalues
eigenvals(pca)
summary(eigenvals(pca))


data(varespec)
head(varespec)
pca <- rda(decostand(varespec, method = "hellinger", scale = F))
pca

#Fitting environmental data
set.seed(42)#use this to make this permutation analysis reproducibe, otherwise it will be different
#every time, especially when number of permutations is low.
ev <- envfit(pca ~ ., data = siteEnv3, choices = 1:2, scaling = 'symmetric', permutations = 1000)
#options for scaling: "species", "sites", "symmetric" (scales both for species and sites), "none" (raw scores).
#cca() hill = T, rda(): correlation = T (standardizes scores among species).
ev

plot(pca, display = "sites", type = "n", scaling = "symmetric")
points(pca, display = "sites", scaling = "symmetric")
plot(ev, add = T)

surf <- ordisurf(pca ~ FineWD, data = siteEnv3, knots = 10, isotropic = T, main = NULL)
surf <- ordisurf(pca ~ CoarseWD, data = siteEnv3, knots = 10, isotropic = T, main = NULL)
surf <- ordisurf(pca ~ RxFireMgmt, data = siteEnv3, knots = 10, isotropic = T, main = NULL)
surf <- ordisurf(pca ~ Litter, data = siteEnv3, knots = 10, isotropic = T, main = NULL)
surf <- ordisurf(pca ~ Canopy, data = siteEnv3, knots = 10, isotropic = T, main = NULL)
surf <- ordisurf(pca ~ mfri_20yr, data = siteEnv3, knots = 10, isotropic = T, main = NULL)

summary(surf)
colnames(siteEnv3)

###################################################################################################
###################################################################################################
#STEP #XX: CA
ca <- cca(bo)#default scaling is 1.
plot(ca)

ca_test <- cca(varespec)
plot(ca_test)

###################################################################################################
###################################################################################################
#STEP #XX: CCA
m1 <- cca(bo ~ ., data = siteEnv3)
m1
summary(m1)
set.seed(45)#use this so results are reproducible
pstat <- permustats(anova(m1))#anova is testing the entire model for significance.
summary(pstat)
densityplot(pstat)
#about 98% of the permuted F values were lower than the observed F statistic
#this is good, we want the Observed F to be in the extreme upper tail, i.e., we want it to be big.

set.seed(45)
perm <- anova(m1)
perm
#P-value is 0.017, about 98% of the permuted values had an F-statistics that were smaller than observed.

#Arguments for constrained ordination with cca
args(anova.cca)
#by arg specifies what kind of test to use.
#there are four. Be careful with sequential testing because order of terms added can change if 
#terms are not orthogonal.
#by axis -- tests by first axis, then second axis is tested.
#how much extra info can be explained by additional axes.

upr <- cca(bo ~ ., data = siteEnv3)
lwr <- cca(bo ~ 1, data = siteEnv3)
set.seed(1)
mods <- ordistep(lwr, scope = formula(upr), trace = 0)

set.seed(45)
anova(mods, by = "axis")
#by "terms" can result in variables having different significance depending on the order they are added. Be careful with this unless you have a robust experiment.
#by "margin" is more helpful
#effect of variables given other variables.
anova(mods, by = "margin")

#use decorana() to check gradient length.
#cca is good for a dataset with a long gradient.
plot(m1)
#no arch in site distribution across CCA1 -- indicates shorter gradient length.

#Forward selection model building
m2 <- ordistep(lwr, scope = formula(m1), trace = F)
m2
#With two of the environmental variables (RxFireMgmt and FineWD we can explain about 28% of the data)
#The first axis (CCA1) explains about 65% of the variation in the model (CCA1 Eigenvalue: 0.12763 / Constrained Inertia: 0.1956)
0.12763/0.1956

#Still a lot of structure in the data not explained by the constrained model. 
#First unconstrained axis explains more of the variation than the first constrained axis
#The first 3 unconstrained axes esplain more variation than any of the two constrained axes.

plot(m2)
#RxFireMgmt and FineWD are nearly at 45 degs from the CCA axes so difficult to say which they best line up with.
m2$anova
#RxFireMgmt was the first term to be selected
#FineWD was the second, and last term.

#Can also use rda() instead of cca() if you are worried about the arch effect.
#Do a Hellinger transformation on species data
spph <- decostand(bo, method = "hellinger")
m3 <- rda(spph ~ ., data = siteEnv3)
lwr <- rda(spph ~ 1., data = siteEnv3)
m4 <- ordistep(lwr, scope = formula(m3), trace = F)
m4#similar output to cca
plot(m4)#similar to cca
#species are sort of squished into the center of the plot. This is not an absolute product of the model,
#it's their relative positions to each other, the sites, and the axes that matter.
#You can spread them out by changing your scaling.

#Stepwise with adjusted R2
m5 <- ordiR2step(lwr, scope = formula(m3), trace = F)
m5$anova

m6 <- ordiR2step(lwr, scope = formula(m1), trace = F)
m6$anova

#Restricted permutation tests
#need to do this to account for dependencies among observations
#in this case you have spatial correlation because half of the sites
#are at Eglin and half are at Appalachicola/St. Marks
#permute library is a dependency on vegan so you don't need to load it seperately from vegan
#hierarchy
#sample level (plots for your study)
#plot level (sites for your study)
#block level (region for your study)
#blocks are not permuted
#sites are not exchanged between plots.
#variation between blocks must be excluded from tests
#Use + Condition(blocks)
#function for changing data type
#example: env <- transform(env, year = as.numeric(as.character(year)))
#useful when R imports data in the wrong format.
env <- transform(siteEnv3, Region = as.factor(Region))
c1 <- rda(bo ~ Litter + FineWD + mfri_20yr + StDevFRI_20YR + Season_20yr + RxFireMgmt + Condition(Region), data = env)
h <- how(within = Within(type = "free"), plots = Plots(strata = env$Region))
set.seed(42)
tt <- anova(c1, permutations = h, model = "reduced")
tt
summary(tt)
summary(c1)
c1$anova


length(env$Region)
length(bo[,1])

###################################################################################################
###################################################################################################
#STEP #XX: RDA
rda1 <- rda(bo ~ ., data = siteEnv3)
rda1
#Inertia is variance, this is because you did not standardize species variables (from tutorial).
#shouldn't be the env. variables that get standardized? If it was standardized Total = 1.00.



plot(ca)

ca_test <- cca(varespec)
plot(ca_test)
