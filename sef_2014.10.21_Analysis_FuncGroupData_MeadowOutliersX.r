###################################################################################################
############################-----start-----########################################################
#The purpose of this script is to explore the environmental data for the community ecology
#multivariate analysis.

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
library(simba)
library(fields)#for set.panel()
library(labdsv)#Brooke's recommendation


###################################################################################################
###################################################################################################
#STEP #2: OPEN AND ADJUST BIOMASS DATA


#########################################################
#2A: Open biomass data

#Plot-level biomass data. This data has been modified so that vine species are consistent
#across all sites and outlier plots located in uncharacteristics areas of sites (wetlands or meadows) 
#have been removed.
plotBiomass <- read.table(
  "C:/usfs_sef_data_output/sef_Ecology_BiomassPlotMatrix_Ep1_FuncGroupOulierX_2014-10-21_10.13.08.csv", 
  header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE,
  stringsAsFactors = F)

#########################################################
#2B The steps below create a more appropriate dataset by removing and combining some categories.

#For cover, remove non-plant categories, these are unecessary for our analysis.
#dead woody, 
plotBiomass2 <- plotBiomass[,!colnames(plotBiomass) == "dead.woody"]


#########################################################
#2C: calculate site means for untransformed plot-level data.
#USED FOR DATA SCREENING ONLY - DATA NOT USED IN ANALYSIS.

#Calculate untransformed site means from plot-level data 
siteBiomass <- summarize(X = plotBiomass2[,3:length(plotBiomass2[1,])], by = plotBiomass2$siteName, 
                         colMeans, stat.name = colnames(plotBiomass2)[3])
colnames(siteBiomass)[1] <- "siteName"

#Remove site names from matrix columns and assign as row names.
#All biomass >> UNTRANSFORMED
siteBiomass2 <- siteBiomass[,2:(length(siteBiomass[1,]))]
rownames(siteBiomass2) <- siteBiomass[,1]

#########################################################
#2d: Conduct log transformation on plot-level data and then calculate site means.

#Conduct log-transformation on plot-level data
#>>>>Justification for log-transformation?
plotBiomassLogTrans <- data.trans(plotBiomass2[,3:length(plotBiomass2[1,])], method = 'log', 
                                  plot = F)

#Calculate site means from log-transformed plot-level data (USED IN ANALYSIS)
siteBiomassLogTrans <- summarize(X = plotBiomassLogTrans, by = plotBiomass$siteName, colMeans, 
                                 stat.name = colnames(plotBiomass)[3])
colnames(siteBiomassLogTrans)[1] <- "siteName"

#Remove site names from matrix columns and assign as row names.
#All biomass >> LOG-TRANSFORMED
siteBiomassLogTrans2 <- siteBiomassLogTrans[,2:(length(siteBiomassLogTrans[1,]))]
rownames(siteBiomassLogTrans2) <- siteBiomassLogTrans[,1]

###################################################################################################
###################################################################################################
#STEP #3: OPEN AND ADJUST COVER DATA

#########################################################
#3A: Open cover data

#Plot-level cover data. This data has been modified so that vine species are consistent
#across all sites and outlier plots located in uncharacteristics areas of sites (wetlands or meadows) 
#have been removed.
plotCover <- read.table(
  "C:/usfs_sef_data_output/sef_Ecology_CoverPlotMatrix_Ep1_FuncGroupOulierX_2014-10-21_10.13.08.csv", 
  header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE,
  stringsAsFactors = F)

#########################################################
#3B The steps below create a more appropriate dataset by removing and combining some categories.

#For cover, remove non-plant categories, these are unecessary for our analysis.
#litter, bare.soil, and dead woody. 
plotCover2 <- plotCover[,!colnames(plotCover) %in% c("litter", "soil", "dead.woody")]

#########################################################
#3C: calculate site means for untransformed plot-level data.
#USED FOR DATA SCREENING ONLY - DATA NOT USED IN ANALYSIS.

#Calculate untransformed site means from plot-level data 
siteCover <- summarize(X = plotCover2[,3:length(plotCover2[1,])], by = plotCover2$siteName, 
                       colMeans, stat.name = colnames(plotCover2)[3])
colnames(siteCover)[1] <- "siteName"

#Remove site names from matrix columns and assign as row names.
#All cover >> UNTRANSFORMED
siteCover2 <- siteCover[,2:(length(siteCover[1,]))]
rownames(siteCover2) <- siteCover[,1]

#########################################################
#3d: Conduct transformation on plot-level data and then calculate site means?
#Is this necessary for cover data. This data is bounded (0-100) so a log-trans
#may not be appropriate if data is not normal.

###################################################################################################
###################################################################################################
#STEP #4: OPEN AND ADJUST ENVIRONMENTAL DATA

#########################################################
#4a: Open environmental matrix
siteEnv <- read.table(
  "C:/usfs_sef_data_output/2014.03.13_EnvironmentalMatrix.csv", 
  header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE,
  stringsAsFactors = F)

#########################################################
#4b: Re-assign column 1 (site names) to  a row name.
siteEnv2 <- siteEnv[,2:(length(siteEnv[1,]))]
rownames(siteEnv2) <- siteEnv[,1]

#########################################################
#4c: Environmental data, remove data you will not analyze
siteEnv3 <- siteEnv2[,-c(12,13)]
#Column 12 is soil drainage. This data is from USDA Soil Survey and not accurate
#Column 13 is region. Only two regions, at this point I don't believe this adds 
#substance to the analysis because there are no physical properties for each region
#that can confound how fire characteristics drive understory plant composition

###################################################################################################
###################################################################################################
#STEP #5: STANDARDIZE ROW ORDER OF DATA SETS (I.E. ORDER SITE NAMES (ROW NAMES) OF BIOLOGICAL DATA 
#ACCORDING TO ORDER IN ENVIRONMENTAL DATA SET) AND REMOVE COLUMNS WITH NO VALUES > 0. 

#Reorder biological data (biomass) so sites are in same order as environmental data
#Step not necessary for cover data... sites are already in same order as environmental data.
ro <- data.frame(siteNo = 1:length(rownames(siteEnv3)), siteName = I(rownames(siteEnv3)))
ro2 <- ro[order(ro[,2]),]

#Untransformed biomass data
siteBiomass3 <- data.frame(so = ro2[,1], siteBiomass2)
siteBiomass4 <- siteBiomass3[order(siteBiomass3$so),]
siteBiomass5 <- siteBiomass4[,!colnames(siteBiomass4) %in% c("so")]

#Log-tranformed biomass data
siteBiomassLogTrans3 <- data.frame(so = ro2[,1], siteBiomassLogTrans2)
siteBiomassLogTrans4 <- siteBiomassLogTrans3[order(siteBiomassLogTrans3$so),]
siteBiomassLogTrans5 <- siteBiomassLogTrans4[,!colnames(siteBiomassLogTrans4) %in% c("so")]

#Remove species with no occurences, this is necessary before you can use foa.plots below
biomassOrig <- drop.var(siteBiomass5, min.fo = 1)
biomassLog <- drop.var(siteBiomassLogTrans5, min.fo = 1)
coverOrig <- drop.var(siteCover2, min.fo = 1)


###################################################################################################
###################################################################################################
#STEP #6: DATA SCREENING

#########################################################
#6a: Summary stats for biomass data

#Summary stats
round(stat.desc(biomassOrig),2)

#Look at how data is distributed
foa.plots(biomassOrig)

#Drop rare species < 10% of sites (this removes any species that occur at 2 sites or less).
#bo <- drop.var(biomassOrig, min.po=10)
#bl <- drop.var(biomassLog, min.po=10)
#Do not drop functional groups
bo <- biomassOrig
bl <- biomassLog

#Check to see how many species were dropped for untransformed biomass data

#length(biomassOrig[1,])#
#length(bo[1,])#
#length(biomassLog[1,])#
#length(bl[1,])
#40 to 33 variables (7 drops)

str(bo)#22 obs and 8 variables
str(bl)#22 obs and 8 variables

#Convert data.frame to a matrix (needed for uv.plots)
bo2 <- data.matrix(frame = bo, rownames.force = NA)
bl2 <- data.matrix(frame = bl, rownames.force = NA)


#uv.plots displays histogram, box and whisker, cumulative distribution and normal q-q plots
#in one pane for each variable.
uv.plots(bo2)
uv.plots(bl2)

#What percent of values are zero, if it is over 50% you should consider changing the data
#to presence/absence
length(bo2[bo2 == 0])/(length(bo2[1,])*length(bo2[,1]))
#7% zero values >>> do not do binary transformation
length(bl2[bl2 == 0])/(length(bl2[1,])*length(bl2[,1]))
#7% zero values >>> do not do binary transformation

#Look at correlation between species variables.
chart.Correlation(bo2, method = "pearson")
chart.Correlation(bl2, method = "pearson")

#########################################################
#6b: Summary stats for cover data

#Summary stats
round(stat.desc(coverOrig),2)

#Look at how data is distributed
foa.plots(coverOrig)

#Drop rare species < 10% of sites (this removes any species that occur at 2 sites or less).
#co <- drop.var(coverOrig, min.po=10)
#Do not drop functional groups
co <- coverOrig

#Check to see how many species were dropped for untransformed biomass data
#length(coverOrig[1,])#
#length(co[1,])#
#46 to 35 variables (11 drops)

str(co)#22 obs and 8 variables

#Convert data.frame to a matrix (needed for uv.plots)
co2 <- data.matrix(frame = co, rownames.force = NA)

#uv.plots displays histogram, box and whisker, cumulative distribution and normal q-q plots
#in one pane for each variable.
uv.plots(co2)

#What percent of values are zero, if it is over 50% you should consider changing the data
#to presence/absence
length(co2[co2 == 0])/(length(co2[1,])*length(co2[,1]))
#45% zero values >>> do not do binary transformation

#Look at correlation between species variables.
chart.Correlation(co2, method = "pearson")

#########################################################
#6c: Summary stats for environmental data

#Show how environmental variables are correlated.
chart.Correlation(siteEnv3)

#Effect of column-standardization on untransformed data.
siteEnv4 <- data.stand(siteEnv3, method = 'total', margin = 'column', plot = F)

###################################################################################################
###################################################################################################
#STEP #7: DATA STANDARDIZATION

#Not necessary to examine effect of column-standardizations on data. Units are the same. 

#Don't want to do a standardization by row but here's the effects on data.


###################################################################################################
###################################################################################################
#STEP #8: OUTLIERS

#Find univariate outliers
uv.outliers(bo2, id = names(bo2)[1]:names(bo2)[length(names(bo2))], sd.limit = 1)
#did not figure out this function because this is not so important in multivariate space.
#I don't think Euclidean distance is appropriate for biological data, should I be using
#a different method? At any rate, 3 is the recommended sd.limit. None of your sites
#exceeds 3 for the Euclidean method.

mv.outliers(bo2, method = 'mahalanobis', sd.limit=1)#E103B_S3
mv.outliers(bl2, method = 'mahalanobis', sd.limit=1)#A18
mv.outliers(co2, method = 'mahalanobis', sd.limit=1)#E403B

###################################################################################################
###################################################################################################
#STEP #9: PCoA


###################################################################################################
#STEP #9a: Principal Coordinate Analysis on biomass data.

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
#STEP #9b: Principal Coordinate Analysis on original biomass data for log-transformed data.

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
#STEP #10: Gradient length diagnostics with DCA for untransformed biomass data, log-transformed
#biomass data, and untransformed cover data.

#Create a bindary dataset for species presence absence
#Biomass - untransformed
speocc <- data.trans(bo2, method = "power", exp = 0, plot = F)

decorana(speocc, ira = 0)
#Axis length is 0.74344

#Biomass - log-transformed
speocc <- data.trans(bl2, method = "power", exp = 0, plot = F)

decorana(speocc, ira = 0)
#Axis length is 0.74344

#Cover - intransformed
speocc <- data.trans(co2, method = "power", exp = 0, plot = F)

decorana(speocc, ira = 0)
#Axis length is 0.59333

#Short < 1.5 gradient lengths

###################################################################################################
###################################################################################################
#STEP #11: Detect and remove(?) outliers

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
#STEP #12: Direct gradient analysis with CAP


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
#STEP #13: #Distance-based RDA

###################################################################################################
#STEP #13a: untransformed biomass (manhattan distance)
y <- vector('numeric', length = 1000)
for (i in 1:1000)
{
gcap1 <- capscale(bo2 ~ FireRotation + BurnFreeInterval + Season + RxFireMgmt, 
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
#STEP #14b: untransformed biomass (bray-curtis distance)
y <- vector('numeric', length = 1000)
for (i in 1:1000)
{
  
gcap2 <- capscale(bo2 ~ FireRotation + BurnFreeInterval + Season + RxFireMgmt, 
                  data = siteEnv4, distance = "bray")

gcap2s <- summary(gcap2)
x <- anova(gcap1)#
y[i] <- x[1,5]
}
hist(y)
median(y)
range(y)

anova(gcap2, by = "axis")#
anova(gcap2, by = "term")#
plot(gcap2, choices = c(1,2), display = c("wa", "sp", "bp"), scaling = 3, 
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
###################################################################################################
###################################################################################################
###################################################################################################



###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################



#ANALYSIS WITH FUNCTIONAL GROUPS STOPS HERE...


###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################


###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

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
blr(log, emat$FireRotation, genus, Genus)

#3-11
#Grass
#Plot with sites
genus <- ab2mat$grass + ab2mat$bluestem
Genus <- "Graminoid (ex. Aristida)"

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

#10
#Wiregrass
#Plot with sites
genus <- ab2mat$wiregrass
Genus <- "Wiregrass"

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)


#12
#Yucca
#Plot with sites
genus <- ab2mat$yucca
Genus <- "Yucca"

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

###################################################################################################
#1-2
#Forb
#Plot with sites
genus <- ab2mat$dead.forb + ab2mat$live.forb
Genus <- "Forb"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)


#3-11
#Grass
#Plot with sites
genus <- ab2mat$grass + ab2mat$bluestem
Genus <- "Graminoid (ex. Aristida)"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)


#10
#Wiregrass
#Plot with sites
genus <- ab2mat$wiregrass
Genus <- "Wiregrass"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)


#12
#Yucca
#Plot with sites
genus <- ab2mat$yucca
Genus <- "Yucca"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

###################################################################################################
###################################################################################################
#PANEL PLOT BY ORIGINAL DATA __2
dev.off()
set.panel(2,4)
cf <- 0.5

#4-5
#All Palmetto
#Plot with sites
genus <- ab2mat$dead.palmetto + ab2mat$live.palmetto
Genus <- "Palmetto"

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

#23
#Carpet Oak
#Plot with sites
genus <- ab2mat$oak.groundcover
Genus <- "Grouncover Oak"

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

#25
#Gallberry
#Plot with sites
genus <- ab2mat$gallberry
Genus <- "Gallberry"

#Linear regression
logmodel <- lm(genus ~ emat$FireRotation)
y <- genus
plot(emat$FireRotation, genus, 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)", main = "Gallberry")
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
blr(log, emat$FireRotation, genus, Genus)

###################################################################################################
#4-5
#All Palmetto
#Plot with sites
genus <- ab2mat$dead.palmetto + ab2mat$live.palmetto
Genus <- "Palmetto"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#23
#Carpet Oak
#Plot with sites
genus <- ab2mat$oak.groundcover
Genus <- "Grouncover Oak"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#25
#Gallberry
#Plot with sites
genus <- ab2mat$gallberry
Genus <- "Gallberry"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#29
#Mound Oak
#Plot with sites
genus <- ab2mat$oak.shrub
Genus <- "Shrub Oak"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

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
blr(log, emat$FireRotation, genus, Genus)

#19
#Darrow's blueberry
#Plot with sites
genus <- ab2mat$Darrows.blueberry
Genus <- "Darrow's blueberry"

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

#35
#Highbush blueberry
#Plot with sites
genus <- ab2mat$highbush.blueberry
Genus <- "Highbush Blueberry"

#Boundary layer regression
blr(exp, emat$FireRotation, genus, Genus)

#43
#Sweetbay
#Plot with sites
genus <- ab2mat$sweetbay
Genus <- "Sweetbay"

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

###################################################################################################
#30
#Rough Fetterbush
#Plot with sites
genus <- ab2mat$rough.fetterbush
Genus <- "Rough Fetterbush"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#19
#Darrow's blueberry
#Plot with sites
genus <- ab2mat$Darrows.blueberry
Genus <- "Darrow's blueberry"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#35
#Highbush blueberry
#Plot with sites
genus <- ab2mat$highbush.blueberry
Genus <- "Highbush Blueberry"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#43
#Sweetbay
#Plot with sites
genus <- ab2mat$sweetbay
Genus <- "Sweetbay"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

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
blr(exp, emat$FireRotation, genus, Genus)

#38
#Titi
#Plot with sites
genus <- ab2mat$titi
Genus <- "Titi"

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

#41
#Oak Tree
#Plot with sites
genus <- ab2mat$oak.tree
Genus <- "Oak Tree"

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

#15
#Greenbriar
#Plot with sites
genus <- ab2mat$greenbriar
Genus <- "Greenbriar"

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

###################################################################################################
#18
#Huckleberry
#Plot with sites
genus <- ab2mat$huckleberry
Genus <- "Huckleberry"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#38
#Titi
#Plot with sites
genus <- ab2mat$titi
Genus <- "Titi"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#41
#Oak Tree
#Plot with sites
genus <- ab2mat$oak.tree
Genus <- "Oak Tree"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#15
#Greenbriar
#Plot with sites
genus <- ab2mat$greenbriar
Genus <- "Greenbriar"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

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
blr(log, emat$FireRotation, genus, Genus)

#37
#Sparkleberry
#Plot with sites
genus <- ab2mat$sparkleberry
Genus <- "Sparkleberry"

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

#32
#Smooth Fetterbush
#Plot with sites
genus <- ab2mat$smooth.fetterbush
Genus <- "Smooth Fetterbush"

#Boundary layer regression
blr(exp, emat$FireRotation, genus, Genus)

#33
#Large Gallberry
#Plot with sites
genus <- ab2mat$large.gallberry
Genus <- "Large Gallberry"

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

###################################################################################################
#14
#Yellow Jessamine
#Plot with sites
genus <- ab2mat$yellow.jessamine
Genus <- "Yellow Jessamine"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#37
#Sparkleberry
#Plot with sites
genus <- ab2mat$sparkleberry
Genus <- "Sparkleberry"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#32
#Smooth Fetterbush
#Plot with sites
genus <- ab2mat$smooth.fetterbush
Genus <- "Smooth Fetterbush"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#33
#Large Gallberry
#Plot with sites
genus <- ab2mat$large.gallberry
Genus <- "Large Gallberry"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

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
blr(log, emat$FireRotation, genus, Genus)

#36
#Yaupon
#Plot with sites
genus <- ab2mat$yaupon
Genus <- "Yaupon"

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

#42
#Sweetleaf
#Plot with sites
genus <- ab2mat$sweetleaf
Genus <- "Sweetleaf"

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

#21
#Wax Myrtle
#Plot with sites
genus <- ab2mat$wax.myrtle
Genus <- "Wax Myrtle"

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

###################################################################################################
#13
#Muscadine Grape
#Plot with sites
genus <- ab2mat$muscadine.grape
Genus <- "Muscadine Grape"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#36
#Yaupon
#Plot with sites
genus <- ab2mat$yaupon
Genus <- "Yaupon"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#42
#Sweetleaf
#Plot with sites
genus <- ab2mat$sweetleaf
Genus <- "Sweetleaf"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#21
#Wax Myrtle
#Plot with sites
genus <- ab2mat$wax.myrtle
Genus <- "Wax Myrtle"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

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
blr(exp, emat$FireRotation, genus, Genus)

#20
#St. Andrew's Cross
#Plot with sites
genus <- ab2mat$St..Andrews.cross
Genus <- "St. Andrew's Cross"

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

#16
#Gopher Apple
#Plot with sites
genus <- ab2mat$gopher.apple
Genus <- "Gopher Apple"

#Boundary layer regression
blr(exp, emat$FireRotation, genus, Genus)

#31
#Sweet Pepperbush
#Plot with sites
genus <- ab2mat$sweet.pepperbush
Genus <- "Sweet Pepperbush"

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

###################################################################################################
#17
#Wicky
#Plot with sites
genus <- ab2mat$wicky
Genus <- "Wicky"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#20
#St. Andrew's Cross
#Plot with sites
genus <- ab2mat$St..Andrews.cross
Genus <- "St. Andrew's Cross"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#16
#Gopher Apple
#Plot with sites
genus <- ab2mat$gopher.apple
Genus <- "Gopher Apple"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#31
#Sweet Pepperbush
#Plot with sites
genus <- ab2mat$sweet.pepperbush
Genus <- "Sweet Pepperbush"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

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
blr(log, emat$FireRotation, genus, Genus)

#24
#Chokeberry
#Plot with sites
genus <- ab2mat$chokeberry
Genus <- "Chokeberry"

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

#22
#Conradina
#Plot with sites
genus <- ab2mat$conradina
Genus <- "Conradina"

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

#40
#American Holly
#Plot with sites
genus <- ab2mat$American.holly
Genus <- "American Holly"

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

###################################################################################################
#39
#Fringe Tree
#Plot with sites
genus <- ab2mat$fringe.tree
Genus <- "Fringe Tree"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#24
#Chokeberry
#Plot with sites
genus <- ab2mat$chokeberry
Genus <- "Chokeberry"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#22
#Conradina
#Plot with sites
genus <- ab2mat$conradina
Genus <- "Conradina"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#40
#American Holly
#Plot with sites
genus <- ab2mat$American.holly
Genus <- "American Holly"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

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
blr(log, emat$FireRotation, genus, Genus)

#27
#St. Johnswort
#Plot with sites
genus <- ab2mat$St..Johnswort
Genus <- "St. Johnswort"

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

#28
#Narrowleaf Paw Paw
#Plot with sites
genus <- ab2mat$narrowleaf.pawpaw
Genus <- "Narrowleaf Paw Paw"

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

#34
#Winged Sumac
#Plot with sites
genus <- ab2mat$winged.sumac
Genus <- "Winged Sumac"

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

###################################################################################################
#26
#Dewberry
#Plot with sites
genus <- ab2mat$dewberry
Genus <- "Dewberry"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#27
#St. Johnswort
#Plot with sites
genus <- ab2mat$St..Johnswort
Genus <- "St. Johnswort"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#28
#Narrowleaf Paw Paw
#Plot with sites
genus <- ab2mat$narrowleaf.pawpaw
Genus <- "Narrowleaf Paw Paw"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#34
#Winged Sumac
#Plot with sites
genus <- ab2mat$winged.sumac
Genus <- "Winged Sumac"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

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
blr(log, emat$FireRotation, genus, Genus)

#7
#Dead Shrub
#Plot with sites
genus <- ab2mat$dead.shrub
Genus <- "Dead shrub"

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

#8
#Live Tree
#Plot with sites
genus <- ab2mat$live.tree
Genus <- "Live tree"

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

#9
#Dead Tree
#Plot with sites
genus <- ab2mat$dead.tree
Genus <- "Dead tree"

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

###################################################################################################
#6
#Live Shrub
#Plot with sites
genus <- ab2mat$live.shrub
Genus <- "Live shrub"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#7
#Dead Shrub
#Plot with sites
genus <- ab2mat$dead.shrub
Genus <- "Dead shrub"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#8
#Live Tree
#Plot with sites
genus <- ab2mat$live.tree
Genus <- "Live tree"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#9
#Dead Tree
#Plot with sites
genus <- ab2mat$dead.tree
Genus <- "Dead tree"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)


















###################################################################################################
#1
#live.forb
#Plot with sites
genus <- ab2mat$live.forb
Genus <- "Live forb"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

###################################################################################################
#2
#dead.forb
#Plot with sites
genus <- ab2mat$dead.forb
Genus <- "Dead forb"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

###################################################################################################
#3
#grass
#Plot with sites
genus <- ab2mat$grass
Genus <- "Grass (ex. bunchgrass)"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

###################################################################################################
#4
#Live Palmetto
#Plot with sites
genus <- ab2mat$live.palmetto
Genus <- "Live palmetto"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

###################################################################################################
#5
#Dead Palmetto
#Plot with sites
genus <- ab2mat$dead.palmetto
Genus <- "Dead palmetto"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

###################################################################################################
#5-7-9
#Dead Woody
#Plot with sites
genus <- ab2mat$dead.palmetto + ab2mat$dead.shrub + ab2mat$dead.tree
Genus <- "Dead Woody"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

###################################################################################################
#11
#Bluestem
#Plot with sites
genus <- ab2mat$bluestem.grass
Genus <- "Bluestem"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)



##################################################################################################
#23-29-41
#Oak
#Plot with sites
genus <- ab2mat$oak.groundcover + ab2mat$oak.shrub + ab2mat$oak.tree
Genus <- "Oak"
plot(emat$FireRotation, genus, type = "n", 
     xlab = "Fire Rotation (years)", ylab = "Biomass (Mg/ha)")
text(emat$FireRotation, genus, labels = row.names(ab2mat), cex = cf)

#Boundary layer regression
blr(log, emat$FireRotation, genus, Genus)

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
