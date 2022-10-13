###################################################################################################
############################-----start-----########################################################
#The purpose of this script is to identify and remove outlier plots from the biomass plot data 

#This is in response to outliers within the boundary layer regression and the realization
#that some can be traced to individual plots that were located outside or in the ecotone area
#of stand boundaries.

#The point here is not to identify any old outlier but to use anamolous grtass cover data
#to help locate plots that fall outside the site selection criteria for mesic pine flatwoods.

#This script was re-habed/modified on 13-October-2022.
#Unchanged script is first version on GitHub repo (push to remote on 13-Oct-2022 0930 EST.
#Original script name: sef_2014.12.19_PlotBiomassData_OutlierDetection
#Original script location: C:\Users\james\Box\01. james.cronan Workspace\Research\2009_01_SEF\r_scripts

###################################################################################################
###################################################################################################
#STEP #1: ADMINISTRATIVE TASKS

#Reset functions
rm(list=ls())
dev.off()

#Libraries
library(stats)
library(mvoutlier)#for outlier detection
library(fields)#for set.panel function

###################################################################################################
###################################################################################################
#STEP #2: Open plot-level vine-corrected dataset

#Set working directory
setwd(paste("C:/Users/james/Box/01. james.cronan Workspace/Research/UW_PHD/Dissertation/", 
"4_Chapter_4/Data/Understory_Vegetation_FlatFiles/stage_3_vine_class_corrected_data/outputs", sep = ""))


#Open plot biomass data with field measurements standardized to a 2-year rough. This means
# 2-yr post-fire (2011-2012; sampling episode 4) measurements for sites that were burned in 2009-2010 and pre-fire
#(2009-2010; sampling episode 1) measurements for sites that were not burned.
pbd1 <- read.table(
  "20221011_Biomass_OriginalVegClasses_2yrRough_Episodes_1_andf_4.csv", 
  header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE,
  stringsAsFactors = F)

###STOPPED HERE, NOT SURE I WANT TO DO THIS QUITE YET. THIS ISN'T REALLY NEEDED FOR THE MULTIVARIATE 
#ANALYSIS SINCE ANY EFFECT IS LIKELY TO BE DROWNED OUT. MORE IMPORTANT FOR BOUNDARY LAYER REGRESSION
#WHERE I AM LOOKING AT INDIVIDUAL SPECIES.

#Above note was made prior to 2022 analysis. Not sure when. Could have been for 2019 conference presentation
#or during earlier (2014) efforts on this analysis.











###################################################################################################
###################################################################################################
#STEP #3: Examine outliers within sites

#Examine outliers for entire dataset


#Examine outliers for specific cases you are interested in.

#Non bunchgrasses, A18 is an outlier, loading is high relative to long FRI
xx <- pbd1[pbd1$siteName == "A18",]

dev.off()
set.panel(1,2)
cf <- 0.7
title <- "Non-bunchgrasses, Site A18"
plant <- round((xx$grass/120)*100,1)
hist(plant, main = title, xlab = "Percent Cover (%)")
plot(xx$PlotNo, plant, type = "n", xlab = "Plot Number", ylab = "Percent Cover (%)")
text(xx$PlotNo, plant, labels = xx$PlotNo, cex = cf)
#There are no outliers

#what about wiregrasses in A18, is it possible these samples were mislabeled?
xx <- pbd[pbd$siteName == "A18",]

dev.off()
set.panel(1,2)
cf <- 0.7
title <- "Wiregrass, Site A18"
plant <- round((xx$wiregrass/120)*100,1)
hist(plant, main = title, xlab = "Percent Cover (%)")
plot(xx$PlotNo, plant, type = "n", xlab = "Plot Number", ylab = "Percent Cover (%)")
text(xx$PlotNo, plant, labels = xx$PlotNo, cex = cf)
#There are two minor outliers but they are low percentages and not due to plots in an open meadow.

#Wiregrass, E501B, E505, E503C_S4, and E508A are an outliers, loading is high relative to long FRI
#What does cover look like?
#E501B
xx <- pbd[pbd$siteName == "E501B",]

dev.off()
set.panel(1,2)
cf <- 0.7
title <- "Wiregrass, Site E501B"
plant <- round((xx$wiregrass/120)*100,1)
hist(plant, main = title, xlab = "Percent Cover (%)")
plot(xx$PlotNo, plant, type = "n", xlab = "Plot Number", ylab = "Percent Cover (%)")
text(xx$PlotNo, plant, labels = xx$PlotNo, cex = cf)
#There is a very wide range of cover values (~5-85%) but none areoutliers

#E505
xx <- pbd[pbd$siteName == "E505",]

dev.off()
set.panel(1,2)
cf <- 0.7
title <- "Wiregrass, Site E505"
plant <- round((xx$wiregrass/120)*100,1)
hist(plant, main = title, xlab = "Percent Cover (%)")
plot(xx$PlotNo, plant, type = "n", xlab = "Plot Number", ylab = "Percent Cover (%)")
text(xx$PlotNo, plant, labels = xx$PlotNo, cex = cf)
#684, and to some extent, 685 are outliers. Both these sites are in the transition area between the pine
#flatwoods and an adjacent meadow. Looking at the plot photos it is clear that cover lines
#in 684 (wiregrass cover is near 100% while most plots are below 20%) are within the meadow and
#this plot should be dropped. It is less clear that 685 should be dropped but lines do
#appear to extend into meadow in phots and wiregrass cover is ~55%. I would drop them both.


#E503_S4
xx <- pbd[pbd$siteName == "E503-S4",]

dev.off()
set.panel(1,2)
cf <- 0.7
title <- "Wiregrass, Site E503_S4"
plant <- round((xx$wiregrass/120)*100,1)
hist(plant, main = title, xlab = "Percent Cover (%)")
plot(xx$PlotNo, plant, type = "n", xlab = "Plot Number", ylab = "Percent Cover (%)")
text(xx$PlotNo, plant, labels = xx$PlotNo, cex = cf)
#One site (874) is probably an outlier but I don't recall that this plot is outside of the
#flatwoods unit. Do not remove it.


#E508A
xx <- pbd[pbd$siteName == "E508A",]

dev.off()
set.panel(1,2)
cf <- 0.7
title <- "Wiregrass, Site E508A"
plant <- round((xx$wiregrass/120)*100,1)
hist(plant, main = title, xlab = "Percent Cover (%)")
plot(xx$PlotNo, plant, type = "n", xlab = "Plot Number", ylab = "Percent Cover (%)")
text(xx$PlotNo, plant, labels = xx$PlotNo, cex = cf)
#Five plots could be considered outliers: 910, 911, 913, 914 (60-80%), and especially 915 (100%).
#I checked the map. 915 is right on the edge of the meadow and I remember both cover lines
#extending into the meadow. I think the rest of the sites are fully in the forest, its
#just wetter than the rest of the site so grass cover is higher in the understory.
#Remove 915

#Bluestem, S209 and S411 are there outliers, loading is high relative to other plots
#S209
xx <- pbd[pbd$siteName == "S209",]

dev.off()
set.panel(1,2)
cf <- 0.7
title <- "Bluestem, Site S209"
plant <- round((xx$bluestem/120)*100,1)
hist(plant, main = title, xlab = "Percent Cover (%)")
plot(xx$PlotNo, plant, type = "n", xlab = "Plot Number", ylab = "Percent Cover (%)")
text(xx$PlotNo, plant, labels = xx$PlotNo, cex = cf)
#Three plots are outliers: 346, 355, and 356. These all likely have cover lines in the
#wetland and are not characteristic of pine flatwoods.
#I checked the map and these plots definitely fall along the wetland that bisects the stand.
#Remove all three

#S411
xx <- pbd[pbd$siteName == "S411",]

dev.off()
set.panel(1,2)
cf <- 0.7
title <- "Bluestem, Site S411"
plant <- round((xx$bluestem/120)*100,1)
hist(plant, main = title, xlab = "Percent Cover (%)")
plot(xx$PlotNo, plant, type = "n", xlab = "Plot Number", ylab = "Percent Cover (%)")
text(xx$PlotNo, plant, labels = xx$PlotNo, cex = cf)
#There are no outliers.

#S330
xx <- pbd[pbd$siteName == "S330",]

dev.off()
set.panel(1,2)
cf <- 0.7
title <- "Grass, Site S330"
plant <- round((xx$grass/120)*100,1)
hist(plant, main = title, xlab = "Percent Cover (%)")
plot(xx$PlotNo, plant, type = "n", xlab = "Plot Number", ylab = "Percent Cover (%)")
text(xx$PlotNo, plant, labels = xx$PlotNo, cex = cf)
#Plot 107 is an outlier and this plot is, as with plots in S209, in a wetland and is not
#characteristic of pine flatwoods.
#There is no map of this site but photos show this plot has heavy grass cover because it
#is situated in a small wet depression.

#For comparison, these are the clip plots I removed.
#pbd1 <- pbd[!(pbd$PlotNo %in% c(683, 684, 909, 913, 914, 345, 356, 105, 106)),]
pbd1 <- pbd[!(pbd$PlotNo %in% c(684, 685, 915, 345, 355, 356, 107)),]

###################################################################################################
###################################################################################################
#STEP #4: Summarize data by site.

sbd1 <- summarize(X = pbd1[,3:58], by = pbd1$siteName, colMeans, stat.name = "live.forb")
colnames(sbd1)[1] <- "siteName"


###################################################################################################
###################################################################################################
#STEP #5: Save relevant files

#set a date used to name the file.
dt <- Sys.Date()
tm <- format(Sys.time(), format = "%H.%M.%S", 
             tz = "", usetz = FALSE)


#Save the plot-level data in matrix (rows: sites, cols: species) form (Episode 1).
#Add row names as a column
write.table(pbd1, file = paste("c:\\usfs_sef_data_output\\sef_Ecology_CoverPlotMatrix_Ep1_OriginalOulierX_",
                               dt,"_",tm,".csv",sep = ""), append = F, quote = T, sep = ",", eol = "\n", na = "NA", 
            dec = ".", row.names = F, col.names = T, qmethod = c("escape", "double"))#

#Save the site-level data in matrix (rows: sites, cols: species) form (Episode 1).
#Add row names as a column
write.table(sbd1, file = paste("c:\\usfs_sef_data_output\\sef_Ecology_CoverSiteMatrix_Ep1_OriginalOulierX_",
                              dt,"_",tm,".csv",sep = ""), append = F, quote = T, sep = ",", eol = "\n", na = "NA", 
            dec = ".", row.names = F, col.names = T, qmethod = c("escape", "double"))#


