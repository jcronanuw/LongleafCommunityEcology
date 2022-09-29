

###################################################################################################
############################-----start-----########################################################
#Data Processing script for Chapter 4 (Longleaf community ecology) of dissertation.
#Purpose of this script is to access data in Access database archive for the season-of-burn
#fuels study (Cronan et al. 2015 - Forest Ecology & Management/aka Chapter 3 of dissertation)
#and generate flat files (csv files) as an input for multivariate/boundary layer regression analysis
#of understory biomass data.

#This script is an edited version of:
#C:\Users\james\Box\01. james.cronan Workspace\Research\2009_01_SEF\r_scripts\sef_2012.08.24_BiomassDataFINAL.r

#PREVIOUS NOTES
#Standing biomass data analysis.

#NOTES BELOW DESCRIBE DATA USED IN THIS SCRIPT THAT HAS BEEN CHANGED IN ACCESS DB.

#NOTE: Plot numbers for E103B_S3 have been changed to a unique set of numbers. They
#are not the actual plot numbers for this site.

#NOTE: Site A319, plot 324C, episode 2(post) was covered by water and not sampled. 
#This data has been imputed from site(post) averages and is now in Access DB.

#NOTE: Site E100B-East, plot 141D, episode 1(pre) is missing data. Field sheet 
#exists but there are no weights. Also no weights for 142D and 143D but they were 
#entered into Access DB Hmm. A copy with oven weights must be around somewhere. 
#This data has been imputed from site(pre) averages and is now in Access DB.

#SYSTEM TIME FUNCTIONS LOCATED AT LINES:
#The second two system times are disabled so they do not overlap with the gates. The first system
#has two endings: one just prior to the second (line 1792; disabled) and on after the gates that
#times the entire script.
#38     #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>SYSTEM TIME>>>>>>>>>>>>>>
#1804   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>SYSTEM TIME>>>>>>>>>>>>>>
#1918   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>SYSTEM TIME>>>>>>>>>>>>>>

###################################################################################################
###################################################################################################
#STEP #1: ADMINISTRATIVE TASKS

#Reset functions
rm(list=ls())
dev.off()

#Libraries
library(Hmisc)
library(e1071)#used for skewness() and kurtosis() functions
library(RODBC)#used for direct import from Access

system.time({#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>BEGIN SYSTEM TIME>>>>>>>>
#Administrative parameters:

#List sampling episodes you want to look at:
sampleEpisode <- c(1,2,3,4)
#1 = prefire
#2 = postfire
#3 = 1 year postfire
#4 = 2 postfire

#To remove any of these sites type their site ID into the SiteRemoval object
SiteRemoval <- c(28,29,30,31)
#Why are you removing these sites:
#Primary reasons is to remove sites with potential issues.
#Secondary reason is to create an equal sample size for each block.
#6(A18). This site was burned on the same day as A213 (pseudoreplication)
#15(E807B). This site was burned before the second year sampling date.
#26(S411). There is no temperature data for this site.

#2022
#6 - keep this site. Since you are analysing the affect of a 20-yr fire history and not the current 
#fire it does not matter that these two sites were burned on the same day.
#15 - drop this site. Unless one of the sampling episodes was collected 2 years after the last burn
#26 - keep this site. You will not use temperature data in this analysis.



#Create a formal species list without excess listings.
#Inputs here are species material codes.
SpeciesRemoval <- c(1,6,7,15,16,19,20,21,23,24,25,35,36,37,55,66,67,68,70,73,75,79,80,
85,86,89,90)

#Set an upper ceiling for subsample gross weights. This is used in the error checking section
#to signal a possible mix up of samples weights. For example, if someone accidently enetered
#one of the other categories (measured in grams) into the gross category this ceiling should
#help identify those errors since subsample weights are typically in the hundreds of grams.
GrsWt_Ceiling <- 30#Hard to imagine a subsample weighing more than 30 kilograms.
#This object is only used in step 2i

#Number of plots per site. This code will not handle variable plot number per site. If plots
#are missing data then the data must be imputed. Code will correct degrees of freedom for 
#statistical analysis by subtracting one degree of freedom for site and species material with 
#imputed data
#This object is not invoked until step 6.
n <- 20

#Alpha level
a.level <- 0.05

#List sampling areas in units that correspond with database
a.sizes <- c(1,4)

#Customized variance function
#x is the data vector
#y is the degrees of freedom which may vary when data is imputed.
ivar <- function(x,y)
{
(sum((x - mean(x))^2))/(y)
}

#Confidence Interval
CI.xbar <- function(xbar,alpha,tailed,s,n)
{
ll <- (xbar-(qt((1-(alpha/tailed)),n-1)*(s/sqrt(n))))
ul <- (xbar+(qt((1-(alpha/tailed)),n-1)*(s/sqrt(n))))
cbind(ll,ul)
}

###################################################################################################
###################################################################################################
#STEP #2: OPEN, ADJUST, AND REVIEW DATA FILES

###################################################################################################
#STEP #2a: Open data files





#Import Access database into working environment.
dbStanding <- file.path("C:/Users/james/Box/01. james.cronan Workspace/Research/2009_01_SEF/SEF_DataArchive/20111122_sef_db_archive/jim_floridaDB_BackEnd.mdb")
channel <- odbcConnectAccess(dbStanding)

#Extract standing biomass data:
db.1 <- sqlFetch(channel,"tbl_clipPlot")
sdc.clipPlot <- db.1[,1:13]

#Adjust factors to character. Old code (prior to RODBC function imported clip plot
#alpha numeric data as characters, not factors as RODBC does. Ensuing code is better suited to
#character analysis.
i <- sapply(sdc.clipPlot, is.factor)
sdc.clipPlot[i] <- lapply(sdc.clipPlot[i], as.character)

#Extract plot metadata
db.2 <- sqlFetch(channel,"tbl_plot")
sdc.plot <- db.2[,1:6]

#Extract sampling episode ID metadata.
db.3 <- sqlFetch(channel,"tbl_plotSamplingEpisode")
sdc.plotSamplingEpisode <- db.3[,1:3]

#Extract site metadata:
db.4 <- sqlFetch(channel,"tbl_site")
sdc.site <- db.4[,c(1,2,3,4,12)]

#Extract species look-up table:
db.5 <- sqlFetch(channel,"lut_speciesMaterial")
sdc.speciesMaterial <- db.5[,c(1,2,6)]

#Adjust factors to character. Old code (prior to RODBC function imported clip plot
#alpha numeric data as characters, not factors as RODBC does. Ensuing code is better suited to
#character analysis.
i <- sapply(sdc.speciesMaterial, is.factor)
sdc.speciesMaterial[i] <- lapply(sdc.speciesMaterial[i], as.character)

#Extract clip plot sampling area data.
db.6 <- sqlFetch(channel,"tbl_siteClipPlotSizes")
sdc.siteClipPlotSizes <- db.6

#Adjust factors to character. Old code (prior to RODBC function imported clip plot
#alpha numeric data as characters, not factors as RODBC does). Ensuing code is better suited to
#character analysis.
i <- sapply(sdc.siteClipPlotSizes, is.factor)
sdc.siteClipPlotSizes[i] <- lapply(sdc.siteClipPlotSizes[i], as.character)

#Clean up objects
odbcClose(channel)
rm(db.1, db.2, db.3, db.4, db.5, db.6)

###################################################################################################
#STEP #2b: Remove unwanted sites and sampling episodes, then prereview site level data.
#NOTE:
#1: Site IDs to remove are input in SiteRemoval object in step 1.
#2: Sampling Episodes to retain are input in sampleEpisode object in step 1.

#Set up a log to record data properties including errors.
EL <- list()

#Make sure the start and end plot numbers (TagNo) are correct. These will be used to 
#test for errors in all other files with plot number information.

#Remove unwanted sites
a.site <- sdc.site[sdc.site[,1] %in% SiteRemoval == F,]

#Remove plots associated with unwanted sites.
a.plot <- sdc.plot[sdc.plot[,2] %in% SiteRemoval == F,]

#Remove sampling episode IDs associated with unwanted sites.
a.episode <- sdc.plotSamplingEpisode[sdc.plotSamplingEpisode[,
2] %in% unique(a.plot[,1]) == T & sdc.plotSamplingEpisode[,
3] %in% sampleEpisode == T,]

#Remove samples associated with unwanted sites
a.clip <- sdc.clipPlot[sdc.clipPlot[,2] %in% unique(a.episode[,1]) == T,]

#Remove clip plot sampling areas associated with unwanted sites
a.area <- sdc.siteClipPlotSizes[sdc.siteClipPlotSizes[,2] %in% SiteRemoval == F,]

#Review site data
EL[[1]] <- a.site[order(a.site$siteName),]
names(EL)[[1]] <- paste("1. Review site data: There are", length(a.site[,1]), 
"sites being imported into this analysis")

###################################################################################################
#STEP #2c: Update the species menu to include specific entries that eliminates
#ambiguous species materials.

#First you need to create a formal species list
SpeciesMenu <- sdc.speciesMaterial[sdc.speciesMaterial[,1] %in% SpeciesRemoval == F,]

#Second, add duplicates to species materials with ambiguous status (L/D)
SpeciesMenuA <- SpeciesMenu[is.na(SpeciesMenu[,3]) == T,]

#Third, add species with multiple plant parts
pps <- unique(a.clip[,3][is.na(a.clip[,5]) == F])
npp <- mapply(function(y){length(unique(a.clip[,5][is.na(a.clip[,5]) == F & a.clip[,3] == y]))-1},
pps)
tpp <- rep(pps,npp)

#Fourth, ammend ambiguous species list.
for(i in 1:length(tpp))
SpeciesMenuA <- rbind(SpeciesMenuA,SpeciesMenu[SpeciesMenu[,1] == tpp[i],])

#Fifth, add new species ID
SpeciesMenuA[,1] <- seq(max(sdc.speciesMaterial[,1])+1,
max(sdc.speciesMaterial[,1])+length(SpeciesMenuA[,1]),1)
#Use sdc.speciesMaterial in case an old species code was the max. Ensures no duplication.

#Sixth, add new common names
SpeciesMenuA[,2] <- c("dead palmetto", "dead shrub", "dead forb", "dead tree", "live palm rachis")

#Seventh, add a new default L/D status
SpeciesMenuA[,3] <- c(
rep("D", length(SpeciesMenuA[,1][grep("dead", SpeciesMenuA[,2])])),
rep("L", length(SpeciesMenuA[,1][grep("live", SpeciesMenuA[,2])])))

#Eighth, add these new species materials to the subset.
SpeciesMenu <- rbind(SpeciesMenu,SpeciesMenuA)

#Ninth, now turn the remaining ambiguous species materials into live plants.
SpeciesMenu[,2][SpeciesMenu[,1] == 18] <- "live palm frond"
SpeciesMenu[,2][SpeciesMenu[,1] == 29] <- "live shrub"
SpeciesMenu[,2][SpeciesMenu[,1] == 32] <- "live forb"
SpeciesMenu[,2][SpeciesMenu[,1] == 65] <- "live tree"
SpeciesMenu[,3][is.na(SpeciesMenu[,3]) == T] <- "L"

###################################################################################################
#STEP #2d: Apply the changes made to the species material list to the clip plot
#dataset (a.clip)

#First, assign new species codes to species materials with ambiguous status.

#Convert "palmetto" into "dead palmetto". Need to change species code.
a.clip[,3][a.clip[,3] == 18 & a.clip[,4] == "D"] <- SpeciesMenuA[1,1]

#Convert "shrub" into "dead shrub". Need to change species code.
a.clip[,3][a.clip[,3] == 29 & a.clip[,4] == "D"] <- SpeciesMenuA[2,1]

#Convert "forb" into "dead forb". Need to change species code.
a.clip[,3][a.clip[,3] == 32 & a.clip[,4] == "D"] <- SpeciesMenuA[3,1]

#Convert "tree" into "dead tree". Need to change species code.
a.clip[,3][a.clip[,3] == 65 & a.clip[,4] == "D"] <- SpeciesMenuA[4,1]

#Convert "palmetto" into "live palm rachis". Need to change species code.
a.clip[,3][a.clip[,3] == 18 & a.clip[,5] == "rachis"] <- SpeciesMenuA[5,1]

#Second, reassign species names in sc.sn to reflect changes made in step 4a.
sc.sn <- as.character(mapply(function(y) SpeciesMenu[,2][SpeciesMenu[,1]==y],a.clip[,3]))

###################################################################################################
#STEP #2e: Apply the changes made to the species material list to the clip plot
#sampling area dataset (a.area).
#NOTE: revised a.area will not be needed until step 5.

#Convert "palmetto" into "dead palmetto". Need to change species code.
a.area[,3][a.area[,3] == 18 & a.area[,4] == "D"] <- SpeciesMenuA[1,1]

#Convert "shrub" into "dead shrub". Need to change species code.
a.area[,3][a.area[,3] == 29 & a.area[,4] == "D"] <- SpeciesMenuA[2,1]

#Convert "forb" into "dead forb". Need to change species code.
a.area[,3][a.area[,3] == 32 & a.area[,4] == "D"] <- SpeciesMenuA[3,1]

#Convert "tree" into "dead tree". Need to change species code.
a.area[,3][a.area[,3] == 65 & a.area[,4] == "D"] <- SpeciesMenuA[4,1]

#Add "live palm rachis" by duplicating all palmetto. Need to change species code.
lpr.area <- a.area[a.area[,3] == 18,]
lpr.area[,3] <- SpeciesMenuA[5,1]
a.area <- rbind(a.area,lpr.area)

###################################################################################################
#Step #2f: Add additional sample info to the standing biomass dataset 
#(a.clip).
#These objects make the data frame more intuitive by listing plot number (pID), 
#site ID (sID), site name (sNo), sampling episode (seID; 1 = pre, 2 = post, and 
#3 = year 1, 4 = year 2), and common names for the species (sc.sn)

#Lines up plot IDs with other data in clip plot data table.
pID <- mapply(function(y) a.episode[,2][a.episode[,1]==y], a.clip[,2])

#Lines up plot numbers with other data in clip plot data table.
pNo <- mapply(function(y) a.plot[,3][a.plot[,1]==y],pID)

#Lines up site IDs with other data in clip plot data table.
sID <- mapply(function(y) a.plot[,2][a.plot[,1]==y],pID)

#Lines up site numbers with other data in clip plot data table.
sNo <- mapply(function(y) a.site[,2][a.site[,1]==y],sID)

#Lines up sampling episode with other data in clip plot data table.
seID <- mapply(function(y) a.episode[,3][a.episode[,1]==y], a.clip[,2])

#Lines up species common names with other data in clip plot data table.
sc.sn <- as.character(mapply(function(y) SpeciesMenu[,2]
[SpeciesMenu[,1]==y], a.clip[,3]))

###################################################################################################
###################################################################################################
#STEP #3: ERROR CHECKING

###################################################################################################
#STEP #3a: Check structure of the plot-level metadata against the site-level 
#metadata (i.e. a.site vs. a.plot).

#Create an object to measure expected plot metadata
test.eppm <- sort(as.vector(mapply(function(x) seq(a.site[x,3],a.site[x,4],1),
1:length(a.site[,1]))))
test.espm <- mapply(function(y)
a.site[,2][a.site[,3] <= test.eppm[y] & a.site[,4] >= test.eppm[y]], 1:length(test.eppm))

#Show all plots in plot metadata.
test.appm <- sort(a.plot[,3])

#Bin actual plots against expected plots
#This function takes each of the expected plot numbers and measures how many times
#that plot occurs in the actual plot metadata. It would be expected that there 
#would be one of each expected plot. A length of zero means that the plot
#was not entered into the database and a length greater than one means the plot
#was entered more than once.
#The test.na vector will be used to identify 
#actual plots not represented in the expected plots (i.e. incorrectly enetered 
#plot numbers).
test.lapm <- mapply(function(y) length(test.appm[test.appm == y]),test.eppm)
test.lapm[length(test.lapm)] <- ifelse(length(test.appm[test.appm %in% test.eppm == F])==0,1,0)

#Record actual plots that do not match up with any expected plots
test.napm <- vector() 
test.napm <- test.appm[is.na(match(test.appm,test.eppm)) == T]

#Create a table that describes details of errors.
EL.2 <- data.frame(Site = test.espm, Plot = test.eppm, Occurences = test.lapm)
EL[[2]] <- if(length(EL.2[EL.2$Occurences != 1,1])==0)
{
EL[[2]] <- "No errors: actual plot numbers match expected plot numbers"
} else
{
EL[[2]] <- EL.2[EL.2$Occurences != 1,]
}

names(EL)[[2]] <- paste("2. Test plot metadata: table will show instances were expected plots", 
"are missing/duplicated.")

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 1 >>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if(is.list(EL[[2]]) == F)
{

#Add vector of unexpected plot numbers
EL[[3]] <- ifelse(length(test.napm)==0,"No errors: no unexpected plot numbers were detected",
test.napm)
names(EL)[[3]] <- "3. Test plot metadata: listed plot numbers are unexpected"

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 2 >>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if(is.list(EL[[3]]) == F)
{

###################################################################################################
#STEP #3b: Check structure of the sampling episode-level metadata against the site-level 
#metadata (i.e. a.site vs. a.episode).

#Create objects to measure expected sampling episode metadata
test.epem <- rep(sort(as.vector(mapply(function(x) seq(a.site[x,3],a.site[x,4],1),
1:length(a.site[,1])))),length(sampleEpisode))
test.esem <- mapply(function(y)
a.site[,2][a.site[,3] <= test.eppm[y] & a.site[,4] >= test.eppm[y]], 1:length(test.eppm))
test.eeem <- as.vector(mapply(function(y) rep(y,length(test.eppm)),sampleEpisode))

#Bin actual sampling episode IDs against those expected for each plot and sampling episode.
test.laem <- mapply(function(x,y) 
{
length(a.episode[,1][a.episode[,2] == a.plot[,1][a.plot[,3] == x] & a.episode[,3] == y])
},
test.epem,test.eeem)

#Record actual sampling episode IDs that do not match up with those to be expected.
#This object shows the expected sampling episode IDs.
test.eiem <- vector() 
test.eiem <- mapply(function(x,y)
{
a.episode[,1][a.episode[,2] == a.plot[,1][a.plot[,3] == x] & a.episode[,3] == y]
},
test.epem,test.eeem)
#This object compares the expected sampling episode IDs calculated above with the actual
#sampling episode IDs and returns and actual IDs that do not match any of the expected IDs.
test.naem <- vector() 
test.naem <- a.episode[,1][is.na(match(a.episode[,1],test.eiem)) == T]

#Create a table that describes details of errors.
EL.4 <- data.frame(Site = test.esem, Plot = test.epem, Occurences = test.laem)
EL[[4]] <- if(length(EL.4[EL.4$Occurences != 1,1])==0)
{
EL[[4]] <- "No errors: actual sampling episode IDs match expected sampling episode IDs"
} else
{
EL[[4]] <- EL.4[EL.4$Occurences != 1,]
}

names(EL)[[4]] <- paste("4. Test sample episode metadata: Table will show instances",
"were expected sampling episode IDs are missing/duplicated.")

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 3 >>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if(is.list(EL[[4]]) == F)
{


#Add vector of unexpected plot numbers
EL[[5]] <- ifelse(length(test.naem)==0, 
"No errors: no unexpected sampling episode IDs were detected", test.naem)
names(EL)[[5]] <- "5. Test sample episode metadata: listed sample episode IDs are unexpected"

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 4 >>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if(is.list(EL[[5]]) == F)
{

###################################################################################################
#STEP #3c: Check structure of the sample-level data against the site-level 
#metadata (i.e. a.site vs. a.clip).

#Set up a data.frame of expected plots in the sample data.base (a.clip)
eecp <- as.vector(mapply(function(y) rep(y, length(test.eppm)),sampleEpisode))
evcp <- data.frame(Site = rep(test.espm,length(sampleEpisode)), 
Plot = rep(test.eppm,length(sampleEpisode)), Episode = eecp, 
EpisodeId = mapply(function(x,y)
a.episode[,1][a.episode[,2] == a.plot[,1][a.plot[,3] == x] & a.episode[,3] == y], 
rep(test.eppm,length(sampleEpisode)), eecp),
Occurences = rep(0,length(eecp)))

#Measure the number of samples in each plot and sampling episode.
evcp$Occurences <- mapply(function(x,y)
{
length(a.clip[,2][a.clip[,2] == a.episode[,1][a.episode[,2] == 
a.plot[,1][a.plot[,3] == x] & a.episode[,3] == y]])
},
evcp[,2], evcp[,3])

#Measure unexpected sampling episode IDs in the sample dataset
naEpisodes <- vector()
expectedEpisodeIds <- sort(unique(evcp[,4]))
actualEpisodeIds <- sort(unique(a.episode[,1]))
test.aecp <- a.episode[,1][is.na(match(a.episode[,1],test.eiem)) == T]
#The result for this object should be zero. Any sampling episode IDs that
#do not fit into the expected data matrix indicate an error.

#Isolate plots/sampling episodes where the number of samples is zero.
#NOTE: this is not necessarily a bad thing. Some clip plots may not have
#any standing biomass.
if(length(evcp[evcp$Occurences == 0,1])==0)
{
EL[[6]] <- "There were no plots without standing biomass for all sampling episodes"
} else 
{
EL[[6]] <- evcp[evcp$Occurences == 0,]
}
names(EL)[[6]] <- paste("6. For each sampling episode, review plots were no standing", 
"biomass was sampled.")

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE >>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#There is no need for a gate at error #6, possible that a dataset with no errors may register here.

#Add vector of unexpected sampling episode IDs
EL[[7]] <- ifelse(length(naEpisodes)==0,
"No errors: actual sampling episode IDs match expected sampling episode IDs", 
naEpisodes)
names(EL)[[7]] <- paste("7. Test standing clip plot sample data:", 
"listed sample episode IDs are unexpected")

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 5 >>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if(is.list(EL[[7]]) == F)
{
###################################################################################################
#Step #3d
#Check for errors on the sampling area menu (a.area)

#In looking at this dataset, one of the most obvious problems is a mismatch between
#species listed in the clip plot sizes data (sdc.siteClipPlotSizes) and the clip 
#plot sample data (sdc.clipPlot). This subsection checks for the following errors.
#1. Are all species materials in the clip plot dataset (a.clip) represented in 
#the sampling area menu (a.area).
#2. For species represented has a sampling area been entered (i.e. value has not
#been left blank)
#3. And is area acceptable (i.e. 1m2 or 4m2).
#4. Check for duplicate entries in the sampling area menu. This can cause problems
#if there are two entries with different sampling areas.
#5. Incorrect status. This is not checked directly but since species material are 
#segregated in this subsection by status and check against the sample data the
#species material in the sampling material menu will show up as the following error:
#"Species missing from sampling area menu".

#What this subsection does NOT check for:
#1. Errors for entries in the sampling area menu that were not sampled. For example. If there
#is an entry for grape in the sampling area menu, the status is listed as dead and the
#area is listed as 10m2 but grape was not sampled at that site no error will be returned.
#This should be fine since that information will never be used in the loading calculations.
#2. Errors in the clip plot dataset. For example, duplicate sample entries in the clip plot
#dataset (a.clip) not due to multiple bags or samples with multiple plant parts (i.e. palm). 
#These errors are scanned in further substeps.

#The goal here is to create a list of unique species for each site for the clip plot
#dataset (a.clip) and then compare it against the sampling area list in the 
#corresponding clip plot size data table (a.area).

#Any mismatches should be identified and fixed in the Access database. 
#This creates a list, each element contains a vector listing the ordered and unique
#speciesMaterialID numbers by site from the clip plot table (sl).
#This creates a list of all species entries, each level representing a site.
sl <- list()
for(i in 1:length(a.site$siteID))
sl[[i]] <- 
as.numeric(a.clip$speciesMaterialID[a.clip$plotSamplingEpisodeID %in% 
a.episode$plotSamplingEpisodeID[a.episode$plotID %in% 
a.plot$plotID[a.plot$siteID == a.site$siteID[i]]]])

#This creates a list that shows unique species codes by site.
#Status or plant part considerations have been removed because each species material
#has a unique code added in step 2.
slA <- list()
for(i in 1:length(a.site$siteID))
{
slA[[i]] <- sort(unique(sl[[i]]))
}

#Now do the same thing you did above, but for the siteClipPlotSizes table

#This creates a list, each element contains a vector listing the ordered and unique
#speciesMaterialID numbers by site from the clip plot sizes table (szl).
szA <- list()
for(i in 1:length(a.site$siteID))
szA[[i]] <- 
a.area$speciesMaterialID[a.area$siteID == a.site$siteID[i]]

#Set up lists to describe sampling area menu test.
setA <- list()

for(i in 1:length(a.site$siteID))
{
Site <- as.character(a.site$siteName[i])

#Create a common name species list that correlates with speciesMaterialIDs for all materials.
scnA <- mapply(function(y)
{
ifelse(y %in% SpeciesMenu$speciesMaterialID,
SpeciesMenu$commonName[SpeciesMenu$speciesMaterialID %in% y],
paste("Excluded Species:",
sdc.speciesMaterial$commonName[sdc.speciesMaterial$speciesMaterialID %in% y], sep = " "))
}, slA[[i]])

#Test for species material presence in sampling area menu.
#T = species present in clip plot dataset is present in sampling area menu.
#F = species present in clip plot dataset is not present in sampling area menu. 
mlA <- slA[[i]] %in% szA[[i]]#numbers generated are not species material IDs

#Subset the sampling area table for each site.
ssa1A <- a.area[a.area$siteID == a.site$siteID[i],]

#This object shows the sampling area entered in the sampling area menu for each species 
#material sampled for a given site. If a species is present in the clip plot dataset
#but not the species material sampling area menu a NA is produced.
ssa2A <- mapply(function(y)
{
ifelse(length(ssa1A$clipPlotSize[ssa1A$speciesMaterialID == y]) > 1, 9999, 
ssa1A$clipPlotSize[ssa1A$speciesMaterialID == y])
},slA[[i]])

#Create a vector that Explains any errors detected. There are three and they are
#mutually exclusive.
#1. No corresponding species in species materal sampling area menu (i.e. a "No" in
#mlL or mlD).
#2. Sampling area in species material sampling area menu is blank.
#3. Sampling area in species material sampling area menu is incorrect (i.e. not
#1 or 4).
#Create a list of errors:
error <- c("None", "Species missing from sampling area menu", 
"Sampling area not enetered", "Incorrect sampling area entered", 
"Duplicate species material")

#Create a vector that flags and describes errors.
#... for standing live species materials.
#This dataframe shows a row without errors.
exA <- mapply(function(y)
{
#Select error message.
ifelse(mlA[y] == F, error[2], ifelse(is.na(ssa2A[y]) == T, error[3], 
ifelse(ssa2A[y] == 9999, error[5], ifelse(ssa2A[y] %in% a.sizes, 
error[1],error[4])))) 
}, 1:length(mlA))

#Create a table for standing species materials that describes how clip plot sample-level
#data sorresponds to species material data in the sample area menu.
salA <- data.frame(Species_Material_ID = slA[[i]],
Present_In_Menu = mlA, Sampling_Area_m2 = ssa2A, commonName = I(scnA), Error = I(exA))

#Place data generated in this loop into a list
setA[[i]] <- salA

names(setA)[[i]] <- as.character(a.site$siteName[i])
}

#These are list files, they will display how the clip plot samples
#match up with the sampling areas. NAs indicate where a clip plot sample
#(a.clip) for a site does not have a correpsonding sample area in sdc.siteClipPlotSizes.
#This is the file for standing dead materials: setD
#This is the file fpr standing live materials: setL

#Scan the sampling error file (setA) for errors and report findings in the log.
t1 <- vector()
r1 <- list()
for(i in 1:length(setA))
{
t1[i] <- all(setA[[i]][,5]=="None")==T
if(t1[i] == F)
{
r1[[length(r1)+1]] <- setA[[i]]
names(r1)[[length(r1)]] <- names(setA[i])
}
}

#Register errors in the log.
if(length(r1)==0)
{
EL[[8]] <- paste("No errors: sampling area data corresponds with clip plot dataset and contains",
"acceptable values")
} else 
{
EL[[8]] <- r1
}
names(EL)[[8]] <- "8. Test for errors in sampling area menu."

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 6 >>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if(is.list(EL[[8]]) == F)
{

###################################################################################################
#Step #3e: Find instances within the clip plot dataset where numeric species codes 
#are unattached to the master species list (SpeciesMenu).
#NOTE: the "Species not listed" error should never occur because the Access DB prevents users 
#from entering species not present in the original master species list (sdc.speciesMaterial)
#However, this subsection will detect instances when undesired species are in the clip plot
#database because it uses the revised master species list (SpeciesMenu) as a filte
ov <- sort(unique(a.clip$speciesMaterialID)) %in% 
sort(unique(SpeciesMenu$speciesMaterialID))
av <- sort(unique(a.clip$speciesMaterialID))
zv <- av[ov == F]

#Register errors in the log.
if(length(zv)==0)
{
EL[[9]] <- paste("No errors: species material IDs in the clip plot dataset are represented",
"in the species material master list (sdc.speciesMaterial)")
} else 
{
EL[[9]] <- as.data.frame(cbind(Species_Material_ID = zv, 
Common_Name = ifelse(length(SpeciesMenu[,2][SpeciesMenu[,1] %in% zv]) > 0, 
SpeciesMenu[,2][SpeciesMenu[,1] %in% zv],"Species not listed")))
}
names(EL)[[9]] <- "9. Test for errors in species material master list."

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 7 >>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if(is.list(EL[[9]]) == F)
{

###################################################################################################
#Step 3f: check to make sure species status is correct in the clip plot dataset 
#(a.clip) and the sampling area menu (a.area).

#Subset species material codes for live anddead.
liveSpecies <- SpeciesMenu[,1][SpeciesMenu[,3] == "L"]
deadSpecies <- SpeciesMenu[,1][SpeciesMenu[,3] == "D"]

#Check sampling area menu:

#Isolate all live species codes from the sampling area menu.
#This isolates the table with rows where status = L
areaL <- a.area[a.area$liveORdead == "L",]

#Matches species material IDs in areaL against segregated species material IDs that
#should have live status in the revised master species list (SpeciesMenu). A FALSE 
#reading indicates location in the sample area menu without a corresponding entry in the 
#revised master list (i.e. the only way this can occur is if status is entered incorrectly).
areaLt <- areaL[,3] %in% liveSpecies
#If there are no FALSE readings, there are no errors.

#Isolate all dead species codes from the sampling area menu.
#This isolates the table with rows where status = D
areaD <- a.area[a.area$liveORdead == "D",]

#Matches species material IDs in areaD against segregated species material IDs that
#should have dead status in the revised master species list (SpeciesMenu). A FALSE 
#reading indicates location in the sample area menu without a corresponding entry in the 
#revised master list (i.e. the only way this can occur is if status is entered incorrectly).
areaDt <- areaD[,3] %in% deadSpecies
#If there are no FALSE readings, there are no errors.

#Combine test results
area.s <- rbind(areaL,areaD)
area.t <- c(areaLt, areaDt)
area.e <- c(rep("Status should be LIVE",length(areaL[,3])),
rep("Status should be DEAD",length(areaD[,3])))

#Register errors in the log.
if(all(area.t==T))
{
EL[[10]] <- paste("No errors: species material in the sampling area menu have no",
"detectable status errors")
} else 
{
EL[[10]] <- data.frame(
area.s[,1:2][area.t == F,], 
Site_Name = I(mapply(function(y){a.site[,2][a.site[,1] == y]},area.s[,2][area.t == F])), 
speciesMaterialID = area.s[,3][area.t == F], 
Common_Name = I(mapply(function(y){SpeciesMenu[,2][SpeciesMenu[,1] == y]},
area.s[,3][area.t == F])), 
Status = area.s[,4][area.t == F], 
Error = I(area.e[area.t == F]))
}
names(EL)[[10]] <- "10. Test for status errors in sampling area menu."

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 8 >>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if(is.list(EL[[10]]) == F)
{

#Next
#Isolate all live species codes from the clip plot dataset (a.clip).
#This isolates the table with rows where status = L
clipL <- a.clip[a.clip$liveORdead == "L",]

#Matches species material IDs in clipL against segregated species material IDs that
#should have live status in the revised master species list (SpeciesMenu). A FALSE 
#reading indicates a sample in the clip plot dataset without a corresponding entry in the 
#revised master list (i.e. the only way this can occur is if status is entered incorrectly).
clipLt <- clipL[,3] %in% liveSpecies

#Next
#Isolate all dead species codes from the clip plot data (sdc.clipPlot).
#This isolates the table with rows where status = D
clipD <- a.clip[a.clip$liveORdead == "D",]

#Matches species material IDs in clipD against segregated species material IDs that
#should have dead status in the revised master species list (SpeciesMenu). A FALSE 
#reading indicates a sample in the clip plot dataset without a corresponding entry in the 
#revised master list (i.e. the only way this can occur is if status is entered incorrectly).
clipDt <- clipD[,3] %in% deadSpecies

#Next, identify instances where sample status is not recognized (i.e. it is not L and D).
clipA <- a.clip[a.clip[,4] %in% c("L","D") == F,]
clipAt <- a.clip[,4] %in% c("L","D")

#Combine test results
clip.s <- rbind(clipL, clipD)
clip.t <- c(clipLt, clipDt)
clip.e <- c(rep("Status should be LIVE",length(clipL[,3])),
rep("Status should be DEAD",length(clipD[,3])))
 
#Error message for when status is not recognized.
clip.a <- rep("Status not recognized",length(clipA[,3]))

#Register errors in the log.
if(all(c(clip.t, clipAt)==T))
{
EL[[11]] <- paste("No errors: species material in the clip plot dataset have no",
"detectable status errors")
} else 
{
EL[[11]] <- 
rbind(data.frame(
clipPlotID = clip.s[,1][clip.t == F],
Site_ID = mapply(function(y){sID[a.clip[,1] == y]},clip.s[,1][clip.t == F]), 
Site_Name = mapply(function(y){sNo[a.clip[,1] == y]},clip.s[,1][clip.t == F]), 
Plot_ID = mapply(function(y){pID[a.clip[,1] == y]},clip.s[,1][clip.t == F]), 
Plot_No = mapply(function(y){pNo[a.clip[,1] == y]},clip.s[,1][clip.t == F]), 
Sampling_Episode_ID = clip.s[,2][clip.t == F], 
Sampling_Episode = mapply(function(y){seID[a.clip[,1] == y]},clip.s[,1][clip.t == F]), 
Species_Material_ID = clip.s[,3][clip.t == F],
Common_Name = mapply(function(y){sc.sn[a.clip[,1] == y]},clip.s[,1][clip.t == F]), 
clip.s[,4:7][clip.t == F,], 
Error = I(clip.e[clip.t == F])), 
data.frame(
clipPlotID = clipA[,1], 
Site_ID = sID[a.clip[,1] == clipA[,1]], 
Site_Name = sNo[a.clip[,1] == clipA[,1]], 
Plot_ID = pID[a.clip[,1] == clipA[,1]], 
Plot_No = sNo[a.clip[,1] == clipA[,1]], 
Sampling_Episode_ID =  clipA[,2], 
Sampling_Episode = seID[a.clip[,1] == clipA[,1]], 
Species_Material_ID = clipA[,3], 
Common_Name = sc.sn[a.clip[,1] == clipA[,1]], 
clipA[,4:7], 
Error = I(clip.a)))
}
names(EL)[[11]] <- "11. Test for status errors in clip plot dataset."

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 9 >>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if(is.list(EL[[11]]) == F)
{

###################################################################################################
#Step 3g: Check and make sure that plant part (rachis or frond) status is correct. 

#Things to check:
#1: Plant parts should only be associated with live palm.
#**************************************************************************************************#
#**************************************************************************************************#
#2: There is a frond and fronds plant part, there should only
#be frond. Have Paige fix this.
#The standard is frond. There should be no samples labeled as "fronds"
#**************************************************************************************************#
#**************************************************************************************************#

#Make sure frond or rachis is only associated with live palmetto (ID = 18).

#Make sure fronds, frond, or rachis is only associated with live palmetto (ID = 18).
#First isolate species codes for live palm.
lp.code <- SpeciesMenu[,1][grep("live palm",SpeciesMenu[,2])]

#Subset data where speciesMaterialID is not palm
np <- a.clip[a.clip[,3] %in% lp.code == F,]
#Subset any of these samples that have a plant part listed.
npwpp <- np[is.na(np[,5]) == F,]
error.np <- rep("Non-palm species material in rachis or frond category",length(npwpp[,1]))

#Are there any live palms labeled with a "blank" plant part.
#First isolate all palmetto samples
ap <- a.clip[a.clip[,3] %in% lp.code == T,]

#Subset any of these samples that have a blank or incoreectly named plant part.
lpwopp <- ap[ap[,5] %in% c("frond", "rachis") == F,]
error.pp <- rep("Live palm plant part is blank or has unaccepted name",length(lpwopp[,1]))

#Combine test results.
pp.s <- rbind(npwpp, lpwopp)
pp.e <- c(error.np, error.pp)

#Register errors in the log.
if(length(pp.s[,1]) == 0)
{
EL[[12]] <- paste("No errors: no samples in the clip plot dataset have errors",
"detectable in the plant part category")
} else 
{
EL[[12]] <- 
data.frame(
clipPlotID = pp.s[,1], 
Site_ID = mapply(function(y){sID[a.clip[,1] == y]},pp.s[,1]), 
Site_Name = mapply(function(y){sNo[a.clip[,1] == y]},pp.s[,1]), 
Plot_ID = mapply(function(y){pID[a.clip[,1] == y]},pp.s[,1]), 
Plot_No = mapply(function(y){pNo[a.clip[,1] == y]},pp.s[,1]), 
Sampling_Episode_ID = pp.s[,2], 
Sampling_Episode = mapply(function(y){seID[a.clip[,1] == y]},pp.s[,1]), 
Species_Material_ID = pp.s[,3], 
Common_Name = mapply(function(y){sc.sn[a.clip[,1] == y]},pp.s[,1]), 
pp.s[,4:7], 
Error = I(pp.e))
}
names(EL)[[12]] <- "12. Test for plant part errors in clip plot dataset."


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 10 >>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if(is.list(EL[[12]]) == F)
{

###################################################################################################
#Step 3h: Check bag numbers. There are two potential errors. The first is that a bag
#was labeled as a number other than one when in fact there was only one bag. The other
#is that there were multiple bags but they were not labeled appropriately.

#The loop below analyzes both situations and errors are isolated.

#Creates an object that lists all unique sampling episodes
seid <- sort(unique(a.clip$plotSamplingEpisodeID))

#Create files that will store information on errors that have been detected.
BagLog <- list()
SaEpLog <- vector(mode = "numeric")

#Error checking loop.
for(i in 1:length(seid))
{
#isolate samples for sampling episide ID[i]
a1 <- a.clip[a.clip[,2] == seid[i],]
#isolate unique species material combinations.
a3 <- unique(a1[,3])
#count the number of samples for each unique combination.
a4 <- mapply(function(y) length(a1[,3][a1[,3] == y]), a3)
#make a list showing bag numbers per sample. Conditional statement handles
#instances where there is just one sample per plot
if(length(a4) ==  1)
a5 <- sort(a1[,6]) else
a5 <- mapply(function(y) sort(a1[,6][a1[,3] == y]), a3)
#make a list showing what the bag numbers per sample shuold be.
if(length(a4) ==  1)
a6 <- 1:length(a5) else
a6 <- mapply(function(y) 1:length(a5[[y]]), 1:length(a5))
#On occassion a plot may have multiple bags and only one sample material. In this 
#case a5 will be a matrix and a6 will be a vector. This code just puts them on the
#same footing so they can be compared.
#Test expected vs actual bag numbers. NA indicates that bags have not been
#numbered properly
if(length(a4) ==  1)
a7 <- all((a6 %in% a5)==T) else
a7 <- mapply(function(y) all((a6[[y]] %in% a5[[y]]) == T), 1:length(a5))
#Add the common name to the table.
a8 <- mapply(function(y) SpeciesMenu$commonName[
SpeciesMenu$speciesMaterialID == a3[y]],1:length(a3))

#If any errors are detected read the values into the
#the list[[i]] corresponding with the sampling episode, if not then a zero is entered
#into the list for that sampling episode.
if(any(a7 == F))
{
BagLog[[i]] <- data.frame(CommonName = I(a8[a7 == F]), SpeciesCode = a3[a7 == F], 
NumberOfSampels = a4[a7 == F], NumberedCorrectly = a7[a7 == F])
SaEpLog[[i]] <- seid[i]
} else 
{
BagLog[[i]] <- 0
SaEpLog[[i]] <- 0
}
#Set up Header so each object in the list shows information that can be used to locate
#the samples
z <- a.episode$plotID[
a.episode$plotSamplingEpisode == seid[i]]#plotID
a <- a.episode$samplingEpisodeID[
a.episode$plotSamplingEpisodeID==seid[i]]#episodeID
b <- a.plot$plotNo[a.plot$plotID==z]#plot number
cc <- a.plot$siteID[a.plot$plotID==z]#site ID
d <- as.character(a.site$siteName[a.site$siteID==cc])#site number
names(BagLog)[[i]] <- paste(c("SamplingEpisode: ",
seid[i],"; Site: ",d,"; Plot: ",b,"; Episode: ",a),collapse="")
}

#Subset the list so that only sampling episodes with multiple samples are shown.
cbl <- lapply(BagLog, function(y) length(y) > 1)
condensedBL <- BagLog[unlist(cbl)]#comprehensive list of species with sample bag numbering 
#errors. List provides details about species and location (i.e. site, plot, and episode)
condSL <- SaEpLog[SaEpLog != 0]#summary that only shows sampling episode IDs with errors.

#If condSL has a length of one or greater then errors are present and are recorded in the 
#ErrorLog

#Register errors in the log.
if(length(condSL) == 0)
{
EL[[13]] <- paste("No errors: no single samples have a bag number other than 1 and",
"samples with multiple bags are sequentially numbered starting with one.")
} else 
{
EL[[13]] <- condensedBL
}
names(EL)[[13]] <- "13. Test for correct bag numbers in clip plot dataset."
#3/13/2012: one error comes up at site E508A, plot 913(year 1). There are
#two grass samples and both are labelled as bag 1. This has been corrected
#in Access DB (12/05/2012)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 11 >>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if(is.list(EL[[13]]) == F)
{

###################################################################################################
#Step 3i: Check to see if any sample type designations (Sub and All) are mislabeled. 

#First we need an alternative method of isolating All and Sub samples. Once this is
#done we will verify that samples are labeled correctly in the sampleType column.

#Isolate data where a subsample was collected based on data in subsample columns.
ssA <- a.clip[mapply(function(y) any(is.na(a.clip[y,8:11]) == F), 1:length(a.clip[,1])),]

#Next isolate all data that fits requirements of All samples (i.e. all NA values in the
#subsample columns (#8-11). 
ssB <- a.clip[mapply(function(y) all(is.na(a.clip[y,8:11]) == T), 1:length(a.clip[,1])),]

#Next test sample type label against what the data suggest the sample type is.
ssC <- rbind(ssA[ssA[,7] != "Sub",], ssB[ssB[,7] != "All",])

#Attach pertinent data to error file.
ssCep <- mapply(function(y) a.episode[,3][a.episode[,1] %in% y], ssC[,2])
ssCpI <- mapply(function(y) a.episode[,2][a.episode[,1] %in% y], ssC[,2])
ssCpN <- mapply(function(y) a.plot[,3][a.plot[,1] %in% y], ssCpI)
ssCsI <- mapply(function(y) a.plot[,2][a.plot[,1] %in% y], ssCpI)
ssCsN <- mapply(function(y) a.site[,2][a.site[,1] %in% y], ssCsI)
ssCsp <- mapply(function(y) SpeciesMenu[,2][SpeciesMenu[,1] %in% y], ssC[,3])

#Log any errors
if(length(ssCpN) > 0)
{
ssCC <- cbind(SiteName = ssCsN, PlotNo = ssCpN, Episode = ssCep, CommonName = ssCsp, ssC[,4:13])
} else
{
ssCC <- vector()
}

#Register errors in the log.
if(length(ssCC) == 0)
{
EL[[14]] <- paste("No errors: no samples have incorrect sample type (i.e. All and Sub)",
"labels.")
} else 
{
EL[[14]] <- ssCC
}
names(EL)[[14]] <- "14. Test for correct sample type labels."

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 12 >>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if(is.list(EL[[14]]) == F)
{

###################################################################################################
#Step 3j: Check to see if any sub-sample weights in individual categories (i.e. total, 
#subsample wet-weight, and subsample dry-weight) have tare and gross weights that were entered in
#reverse order. In this step we ignore any subsamples with NA (missing values).

#Okay, check for weights that were entered in reverse order for the subsamples.
#First, check the total wet weights.
ssA1 <- ssA[mapply(function(y) all(is.na(ssA[y,8:9]) == F), 1:length(ssA[,1])),] 
wsa <- ssA1[(ssA1[,8] > ssA1[,9]*1000) == T,]#col 9 is in Kg, need to multiply by 1000.

#Second, check the sub-sample wet weights.
ssA2 <- ssA[mapply(function(y) all(is.na(ssA[y,10:11]) == F), 1:length(ssA[,1])),] 
wss <- ssA2[ssA2[,10] > ssA2[,11],] 

#Third, check the sub-sample dry weights.
ssA3 <- ssA[mapply(function(y) all(is.na(ssA[y,12:13]) == F), 1:length(ssA[,1])),] 
dss <- ssA3[ssA3[,12] > ssA3[,13],] 

#Combine recorded errors
css <- rbind(wsa, wss, dss)

#Attach pertinent data to error file.
cssty <- c(rep("Total weight reversed", length(wsa[,1])), 
rep("Subsample wet-weight reversed", length(wss[,1])), 
rep("Subsample dry-weight reversed", length(dss[,1])))
cssep <- mapply(function(y) a.episode[,3][a.episode[,1] %in% y], css[,2])
csspI <- mapply(function(y) a.episode[,2][a.episode[,1] %in% y], css[,2])
csspN <- mapply(function(y) a.plot[,3][a.plot[,1] %in% y], csspI)
csssI <- mapply(function(y) a.plot[,2][a.plot[,1] %in% y], csspI)
csssN <- mapply(function(y) a.site[,2][a.site[,1] %in% y], csssI)
csssp <- mapply(function(y) SpeciesMenu[,2][SpeciesMenu[,1] %in% y], css[,3])

#Log any errors
if(length(css[,1]) > 0)
{
cssC <- cbind(ErrorType = cssty, SiteName = csssN, PlotNo = csspN, Episode = cssep, 
CommonName = csssp, css[,4:13])
} else
{
cssC <- vector()
}

#Register errors in the log.
if(length(cssC) == 0)
{
EL[[15]] <- paste("No errors: no subsamples have mass data that appears to be reversed")
} else 
{
EL[[15]] <- cssC
}
names(EL)[[15]] <- paste("15. Test for subsample mass tare and gross values that have been", 
"enetered in reverse order.")

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 13 >>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if(is.list(EL[[15]]) == F)
{

###################################################################################################
#Step 3k: Check to see if any sub-sample weights among individual categories (i.e. total, 
#subsample wet-weight, and subsample dry-weight) have were entered in the wrong category In this 
#step we ignore any subsamples with NA (missing values).

#Difficult to look for mix ups between the total sample and subsample wet weight
#since the gross total is measured in Kgs. Instead a ceiling value will weed out
#potential mix ups since total gross weights are typically < 5 Kg and gross subsample
#wet/dry-weights are typically > 100g. If either subsample weights are accidently
#entered into this category they should be flagged by this process. The ceiling 
#weight is entered in step 1 (administrative tasks)
ssn1 <- ssA[mapply(function(y) all(is.na(ssA[y,8:9]) == F), 1:length(ssA[,1])),] 
muf <- ssn1[ssn1[,9] > GrsWt_Ceiling,]

#Second, check the mix ups between the wet-weight and dry-weight subsamples. This is
#straightforward since the net of the dry-weight must be less than the net of the
#wet-weight.
ssn2 <- ssA[mapply(function(y) all(is.na(ssA[y,10:13]) == F), 1:length(ssA[,1])),] 
mus <- ssn2[(ssn2[,11] - ssn2[,10]) < (ssn2[,13] - ssn2[,12]),]

#Combine recorded errors
mua <- rbind(muf, mus)

#Attach pertinent data to error file.
muaty <- c(rep("Potential mix up between total weight and another category", length(muf[,1])), 
rep("Subsample categories wet-weight and dry-weight appear reversed", length(mus[,1])))
muaep <- mapply(function(y) a.episode[,3][a.episode[,1] %in% y], mua[,2])
muapI <- mapply(function(y) a.episode[,2][a.episode[,1] %in% y], mua[,2])
muapN <- mapply(function(y) a.plot[,3][a.plot[,1] %in% y], muapI)
muasI <- mapply(function(y) a.plot[,2][a.plot[,1] %in% y], muapI)
muasN <- mapply(function(y) a.site[,2][a.site[,1] %in% y], muasI)
muasp <- mapply(function(y) SpeciesMenu[,2][SpeciesMenu[,1] %in% y], mua[,3])

#Log any errors
if(length(mua[,1]) > 0)
{
muaC <- cbind(ErrorType = muaty, SiteName = muasN, PlotNo = muapN, Episode = muaep, 
CommonName = muasp, mua[,4:13])
} else
{
muaC <- vector()
}

#Register errors in the log.
if(length(muaC) == 0)
{
EL[[16]] <- paste("No errors: no subsample types appear to be miscategorized")
} else 
{
EL[[16]] <- muaC
}
names(EL)[[16]] <- paste("16. Test for subsample types that have been incorrectly categorized.", 
"For example, the subsample wet-weight gross and tare weights were accidently enetered into.", 
"the total weight category.")

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 14 >>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if(is.list(EL[[16]]) == F)
{

###################################################################################################
#Step 3l: Check "All" samples for problems with total oven dry weights including:
#1: weights that are blank (i.e. not entered).
#2: reverse order weights

#This subsets "All" clip plot data sample types without NA values in the dry oven-
#weight columns. 
tod1 <- ssB[mapply(function(y) all(is.na(ssB[y,12:13]) == F), 1:length(ssB[,1])),]

#This subsets "All" clip plot data sample types with NA values (reverse of above) 
#in the dry oven-weight columns. 
tod2 <- ssB[mapply(function(y) any(is.na(ssB[y,12:13]) == T), 1:length(ssB[,1])),]

#Subset the data again to show all entries where gross and tare 
#oven-dry weights were reversed.
sc <- tod1[tod1[,12] > tod1[,13],]

#Combine recorded errors
ase <- rbind(tod2, sc)

#Attach pertinent data to error file.
asety <- c(rep("Sample is missing oven-dry weight value(s)", length(tod2[,1])), 
rep("Sample dry-weight gross and tare appear reversed", length(sc[,1])))
aseep <- mapply(function(y) a.episode[,3][a.episode[,1] %in% y], ase[,2])
asepI <- mapply(function(y) a.episode[,2][a.episode[,1] %in% y], ase[,2])
asepN <- mapply(function(y) a.plot[,3][a.plot[,1] %in% y], asepI)
asesI <- mapply(function(y) a.plot[,2][a.plot[,1] %in% y], asepI)
asesN <- mapply(function(y) a.site[,2][a.site[,1] %in% y], asesI)
asesp <- mapply(function(y) SpeciesMenu[,2][SpeciesMenu[,1] %in% y], ase[,3])

#Log any errors
if(length(ase[,1]) > 0)
{
aseC <- cbind(ErrorType = asety, SiteName = asesN, PlotNo = asepN, Episode = aseep, 
CommonName = asesp, ase[,4:13])
} else
{
aseC <- vector()
}

#Register errors in the log.
if(length(aseC) == 0)
{
EL[[17]] <- paste("No errors: no All sample types appear to be weight value errors")
} else 
{
EL[[17]] <- aseC
}
names(EL)[[17]] <- paste("17. Test for the following errors in the weight columns of All", 
"sample types: missing values and reverse order of weights.")

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE >>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#No gate needed. Missing values* are handled in the impute section.
#*Except missing tare values, those can be determined manually.

#*************************************************************************************************#
#*************************************************************************************************#
#2
#There are three locations where values were missing.
#A319 (post), plot 325, burned palm. ClipPlot ID 1098
#E403B (pre), plot 819, live forb. ClipPlot ID 3240
#S411 (1 year), plot 929, dead shrub. ClipPlot ID 5811
#For the first (clip plot ID 1098) I am using values found on an earlier version 
#of the Access DB For each of the last two (clip plot ID 3240 and 5811) calculate 
#the average for the site and use that value. These issues are dealt with in step 
#four by imputing mean site data. These actions are taken later in theis code
#(step 4a).
#*************************************************************************************************#
#*************************************************************************************************#

###################################################################################################
#Step 3m: Check sub-sample data for outliers:

#Subset sub-samples without missing values.
ssD <- ssA[mapply(function(y) all(is.na(ssA[y,8:13]) == F), 1:length(ssA[,1])),] 

#Moderate and extreme outlier multiplier.
mom <- 1.5
eom <- 3.0

#Subset live and dead subsamples from ssD.
ssLive <- ssD[ssD$liveORdead == "L",]
ssDead <- ssD[ssD$liveORdead == "D",]

#Calculate net subsample wet-weights.
netwetLive <- ssLive[,11] - ssLive[,10]
netwetDead <- ssDead[,11] - ssDead[,10]

#Calculate net subsample dry-weights.
netdryLive <- ssLive[,13] - ssLive[,12]
netdryDead <- ssDead[,13] - ssDead[,12]

#Calculate dry;total weight ratio
dwrLive <- netdryLive/netwetLive
dwrDead <- netdryDead/netwetDead

#Calculate quantiles and interquartile range for live subsamples
qLive <- quantile(dwrLive)
iqrLive <- IQR(dwrLive)

#Calculate lower and upper limits of data considered not to be an outlier for live subsamples
MLLLive <- qLive[2] - iqrLive*mom
MULLive <- qLive[4] + iqrLive*mom
ELLLive <- qLive[2] - iqrLive*eom
EULLive <- qLive[4] + iqrLive*eom

#Identify outliers for live subsamples
outlierLive <- ifelse((dwrLive < ELLLive) == T | (dwrLive > EULLive) == T, "Extreme", 
ifelse((dwrLive < MLLLive) == T | (dwrLive > MULLive) == T, "Moderate", "No"))

#Calculate quantiles and interquartile range for dead subsamples
qDead <- quantile(dwrDead)
iqrDead <- IQR(dwrDead)

#Calculate lower and upper limits of data considered not to be an outlier for dead subsamples
MLLDead <- qDead[2] - iqrDead*mom
MULDead <- qDead[4] + iqrDead*mom
ELLDead <- qDead[2] - iqrDead*eom
EULDead <- qDead[4] + iqrDead*eom

#Identify outliers for dead subsamples
outlierDead <- ifelse((dwrDead < ELLDead) == T | (dwrDead > EULDead) == T, "Extreme", 
ifelse((dwrDead < MLLDead) == T | (dwrDead > MULDead) == T, "Moderate", "No"))

#Create an informative data table for live subsamples. 
ldw <- data.frame("ClipPlotID" = ssLive[,1], 
"siteName" = I(sNo[a.clip[,1] %in% ssLive[,1]]), 
"plotNo" = pNo[a.clip[,1] %in% ssLive[,1]], 
"Episode" = seID[a.clip[,1] %in% ssLive[,1]],
"Species"= I(sc.sn[a.clip[,1] %in% ssLive[,1]]),
ssLive[,4:7], ssLive[,10:11], "netWet" = netwetLive, ssLive[,12:13], "netDry" = netdryLive, 
"DryWetRatio" = round(dwrLive,2), "Outlier" = outlierLive)

#Create an informative data table for dead subsamples. 
ddw <- data.frame("ClipPlotID" = ssDead[,1], 
"siteName" = I(sNo[a.clip[,1] %in% ssDead[,1]]), 
"plotNo" = pNo[a.clip[,1] %in% ssDead[,1]], 
"Episode" = seID[a.clip[,1] %in% ssDead[,1]],
"Species"= I(sc.sn[a.clip[,1] %in% ssDead[,1]]),
ssDead[,4:7], ssDead[,10:11], "netWet" = netwetDead, ssDead[,12:13], "netDry" = netdryDead, 
"DryWetRatio" = round(dwrDead,2), "Outlier" = outlierDead)

#Combine and data tables and put back in original order (by clip plot ID)
dw <- rbind(ldw,ddw)
dw <- dw[order(dw[,1]),]

#Register errors in the log.
if(length(dw[,17][dw[,17] != "No"]) < 1)
{
EL[[18]] <- paste("No errors: no sub sample dry:wet weight ratios are outliers")
} else 
{
EL[[18]] <- dw[dw[,17] != "No",]
}
names(EL)[[18]] <- paste("18. Test for outliers in subsample dry:wet weight ratios")

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE >>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#No gate is needed here.

###################################################################################################
###################################################################################################
#Note before starting this step all above errors should be corrected in the 
#Datasheets and Access DB. Not doing so can cause cascading errors in the steps
#below.

###################################################################################################
###################################################################################################
#STEP #4: Calculate oven-dry net weights for subsamples
#Includes corrections to errors identified in step that must be approximated.

###################################################################################################
#Step 4a: calculate total dry weights for subsamples:

#For sub-samples with NA values in sub-sample weights.
#Calculate mean dry:wet ratios for live and dead sub-samples, these values
#will be applied to sub-samples that have missing weights.
meanLiveRatio <- mean(dwrLive)
meanDeadRatio <- mean(dwrDead)

#Subset subsamples with NA values in any of the weight columns (8 thru 13).
ssE <- ssA[mapply(function(y) any(is.na(ssA[y,8:13]) == T), 1:length(ssA[,1])),]

#Calculate total dry weight for these samples based on total wet weights
#and mean ratios
totalDryWeightNA <- vector(mode = "numeric")
for(i in 1:length(ssE[,1]))
{
if(ssE$liveORdead[i] == "L")
totalDryWeightNA[i] <- (ssE[i,9]*1000 - ssE[i,8])*meanLiveRatio else
totalDryWeightNA[i] <- (ssE[i,9]*1000 - ssE[i,8])*meanDeadRatio
}

#For live sub-samples.
totalDryWeightLive <- (ssLive[,9]*1000 - ssLive[,8])*dwrLive

#For dead sub-samples.
totalDryWeightDead <- (ssDead[,9]*1000 - ssDead[,8])*dwrDead

totalDryWeightSub <- c(totalDryWeightDead,totalDryWeightLive,totalDryWeightNA)
SubID <- c(ssDead[,1],ssLive[,1],ssE[,1])

SubWeight <- cbind("Total" = round(totalDryWeightSub,2), "clipPlotID" = SubID)

###################################################################################################
###################################################################################################
#STEP #5: Calculate total dry weights for all samples.

###################################################################################################
#STEP #5a: Add two more columns to a.clip and partition total dry weights from
#sub-sample dry weights which are currently in the same column.

#First rename cols 10 and 13 to reflect nature of data
#Also, units in row 13 are incorrect, this data is in grams not kg. 
names(a.clip[,10]) <- "wet_subsample_bagged_tare_g"
names(a.clip[,11]) <- "wet_subsample_bagged_gross_g"
names(a.clip[,12]) <- "dry_subsample_bagged_tare_g"
names(a.clip[,13]) <- "dry_subsample_bagged_gross_g"

#Next populate vectors with total dry weight data.
tdt <- ifelse(a.clip[,7] == "All", a.clip[,12], NA)
tdg <- ifelse(a.clip[,7] == "All", a.clip[,13], NA)
a.clip[,12] <- ifelse(a.clip[,7] == "All", NA, a.clip[,12])
a.clip[,13] <- ifelse(a.clip[,7] == "All", NA, a.clip[,13])

#Add partitioned data back into a.clip as new columns
a.clip <- cbind(a.clip, dry_total_tare_g = tdt, dry_total_gross_g = tdg)

#Add calculated dry total weights from subsamples into the new total dry weight column.
#To do this just use the add the net total oven-dry weight calculated in step
#3a (SubWeight object) into the column 15 of a.clip and use 
#a tare weight of zero (column 14). Now the sub and all samples have equivalent values
#in cols 14, and 15 and they can both be operated on in the same manner
#for loading by plot and site.
for (i in SubWeight[,2])
{
a.clip[,14][a.clip[,1] == i] <- 0
a.clip[,15][a.clip[,1] == i] <- SubWeight[,1][SubWeight[,2] == i]
}

###################################################################################################
#STEP #5b: Take blank All sample weights identified in Error Check #17 and impute values
#based on mean weight of other plots at the site. 

#Restructure error log from step 3l to display a larger amount of sample information.
aseI <- cbind(ErrorType = asety, SiteID = asesI, SiteName = asesN, 
PlotID = asepI, PlotNo = asepN, Episode = aseep, CommonName = asesp, 
ase)


#Isolate errors in EL17 caused by missing sample weights (i.e. not reversed tare/gross weights)
#These are errors that can be repaired with the impute function.
aseIb <- aseI[asety == "Sample is missing oven-dry weight value(s)",]

#Remove error message to make table easier to display.
imp <- aseIb[,-1]

#Determine the number of samples that can be used for imputing data at site/episode,
#site/all and all sites/episodes.

#First step is to isolate all samples with NA values for weights so they are not
#counted as possible data that can be used for imputing calculations.
no.clip <- mapply(function(y) a.clip[,1][a.clip[,2] == imp[y,8] & a.clip[,3] == imp[y,9]],
1:length(imp[,1]))

#Convert list into vector. This shows all clip plot IDs that cannot be incorporated
#into imputed data. This includes those with NA in the the sample weights and additional
#bags attached to those samples.
no.clip <- unlist(no.clip)

#Remove these samples from total clip plot table
ok.imp <- a.clip[(a.clip[,1] %in% no.clip) == F,]

#Second, determine the different types of missing data
#1 = tare value missing/gross value present
#2 = tare value present/gross value missing
#3 = tare and gross values missing
e.type <- vector(mode = "numeric")
for(i in 1:length(imp[,1]))
{
if(all(is.na(imp[i,18:19]) == T))
e.type[i] <- 3 else
if(is.na(imp[i,18]) == T)
e.type[i] <- 2 else
e.type[i] <- 1
}

#Third, count the number of imputable samples within the site and sampling episode.
n.se <- mapply(function(y) 
{
length(ok.imp[,1][
ok.imp[,2] %in% a.episode[,1][
a.episode[,2] %in% a.plot[,1][a.plot[,2] == imp[y,1]] & a.episode[,3] == imp[y,5]] & 
ok.imp[,3] == imp[y,9] & ok.imp[,6] == 1])
}, 
1:length(imp[,1]))

#Fourth, count the number of imputable samples within the site and all sampling episodes.
n.sa <- mapply(function(y) 
{
length(ok.imp[,1][
ok.imp[,2] %in% a.episode[,1][a.episode[,2] %in% a.plot[,1][a.plot[,2] == imp[y,1]]] & 
ok.imp[,3] == imp[y,9] & ok.imp[,5] == imp[y,11] & ok.imp[,6] == 1])
}, 
1:length(imp[,1]))

#Fifth, count the number of imputable samples among all sites and sampling episodes.
n.aa <- mapply(function(y) 
{
length(ok.imp[,1][ok.imp[,3] == imp[y,9] & ok.imp[,5] == imp[y,11] & ok.imp[,6] == 1])
}, 
1:length(imp[,1]))

#Sixth, input values for tare weights
#Error type 1 (tare present) remains unchanged
#Error type 2 (only tare absent) remains NA. This is an error that can be corrected in
#the database. Looking at the data can often provde a good idea of what the tare value
#should be.
#Error type 3 gets a zero value
a.clip[,14][a.clip[,1] %in% imp[,7][e.type == 3]] <- 0

#Seventh, impute values for samples that can be imputed from data from the same site and 
#episode.

#This loop runs for each missing gross weight value and imputes a value.
for(i in 1:length(n.aa))
{
#Imputes value using samples of the species material from site and sampling episode.
if(n.se[i] > 0)
{
#Subset data for site, sampling episode, and species material.
ac.temp <- a.clip[,c(1,2,3,6,14,15)][
a.clip[,2] %in% a.episode[,1][a.episode[,2] %in% a.plot[,1][a.plot[,2] == imp[i,1]] & 
a.episode[,3] == imp[i,5]] & a.clip[,3] == imp[i,9],]
#Calculate the net weights, this will be used as a gross weight since tare is set to zero.
ac.temp[,7] <- ac.temp[,6] - ac.temp[,5]
#Determine the sum of net weights. This is for instances when there are multiple bags per subsample.
ac.7 <- mapply(function(y) sum(ac.temp[,7][ac.temp[,2] == y]), unique(ac.temp[,2]))
#Remove extra columns (tare wt and gross wt) from subset data and condense rows to reflect unique samples.
ac.temp <- ac.temp[,1:4][ac.temp[,4] == 1,]
#Add net weights.
ac.temp <- data.frame(ac.temp, NetWt = ac.7)
#Impute missing data values.
ac.temp <- impute(ac.temp, what = "mean")
} else
{
#Imputes value using samples of the species material from site and all sampling episodes. 
#See explanations above for each action.
if(n.sa[i] > 0)
{
ac.temp <- a.clip[,c(1,2,3,6,14,15)][
a.clip[,2] %in% a.episode[,1][a.episode[,2] %in% a.plot[,1][a.plot[,2] == imp[i,1]] & 
a.episode[,3] == imp[i,5]] & a.clip[,3] == imp[i,9],]
ac.temp[,7] <- ac.temp[,6] - ac.temp[,5]
ac.7 <- mapply(function(y) sum(ac.temp[,7][ac.temp[,2] == y]), unique(ac.temp[,2]))
ac.temp <- ac.temp[,1:4][ac.temp[,4] == 1,]
ac.temp <- data.frame(ac.temp, NetWt = ac.7)
ac.temp <- impute(ac.temp, what = "mean")
} else
#Imputes value using samples of the species material from all sites and all sampling episodes. 
#See explanations above for each action.
if(n.aa[i] > 0)
{
ac.temp <- a.clip[,c(1,2,3,6,14,15)][
a.clip[,2] %in% a.episode[,1][a.episode[,2] %in% a.plot[,1][a.plot[,2] == imp[i,1]] & 
a.episode[,3] == imp[i,5]] & a.clip[,3] == imp[i,9],]
ac.temp[,7] <- ac.temp[,6] - ac.temp[,5]
ac.7 <- mapply(function(y) sum(ac.temp[,7][ac.temp[,2] == y]), unique(ac.temp[,2]))
ac.temp <- ac.temp[,1:4][ac.temp[,4] == 1,]
ac.temp <- data.frame(ac.temp, NetWt = ac.7)
ac.temp <- impute(ac.temp, what = "mean")
} else
{
#If for some reason no values were imputed the data will remain as is and will be logged in the
#error file (EL12) below.
ac.temp <- a.clip[,c(1,2,3,6,14,15)][
a.clip[,2] %in% a.episode[,1][a.episode[,2] %in% a.plot[,1][a.plot[,2] == imp[i,1]] & 
a.episode[,3] == imp[i,5]] & a.clip[,3] == imp[i,9],]
ac.temp[,7] <- ac.temp[,6] - ac.temp[,5]
ac.temp <- ac.temp[,c(1,2,3,4,7),]
colnames(ac.temp) <- c(names(ac.temp[,1:4]),"NetWt")
}
}
#Update a.ground with imputed value
a.clip[,15][a.clip[,1] == imp[i,7]] <- round(ac.temp[,5][ac.temp[,1] == imp[i,7]],2)
}

wtr <- cbind(imp[,1:13], e.type, n.se, n.sa, n.aa,
a.clip[,14:15][a.clip[,1] %in% imp[,7],])

#Register errors in the log.
if(length(wtr[,1][is.na(wtr[,18]) == T | is.na(wtr[,19]) == T]) < 1)
{
EL19 <- list()
EL19[[1]] <- paste("19. Test for samples with gross weight but no tare weight", 
"or missing weights that could not be imputed.")
EL19[[2]] <- paste("No errors: all missing values have been imputed.")
EL19[[3]] <- wtr[mapply(function(y) all(c(is.na(wtr[y,18]), is.na(wtr[y,19])) == F), 
1:length(wtr[,1])),]
EL[[19]] <- EL19
} else 
{
EL19 <- list()
EL19[[1]] <- paste("19. Test for samples with gross weight but no tare weight", 
"or missing weights that could not be imputed.")
EL19[[2]] <- paste("Values that were not successfully imputed:")
EL19[[3]] <- wtr[is.na(wtr[,18]) == T | is.na(wtr[,19]) == T,]
EL19[[4]] <- paste("Values that were successfully imputed:")
EL19[[5]] <- wtr[mapply(function(y) all(c(is.na(wtr[y,18]), is.na(wtr[y,19])) == F), 
1:length(wtr[,1])),]
EL[[19]] <- EL19
}

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE >>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#No gate is needed here.

###################################################################################################
###################################################################################################
#STEP #6: Calculate weights/loading at the plot level.

###################################################################################################
#STEP #6a: Condense clip plot data to the plot level. Essentially this means summing all
#samples with multiple bags which currently exist as one entry per bag in the clip plot
#dataset (a.clip). Each important information is condensed here:
#1: Clip Plot ID (zz1)
#2: Sampling Episode ID (zz2_5[1])
#3: Species Code (zz2_5[2])
#4: Live/Dead Staus (zz2_5[3])
#5: Plant Part (zz2_5[4])
#6: Tare weight (zz12)
#7: Gross Weight (zz13)
#8: Plot ID (pIDx)
#9: Plot No. (pNox)
#10: Site ID (sIDx)
#11: Site Name (sNox)
#12: Sampling Episode (seIDx)
#13: Common Species Name (sc.snx)
#14: Has data been imputed? (impx)

#Condense clip plot sample episode IDs, and species codes
#to the sample level.
zz2_5 <- data.frame(EpisodeID = rep(0,length(a.clip[,1][a.clip[,6] ==1])),
SpeciesCode = rep(0,length(a.clip[,1][a.clip[,6] ==1])))
for(e in 1:2)
{
zzz <- mapply(function(y) 
{
mapply(function(x)
{
unique(a.clip[,e+1][a.clip[,2] == y & a.clip[,3] == x])
},unique(a.clip[,3][a.clip[,2] == y]))
},unique(a.clip[,2])
)
zz2_5[,e] <- as.vector(unlist(zzz))
}

#Calculate tare weights by sample (i.e. sum up samples with multiple bags).
zz12 <- mapply(function(y) 
{
mapply(function(x)
{
sum(a.clip[,14][a.clip[,2] == y & a.clip[,3] == x])
},unique(a.clip[,3][a.clip[,2] == y]))
},unique(a.clip[,2])
)
zz12 <- as.vector(unlist(zz12))

#Calculate gross weights by sample (i.e. sum up samples with multiple bags).
zz13 <- mapply(function(y) 
{
mapply(function(x)
{
sum(a.clip[,15][a.clip[,2] == y & a.clip[,3] == x])
},unique(a.clip[,3][a.clip[,2] == y]))
},unique(a.clip[,2])
)
zz13 <- as.vector(unlist(zz13))

#Calculate net weights by sample.
zz14 <- zz13 - zz12

#Condense plot IDs to the sample level (in cases of multiple bags take unique())
pIDx <- mapply(function(y) 
{
mapply(function(x)
{
unique(pID[a.clip[,2] == y & a.clip[,3] == x])
},unique(a.clip[,3][a.clip[,2] == y]))
},unique(a.clip[,2])
)
pIDx <- as.vector(unlist(pIDx))

#Condense plot numbers to the sample level (in cases of multiple bags take unique())
pNox <- mapply(function(y) 
{
mapply(function(x)
{
unique(pNo[a.clip[,2] == y & a.clip[,3] == x])
},unique(a.clip[,3][a.clip[,2] == y]))
},unique(a.clip[,2])
)
pNox <- as.vector(unlist(pNox))

#Condense site IDs to the sample level (in cases of multiple bags take unique())
sIDx <- mapply(function(y) 
{
mapply(function(x)
{
unique(sID[a.clip[,2] == y & a.clip[,3] == x])
},unique(a.clip[,3][a.clip[,2] == y]))
},unique(a.clip[,2])
)
sIDx <- as.vector(unlist(sIDx))

#Condense site names to the sample level (in cases of multiple bags take unique())
sNoC <- as.character(sNo)#This needs to be done b/c unlist() handles factors in a manner
#inconsistent with the goal of this substep.
sNox <- mapply(function(y) 
{
mapply(function(x)
{
unique(sNoC[a.clip[,2] == y & a.clip[,3] == x])
},unique(a.clip[,3][a.clip[,2] == y]))
},unique(a.clip[,2])
)
sNox <- as.vector(unlist(sNox))

#Condense sampling episode to the sample level (in cases of multiple bags take unique())
seIDx <- mapply(function(y) 
{
mapply(function(x)
{
unique(seID[a.clip[,2] == y & a.clip[,3] == x])
},unique(a.clip[,3][a.clip[,2] == y]))
},unique(a.clip[,2])
)
seIDx <- as.vector(unlist(seIDx))

#Condense common species name to the sample level (in cases of multiple bags take unique())
sc.snx <- mapply(function(y) 
{
mapply(function(x)
{
unique(sc.sn[a.clip[,2] == y & a.clip[,3] == x])
},unique(a.clip[,3][a.clip[,2] == y]))
},unique(a.clip[,2])
)
sc.snx <- as.vector(unlist(sc.snx))

#Identifies data points with one or more imputed sample weights
impx <- mapply(function(y) 
{
mapply(function(x)
{
ifelse(any(a.clip[,1][a.clip[,2] == y & a.clip[,3] == x] %in% wtr[,7]),1,2)
},unique(a.clip[,3][a.clip[,2] == y]))
},unique(a.clip[,2])
)
impx <- as.vector(unlist(impx))

#Create a new data frame that lists net weights.
a.sample <- data.frame(SiteName = I(sNox), PlotNo = pNox, SampleEpisode = seIDx, 
CommonName = I(sc.snx), TareWt = zz12, 
GrossWt = zz13, NetWt = zz14, SamplingEpisodeID = zz2_5[,1], PlotID = pIDx, 
SiteID = sIDx, SpeciesCode = zz2_5[,2], Imputed = impx)

#Quick way to test for errors in step 5a. Does the number of samples equal the number of
#entries in the clip plot dataset where bagNo(or a.clip[,6]) = 1?
#Register errors in the log.
EL[[20]] <- ifelse(length(a.sample[,1]) == length(a.clip[,1][a.clip[,6] == 1]),
paste("NO ERROR: Number of samples equals the number of entries where bag number = 1"), 
paste("ERROR: Number of samples is incongruent with expected number of samples based on bag", 
"numbers in the clip plot dataset set to 1"))

names(EL)[[20]] <- "20. Test incongruency for total number of samples in new clip plot dataset."

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE >>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#No gate is needed here.

###################################################################################################
#STEP #6b: Calculate loading in tonnes per ha.

#Calculates loading in tonnes per hectare at the plot level
loadingHa <- mapply(function(y)
{
(a.sample[y,7]/1000000) * (10000/a.area[,5][a.area[,2] == a.sample[y,10] & 
a.area[,3] == a.sample[y,11]])
},1:length(a.sample[,1]))

#Combine loading values into the a.sample object
a.sample <- cbind(a.sample,loadingHa)

#Must be disabled when gates are active
#})# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> END SYSTEM TIME >>>>>>>>

###################################################################################################
###################################################################################################
#STEP #7: Calculate average net weights and summary statistics at the site level.
#This requires populating 20 data points per site for each species to represent zero values.

###################################################################################################
#STEP #7a: Incorporate zero measurements into the plot level data. This requires
#the construction of a symmetrical plot level matrix that contains all possible species
#for each plot. If no sample weights are available for a given species at a given plot then
#that point will be assigned a sample weight of zero.

#Measure the system time. You are creating a data frame that is about ten times larger
#than a.sample.

#This will take about 300 seconds. Or about 210 seconds longer than all preceeding code.
#Must be disabled when gates are active
#system.time({#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>BEGIN SYSTEM TIME>>>>>>>>

#Sampling episode.
vSamplingEpisode <- expand.grid(SpeciesMenu[,1],1:20,a.site[,1], sampleEpisode)[,4]

#Site ID
vSiteID <- expand.grid(SpeciesMenu[,1],1:20,a.site[,1], sampleEpisode)[,3]

#Site Name
vSiteName <- expand.grid(SpeciesMenu[,1],1:20,a.site[,2], sampleEpisode)[,3]

#Plot ID
a.plot <- a.plot[order(a.plot[,2],a.plot[,3]),]#first you need to order plot numbers to sites.
vPlotID <- expand.grid(SpeciesMenu[,1],a.plot[,1],sampleEpisode)[,2]

#Plot Number
vPlotNo <- expand.grid(SpeciesMenu[,1],a.plot[,3],sampleEpisode)[,2]

#Sampling Episode ID
#Need to order episode IDs to plot numbers and site IDs

objects()

#Add a col for plot numbers.
pnus <- mapply(function(y)
{
a.plot[,3][a.plot[,1] == y]
},
a.episode[,2])

#Add a col for site IDs
sids <- mapply(function(y)
{
a.plot[,2][a.plot[,1] == y]
},
a.episode[,2])

b.episode <- cbind(a.episode,pnus,sids)
c.episode <- b.episode[order(b.episode[,3],b.episode[,5],b.episode[,4]),]

vSamplingEpisodeID <- expand.grid(SpeciesMenu[,1],c.episode[,1])[,2]

#Species Code
vSpeciesCode <- expand.grid(SpeciesMenu[,1],1:20,a.site[,1], sampleEpisode)[,1]

#Common Name
vCommonName <- expand.grid(SpeciesMenu[,2],1:20,a.site[,1], sampleEpisode)[,1]

#Tare Weight
vTareWeight <- mapply(function(y)
{
ifelse(length(a.sample[,5][a.sample[,8] == vSamplingEpisodeID[y] & 
a.sample[,11] == vSpeciesCode[y]]) == 1, 
a.sample[,5][a.sample[,8] == vSamplingEpisodeID[y] & a.sample[,11] == vSpeciesCode[y]],
ifelse(length(a.sample[,5][a.sample[,8] == vSamplingEpisodeID[y] & 
a.sample[,11] == vSpeciesCode[y]]) == 0,0,NA))
},
1:length(vSamplingEpisode))

#Gross Weight
vGrossWeight <- mapply(function(y)
{
ifelse(length(a.sample[,6][a.sample[,8] == vSamplingEpisodeID[y] & 
a.sample[,11] == vSpeciesCode[y]]) == 1, 
a.sample[,6][a.sample[,8] == vSamplingEpisodeID[y] & a.sample[,11] == vSpeciesCode[y]],
ifelse(length(a.sample[,6][a.sample[,8] == vSamplingEpisodeID[y] & 
a.sample[,11] == vSpeciesCode[y]]) == 0,0,NA))
},
1:length(vSamplingEpisode))

#Net Weight
vNetWeight <- mapply(function(y)
{
ifelse(length(a.sample[,7][a.sample[,8] == vSamplingEpisodeID[y] & 
a.sample[,11] == vSpeciesCode[y]]) == 1, 
a.sample[,7][a.sample[,8] == vSamplingEpisodeID[y] & a.sample[,11] == vSpeciesCode[y]],
ifelse(length(a.sample[,7][a.sample[,8] == vSamplingEpisodeID[y] & 
a.sample[,11] == vSpeciesCode[y]]) == 0,0,NA))
},
1:length(vSamplingEpisode))

#Imputed (2 = no, 1 = yes)
vImputed <- mapply(function(y)
{
ifelse(length(a.sample[,12][a.sample[,8] == vSamplingEpisodeID[y] & 
a.sample[,11] == vSpeciesCode[y]]) == 1, 
a.sample[,12][a.sample[,8] == vSamplingEpisodeID[y] & a.sample[,11] == vSpeciesCode[y]],
ifelse(length(a.sample[,12][a.sample[,8] == vSamplingEpisodeID[y] & 
a.sample[,11] == vSpeciesCode[y]]) == 0,2,NA))
},
1:length(vSamplingEpisode))

#Metric Loading (ton/ha)
vLoading <- mapply(function(y)
{
ifelse(length(a.sample[,13][a.sample[,8] == vSamplingEpisodeID[y] & 
a.sample[,11] == vSpeciesCode[y]]) == 1, 
a.sample[,13][a.sample[,8] == vSamplingEpisodeID[y] & a.sample[,11] == vSpeciesCode[y]],
ifelse(length(a.sample[,13][a.sample[,8] == vSamplingEpisodeID[y] & 
a.sample[,11] == vSpeciesCode[y]]) == 0,0,NA))
},
1:length(vSamplingEpisode))


SamMat <- data.frame(SiteName = I(vSiteName), PlotNo = vPlotNo, 
SampleEpisode = vSamplingEpisode, CommonName = I(vCommonName), 
TareWt = vTareWeight, GrossWt = vGrossWeight, NetWt = vNetWeight, 
Loading_TonsHa = vLoading, SamplingEpisodeID = vSamplingEpisodeID, 
PlotID = vPlotID, SiteID = vSiteID, SpeciesCode = vSpeciesCode, 
Imputed = vImputed)

#Must be disabled when gates are active
#})# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> END SYSTEM TIME >>>>>>>>

###################################################################################################
#STEP #7b: Calculate loading by site and sampling episode across all possible samples.

#Must be disabled when gates are active
#system.time({#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>BEGIN SYSTEM TIME>>>>>>>>

#Sampling episode
sSamplingEpisode <- expand.grid(SpeciesMenu[,1],a.site[,1], sampleEpisode)[,3]

#Site ID
sSiteID <- expand.grid(SpeciesMenu[,1],a.site[,1], sampleEpisode)[,2]

#Site Name
sSiteName <- expand.grid(SpeciesMenu[,1],a.site[,2], sampleEpisode)[,2]

#Species Code
sSpeciesCode <- expand.grid(SpeciesMenu[,1],a.site[,1], sampleEpisode)[,1]

#Common Name
sCommonName <- expand.grid(SpeciesMenu[,2],a.site[,1], sampleEpisode)[,1]

#Mean Net Weight
sMeanNetWeight <- as.vector(mapply(function(y)
{
mapply(function(x)
{
mapply(function(w)
{
mean(SamMat[,7][SamMat[,12] == w & SamMat[,11] == x & SamMat[,3] == y])
}, unique(SamMat[,12]))
}, unique(SamMat[,11]))
}, unique(SamMat[,3])))

#Mean Loading
sMeanLoading <- as.vector(mapply(function(y)
{
mapply(function(x)
{
mapply(function(w)
{
mean(SamMat[,8][SamMat[,12] == w & SamMat[,11] == x & SamMat[,3] == y])
}, unique(SamMat[,12]))
}, unique(SamMat[,11]))
}, unique(SamMat[,3])))

#Imputed (2 = no, 1 = yes)
sn <- as.vector(mapply(function(y)
{
mapply(function(x)
{
mapply(function(w)
{
length(SamMat[,13][SamMat[,13] == 2 & SamMat[,12] == w & SamMat[,11] == x & SamMat[,3] == y])
}, unique(SamMat[,12]))
}, unique(SamMat[,11]))
}, unique(SamMat[,3])))

#Sample Variance
sVariance <- as.vector(mapply(function(y)
{
mapply(function(x)
{
mapply(function(w)
{
ivar(SamMat[,8][SamMat[,12] == w & SamMat[,11] == x & SamMat[,3] == y], 
length(SamMat[,13][SamMat[,13] == 2 & SamMat[,12] == w & SamMat[,11] == x & SamMat[,3] == y]))
}, unique(SamMat[,12]))
}, unique(SamMat[,11]))
}, unique(SamMat[,3])))

#Sample Standard Deviation
sStdDev <- sqrt(sVariance)

#Sample Coeffiecient of Variation
sCV <- (sStdDev/sMeanLoading)*100

#Sample Standard Error of the Mean
sSSEM <- sStdDev/sqrt(n)

#Confidence Interval
sCI <- CI.xbar(sMeanLoading,a.level,2,sStdDev,sn)

#Symmetry
sG1 <- as.vector(mapply(function(y)
{
mapply(function(x)
{
mapply(function(w)
{
skewness(SamMat[,8][SamMat[,12] == w & SamMat[,11] == x & SamMat[,3] == y])
}, unique(SamMat[,12]))
}, unique(SamMat[,11]))
}, unique(SamMat[,3])))

#Kurtosis
sG2 <- as.vector(mapply(function(y)
{
mapply(function(x)
{
mapply(function(w)
{
kurtosis(SamMat[,8][SamMat[,12] == w & SamMat[,11] == x & SamMat[,3] == y])
}, unique(SamMat[,12]))
}, unique(SamMat[,11]))
}, unique(SamMat[,3])))


#Create a site matrix
SitMat <- data.frame(SiteName = I(sSiteName), SampleEpisode = sSamplingEpisode, 
CommonName = I(sCommonName), NetWt = round(sMeanNetWeight,2), 
Loading_TonsHa = round(sMeanLoading,2),
n = sn, StdDev = round(sStdDev,2), CoefVar = round(sCV,2), LoCI = round(sCI[,1],2), 
UpCI = round(sCI[,2],2), Symmetry = round(sG1,2), Kurtosis = round(sG2,2), 
SiteID = sSiteID, SpeciesCode = sSpeciesCode)

###################################################################################################
###################################################################################################
#STEP #8: Save relevant files

#set a date used to name the file.
dt <- Sys.Date()
tm <- format(Sys.time(), format = "%H.%M.%S", 
tz = "", usetz = FALSE)

#Save the plot-level data.
write.table(SamMat, file = paste("c:\\usfs_sef_data_output\\sef_SampleMatrix_",
dt,"_",tm,".csv",sep = ""), append = F, quote = T, sep = ",", eol = "\n", na = "NA", 
dec = ".", row.names = F, col.names = T, qmethod = c("escape", "double"))#

#Save the site-level data.
write.table(SitMat, file = paste("c:\\usfs_sef_data_output\\sef_SiteMatrix_",
dt,"_",tm,".csv",sep = ""), append = F, quote = T, sep = ",", eol = "\n", na = "NA", 
dec = ".", row.names = F, col.names = T, qmethod = c("escape", "double"))#

#Must be disabled when gates are active
#})# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> END SYSTEM TIME >>>>>>>>

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 14 >>>>>>>>>>>>>>>>
} else
{
print("Gate 14 >> ERROR: Model stopped running because of an error in step 3k: Error[16].")
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 14 >>>>>>>>>>>>>>>>

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 13 >>>>>>>>>>>>>>>>
} else
{
print("Gate 13 >> ERROR: Model stopped running because of an error in step 3j: Error[15].")
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 13 >>>>>>>>>>>>>>>>

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 12 >>>>>>>>>>>>>>>>
} else
{
print("Gate 12 >> ERROR: Model stopped running because of an error in step 3i: Error[14].")
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 12 >>>>>>>>>>>>>>>>

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 11 >>>>>>>>>>>>>>>>
} else
{
print("Gate 11 >> ERROR: Model stopped running because of an error in step 3h: Error[13].")
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 11 >>>>>>>>>>>>>>>>

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 10 >>>>>>>>>>>>>>>>
} else
{
print("Gate 10 >> ERROR: Model stopped running because of an error in step 3g: Error[12].")
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 10 >>>>>>>>>>>>>>>>

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 9 >>>>>>>>>>>>>>>>>
} else
{
print("Gate 9 >> ERROR: Model stopped running because of an error in step 3f: Error[11].")
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 9 >>>>>>>>>>>>>>>>>

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 8 >>>>>>>>>>>>>>>>>
} else
{
print("Gate 8 >> ERROR: Model stopped running because of an error in step 3f: Error[10].")
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 8 >>>>>>>>>>>>>>>>>

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 7 >>>>>>>>>>>>>>>>>
} else
{
print("Gate 7 >> ERROR: Model stopped running because of an error in step 3e: Error[9].")
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 7 >>>>>>>>>>>>>>>>>

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 6 >>>>>>>>>>>>>>>>>
} else
{
print("Gate 6 >> ERROR: Model stopped running because of an error in step 3d: Error[8].")
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 6 >>>>>>>>>>>>>>>>>

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 5 >>>>>>>>>>>>>>>>>
} else
{
print("Gate 5 >> ERROR: Model stopped running because of an error in step 3c: Error[7].")
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 5 >>>>>>>>>>>>>>>>>

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 4 >>>>>>>>>>>>>>>>>
} else
{
print("Gate 4 >> ERROR: Model stopped running because of an error in step 3b: Error[5].")
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 4 >>>>>>>>>>>>>>>>>

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 3 >>>>>>>>>>>>>>>>>
} else
{
print("Gate 3 >> ERROR: Model stopped running because of an error in step 3b: Error[4].")
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 3 >>>>>>>>>>>>>>>>>

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 2 >>>>>>>>>>>>>>>>>
} else
{
print("Gate 2 >> ERROR: Model stopped running because of an error in step 3a: Error[3].")
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 2 >>>>>>>>>>>>>>>>>

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 1 >>>>>>>>>>>>>>>>>
} else
{
print("Gate 1 >> ERROR: Model stopped running because of an error in step 3a: Error[2].")
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GATE 1 >>>>>>>>>>>>>>>>>

})# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> END SYSTEM TIME >>>>>>>>

#NOTE: Errors[1,6] do not have gates because a dataset with no errors may prodice results 
#within these checks. Results point out potential areas of trouble.

#This is an alarm function, it doesn't seem to work though in JGR (works in R program's GUI)
#Located in "utils" library.
alarm()

