

###################################################################################################
############################-----start-----########################################################
#Convert multivariate plot cover and biomass data from original plant categories to functional 
#groups.

###################################################################################################
###################################################################################################
#
#
#
#October 2022 NOTES
#This script was re-habed/modified on 13-October-2022.
#Unchanged script is first version on GitHub repo (push to remote on 13-Oct-2022 0945 EST.
#Original script name: sef_2014.10.20_VineCor_OutlierX_FuncGroup
#Original script location: C:\Users\james\Box\01. james.cronan Workspace\Research\2009_01_SEF\r_scripts

#Script was modified to run on new input files and to produce relassification by genus in addition
#to plant functional groups.
#I also removed reclass of cover data because I don't think I will be using it for this analysis.
#
#
#
###################################################################################################
###################################################################################################

###################################################################################################
###################################################################################################
#STEP #1: ADMINISTRATIVE TASKS

#Reset functions
rm(list=ls())
dev.off()

#Libraries
#none

###################################################################################################
###################################################################################################
#
#
#                                      PLANT FUNCTIONAL GROUPS
#
#
###################################################################################################
###################################################################################################

###################################################################################################
###################################################################################################
#STEP #2: OPEN, ADJUST, AND REVIEW DATA FILES

###################################################################################################
#STEP #2a: Open biomass data


#Set working directory
setwd(paste("C:/Users/james/Box/01. james.cronan Workspace/Research/UW_PHD/Dissertation/4_Chapter_4", 
"/Data/Understory_Vegetation_FlatFiles/stage_4_aggregate/inputs", sep = ""))


#Open plot biomass data with field measurements standardized to a 2-year rough. This means
# 2-yr post-fire (2011-2012; sampling episode 4) measurements for sites that were burned in 2009-2010 and pre-fire
#(2009-2010; sampling episode 1) measurements for sites that were not burned.
plotBiomass <- read.table(
  "20221011_Biomass_OriginalVegClasses_2yrRough_Episodes_1_and_4.csv", 
  header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE,
  stringsAsFactors = F)

###################################################################################################
#STEP #2c: Open fuctional group crosswalk for biomass data

#Open species crosswalk file to aggregate original loading data:
biomassCross <- read.table("20221012_VegClass_lookupTable_simple.csv", 
header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)

###################################################################################################
###################################################################################################
#STEP #3: CONVERT ORGINAL DATA TO PLANT FUNCTIONAL GROUPS

###################################################################################################
#STEP #3a: Biomass data

#Create a vector of unique functional groups
bpfg <- sort(unique(biomassCross$funcGroup))

#Create a data.frame to accept the functional group values.
bfg <- data.frame(matrix(rep(0,length(plotBiomass[,1]) * length(bpfg)), 
                        nrow = length(plotBiomass[,1]), 
                        ncol = length(bpfg)))


for(i in 1:length(bpfg))
{
  bfg[,i] <- apply(plotBiomass[colnames(plotBiomass) %in% 
                    biomassCross$original_name[biomassCross$funcGroup == bpfg[i]]],1,sum)
}


colnames(bfg) <- bpfg

biomassFuncGroup <- data.frame(plotBiomass[,1:3], bfg)

###################################################################################################
###################################################################################################
#STEP #4: Save relevant files

#set a date used to name the file.
dt <- Sys.Date()
tm <- format(Sys.time(), format = "%H.%M.%S", 
tz = "", usetz = FALSE)

setwd(paste("C:/Users/james/Box/01. james.cronan Workspace/Research/UW_PHD/Dissertation/4_Chapter_4", 
"/Data/Understory_Vegetation_FlatFiles/stage_4_aggregate/outputs", sep = ""))

#Save the plot-level biomass data.
write.table(biomassFuncGroup, file = paste("phd_chapter4_biomassPlotMatrix_2yr_Ep_1_4_FuncGroup_", 
                                           dt,"_",tm,".csv", sep = ""), append = F, quote = T, 
            sep = ",", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T, 
            qmethod = c("escape", "double"))#

###################################################################################################
###################################################################################################
#
#
#                                               GENUS
#
#
###################################################################################################
###################################################################################################

###################################################################################################
###################################################################################################
#STEP #5: CONVERT ORGINAL DATA TO GENUS 

###################################################################################################
#STEP #5a: Biomass data

#Create a vector of unique genus
bpge <- sort(unique(biomassCross$genus))

#Create a data.frame to accept the functional group values.
bge <- data.frame(matrix(rep(0,length(plotBiomass[,1]) * length(bpge)), 
                         nrow = length(plotBiomass[,1]), 
                         ncol = length(bpge)))


for(i in 1:length(bpge))
{
  bge[,i] <- apply(plotBiomass[colnames(plotBiomass) %in% 
                                 biomassCross$original_name[biomassCross$genus == bpge[i]]],1,sum)
}


colnames(bge) <- bpge

biomassGenus <- data.frame(plotBiomass[,1:3], bge)

###################################################################################################
###################################################################################################
#STEP #6: Save relevant files

#set a date used to name the file.
dt <- Sys.Date()
tm <- format(Sys.time(), format = "%H.%M.%S", 
             tz = "", usetz = FALSE)

setwd(paste("C:/Users/james/Box/01. james.cronan Workspace/Research/UW_PHD/Dissertation/4_Chapter_4", 
            "/Data/Understory_Vegetation_FlatFiles/stage_4_aggregate/outputs", sep = ""))

#Save the plot-level biomass data.
write.table(biomassGenus, file = paste("phd_chapter4_biomassPlotMatrix_2yr_Ep_1_4_Genus_", 
                                           dt,"_",tm,".csv", sep = ""), append = F, quote = T, 
            sep = ",", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T, 
            qmethod = c("escape", "double"))#




