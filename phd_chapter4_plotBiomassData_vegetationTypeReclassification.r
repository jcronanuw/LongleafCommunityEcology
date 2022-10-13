

###################################################################################################
############################-----start-----########################################################
#Convert multivariate plot cover and biomass data from original plant categories to functional 
#groups.

#This script was re-habed/modified on 13-October-2022.
#Unchanged script is first version on GitHub repo (push to remote on 13-Oct-2022 0945 EST.
#Original script name: sef_2014.10.20_VineCor_OutlierX_FuncGroup
#Original script location: C:\Users\james\Box\01. james.cronan Workspace\Research\2009_01_SEF\r_scripts

###################################################################################################
###################################################################################################
#STEP #1: ADMINISTRATIVE TASKS

#Reset functions
rm(list=ls())
dev.off()

#Libraries




###################################################################################################
###################################################################################################
#STEP #2: OPEN, ADJUST, AND REVIEW DATA FILES

###################################################################################################
#STEP #2a: Open biomass data

#Plot-level biomass data. This data has been modified so that vine species are consistent
#across all sites and outlier plots located in uncharacteristics areas of sites (wetlands or meadows) 
#have been removed.
plotBiomass <- read.table(
  "C:/usfs_sef_data_output/sef_Ecology_BiomassPlotMatrix_Ep1_OriginalOulierX_2014-05-05_15.00.18.csv", 
  header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE,
  stringsAsFactors = F)

###################################################################################################
#STEP #2b: Open cover data
plotCover <- read.table(
  "C:/usfs_sef_data_output/sef_Ecology_CoverPlotMatrix_Ep1_OriginalOulierX_2014-10-15_10.32.12.csv", 
  header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE,
  stringsAsFactors = F)

###################################################################################################
#STEP #2c: Open fuctional group crosswalk for biomass data

#Open species crosswalk file to aggregate origina loading data:
biomassCross <- read.table("C:/usfs_sef_data_output/sef_crosswalkBiomass_Original_to_FuncGroup.csv", 
header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
biomassCross <- cleanup.import(biomassCross)

###################################################################################################
#STEP #2d: Open fuctional group crosswalk for cover data

#Open species crosswalk file to aggregate original cover data:
coverCross <- read.table("C:/usfs_sef_data_output/sef_crosswalkCover_Original_to_FuncGroup.csv", 
                           header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
coverCross <- cleanup.import(coverCross)


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
                    biomassCross$commonName[biomassCross$funcGroup == bpfg[i]]],1,sum)
}


colnames(bfg) <- bpfg

biomassFuncGroup <- data.frame(plotBiomass[,1:2], bfg)

###################################################################################################
#STEP #3b: Cover data

#Create a vector of unique functional groups
cpfg <- sort(unique(coverCross$funcGroup))

#Create a data.frame to accept the functional group values.
cfg <- data.frame(matrix(rep(0,length(plotCover[,1]) * length(cpfg)), 
                         nrow = length(plotCover[,1]), 
                         ncol = length(cpfg)))


for(i in 1:length(cpfg))
{
  cfg[,i] <- apply(plotCover[colnames(plotCover) %in% 
                               coverCross$commonName[coverCross$funcGroup == cpfg[i]]],1,sum)
}


colnames(cfg) <- cpfg

coverFuncGroup <- data.frame(plotCover[,1:2], cfg)

###################################################################################################
###################################################################################################
#STEP #4: Save relevant files

#set a date used to name the file.
dt <- Sys.Date()
tm <- format(Sys.time(), format = "%H.%M.%S", 
tz = "", usetz = FALSE)

#Save the plot-level biomass data.
write.table(biomassFuncGroup, file = paste("c:\\usfs_sef_data_output\\sef_Ecology_BiomassPlotMatrix_Ep1_FuncGroupOulierX_", 
                                           dt,"_",tm,".csv", sep = ""), append = F, quote = T, 
            sep = ",", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T, 
            qmethod = c("escape", "double"))#

#Save the plot-level cover data.
write.table(coverFuncGroup, file = paste("c:\\usfs_sef_data_output\\sef_Ecology_CoverPlotMatrix_Ep1_FuncGroupOulierX_", 
                                         dt,"_",tm,".csv",sep = ""), append = F, quote = T, 
            sep = ",", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T, 
            qmethod = c("escape", "double"))#


#This is an alarm function, it doesn't seem to work though.
#alarm
#cat("\b")
