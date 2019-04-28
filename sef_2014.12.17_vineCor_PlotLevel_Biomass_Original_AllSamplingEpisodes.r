

###################################################################################################
############################-----start-----########################################################
#Purpose of this script is to elimate inconcistencies with how vine species were categorized
#for all sampling epiosdes. Vines were inconsistently categorized by life form or species primarily 
#during the first field sampling season but also for multiple plots in episode 3 and one in episode
#4.

#This script creates a finished dataset for all sites that were sampled pre-fire, yr1, and yr2.

#To fix I calculated the fraction of biomass for each of the three vine species (yellow
#jessamine, greenbriar, and muscadine grape) in episode 4 when there is only on instance
#of vines not being categorized by species for each site. This shows the relative 
#proportion of vine species by biomass at each site. I then applied these coefficients to
#all instances where biomass was categorized as "vine". This occurred at several sites in
#episodes 1 and 3 and one site in episode 4. The calculated biomass amounts were added to
#any existing values for the vine species within a plot (although I could find no instances
#where this occurred. That is, no plots had both vine and one or more of the vine species).

#This technique has limitations, it assumes that the proprtion of vines does not change over time.
#It means you can only look at variance within each site but you can't look at plot-level 
#correlations among species that includes vines because the plot-level data doesn't reflect actual 
#conditions.


###################################################################################################
###################################################################################################
#STEP #1: ADMINISTRATIVE TASKS

#Reset functions
rm(list=ls())
dev.off()

library(Hmisc)#summarize()


###################################################################################################
###################################################################################################
#STEP #2: Create objects that list vine types for each type of dataset: original and by genera

#Create an object that lists column headings in the vine category
vinesG <- c("vine", "Vitus", "Gelsemium", "Smilax")
vinesCN <- c("vine", "muscadine.grape", "yellow.jessamine", "greenbriar")

###################################################################################################
###################################################################################################
#STEP #3: Open plot-level data (PLOT LEVEL)

#Open episode 1 biomass data matrix for all sites sampled in episode 1
PbiomassA <- read.table(
  "C:/usfs_sef_data_output/sef_Ecology_BiomassPlotMatrix_Ep1_Original2014-02-26_15.06.16.csv", 
  header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE,
  stringsAsFactors = F)

#Remove sites E807B, E503C, E507B, and A314 from biomass1a b/c they are not represented in latter 
#episodes.
Psn <- unique(PbiomassA$SiteName)
Psn18 <- Psn[-c(8, 20,21,22)]

PbiomassBtt <- PbiomassA[PbiomassA$SiteName %in% Psn18,]

#Convert list form to matrix with genera as col headings for biomassBtt.
PbiomassBt <- matrix(PbiomassBtt[,8],length(unique(PbiomassBtt[,2])), length(unique(PbiomassBtt[,4])), byrow = T, 
                     dimnames = list(unique(PbiomassBtt[,2]),unique(PbiomassBtt[,4])))

PbiomassB <- data.frame(siteName = I(mapply(function(y) rep(y,20), unique(PbiomassBtt$SiteName))), 
                        PlotNo = unique(PbiomassBtt$PlotNo), PbiomassBt)

#Open episode 4 biomass data matrix for all sites sampled in episode 4
Pbiomass4att <- read.table(
  "C:/usfs_sef_data_output/sef_Ecology_BiomassPlotMatrix_Ep4_Original2014-02-26_15.00.01.csv", 
  header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE,
  stringsAsFactors = F)

#Convert list form to matrix with genera as col headings for biomassCtt.
Pbiomass4at <- matrix(Pbiomass4att[,8],length(unique(Pbiomass4att[,2])), length(unique(Pbiomass4att[,4])), byrow = T, 
                      dimnames = list(unique(Pbiomass4att[,2]),unique(Pbiomass4att[,4])))

Pbiomass4a <- data.frame(siteName = I(mapply(function(y) rep(y,20), unique(Pbiomass4att$SiteName))), 
                         PlotNo = unique(Pbiomass4att$PlotNo), Pbiomass4at)

#Open episode 1 biomass data matrix for all sites sampled in episode 3
Pbiomass3att <- read.table(
  "C:/usfs_sef_data_output/sef_Ecology_BiomassPlotMatrix_Ep3_Original2014-02-26_15.02.49.csv", 
  header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE,
  stringsAsFactors = F)

#Remove site E807B from Pbiomass3att b/c it is not represented in episode 4 (yr2).
Psn_2 <- unique(Pbiomass3att$SiteName)
Psn18_2 <- Psn_2[-8]

Pbiomass3btt <- Pbiomass3att[Pbiomass3att$SiteName %in% Psn18_2,]

#Convert list form to matrix with genera as col headings for biomass3btt.
Pbiomass3at <- matrix(Pbiomass3btt[,8],length(unique(Pbiomass3btt[,2])), length(unique(Pbiomass3btt[,4])), byrow = T, 
                      dimnames = list(unique(Pbiomass3btt[,2]),unique(Pbiomass3btt[,4])))

Pbiomass3a <- data.frame(siteName = I(mapply(function(y) rep(y,20), unique(Pbiomass3btt$SiteName))), 
                         PlotNo = unique(Pbiomass3btt$PlotNo), Pbiomass3at)


###################################################################################################
###################################################################################################
#STEP 4: Look at the difference between sampling episodes 1, 3, and 4 in terms of the percent of vine
#species not classified by species.

ep1_vine_plot <- PbiomassB$vine
ep1_vine_site <- tapply(ep1_vine_plot, PbiomassB$siteName, mean)


ep1_all_plot <- apply(PbiomassB[,colnames(PbiomassB) %in% vinesCN[-1]],1,sum)
ep1_all_site <- tapply(ep1_all_plot, PbiomassB$siteName, mean)

e1 <- round((ep1_vine_site/(ep1_vine_site + ep1_all_site))*100,1)

ep3_vine_plot <- Pbiomass3a$vine
ep3_vine_site <- tapply(ep3_vine_plot, Pbiomass3a$siteName, mean)


ep3_all_plot <- apply(Pbiomass3a[,colnames(Pbiomass3a) %in% vinesCN[-1]],1,sum)
ep3_all_site <- tapply(ep3_all_plot, Pbiomass3a$siteName, mean)

e3 <- round((ep3_vine_site/(ep3_vine_site + ep3_all_site))*100,1)

ep4_vine_plot <- Pbiomass4a$vine
ep4_vine_site <- tapply(ep4_vine_plot, Pbiomass4a$siteName, mean)


ep4_all_plot <- apply(Pbiomass4a[,colnames(Pbiomass4a) %in% vinesCN[-1]],1,sum)
ep4_all_site <- tapply(ep4_all_plot, Pbiomass4a$siteName, mean)

e4 <- round((ep4_vine_site/(ep4_vine_site + ep4_all_site))*100,1)

#check to make sure fractional coeffiecients can be applied to all sites and sampling episodes.
#That is there should be a value of 0 or greater in the fourth sampling episode if there
#is a positive value in sampling episodes 1 or 3.
cbind(e1,e3,e4)

#Both sampling episodes 1 and 3 have 2-3 sites 40% or more of vine biomass not classified into
#species and several sites with smaller percentages of vine biomass unclassified. Sampling
#episode 4 only has one site with a small percentage of biomass unclassified.

#I feel that the best way to deal with this is to determine the fraction of vine biomass
#in each species for sampling episode 4 and then apply this fraction to the unlcassified
#vine biomass in all sampling episodes.

###################################################################################################
###################################################################################################
#STEP #5: Calculate fractional coefficients from episode 4 for each of the three vine species.

#List of common names for vine species.
vineSpecies <- vinesCN[-1]

t1 <- Pbiomass4a[,colnames(Pbiomass4a) %in% vineSpecies]

#Calculate episode 4 site means for each vine species
vineSiteMean_ep4 <- aggregate(.~Pbiomass4a$siteName, data = t1, mean)

#Calculate episode 4 sum of biomass for vine species
vineTotal_ep4 <- apply(t1, 1, sum)
totalSiteMean_ep4 <- tapply(vineTotal_ep4, Pbiomass4a$siteName, mean)

#Calculate fractional coefficients for each vine species and site in episode 4.
fvc <- vineSiteMean_ep4[,2:4]/totalSiteMean_ep4
fvc[fvc == "NaN"] <- 0#for any sites with no vines convert NaN to zero.

#Same as the check above except look at the actual fractional coefficients.
cbind(e1, e3, e4, fvc)
#There no issues. There are fractional coefficients that total 1 for all instances where
#vines were not classified by species.


###################################################################################################
###################################################################################################
#STEP #6: Apply fractional coefficients to biomass entries in the vine category for all episodes.

#Create copies of datasets for episodes 1,3, and 4.
b1 <- PbiomassB
b3 <- Pbiomass3a
b4 <- Pbiomass4a

#Apply the fractional coefficients to each sampling episode
b1_a <- mapply(function(y)
  {
  mapply(function(x)
    {b1[x,colnames(b1) == y] + 
    b1$vine[x] * fvc[vineSiteMean_ep4[,1] == b1$siteName[x],colnames(fvc) == y]
  }, 1:length(b1$siteName))  
  }, vineSpecies)

b3_a <- mapply(function(y)
{
  mapply(function(x)
  {b3[x,colnames(b3) == y] + 
     b3$vine[x] * fvc[vineSiteMean_ep4[,1] == b3$siteName[x],colnames(fvc) == y]
  }, 1:length(b3$siteName))  
}, vineSpecies)

b4_a <- mapply(function(y)
{
  mapply(function(x)
  {b4[x,colnames(b4) == y] + 
     b4$vine[x] * fvc[vineSiteMean_ep4[,1] == b4$siteName[x],colnames(fvc) == y]
  }, 1:length(b4$siteName))  
}, vineSpecies)


#Remove the vine column
b1_b <- b1[,-6]
b3_b <- b3[,-6]
b4_b <- b4[,-6]

#Add the new values for vine species.
b1_b[,15] <- b1_a[,1]#sampling episode 1; muscadine.grape
b3_b[,15] <- b3_a[,1]#sampling episode 3; muscadine.grape
b4_b[,15] <- b4_a[,1]#sampling episode 4; muscadine.grape

b1_b[,16] <- b1_a[,2]#sampling episode 1; yellow.jessamine
b3_b[,16] <- b3_a[,2]#sampling episode 3; yellow.jessamine
b4_b[,16] <- b4_a[,2]#sampling episode 4; yellow.jessamine

b1_b[,17] <- b1_a[,3]#sampling episode 1; greenbriar
b3_b[,17] <- b3_a[,3]#sampling episode 3; greenbriar
b4_b[,17] <- b4_a[,3]#sampling episode 4; greenbriar

###################################################################################################
###################################################################################################
#STEP #7: Clean up dataset

#Check to make sure all rows are in the same order for each samping episode.
cbind(b1_b$PlotNo,b3_b$PlotNo,b4_b$PlotNo)
#Okay, they all look good.

###################################################################################################
###################################################################################################
#STEP #8: Save dataset.

#set a date used to name the file.
dt <- Sys.Date()
tm <- format(Sys.time(), format = "%H.%M.%S", 
             tz = "", usetz = FALSE)

#Save the plot-level data.
#Sampling episode 1
write.table(b1_b, 
            file = paste("c:\\usfs_sef_data_output\\sef_Ecology_BiomassPlotMatrix_Ep1_Original_", 
                         dt,"_",tm,"_vineCor.csv",sep = ""), append = F, quote = T, sep = ",", 
            eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T, 
            qmethod = c("escape", "double"))

#Sampling episode 3
write.table(b3_b, 
            file = paste("c:\\usfs_sef_data_output\\sef_Ecology_BiomassPlotMatrix_Ep3_Original_", 
                         dt,"_",tm,"_vineCor.csv",sep = ""), append = F, quote = T, sep = ",", 
            eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T, 
            qmethod = c("escape", "double"))

#Sampling episode 4
write.table(b4_b, 
            file = paste("c:\\usfs_sef_data_output\\sef_Ecology_BiomassPlotMatrix_Ep4_Original_", 
                         dt,"_",tm,"_vineCor.csv",sep = ""), append = F, quote = T, sep = ",", 
            eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T, 
            qmethod = c("escape", "double"))


