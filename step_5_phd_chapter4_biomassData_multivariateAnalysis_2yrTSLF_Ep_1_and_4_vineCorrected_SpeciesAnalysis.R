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
#STEP #4: OPEN AND ADJUST BIOMASS DATA AT SPECIES LEVEL

#########################################################
#4A: Open biomass data

#Plot-level biomass data. This data has been modified so that vine species are consistent
#across all sites.
#Outlier sites have not been removed.
plotBiomass <- read.table(paste("Understory_Vegetation_FlatFiles/stage_4_aggregate/outputs/", 
                                "phd_chapter4_biomassPlotMatrix_2yr_Ep_1_4_Species_2022-10-19_13.26.07.csv",
                                sep = ""), header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE,
                          stringsAsFactors = F)

#########################################################
#4B: calculate site means for untransformed plot-level data.
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

###################################################################################################
###################################################################################################
#STEP #5: DATA FORMATTING/SCREENING

#hand off object name.
spo <- siteBiomass3

#Convert data.frame to a matrix (needed for uv.plots)
so2 <- data.matrix(frame = spo, rownames.force = NA)

#uv.plots displays histogram, box and whisker, cumulative distribution and normal q-q plots
#in one pane for each variable.
uv.plots(so2)#Many species are right skewed with a long right tail. To reduce the effect of these outliers use log-transformed data.

#What percent of values are zero, if it is over 50% you should consider changing the data
#to presence/absence
length(so2[so2 == 0])/(length(so2[1,])*length(so2[,1]))
#64% zero values

#Look at correlation between species variables.
#generate table (you need to run correlation_matrix() function for this to work:
#https://www.r-bloggers.com/2020/07/create-a-publication-ready-correlation-matrix-with-significance-levels-in-r/
scm <- correlation_matrix(as.data.frame(so2), type = "pearson", show_significance = T, 
                         digits = 2, use = "lower", replace_diagonal = T)
chart.Correlation(so2, method = "pearson")#performanceAlanlytics

#########################################################
#6c: Summary stats for environmental data

#Show how environmental variables are correlated.
chart.Correlation(siteEnv3)
ecm <- correlation_matrix(as.data.frame(siteEnv3), type = "pearson", show_significance = T, 
                          digits = 2, use = "lower", replace_diagonal = T)

#Show dominant species

#Remove non-living categories (dead woody) from the biomass data
last_col <- length(spo[1,])
nextLast_col <- length(spo[1,]) - 1

#List of dominant species
prim <- mapply(function(x) {colnames(spo)[order(spo[x,])][last_col]}, 
               1:length(spo[,1]))
seco <- mapply(function(x) {colnames(spo)[order(spo[x,])][nextLast_col]}, 
               1:length(spo[,1]))

dTable <- data.frame(Sites = I(rownames(spo)), Primary = I(prim), Secondary = I(seco))

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

#Show chart of biomass for the 10 most common species

#Calculate untransformed means for site-level data 
#Caluclate mean biomass for each species
mean_biomass <- apply(spo, 2, mean)

#Select the top 10 species
#Top 14 >> includes non-species groups
top10a <- sort(mean_biomass, decreasing = T)[1:14]

#Top 10: manually remove non species
top10b <- top10a[-c(6,8,9,10)]
t10_table <- spo[,colnames(spo) %in% names(top10b)]

#Convert data frame into format that can be fed into boxplot() function.
#Create vectors to hold biomass (numeric) and species data (factor)
t10_biomass <- vector(mode = 'numeric')#biomass data
t10_spp <- vector()#species data

#Create two sets of vectors from existing table.
for(i in 1:length(t10_table[1,]))
{
  t10_biomass <- c(t10_biomass, t10_table[[i]])
  t10_spp <- c(t10_spp, rep(colnames(t10_table[i]), length(t10_table[,1])))
}

t10_spp <- as.factor(t10_spp)#convert species from character to factors.
t10_biomass <- as.numeric(t10_biomass)#convert biomass from characters to numeric.

#Combine vectors into a new table that can be fed into boxplot() function.
t10_vec <- data.frame(species = t10_spp, biomass = t10_biomass)

new_order <- with(t10_vec, reorder(species, biomass, mean , na.rm=T))

#Generate box plot
boxplot(t10_vec$biomass ~ new_order)
###################################################################################################
###################################################################################################
#STEP 13: BOUNDARY LAYER REGRESSION FOR FIRE ROTATION

#Old file reassignments
emat <- siteEnv3
ab2mat <- spo

ev <- emat$mfri_20yr
etext <- "mFRI (20 years)"

###################################################################################################
###################################################################################################
#PANEL PLOT BY ORIGINAL DATA __1
dev.off()
set.panel(2,5)
cf <- 0.5

#1
#Runner Oak
#Plot with sites
species <- ab2mat$QUMI2
Species <- "Quercus minima"

#Boundary layer regression
blr(log, ev, species, Species)

#2
#Saw Palmetto
#Plot with sites
species <- ab2mat$SERE2
Species <- "Serenoa repens"

#Boundary layer regression
blr(log, ev, species, Species)

#3
#Wiregrass
#Plot with sites
species <- ab2mat$ARST5
Species <- "Aristida stricta"

#Boundary layer regression
blr(log, ev, species, Species)

#4
#Gallberry
#Plot with sites
species <- ab2mat$ILGL
Species <- "Ilex glabra"

#Boundary layer regression
blr(log, ev, species, Species)

#5
#Lyonia
#Plot with sites
species <- ab2mat$LYFE
Species <- "Lyonia ferruginea"

#Boundary layer regression
blr(log, ev, species, Species)

#6
#Darrow's Blueberry
#Plot with sites
species <- ab2mat$VADA
Species <- "Vaccinium darrowii"

#Boundary layer regression
blr(log, ev, species, Species)

#7
#Little Bluestem
#Plot with sites
species <- ab2mat$ANGL10
Species <- "Andropogon glaucopsis"

#Boundary layer regression
blr(log, ev, species, Species)

#8
#Huckleberry
#Plot with sites
species <- ab2mat$GADU
Species <- "Gaylussacia dumosa"

#Boundary layer regression
blr(log, ev, species, Species)

#9
#Greenbriar
#Plot with sites
species <- ab2mat$SMILA2
Species <- "Smilax auriculata"

#Boundary layer regression
blr(log, ev, species, Species)

#10
#Yellow Jessamine
#Plot with sites
species <- ab2mat$GESE
Species <- "Gelsemium sempervirens"

#Boundary layer regression
blr(log, ev, species, Species)

###################################################################################################
###################################################################################################
#END