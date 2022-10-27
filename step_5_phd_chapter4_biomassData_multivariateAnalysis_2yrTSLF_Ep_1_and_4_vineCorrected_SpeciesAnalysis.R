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
blr <- function(v1,v2,v3,v4,v5,v6,v7)
{
  cf <- 0.8
  mv <- vector(mode = "numeric")
  fv <- vector(mode = "numeric")
  le <- vector(mode = "numeric")
  lower <- min(floor(v4*v2)/v2)
  upper.1 <- max(ceiling(v4*v2)/v2)
  upper.2 <- upper.1 - 1/v2
  brk.pt <- (upper.1-lower)/5
  bins <- seq(lower, upper.2, brk.pt)
  
  for(i in 1:length(bins))
  {
    if(i == 1)
    {
      mv[i] <- max(v5[v4 >= bins[i] & v4 <= bins[i] + brk.pt])
      fv[i] <- min(v4[v4 >= bins[i] & v4 <= bins[i] + brk.pt & v5 == mv[i]])
      le[i] <- length(v5[v4 >= bins[i] & v4 <= bins[i] + brk.pt])
      mv[i] <- mv[i] + 0.0001
    } else
    {
      mv[i] <- max(v5[v4 > bins[i] & v4 <= bins[i] + brk.pt])
      fv[i] <- min(v4[v4 > bins[1] & v4 <= bins[i] + brk.pt & v5 == mv[i]])
      le[i] <- length(v5[v4 > bins[i] & v4 <= bins[i] + brk.pt])
      mv[i] <- mv[i] + 0.0001
    }
  }
  
  v42 <- fv+0.00001
  v52 <- mv  
  
  d <- data.frame(v42,v52)
  logmodel <- lm(v52~v1(v42),data=d)
  
  ### fake vector
  v4vec <- seq(0.001,8, length=101)
  logpred <- predict(logmodel, newdata=data.frame(v42=v4vec))
  
  #Plot with boundary layer regression
  plot(v4, v5, pch = 1, xlab = v7, ylab = "biomass (Mg/ha)", main = v6)
  points(fv, mv, pch = 16)
  lines(v4vec,logpred)
  par(cex = 0.9)
  #text(4.0, (max(y) - max(y)/16), paste("Intercept", round(logmodel$coefficients[1],4)), cex = cf)
  #text(5.8, (max(y) - max(y)/16), paste("Slope", round(logmodel$coefficients[2],4)), cex = cf)
  xcoord <- (upper.2 - lower)*v3 + lower
  ms <- summary(logmodel)
  text(xcoord, (max(v5) - max(v5)/22), paste("r-squared", round(ms$r.squared,2)), cex = cf)
  text(xcoord, (max(v5) - max(v5)/8), paste("P-value", round(ms$coefficients[8],3)), cex = cf)
  print(summary(logmodel))
  par(cex = 0.65)
}

###################################################################################################
############################################FUNCTION##############################################

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

#Calculate total understory biomass for each plot
tot_biomass <- apply(plotBiomass[4:length(plotBiomass[1,])], 1, sum)

#Calculate total understory biomass for each plot
site_mean <- summarize(X = tot_biomass, by = plotBiomass$siteName, 
                       mean, stat.name = "mean")
round(mean(site_mean$mean),2)
round(range(site_mean$mean),2)

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

#Show summary stats
sss <- round(stat.desc(spo),2)

#Mean biomass values for all sites for all species ordered from highest to lowest
sort(sss['mean',], decreasing = T)

#Average biomass for dwarf line oak
sp <- "LYFE"
sss['mean',sp]
sss['SE.mean',sp]

#Percent of total understory (living) biomass for QUMI2
round(mean((spo$QUMI2/site_mean$mean)*100),1)


#uv.plots displays histogram, box and whisker, cumulative distribution and normal q-q plots
#in one pane for each variable.
uv.plots(spo)#Many species are right skewed with a long right tail. To reduce the effect of these outliers use log-transformed data.

#What percent of values are zero, if it is over 50% you should consider changing the data
#to presence/absence
length(spo[spo == 0])/(length(spo[1,])*length(spo[,1]))
#64% zero values

#Look at correlation between species variables.
#generate table (you need to run correlation_matrix() function for this to work:
#https://www.r-bloggers.com/2020/07/create-a-publication-ready-correlation-matrix-with-significance-levels-in-r/
scm <- correlation_matrix(as.data.frame(spo), type = "pearson", show_significance = T, 
                         digits = 2, use = "lower", replace_diagonal = T)
chart.Correlation(spo, method = "pearson")#performanceAlanlytics

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

#ARST5
wg <- data.frame(site = rownames(spo), ARST5 = round(spo$ARST5,2), mfri = siteEnv3$mfri_20yr)
wg[order(wg$mfri),]

#ANGL10
bs <- data.frame(site = rownames(spo), ANGL10 = round(spo$ANGL10,2), mfri = siteEnv3$mfri_20yr)
bs[order(bs$mfri),]

#GESE
yj <- data.frame(site = rownames(spo), GESE = round(spo$GESE,2), gd = siteEnv3$Season_20yr)
yj[order(yj$gd),]


###################################################################################################
###################################################################################################
#STEP 13: BOUNDARY LAYER REGRESSION

#Old file reassignments
emat <- siteEnv3
ab2mat <- spo

#ENVIRONMENTAL VARIABLES

#Fire Rotation
#Environmental variable
ev <- emat$mfri_20yr
#Label for X-axis
envLabel <- "mFRI (Years)"
#Decimals for calculating break points for fixed width windows (1 is nearest whole number, 10 is nearest
#tenth and so on)
dec <- 1

#Season
#Environmental variable
ev <- emat$Season_20yr
#Label for X-axis
envLabel <- "Season"
#Decimals for calculating break points for fixed width windows (1 is nearest whole number, 10 is nearest
#tenth and so on)
dec <- 10

#STATIC VARIABLES
#Type of regression model to use
model <- log
#Multiplier for x coordinate of legend (0-1).
#0 = minimum x coordinate
#1 = maximum x coordinate
xcoord <- 0.8

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
speciesLabel <- "Quercus minima"

#Boundary layer regression
blr(model, dec, xcoord, ev, species, speciesLabel, envLabel)

#2
#Saw Palmetto
#Plot with sites
species <- ab2mat$SERE2
speciesLabel <- "Serenoa repens"

#Boundary layer regression
blr(model, dec, xcoord, ev, species, speciesLabel, envLabel)

#3
#Wiregrass
#Plot with sites
species <- ab2mat$ARST5
speciesLabel <- "Aristida stricta"

#Boundary layer regression
blr(model, dec, xcoord, ev, species, speciesLabel, envLabel)

#4
#Gallberry
#Plot with sites
species <- ab2mat$ILGL
speciesLabel <- "Ilex glabra"

#Boundary layer regression
blr(model, dec, xcoord, ev, species, speciesLabel, envLabel)

#5
#Lyonia
#Plot with sites
species <- ab2mat$LYFE
speciesLabel <- "Lyonia ferruginea"

#Boundary layer regression
blr(model, dec, xcoord, ev, species, speciesLabel, envLabel)

#6
#Darrow's Blueberry
#Plot with sites
species <- ab2mat$VADA
speciesLabel <- "Vaccinium darrowii"

#Boundary layer regression
blr(model, dec, xcoord, ev, species, speciesLabel, envLabel)

#7
#Little Bluestem
#Plot with sites
species <- ab2mat$ANGL10
speciesLabel <- "Andropogon glaucopsis"

#Boundary layer regression
blr(model, dec, xcoord, ev, species, speciesLabel, envLabel)

#8
#Huckleberry
#Plot with sites
species <- ab2mat$GADU
speciesLabel <- "Gaylussacia dumosa"

#Boundary layer regression
blr(model, dec, xcoord, ev, species, speciesLabel, envLabel)

#9
#Greenbriar
#Plot with sites
species <- ab2mat$SMILA2
speciesLabel <- "Smilax auriculata"

#Boundary layer regression
blr(model, dec, xcoord, ev, species, speciesLabel, envLabel)

#10
#Yellow Jessamine
#Plot with sites
species <- ab2mat$GESE
speciesLabel <- "Gelsemium sempervirens"

#Boundary layer regression
blr(model, dec, xcoord, ev, species, speciesLabel, envLabel)

###################################################################################################
###################################################################################################
#END



