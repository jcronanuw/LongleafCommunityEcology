###################################################################################################
############################-----start-----########################################################
#                             ANALYSIS -- PLANT FUNCTIONAL GROUPS
#The purpose of this script is to conduct a constrained ordination of genera-level understory
#composition data using fire regime characteristics and forest structure (overstory and forest
#floor) using constraints.

#Run scripts prior to running this analysis
#correlation_matrix
#biostats

#Which computer are you using?
#Forest Service>>>>>>> FS
#Personal >>>>>>>>>>>> JC
computer <- "JC"

###################################################################################################
###################################################################################################
#STEP #1: ADMINISTRATIVE TASKS

#Reset functions
#rm(list=ls())
#do not run this function b/c it will erase biostats and correlation_matrix scripts
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
siteEnv <- read.table(paste("Fire_History/", "phd_chapter4_environmentalData.csv", sep = ""), 
                      header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE, stringsAsFactors = F, fill = T)

#########################################################
#2b: Re-assign column 1 (site names) to  a row name.
#Also remove column 2 (site numbers)
siteEnv2 <- siteEnv[,3:(length(siteEnv[1,]))]
rownames(siteEnv2) <- siteEnv[,1]

#########################################################
#2c: Environmental data, remove data you will not analyze
#siteEnv3 <- siteEnv2[,-c(11,12)]
#Column 11 is soil drainage. This data is from USDA Soil Survey and not accurate
#Column 12 is region. Only two regions, at this point I don't believe this adds 
#substance to the analysis because there are no physical properties for each region
#that can confound how fire characteristics drive understory plant composition
siteEnv3 <- siteEnv2[,-c(11)]
#maintain region column for restricted permutation tests

###################################################################################################
###################################################################################################
#STEP #3: OPEN AND ADJUST BIOMASS DATA

#########################################################
#3a: Open biomass data

#Plot-level biomass data. This data has been modified so that vine species are consistent
#across all sites.
#Outlier sites have not been removed.
plotBiomass <- read.table(paste("Understory_Vegetation_FlatFiles/stage_4_aggregate/outputs/", 
                                "phd_chapter4_biomassPlotMatrix_2yr_Ep_1_4_FuncGroup_2022-10-24_16.00.38.csv",
                                sep = ""), header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE,
                          stringsAsFactors = F)

#########################################################
#3b: calculate site means for untransformed plot-level data.

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
#STEP #4: DATASET ADJUSTMENTS

###################################################################################################
#Standardize row order of data sets (i.e. order site names (row names) of biological data 
#according to order in environmental data set) and remove columns with no values > 0. 

#Reorder biological data (biomass) so sites are in same order as environmental data
#Step not necessary for cover data... sites are already in same order as environmental data.
ro <- data.frame(siteNo = 1:length(rownames(siteEnv3)), siteName = I(rownames(siteEnv3)))
ro2 <- ro[order(ro[,2]),]

#Untranformed biomass data
siteBiomass4 <- data.frame(so = ro2[,1], siteBiomass3)
siteBiomass5 <- siteBiomass4[order(siteBiomass4$so),]
siteBiomass6 <- siteBiomass5[,!colnames(siteBiomass5) %in% c("so")]

#Remove taxa with no occurences, this is necessary before you can use foa.plots below
fugr <- drop.var(siteBiomass6, min.fo = 1)#biostats

#Show taxa that were dropped:

dropped <- match(colnames(siteBiomass6), colnames(fugr))
#cbind(colnames(siteBiomass6), dropped)
#length(dropped)
#length(dropped) - length(dropped[is.na(dropped) == T])
#No drops.

###################################################################################################
###################################################################################################
#STEP #5: DATA SCREENING

#########################################################
#5a: Summary stats for biomass data

#Summary stats
pfg_ss <- round(stat.desc(fugr),2)

#sort(pfg_ss[9,], decreasing = T)

#Look at how data is distributed

#foa.plots(fugr)

#Drop rare genus < 10% of sites.
fco <- drop.var(fugr, min.po=10)

#Check to see how many genus were dropped for untransformed biomass data

#length(fugr[1,])#
#length(fco[1,])#

#No drops

#Double check number of plant functional groups (cols) and sites (rows)

#str(fco)

#21 obs and 7 variables

#Check summary stats for each plant functional group.
#uv.plots displays histogram, box and whisker, cumulative distribution and normal q-q plots
#in one pane for each variable.

#uv.plots(fco)

#About half of functional groups  are right skewed with a long right tail. 
#Conduct a log tranformation

#What percent of values are zero? If it is over 50% you should consider changing the data
#to presence/absence

#length(fco[fco == 0])/(length(fco[1,])*length(fco[,1]))

#~8% zero values >>> no need to convert data to presence/absence

#Look at correlation between genus variables.
#generate table (you need to run correlation_matrix() function for this to work:
#https://www.r-bloggers.com/2020/07/create-a-publication-ready-correlation-matrix-with-significance-levels-in-r/

scm.fco <- correlation_matrix(as.data.frame(fco), type = "pearson", show_significance = T, 
                         digits = 2, use = "lower", replace_diagonal = T)
#chart.Correlation(fco, method = "pearson")#performanceAlanlytics

#Understory trees and vines are highly correlated (0.65). All other are below 0.55.

#log transform genera data.
fco_log <- decostand(fco, method = "log")

#Check distributions

#uv.plots(fco_log)

#Improved for all.

#Hellinger transformation on genera data
fco_hel <- decostand(fco, method = "hellinger")

#Check distributions

#uv.plots(fco_hel)

#Improvement over log transformation distributions, they are more 'normal'.

#########################################################
#6c: Summary stats for environmental data

#Show how environmental variables are correlated.

scm.env <- correlation_matrix(as.data.frame(siteEnv3), type = "pearson", show_significance = T, 
                          digits = 2, use = "lower", replace_diagonal = T)

#chart.Correlation(siteEnv3)

#Use 0.55 as cut-off for removing variables with high colinearity 
#Sort of arbitrary, but chose cutoff described in Adam et al. 2013 t.
#Variables with correlation >= 0.55
#1) 10-yr and 20-yr mFRI
#drop 10-yr
#2) 10-yr and 20-yr season of burn ratio
#drop 10-yr
#3)Litter and 10-yr season of burn ratio
#already dropping 10-yr season of burn ratio
#4)Fine DWD loading and region
#ignore relationship with region
#region is a nuisance variable and I don't want to exclude possible
#effects of loading even though they will likely be washed out by
#region as a conditional constraint.

#Remove env. variables selected above

siteEnv4 <- siteEnv3[,-c(6,9)]#dropping the 10-yr values for mean fire return
#interval and ratio of growing:dormant season burns.

#Create a new vector with two regions (east and west) instead of 3 (SMNWR, APNF, and EAFB).

region <- as.factor(ifelse(siteEnv3$Region == 3, 2, siteEnv3$Region))

#Check distributions of environmental data
#Remove categories

check_env <- siteEnv4[,-c(9,10)]
#uv.plots(check_env)

#Right skewed
#1) Coarse DWD
#2) Std Dev on 20-year fire regime

#Log transform right-skewed variables
cdwd_log <- decostand(siteEnv4$CoarseWD, method = "log")
stde_log <- decostand(siteEnv4$StDevFRI_20YR, method = "log")

#Combine fine DWD with litter.
ff <- siteEnv4$Litter + siteEnv4$FineWD

#Create new environmental variables data frame
env <- data.frame(canopy = siteEnv4$Canopy,
                     ff = ff,
                     coarseWD = cdwd_log,
                     mfri = siteEnv4$mfri_20yr,
                     sd_fri = stde_log,
                     gd_ratio = siteEnv4$Season_20yr,
                     rxfire = siteEnv4$RxFireMgmt,
                     region = as.factor(region))
rownames(env) <- rownames(siteEnv4)

#Recheck distributions for environmental variables.

#uv.plots(env)#looks good

#multivariate outliers
#this function if from biostats and I'm not sure how it selects outliers so I am going
#to maintain these sites for now.
#Although clearly S330 stands out as an unusual site for both resp. and expl. matrices.
#E103B_S3 is a bit of an odd site. Half the plots were in a plantation with closed canopy
#and the other half were in a natural flatwoods stand so there is probably an odd mix of species.
#Interesting because I don't remember it being that unusual, in fact it is pretty typical
#for the flatwoods sites we sampled.
#Worth noting that standard deviations are pretty close to 2 for both these sites.
#That number is a bit arbitrary. If it was 2.5 you wouldn't have any outliers

#Environmental data
#mv.outliers(env[1:7], method = "euclidean", sd.limit =2.0)

#biomass
#mv.outliers(fco, method = "bray", sd.limit =2.0)

#Drop S330 because it is an outier within the environmental dataset. This will also create a balanced
#sample design for regional blocks and allow me to randomize within regions for the RDA with region
#as a conditional term.

env <- env[rownames(env) != "S330",]
fco_log <- fco_log[rownames(fco_log) != "S330",]
fco <- fco[rownames(fco) != "S330",]

###################################################################################################
###################################################################################################
#STEP #6: PRELIMINARY ANALYSIS

#########################################################
#7a: Determine plant functional groups with the highest and second highest biomass at each site
#and then display occurrence of dominant and co-dominant plant functional groups with a histogram.

#Create a histogram of primary and secondary dominance across sites by plant functional groups.

#Identify last and second to last columns to calculate which plant functional groups had the highest and
#second highest biomass values at each site
last_col <- length(fco[1,])
nextLast_col <- length(fco[1,]) - 1

#List of dominant plant functonal group
prim <- mapply(function(x) {colnames(fco)[order(fco[x,])][last_col]}, 
               1:length(fco[,1]))
seco <- mapply(function(x) {colnames(fco)[order(fco[x,])][nextLast_col]}, 
               1:length(fco[,1]))

dTable <- data.frame(Sites = I(rownames(fco)), Primary = I(prim), Secondary = I(seco))

udom <- unique(c(unique(dTable[,2]),unique(dTable[,3])))

hprim <- mapply(function(x) {length(dTable[,2][dTable[,2] == udom[x]])}, 1:length(udom))
hseco <- mapply(function(x) {length(dTable[,3][dTable[,3] == udom[x]])}, 1:length(udom))

hdom <- data.frame(Species = I(udom), Primary = hprim, Secondary = hseco)
h2dom <- hdom[order(hdom[,2], decreasing = T),]

#Generate plot of primary and secondary dominant understory plant functional groups.

#par(mai = c(2.5,1.2,1,1))
#barplot(t(cbind(h2dom$Primary,h2dom$Secondary)), main = "", 
#        xaxt = "n", ylab = "Number of sites", xlab = "", axes = F, beside = T, 
#        col = c("white", "dark grey"))
#axis(2)
#text(matrix(c(seq(2,nrow(h2dom)*3,3), rep(-0.25,nrow(h2dom))), nrow = nrow(h2dom), 
#            ncol = 2, byrow = F), srt = 60, adj = 1, xpd = T, labels = paste(h2dom$Species), 
#     cex = 0.95)
#legend(15,6, c("Highest biomass", "Second highest biomass"), fill = c("white", "dark grey"), bty = "n")


#Generate summary stats for cleaned plant functional group data.
fss <- round(stat.desc(fco),2)

#fss
#fco
#herb_mat <- cbind(fco$graminoids, fco$forb)
#wood_mat <- cbind(fco$shrub, fco$subshrub, fco$understory, fco$vine)
#herb_load <- apply(herb_mat, 1, sum)
#wood_load <- apply(wood_mat, 1, sum)
#load_mat <- cbind(herb_load, wood_load)
#rownames(load_mat) <- rownames(fco)
#colnames(load_mat) <- c("herb", "woody")
#load_total <- apply(load_mat, 1, sum)
#herb_perc <- round((herb_load/load_total)*100,1)
#median(herb_perc)
#plot(sort(herb_perc), decreasing = T)

#Mean percent herb loading
#mean(herb_perc)

#Range for percent herb loading
#range(herb_perc)

#Mean and range of live understory loading
#round(mean(load_total),2)
#round(range(load_total),2)

#Shrub loading
#round(stat.desc(wood_load),2)

#FORB
#fo <- data.frame(site = rownames(fugr), forb = round(fugr$forb,2), mfri = siteEnv4$mfri_20yr)
#fo[order(fo$mfri),]

#Generate more specific data used in discussion section
forb <- data.frame(mfri = env$mfri, forb = fco$forb)
gram <- data.frame(mfri = env$mfri, graminoids = fco$graminoids)

#Show loading for specified mfri window used in boundary layer regressions
#Forbs
mfri.min <- 6.0#low end of mfri (nearest whole number)
mfri.max <- 7.0#high end of mfri (nearest whole number)
sort(forb$forb[forb$mfri > mfri.min & forb$mfri <= mfri.max ], decreasing = T)

#Graminoids
mfri.min <- 2.0#low end of mfri (nearest whole number)
mfri.max <- 7.0#high end of mfri (nearest whole number)
sort(gram$graminoids[gram$mfri > mfri.min & gram$mfri <= mfri.max ], decreasing = T)
gram[order(gram$graminoids),]

###################################################################################################
###################################################################################################
#STEP #7: Gradient length diagnostics with DCA for untransformed biomass data, log-transformed
#biomass data, and untransformed cover data.

#Create a bindary dataset for species presence absence
#Biomass - log-transformed
#Note: output is the same whether you use untransformed or transformed data.
speocc <- data.trans(fco_log, method = "power", exp = 0, plot = F)

decorana(speocc, ira = 0)

#The 3-4 cut-off proposed by ter Braak and Smilauer (2002) used as a threshold 
#for selecting a multivariate analysis technique (< 3 PCA/RDA; data are suited for linear method, and
#> 4 CA/CCA; data are suited for unimodal methods).

#Values less than 3-4 indicate a short gradient length (i.e., use linear analysis). This is consistent 
#with sites selection which focused on a narrow set of environmental conditions. 
#I.e. mesic flatwoods with a history of management with frequent rx fire and no other evidence of 
#major disturbance.

###################################################################################################
###################################################################################################
#STEP #8: PCA of the covariance matrix
pca <- rda(decostand(fco_log, method = "hellinger"), scale = T)#scale = T should be PCA of the correlation matrix
#transformation with decostand() -- use hellinger transformation
#arg: scale = T - this standardizes the data and gives it unit variance (i.e., all variables have a variance = 1)
#the rda() function automatically centers the data around the mean for each variable.

scores_1 <- scores(pca)
sc <- scores_1$sites[,1]
plot(sc, env$mfri)
plot(sc, env$coarseWD)
cor(sc, env$mfri)
cor(sc, env$coarseWD)
#Provide info on PCA
pca
summary(pca)

#Biplot for PCA show labels on sites.
biplot(pca)

#Create a simplified biplot showing sites as dots grouped by region.

#Set up biplot with species loadings
par(cex = 1.75, cex.lab = 0.7, cex.axis = 0.7)
scl <- 2
biplot(pca, type=c("text", "none"), col=c("black", "black"), xlab="", 
       ylab="", scaling = scl)

#Add variance explained to each axis.
title(xlab = paste ("PC1 (",
                   round((eigenvals(pca)[1]/sum(eigenvals(pca)))*100,0), 
                   " percent)", sep = ""), 
      ylab = paste("PC2 (",
                   round((eigenvals(pca)[2]/sum(eigenvals(pca)))*100,0), 
                   " percent)", sep = ""), 
      mgp=c(2.2, 2.2, 0))

#Add sites to biplot categorized by region.
par(cex = 2)

points.sites <- ifelse(region == 1, 16, 1)
color <- ifelse(region == 1,"black", "black")
points(pca, pch = points.sites, col = color, cex = 0.75)

par(cex = 1.5)

legend("topright", legend = c("Eglin AFB", "Appalachicola NF/St. Marks NWR"), bty = "n", 
       col = c("black", "black"), pch = c(16,1), pt.bg = c(16,1), cex = 1)


#Show percent of variation explained by first four PC axes.
round((eigenvals(pca)[1]/sum(eigenvals(pca)))*100,1)
round((eigenvals(pca)[2]/sum(eigenvals(pca)))*100,1)
round((eigenvals(pca)[3]/sum(eigenvals(pca)))*100,1)
round((eigenvals(pca)[4]/sum(eigenvals(pca)))*100,1)

#PCA scores
#scores(pca, choices = 1:3, display = "sites", scaling = "sites", correlation = T)
loadings <- scores(pca, choices = 1:4, display = "species", scaling = "species", correlation = T)
sort(abs(loadings[,1]), decreasing = T)
sort(abs(loadings[,2]), decreasing = T)
sort(abs(loadings[,3]), decreasing = T)

#Show biplot for axis 1 and 3 so you can see how bunchgrass is arranged relative to the first axis.
dev.off()
ordiplot (pca, choices = c(1,3), type = "text")
ordiplot (pca, choices = c(2,3), type = "text")

#Plot pca
#biplot(pca, scaling = "sites")
#biplot(pca, scaling = "symmetric")
#biplot(pca, scaling = "species")
#biplot(pca, scaling = "none")

#Use brokan stick method to determine how many PCA axes to retain
screeplot(pca, bstick = T, type = "l", main = NULL)
#Extract eigenvalues
eigenvals(pca)
summary(eigenvals(pca))

#Fitting environmental data
#set.seed(42)#use this to make this permutation analysis reproducibe, otherwise it will be different
#every time, especially when number of permutations is low.
#ev <- envfit(pca ~ ., data = siteEnv3, choices = 1:2, scaling = 'symmetric', permutations = 1000)
#options for scaling: "species", "sites", "symmetric" (scales both for species and sites), "none" (raw scores).
#cca() hill = T, rda(): correlation = T (standardizes scores among species).
#ev

plot(pca, display = "sites", type = "n", scaling = "symmetric")
points(pca, display = "sites", scaling = "symmetric")
plot(ev, add = T)

par(cex = 0.5)
surf <- ordisurf(pca ~ canopy, data = env, knots = 1, isotropic = T, main = NULL)
surf <- ordisurf(pca ~ CoarseWD, data = siteEnv3, knots = 10, isotropic = T, main = NULL)
surf <- ordisurf(pca ~ Litter, data = siteEnv3, knots = 10, isotropic = T, main = NULL)
surf <- ordisurf(pca ~ mfri_20yr, data = siteEnv3, knots = 5, isotropic = T, main = NULL)
summary(surf)

###################################################################################################
###################################################################################################
#STEP #9: RDA of the covariance matrix

#Which transformation do you want to pass to the rda()
spp <- fco_log

#Generate full RDA model
set.seed(1)
pfg_rda_upr <- rda(spp ~ canopy + ff + coarseWD + mfri + sd_fri + gd_ratio + rxfire + Condition(region), data = env)

#Generate the null RDA model
set.seed(41)
pfg_rda_lwr <- rda(spp ~ 1., data = env)

#Conduct stepwise forward selection model building
set.seed(1)
pfg_rda_1 <- ordistep(pfg_rda_lwr, scope = formula(pfg_rda_upr), trace = T)

#View model details
pfg_rda_1
anova(pfg_rda_1, by = "terms")

#Check variance inflation factor of your terms. Two rules of thumb.
#1) Anything higher than 20 should be dropped from analysis. 
#2) Anything higher than 10 should be considered for removal.
vif.cca(pfg_rda_1)

#Here is your RDA model
#It contains a single term (coarseWD). No other terms are significant.
#How much variance is explained by this model:
8.331/37.188 #(Constrained Inertia/Total Inertia)
#22.4%

#Show terms and order they were added to model
#They will be listed in order (top to bottom > first to lest variable)
pfg_rda_1$anova

scores_1 <- scores(pfg_rda_1)
sc <- scores_1$sites[,2]
plot(sc, env$mfri)

plot(sc, env$coarseWD)




#Stepwise with adjusted R2 -- forward selection
#This is important to do when your model has a lot of explanatory variables
#because each time you add a new variable the R2 gets artificially inflated.
#This method reduces the R2 by looking at the number of terms (expl. vars)
#relative to the sample size. As the number of terms increase and the number
#of samples decrease the adjusted R2 receives a greater reduction relative
#to the unadjusted R2.

#Two step solution for testing models.
#Blanchet et al. 2008
#Global test of all constraints
#If test is significant then proceed to stepwise model building
#Add constraints if term has P < 0.05 and adjusted R2 exceeds
#the global adjusted R2.
#This is how ordiR2step() builds models.
set.seed(1)
pfg_rda_2 <- ordiR2step(pfg_rda_lwr, pfg_rda_upr, trace = T)
pfg_rda_2$anova
anova(pfg_rda_2, by = "terms")

#RDA did not identify significant terms
#Do any of the explantory variables differ significantly between the two regions.
#Canopy
t.test(canopy ~ region, data = env)
plot(canopy ~ region, data = env)
#No significance between regions

#Forest floor
t.test(ff ~ region, data = env)
plot(ff ~ region, data = env)
#No significance between regions

#Coarse DWD
t.test(coarseWD ~ region, data = env)
plot(coarseWD ~ region, data = env)
#No significance between regions

#Mean fire return interval (20 years)
t.test(mfri ~ region, data = env)
plot(mfri ~ region, data = env)
#Yes, there is a significant difference

#Standard deviation of mFRI
t.test(sd_fri ~ region, data = env)
plot(sd_fri ~ region, data = env)
#Yes, there is a significant difference

#Ratio of growing:dormant season burns
t.test(gd_ratio ~ region, data = env)
plot(gd_ratio ~ region, data = env)
#No significant difference between regions

#Length of time units were regularly burned
t.test(rxfire ~ region, data = env)
plot(rxfire ~ region, data = env)
#Yes, there is a difference between regions

###################################################################################################
###################################################################################################
#Conduct a RDA only on contraints you have a hypothesis for.
set.seed(1)
pfg_rda_a <- rda(spp ~ mfri + gd_ratio + Condition(region), data = env)
vif.cca(pfg_rda_a)
pfg_rda_a
pfg_rda_a$anova

###################################################################################################
###################################################################################################
#Boundary layer regression for plant functional groups

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


#Old file reassignments
emat <- siteEnv3
ab2mat <- fco

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
set.panel(2,3)
cf <- 0.5

#1
#Shrubs
#Plot with sites
species <- ab2mat$shrub
speciesLabel <- "Shrubs"

#Boundary layer regression
blr(model, dec, xcoord, ev, species, speciesLabel, envLabel)

#2
#Graminoids
#Plot with sites
species <- ab2mat$graminoids
speciesLabel <- "Graminoids"

#Boundary layer regression
blr(model, dec, xcoord, ev, species, speciesLabel, envLabel)

#3
#Sub-shrubs
#Plot with sites
species <- ab2mat$subshrub
speciesLabel <- "Sub-shrub"

#Boundary layer regression
blr(model, dec, xcoord, ev, species, speciesLabel, envLabel)

#4
#Forbs
#Plot with sites
species <- ab2mat$forb
speciesLabel <- "Forbs"

#Boundary layer regression
blr(model, dec, xcoord, ev, species, speciesLabel, envLabel)

#5
#Vines
#Plot with sites
species <- ab2mat$vine
speciesLabel <- "Vines"

#Boundary layer regression
blr(model, dec, xcoord, ev, species, speciesLabel, envLabel)

#6
#Understory trees
#Plot with sites
species <- ab2mat$understory
speciesLabel <- "Understory trees"

#Boundary layer regression
blr(model, dec, xcoord, ev, species, speciesLabel, envLabel)

#END
