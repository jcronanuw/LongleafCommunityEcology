###################################################################################################
############################-----start-----########################################################
#                                        GENERA-LEVEL ANALYSIS
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
                                "phd_chapter4_biomassPlotMatrix_2yr_Ep_1_4_FuncGroup_2022-10-20_14.35.08.csv",
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
cbind(colnames(siteBiomass6), dropped)
length(dropped)
length(dropped) - length(dropped[is.na(dropped) == T])
#No drops.

###################################################################################################
###################################################################################################
#STEP #5: DATA SCREENING

#########################################################
#5a: Summary stats for biomass data

#Summary stats
round(stat.desc(fugr),2)

#Look at how data is distributed
foa.plots(fugr)

#Drop rare genus < 10% of sites.
fco <- drop.var(fugr, min.po=10)

#Check to see how many genus were dropped for untransformed biomass data
length(fugr[1,])#
length(fco[1,])#
#No drops

str(fco)
#21 obs and 7 variables

#uv.plots displays histogram, box and whisker, cumulative distribution and normal q-q plots
#in one pane for each variable.
uv.plots(fco)#About half of functional groups  are right skewed with a long right tail. 
#Conduct a log tranformation

#What percent of values are zero, if it is over 50% you should consider changing the data
#to presence/absence
length(fco[fco == 0])/(length(fco[1,])*length(fco[,1]))
#~8% zero values >>> no need to convert data to presence/absence

#Look at correlation between genus variables.
#generate table (you need to run correlation_matrix() function for this to work:
#https://www.r-bloggers.com/2020/07/create-a-publication-ready-correlation-matrix-with-significance-levels-in-r/
scm <- correlation_matrix(as.data.frame(fco), type = "pearson", show_significance = T, 
                         digits = 2, use = "lower", replace_diagonal = T)
chart.Correlation(fco, method = "pearson")#performanceAlanlytics
#Understory trees and vines are highly correlated (0.65). All other are below 0.55.

#log transform genera data.
fco_log <- decostand(fco, method = "log")

#Check distributions
uv.plots(fco_log)
#Improved for all.

#Hellinger transformation on genera data
fco_hel <- decostand(fco, method = "hellinger")
#Check distributions
uv.plots(fco_hel)
#Improvement over log transformation distributions, they are more 'normal'.

#########################################################
#6c: Summary stats for environmental data

#Show how environmental variables are correlated.
chart.Correlation(siteEnv3)
scm <- correlation_matrix(as.data.frame(siteEnv3), type = "pearson", show_significance = T, 
                          digits = 2, use = "lower", replace_diagonal = T)
#Use 0.55 as cut-off for removing variables with high colinearity 
#Sort of arbitrary, but chose (Adam et al. 2013).
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

#Check distributions of environmental data
#Remove categories
check_env <- siteEnv4[,-c(9,10)]
uv.plots(check_env)
#Right skewed
#1) Coarse DWD
#2) Std Dev on 20-year fire regime

#Log transform right-skewed variables
cdwd_log <- decostand(siteEnv4$CoarseWD, method = "log")
stde_log <- decostand(siteEnv4$StDevFRI_20YR, method = "log")

#Create new environmental variables data frame
env <- data.frame(canopy = siteEnv4$Canopy,
                     litter = siteEnv4$Litter,
                     coarseWD = cdwd_log,
                     fineWD = siteEnv4$FineWD,
                     mfri = siteEnv4$mfri_20yr,
                     sd_fri = stde_log,
                     gd_ratio = siteEnv4$Season_20yr,
                     rxfire = siteEnv4$RxFireMgmt,
                     region = as.factor(siteEnv4$Region))
rownames(env) <- rownames(siteEnv4)


uv.plots(env)#looks good

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
mv.outliers(env[1:8], method = "euclidean", sd.limit =2)#S330

#biomass
mv.outliers(fco_log, method = "bray", sd.limit =2)#S330

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

par(mai = c(2.5,1.2,1,1))
barplot(t(cbind(h2dom$Primary,h2dom$Secondary)), main = "", 
        xaxt = "n", ylab = "Number of sites", xlab = "", axes = F, beside = T, 
        col = c("white", "dark grey"))
axis(2)
text(matrix(c(seq(2,nrow(h2dom)*3,3), rep(-0.25,nrow(h2dom))), nrow = nrow(h2dom), 
            ncol = 2, byrow = F), srt = 60, adj = 1, xpd = T, labels = paste(h2dom$Species), 
     cex = 0.95)

legend(15,6, c("Highest biomass", "Second highest biomass"), fill = c("white", "dark grey"), bty = "n")

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
pca <- rda(decostand(gco, method = "hellinger"), scale = F)#scale = T should be PCA of the correlation matrix
#transformation with decostand() -- use hellinger transformation
#arg: scale = T - this standardizes the data and gives it unit variance (i.e., all variables have a variance = 1)
#the rda() function automatically centers the data around the mean for each variable.

scores(pca, choices = 1:3, display = "sites", scaling = "sites", correlation = T)
scores(pca, choices = 1:3, display = "species", scaling = "species", correlation = F)

pca
#Plot pca
biplot(pca, scaling = "sites")
biplot(pca, scaling = "symmetric")
biplot(pca, scaling = "species")
biplot(pca, scaling = "none")

#Use brokan stick method to determine how many PCA axes to retain
screeplot(pca, bstick = T, type = "l", main = NULL)
#Extract eigenvalues
eigenvals(pca)
summary(eigenvals(pca))

#Plot data with different symbols for each site location
sitetype <- as.numeric (as.factor (siteEnv3$Region))
ordiplot (pca, display = 'sites', type = 'n')
points (pca, pch = sitetype, col = sitetype)

#Determine which genus are most correlated with PC axes.
loadings <- scores(pca, display = "species", scaling = 0)
sort(abs(loadings[,1]), decreasing = T)
sort(abs(loadings[,2]), decreasing = T)


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
#STEP #9: RDA of the covariance matrix

#Which transformation do you want to pass to the rda()
spp <- fco_log

#Generate full RDA model
gen_rda_upr <- rda(spp ~ ., data = env)
gen_rda_upr
ps_upr <- permustats(anova(gen_rda_upr))#anova is testing the entire model for significance.
summary(ps_upr)
densityplot(ps_upr)

set.seed(45)#use this so results are reproducible
pstat_1 <- permustats(anova(gen_rda_upr))#anova is testing the entire model for significance.
summary(pstat_1)
densityplot(pstat_1)
#((1 - Pr(perm) value) *100)% of the permuted F values were lower than the observed F statistic
#this is good if higher than 95%, we want the Observed F to be in the extreme upper tail, i.e., we want it to be big.
set.seed(1)
anova(gen_rda_upr, by = "terms")#check out upr model

#Generate the null RDA model
set.seed(45)
gen_rda_lwr <- rda(spp ~ 1., data = env)

set.seed(1)
gen_rda_1 <- ordistep(gen_rda_lwr, scope = formula(gen_rda_upr), trace = T)
summary(gen_rda_1)
gen_rda_1$anova
anova(gen_rda_1)

#Stepwise with adjusted R2
set.seed(1)
gen_rda_2 <- ordiR2step(gen_rda_lwr, scope = formula(gen_rda_upr), trace = T)
gen_rda_2$anova
anova(gen_rda_2, by = "terms")

#Test starting model, but with region is a conditional block
c1 <- rda(fco_hel ~ fineWD + coarseWD + mfri + Condition(region), data = env)
h <- how(within = Within(type = "free"), plots = Plots(strata = env$region))
set.seed(42)
tt <- anova(c1, permutations = h, model = "reduced")
tt
summary(tt)



###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
#Extra code -- part of previous analyses - may be useful for this analysis



###################################################################################################
#STEP #10b: Principal Coordinate Analysis on original biomass data for log-transformed data.

lcm01 <- vegdist(gco, "manhattan")#untransformed data with manhattan distance measure.
lcm02 <- vegdist(gco_log, "bray")#untransformed data with bray-curtis distance measure.

pcoa01 <- cmdscale(lcm01, k = 10, eig = T, add = F)
pcoa02 <- cmdscale(lcm02, k = 10, eig = T, add = F)

#Percent variation explained by each principal coordinate
round(pcoa01$eig/sum(pcoa01$eig)*100,2)
round(pcoa02$eig/sum(pcoa02$eig)*100,2)

#Assign pcoa object to functions between  #>>>><<<<<# lines
x <- pcoa01
y <- gco
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<##
#Broken stick plot
plot(x$eig/sum(pcoa01$eig)*100-bstick(21)*100,xlab = "PC", 
     ylab="Actual-random % variation explained")
abline(h=0)
dev.off()


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






###################################################################################################
###################################################################################################
#STEP #14: #Distance-based RDA

###################################################################################################
#STEP #14a: untransformed biomass (manhattan distance)
bo2 <- gco
siteEnv4 <- env
y <- vector('numeric', length = 1000)
for (i in 1:1000)
{
gcap1 <- capscale(bo2 ~ mfri + gd_ratio + rxfire, 
                  data = siteEnv4, distance = "bray")

gcap1s <- summary(gcap1)

x <- anova(gcap1)#
y[i] <- x[1,4]
}
hist(y)
median(y)
range(y)

anova(gcap1, by = "axis")#
anova(gcap1, by = "term")#

plot(gcap1, choices = c(1,2), type = "points", display = "wa", scaling = 2)
plot(gcap1, choices = c(1,2), type = "n", display = "lc", scaling = 2)
text(gcap1, choices = c(1,2), labels = row.names(gco), cex = 0.8)
round(intrasetcor(gcap1), 5)
round(intersetcor(gcap1), 5)

plot(gcap1, choices = c(1,2), display = c("wa", "sp", "bp"), scaling = 3, 
     title = "Biomass by genera (model significance: p =)", xlab = "CAP1 (%)", 
     ylab = "CAP2 (%)")

#How does mFRI vary with CAP2 scores
plot(gcap1s$sites[,2], env$mfri[order(gcap1s$sites[,2], decreasing = T)], type = "n")
text(gcap1s$sites[,2], env$mfri[order(gcap1s$sites[,2], decreasing = T)], labels = rownames(env))

model.1 <- lm(as.vector(env$mfri[order(gcap1s$sites[,2], decreasing = T)] ~ gcap1s$sites[,2]))
summary(model.1)
abline(model.1)
segments(gcap1s$sites[,2], env$mfri[order(gcap1s$sites[,2], decreasing = T)], 
         gcap1s$sites[,2], predict(model.1))

#How does FineWD vary with CAP2 scores
plot(gcap1s$sites[,1], env$fineWD[order(gcap1s$sites[,1], decreasing = T)], type = "n")
text(gcap1s$sites[,1], env$fineWD[order(gcap1s$sites[,1], decreasing = T)], labels = rownames(env))

model.1 <- lm(as.vector(env$fineWD[order(gcap1s$sites[,1], decreasing = T)] ~ gcap1s$sites[,1]))
summary(model.1)
abline(model.1)
segments(gcap1s$sites[,1], env$fineWD[order(gcap1s$sites[,1], decreasing = T)], 
         gcap1s$sites[,1], predict(model.1))

#How does Litter vary with CAP2 scores
plot(gcap1s$sites[,1], env$litter[order(gcap1s$sites[,1], decreasing = T)], type = "n")
text(gcap1s$sites[,1], env$litter[order(gcap1s$sites[,1], decreasing = T)], labels = rownames(env))

model.1 <- lm(as.vector(env$litter[order(gcap1s$sites[,1], decreasing = T)] ~ gcap1s$sites[,1]))
summary(model.1)
abline(model.1)
segments(gcap1s$sites[,1], env$litter[order(gcap1s$sites[,1], decreasing = T)], 
         gcap1s$sites[,1], predict(model.1))



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
pca <- rda(decostand(bl, method = "hellinger"), scale = T)#scale = T should be PCA of the correlation matrix
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
m1 <- cca(gco ~ ., data = env)
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
