###########################################################################################
# FDI_Random_Forest.R
# Created by Matt Norris 7/15/2015
# Modified by Matt Norris 12/11/2015

################# Library Installation, Setting WD #############################

#Package Installation only needs to be run once
#package documentation in working directory's Packages folder

install.packages("caret")
install.packages("mice")
install.packages("readstata13")
install.packages('rattle')
install.packages('rpart.plot')
install.packages('RColorBrewer')
install.packages('randomForest')
install.packages('party')
install.packages('CALIBERrfimpute')
install.packages("corrplot")

#This section needs to be run once per session

setwd()

library(readstata13)
library(foreign)
library(rpart)
library(rattle)
library(rpart.plot)
library(RColorBrewer)
library(randomForest)
library(party)
library(CALIBERrfimpute)
library(missForest)
library(survival)
require(graphics)
library(caret)
require(graphics)
library(corrplot)
library(mice)

######################## Data Prep ################## #################

master_ipe <- read.dta13("IPE_Master_1988-2014 public_v1_2.dta")

# Excel for Inspection
# write.csv(master_ipe, file = "ipe_explore.csv", row.names = TRUE)

# Create string vector for subsetting

myvars       <- c(
#dependent variables
   "fdidollars_WB" 
   ,"fdi_stocks_UTD"
#ID variables
   ,"year"
   ,"gwno"
   ,"country"
   ,"incomegroup_FK"
#independent
   ,"gdp_full" 
   ,"gdppc_full" 
   ,"growth_full"
   ,"trade_WB"
   ,"MFNtariff_MFN"
   ,"exp_WB"
   ,"combinedoil_AE"
   ,"natresource_rents_WB"
   ,"govt_consump_WB"
   ,"surdef_WB"
   ,"democ_P4"
   ,"yr_sch_BL" 
   ,"kai_FK"
   )

###################### Feature Engineering #######################

################################################################

# Cut to variables of interest

#Data trimming for testing purposes

ipe_trim  <- master_ipe[myvars]

# Lower computation times in testing are advised, random subsets
# should be used unless doing final run

#ipe_trim_short <- ipe_trim[sample(nrow(ipe_trim), "50"), ]

# Finally, unclass converts our characters to factors 
# which is necessary for missing variable imputation
# comment out first comment for quick run, second for complete run

final_ipe <-as.data.frame(unclass(ipe_trim))
#final_ipe <-as.data.frame(unclass(ipe_trim_short))

############################ Eliminating Missing Values ##############################

#Random Forests will not work with missing values so imputation is necessary

# Setting seeds so we get reproducible results (works for remainder of program)

set.seed(1)

#Random Forest Multiple Imputation by Chained Equations (MICE)
#is extremely resource intensive, suggest the roughfix for testing
#and quick estimates

ipe_wo_missing <- mice(final_ipe, meth = "rf", ntree = 3)

#Roughfix method

#ipe_wo_missing <- na.roughfix(final_ipe)

#Used to merge MICE imputations with unimputed dataset,
#comment out for roughfix run

ipe_wo_missing <- complete(ipe_wo_missing,5)

################################# Random Forest Model ####################################

# rforest() can only handle 52 factors so we drop them now
# these countries can be matched back in by country code if necessary

ipe_wo_missing$country <- NULL

#Note: reduce ntree for faster run to 501, in final run, ntree=5000 is sufficient

rf.fdi_stocks_UTD <- randomForest(fdi_stocks_UTD ~ ., ipe_wo_missing, proximity=TRUE, importance=TRUE, ntree=501)
rf.fdidollars_WB <- randomForest(fdidollars_WB ~ ., ipe_wo_missing, proximity=TRUE, importance=TRUE, ntree=501)

####################### Random Forest Variable Importance #############################

varImpPlot(rf.fdi_stocks_UTD)
varImpPlot(rf.fdidollars_WB)

############################# Proximity Matrix ########################

par(mfrow=c(1, 1)) #sets formatting for proximity matrix 

#This MDSplot gives an example of plotting a proximity matrix
#There are other ways to plot this, but I've found this the most effective
#to use this properly, the response variable must be a factor
#after dropping countries, we only have incomegroup_FK as factor variable
#I wanted to include this code for future analysis
#I will explain the interpretation in the documentation

MDSplot(rf.fdi_stocks_UTD, ipe_wo_missing$incomegroup_FK, palette=brewer.pal(4, "GnBu"))
legend("topleft", legend=levels(ipe_wo_missing$incomegroup_FK), fill=brewer.pal(4, "GnBu")) 

######################### Partial Dependence Plots #######################

#This is resource intensive so minimize processing time, so again use random subsets
par(mfrow=c(2, 3)) #formats for partial dependence plots

imp_fdi_stocks_UTD <- importance(rf.fdi_stocks_UTD)
imp_fdi_stocks_UTD_list <- rownames(imp_fdi_stocks_UTD)[order(imp_fdi_stocks_UTD[, 1], decreasing=TRUE)]
for (i in seq_along(imp_fdi_stocks_UTD_list) ){
  partialPlot(rf.fdi_stocks_UTD, ipe_wo_missing, imp_fdi_stocks_UTD_list[i], 
              xlab=imp_fdi_stocks_UTD_list[i], main=paste("Partial Dependence on", 
              imp_fdi_stocks_UTD_list[i]), ylab="fdi_stocks_UTD", plot=TRUE, n.pt=5,"test")
}

imp_fdidollars_WB <- importance(rf.fdidollars_WB)
imp_fdidollars_WB_list <- rownames(imp_fdidollars_WB)[order(imp_fdidollars_WB[, 1], decreasing=TRUE)]
for (i in seq_along(imp_fdidollars_WB_list) ){
  partialPlot(rf.fdidollars_WB, ipe_wo_missing, imp_fdidollars_WB_list[i], 
              xlab=imp_fdidollars_WB_list[i], main=paste("Partial Dependence on", 
              imp_fdidollars_WB_list[i]), ylab="fdidollars_WB", plot=TRUE, n.pt=5,"test")
}

##################### Correlations ###########################

ipe_wo_missing$incomegroup_FK <- NULL

par(mfrow=c(1, 1))
cor.ipe_wo_missing <- cor(ipe_wo_missing, use="complete.obs", method="pearson") 
corrplot(cor.ipe_wo_missing, method = "circle")

#Note: This corrplot can cause issues for formatting, so clear environment
#if necessary
######################## End of Program ######################
