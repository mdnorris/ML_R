###########################################################################################
# CI_Scoring_Master_RF.R
# Created by Matt Norris 7/15/2015
# Modified by Matt Norris 12/13/2015

################# Libary Installation, Setting WD #############################

install.packages("corrplot")
install.packages("gmodels")
install.packages("caret")
install.packages('rattle')
install.packages('rpart.plot')
install.packages('RColorBrewer')
install.packages('randomForest')
install.packages('party')
install.packages('CALIBERrfimpute') 
install.packages('dplyr') 

setwd("C:/Drive/Code/Ben Graham")

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
library(dplyr)
library(gmodels)
library(corrplot)

######################## Data Prep ###################################

#Import full scoring method dataset
all_CI_vars  <- read.csv("all_vars_100.csv", header = TRUE, sep = ",")

#I included a bit of code that worked for me to directly import a SPSS dataset and perform the same operations
#dataset = read.spss("C:\\PathToFile\\MyDataFile.sav", to.data.frame=TRUE)

#### Create string vectors to subset data

#Ranking with health score, hazard score, health score, climate vulnerability,
#and water quality score; these strings are variable names or column headers, and 
#the exist to easily subset data further in the program

myvars.CIscoreCVWQ <- c(
  "stHazScore"
  ,"stHealthScore" 
  ,"stSVscore" 
  ,"stWaterScore"
  ,"st_CIscoreCVWQ"
  ,"stCCVscore"
)

#this next section will actually create the subsetted data 
#will be imputed and put through the random forest algorithm
#an example such as this line of code: SVscore <- SVscore[sample(nrow(SVscore), "2000"), ]
#will drop observations and so make the time consuming operations of imputing and 
#bootstrapping perform more quickly, perhaps ideal for testing the program

myvars.CCVscore  <- colnames(subset(all_CI_vars, select = c(tract_novehicle:tract_warm_nights_change,stCCVscore)))
myvars.HazScore  <- colnames(subset(all_CI_vars, 
                                    select = c(tr_sl_popwt_senslu:tr_sl_popwt_distwt_haz_senslu_traffic_01,stHazScore)))
myvars.SVscore  <- colnames(subset(all_CI_vars, select = c(tract_population:tract_turnout,stSVscore)))
myvars.HealthScore  <- colnames(subset(all_CI_vars, select = c(tract_pm25:tract_rhaz,stHealthScore)))

CIscoreCVWQ  <- all_CI_vars[myvars.CIscoreCVWQ]
CCVscore <- all_CI_vars[myvars.CCVscore]
HazScore <- all_CI_vars[myvars.HazScore]
SVscore  <- all_CI_vars[myvars.SVscore]
HealthScore <- all_CI_vars[myvars.HealthScore]

# Setting seeds so we get reproducible results (works for remainder of program)

set.seed(1)

#################################################################################
############################ CIscoreCVWQ ########################################

#Fix missing values, all factor variables with less than 1% missing, omit observations
#and calculate random forest with proximity and important matrices
#This was just for example purposes of how variables may steal importance from each other
#so it will be ommitted in this thorough run

nm.CIscoreCVWQ <- na.omit(CIscoreCVWQ)
rf.CIscoreCVWQ <- randomForest(st_CIscoreCVWQ ~ ., nm.CIscoreCVWQ,
                              proximity=TRUE, importance=TRUE, ntree=500)

#Format and plot variable importance
par(mfrow=c(1, 1))
varImpPlot(rf.CIscoreCVWQ)

#Formats and plot partial dependence
par(mfrow=c(1, 1))
imp <- importance(rf.CIscoreCVWQ)
imp.list <- rownames(imp)[order(imp[, 1], decreasing=TRUE)]
for (i in seq_along(imp.list) ){
 partialPlot(rf.CIscoreCVWQ, nm.CIscoreCVWQ, imp.list[i],
             xlab=imp.list[i], main=paste("Partial Dependence on",
                                          imp.list[i]), ylab="rf.CIscoreCVWQ", plot=TRUE, n.pt=25,"test")
}

#################################################################################
############################ CCVscore ########################################

#Fix missing values, all factor variables with less than 1% missing, omit observations
#and calculate random forest with proximity and important matrices
nm.CCVscore  <- na.roughfix(CCVscore)

#Change to factors for classification random forest, create proximity and importance plots
nm.CCVscore$stCCVscore <- as.factor(nm.CCVscore$stCCVscore)
rf.CCVscore <- randomForest(stCCVscore ~ ., nm.CCVscore, proximity=TRUE, importance=TRUE, ntree=500)

#Plot variable importance
varImpPlot(rf.CCVscore)

#formats for partial dependence plots
par(mfrow=c(1, 1)) 
imp <- importance(rf.CCVscore)
imp.list <- rownames(imp)[order(imp[, 1], decreasing=TRUE)]
for (i in seq_along(imp.list) ){
  partialPlot(rf.CCVscore, nm.CCVscore, imp.list[i], 
              xlab=imp.list[i], main=paste("Partial Dependence on", 
                                           imp.list[i]), ylab="rf.CCVscore", plot=TRUE, n.pt=25)
}

#################################################################################
############################ HazScore ########################################

nm.HazScore <- mice(HazScore, meth = "rf", ntree = 5)
nm.HazScore <- complete(nm.HazScore,5)
#Fill in missings and format for correlations
nm.HazScore  <- na.roughfix(HazScore)
par(mfrow=c(1, 1))

#Convert dependent to factor and run classification random forest
nm.HazScore$stHazScore <- as.factor(nm.HazScore$stHazScore)
rf.HazScore <- randomForest(stHazScore ~ ., nm.HazScore, proximity=TRUE, importance=TRUE, ntree=500)

#Plot variable importance
varImpPlot(rf.HazScore)

#Format and plot partial dependence
par(mfrow=c(1, 2)) 
imp <- importance(rf.HazScore)
imp.list <- rownames(imp)[order(imp[, 1], decreasing=TRUE)]
for (i in seq_along(imp.list) ){
  partialPlot(rf.HazScore, nm.HazScore, imp.list[i], 
              xlab=imp.list[i], main=paste("Partial Dependence on", 
                                           imp.list[i]), ylab="rf.HazScore", plot=TRUE)
}

#################################################################################
############################ SVscore ########################################

#Fill in missing values, format for correlations and variable importance plot
nm.SVscore  <- na.roughfix(SVscore)
par(mfrow=c(1, 1))


#convert dependent to factor and run classification random forest
nm.SVscore$stSVscore <- as.factor(nm.SVscore$stSVscore)
rf.SVscore <- randomForest(stSVscore ~ ., nm.SVscore, proximity=TRUE, importance=TRUE, ntree=500)
par(mfrow=c(1, 1))

#Plot variable importance
varImpPlot(rf.SVscore)

#formats and plot partial dependence 
par(mfrow=c(1, 1)) 
imp <- importance(rf.SVscore)
imp.list <- rownames(imp)[order(imp[, 1], decreasing=TRUE)]
for (i in seq_along(imp.list) ){
  partialPlot(rf.SVscore, nm.SVscore, imp.list[i], 
              xlab=imp.list[i], main=paste("Partial Dependence on", 
                                           imp.list[i]), ylab="rf.SVscore", plot=TRUE)
}

#################################################################################
############################ HealthScore ########################################

#Format and run correlations
nm.HealthScore  <- na.roughfix(HealthScore)
par(mfrow=c(1, 1))

#convert to factor and run classification random forest
nm.HealthScore$stHealthScore <- as.factor(nm.HealthScore$stHealthScore)
rf.HealthScore <- randomForest(stHealthScore ~ ., nm.HealthScore, proximity=TRUE, importance=TRUE, ntree=500)

#plot variable importance
varImpPlot(rf.HealthScore)

#Format and run partial dependence plots
par(mfrow=c(1, 2)) #formats for partial dependence plots
imp <- importance(rf.HealthScore)
imp.list <- rownames(imp)[order(imp[, 1], decreasing=TRUE)]
for (i in seq_along(imp.list) ){
  partialPlot(rf.HealthScore, nm.HealthScore, imp.list[i], 
              xlab=imp.list[i], main=paste("Partial Dependence on", 
                                           imp.list[i]), ylab="rf.HealthScore", plot=TRUE, n.pt=5)
}

############################# Correlation Matrix ########################

nm.CIscoreCVWQ <- na.omit(CIscoreCVWQ)
nm.CCVscore  <- na.roughfix(CCVscore)
nm.HazScore  <- na.roughfix(HazScore)
nm.SVscore  <- na.roughfix(SVscore)
nm.HealthScore  <- na.roughfix(HealthScore)

# format and plot correlation matrix for CCVscore

par(mfrow=c(1, 1))

cor.CCVscore <- cor(nm.CCVscore, use="complete.obs", method="pearson") 
corrplot(cor.CCVscore, method = "circle")

#  plot correlation matrix for Hazscore

cor.HazScore <- cor(nm.HazScore, use="complete.obs", method="pearson") 
corrplot(cor.HazScore, method = "circle")

#plot correlation matrix for SVscore

cor.SVscore <- cor(nm.SVscore, use="complete.obs", method="pearson") 
corrplot(cor.SVscore, method = "circle")

#plot correlation matrix for Healthscore

cor.HealthScore <- cor(nm.HealthScore, use="complete.obs", method="pearson") 
corrplot(cor.HealthScore, method = "circle")

#Not going to work without supercomputer
MDSplot(rf.HazScore, rf.HazScore$stHazScore,fill=brewer.pal(4, "GnBu"))
legend("topleft", legend=levels(rf.HazScore$stHazScore), fill=brewer.pal(4, "GnBu")) 


