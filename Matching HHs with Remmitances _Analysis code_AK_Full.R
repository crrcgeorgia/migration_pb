######################################################
#Matching HHs with Remmitances VS HHs without Remittances
######################################################
title: "Migration Matching analysis code"

library(Matching)
library(MASS)
library(RTextTools)
library(haven)
library(descr)
library(Matching)
library(rgenoud)
library(psych)
library(sandwich)
library(lmtest)
library(multiwayvcov)
library(survey)
library(margins)
library(ggeffects)
library(sjmisc)
library(ggplot2)
library(effects)
library(MASS)
library(ggeffects)
library(emmeans)
library(erer)
library(sjPlot)
library(dplyr)
library(janitor)


rm(list = ls())
getwd()
setwd("D:\\anano\\Policy Briefs\\migration\\migration moneytot DKRA-it")
CB2019<-read_dta("CB2019_Georgia_response_30Jan2020 - UPDATED.dta")
designgeoun <- svydesign(id=~PSU,weights=~HHWT, strat=~SUBSTRATUM, data=CB2019)
names(designgeoun)

table(CB2019$MONYTOT)
sum(is.na(CB2019$WithoutMWithR))
sum(is.na(CB2019$SUBSTRATUM))
sum(is.na(CB2019$ETHNIC))
sum(is.na(CB2019$HHSIZE))
sum(is.na(CB2019$RELCONDrec))
sum(is.na(CB2019$RESPEDU))
sum(is.na(CB2019$shareFemale))
sum(is.na(CB2019$ECONSTNrec))
sum(is.na(CB2019$RATEHAPrec))
cols2=c('SUBSTRATUM', 'ETHNIC', 'HHSIZE','RESPEDU', 'shareFemale')
sum(is.na(CB2019[,cols2]))

#WithoutMWithR - 1 is with remittances, 0 is without remittances
CB2019_cln <- CB2019[!is.na(CB2019$WithoutMWithR),]
dim(CB2019)
dim(CB2019_cln)
CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$RATEHAPrec),]

CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$ETHNIC),]
dim(CB2019_cln)

CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$RELCONDrec),]
dim(CB2019_cln)

CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$RESPEDU),]
dim(CB2019_cln)

#economic situation
CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$ECONSTNrec),]
dim(CB2019_cln)

CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$wealth_index),]
dim(CB2019_cln)


#Estimand --	A character string for the estimand. 
##The default estimand is "ATT", the sample average treatment effect for the treated. 
#"ATE" is the sample average treatment effect (for all), and 
#"ATC" is the sample average treatment effect for the controls.
#M is a scalar for the number of matches which should be found. The default is one-to-one matching. Also see the ties option.

#Tr=CB2019_cln$WithoutMWithR - 1 if "With remittances from emigrant", 0 if "without"
CB2019_cln$SUBSTRATUM<-as_factor(CB2019_cln$SUBSTRATUM)
CB2019_cln$ETHNIC<-as_factor(CB2019_cln$ETHNIC)
CB2019_cln$RESPEDU<-as_factor(CB2019_cln$RESPEDU)
CB2019_cln$ECONSTNrec<-as_factor(CB2019_cln$ECONSTNrec)
CB2019_cln$MONYTOTcontext<-as_factor(CB2019_cln$MONYTOTcontext)
CB2019_cln$SPENDMOrec<-as_factor(CB2019_cln$SPENDMOrec)
CB2019_cln$RATEHAPrec<-as_factor(CB2019_cln$RATEHAPrec1)
CB2019_cln$RATEHAPrec<-as_factor(CB2019_cln$wealth_index)


psmodel<-glm(CB2019_cln$WithoutMWithR~ CB2019_cln$SUBSTRATUM
             + CB2019_cln$ETHNIC + CB2019_cln$HHSIZE 
             + CB2019_cln$RESPEDU + CB2019_cln$shareFemale, family = "binomial")
hist(psmodel$fitted.values)
CB2019_cln$prop.score<-psmodel$fitted.values


X <- cbind(CB2019_cln$SUBSTRATUM, CB2019_cln$ETHNIC,
           CB2019_cln$RESPEDU, CB2019_cln$HHSIZE,
           CB2019_cln$shareFemale,
           CB2019_cln$prop.score) 
BalanceMat <- cbind(CB2019_cln$SUBSTRATUM, CB2019_cln$ETHNIC,
                    CB2019_cln$RESPEDU, CB2019_cln$HHSIZE,
                    CB2019_cln$shareFemale, CB2019_cln$prop.score)
dim(BalanceMat)

#Now we are looking for optimal weights
genout <- GenMatch(Tr=CB2019_cln$WithoutMWithR, X=X,
                   BalanceMatrix=BalanceMat, estimand = "ATE", 
                   pop.size = 1000, wait.generations = 1)


#Matching
mout<- Match(Tr=CB2019_cln$WithoutMWithR, X=X, 
             Weight.matrix = genout)
summary(mout)


#Determining if balance has actually obtained or not on th variables of interst
mb<- MatchBalance(CB2019_cln$WithoutMWithR ~ CB2019_cln$SUBSTRATUM
                  + CB2019_cln$ETHNIC + CB2019_cln$HHSIZE 
                  + CB2019_cln$RESPEDU + CB2019_cln$shareFemale 
                  +CB2019_cln$prop.score, 
                  match.out = mout, nboots = 1000)
summary(mb)
treated<-CB2019_cln[mout$index.treated,]
treated$mweights<-mout$weights
treated$mweights

control<-CB2019_cln[mout$index.control,]
control$mweights<-mout$weights
control$mweights

MergedData <- merge(treated, control, all=TRUE)
glm(MergedData$RELCONDrec~MergedData$WithoutMWithR,w=MergedData$mweights)


#######################################
#######################################

##saving merged data as dta file 
setwd("D:\\anano\\Policy Briefs\\migration\\migration moneytot DKRA-it\\Matching HHs with Remittances VS HHs without Remittances")
merged_data_3 = data.frame(MergedData)
names(merged_data_3)
getwd()
names(merged_data_3) = gsub("[.]", "_", names(merged_data_3))
new_path_3 = file.path(getwd(), "mergedhhrVShhwr_withoutWI.dta")

write_dta(merged_data_3, new_path_3)


#######################################
#######################################
#######################################


#work in merged data
#setwd("D:\\anano\\Policy Briefs\\migration\\migration moneytot DKRA-it\\Matching HHs with Remittances VS HHs without Remittances")

MergedData1<-read_dta("mergedhhrVShhwr_withoutWI.dta")

MergedData1$RELCONDrec<-as.factor(MergedData1$RELCONDrec)
MergedData1$SUBSTRATUM<-as_factor(MergedData1$SUBSTRATUM)
MergedData1$ETHNIC<-as_factor(MergedData1$ETHNIC)
MergedData1$RESPEDU<-as_factor(MergedData1$RESPEDU)
MergedData1$ECONSTNrec<-as_factor(MergedData1$ECONSTNrec)
MergedData1$MONYTOTcontext<-as_factor(MergedData1$MONYTOTcontext)
MergedData1$SPENDMOrec<-as_factor(MergedData1$SPENDMOrec)
MergedData1$RATEHAPrec<-as_factor(MergedData1$RATEHAPrec1)
MergedData1$RELCONDrec<-factor(MergedData1$RELCONDrec, levels = c(1,2,3), labels = c("Relatively poor", "Fair", "Relatively good"))
str(MergedData1$WithoutMWithR)
MergedData1$WithoutMWithR<-as_factor(MergedData1$WithoutMWithR)
MergedData1$RATEHAPrec1<-as_factor(MergedData1$RATEHAPrec1)



###################################
#checkbalance
###################################
#MergedData1$WithoutMWithR<-as.numeric(MergedData1$WithoutMWithR)
#checkbalance <- MatchBalance(MergedData1$WithoutMWithR ~ MergedData1$SUBSTRATUM
   #                          + MergedData1$ETHNIC + MergedData1$HHSIZE 
    #                         + MergedData1$RESPEDU + MergedData1$shareFemale 
     #                        + MergedData1$prop_score, 
      #                       nboots = 1000)
#str(MergedData1$WithoutMWithR)
###################################
#checkbalance
###################################

#model 1 - C14 real economic condition
model1 <- polr(MergedData1$RELCONDrec~MergedData1$WithoutMWithR, w=MergedData1$mweights, data=MergedData1)
table(MergedData1$RELCONDrec)
table(MergedData1$WithoutMWithR)

#RELCONDrec 1 "Relatively poor", 2 "Fair", 3 "Relatively good"
summary(model1)
#Marginal effects and estimated marginal means from regression models
# ggeffects package computes predicted marginal means for the response
plot(ggemmeans(model1, "WithoutMWithR"))
ggemmeans(model1, "WithoutMWithR")



(ci <- confint(model1))
exp(cbind(coef(model1),t(ci)))
control$mweights
mout$weights




#model 2 - hh income
##1 "Up to USD 100", 2 "USD 101-400",3 "More than USD 400", 98 "DKRA"
model2 <- polr(MergedData1$MONYTOTcontext~MergedData1$WithoutMWithR, 
               w=MergedData1$mweights, data=MergedData1)
table(MergedData1$MONYTOTcontext)
table(MergedData1$WithoutMWithR)
summary(model2)
#Marginal effects and estimated marginal means from regression models
# ggeffects package computes predicted marginal means for the response
plot(ggemmeans(model2, "WithoutMWithR"))
ggemmeans(model2, "WithoutMWithR")

(ci <- confint(model2))
exp(cbind(coef(model2),t(ci)))


#model 3 - hh spending
##1 "Spending-Up to USD 100", 2 "Spending-USD 101-400",3 "Spending-More than USD 400"
model3 <- polr(MergedData1$SPENDMOrec~MergedData1$WithoutMWithR, 
               w=MergedData1$mweights, data=MergedData1)
table(MergedData1$SPENDMOrec)
table(MergedData1$WithoutMWithR)
summary(model3)
#Marginal effects and estimated marginal means from regression models
# ggeffects package computes predicted marginal means for the response
plot(ggemmeans(model3, "WithoutMWithR"))
ggemmeans(model3, "WithoutMWithR")

(ci <- confint(model3))
exp(cbind(coef(model3),t(ci)))
control$mweights
mout$weights

#model 4 - C1 economic situation
##4 "Can afford to buy expensive durables",3 "Enough for food/clothes",
##2 "Enough for food only", 1 "Not enough for food"

model4 <- polr(MergedData1$ECONSTNrec~MergedData1$WithoutMWithR, 
               w=MergedData1$mweights, data=MergedData1)
table(MergedData1$ECONSTNrec)
table(MergedData1$WithoutMWithR)
summary(model4)
#Marginal effects and estimated marginal means from regression models
# ggeffects package computes predicted marginal means for the response
plot(ggemmeans(model4, "WithoutMWithR"))
ggemmeans(model4, "WithoutMWithR")

(ci <- confint(model4))
exp(cbind(coef(model4),t(ci)))

#model 5 - C2 happiness
#happiness

model5 <- polr(MergedData1$RATEHAPrec1~MergedData1$WithoutMWithR, 
               w=MergedData1$mweights, data=MergedData1)
table(MergedData1$RATEHAPrec1)
table(MergedData1$WithoutMWithR)
summary(model5)
#Marginal effects and estimated marginal means from regression models
# ggeffects package computes predicted marginal means for the response
plot(ggemmeans(model5, "WithoutMWithR"))
ggemmeans(model5, "WithoutMWithR")


(ci <- confint(model5))
exp(cbind(coef(model5),t(ci)))

#Happiness
MergedData1$RATEHAPrec2 <- 0

MergedData1[MergedData1$RATEHAPrec1 %in% c(1,2), "RATEHAPrec2"] <- 1
MergedData1[MergedData1$RATEHAPrec1 %in% c(3), "RATEHAPrec2"] <- 2
MergedData1[MergedData1$RATEHAPrec1 %in% c(4,5), "RATEHAPrec2"] <- 3

MergedData1$RATEHAPrec2
table(MergedData1$RATEHAPrec2)
MergedData1$RATEHAPrec2<-as_factor(MergedData1$RATEHAPrec2)
###################################



model6 <- polr(MergedData1$RATEHAPrec2~MergedData1$WithoutMWithR, 
               w=MergedData1$mweights, data=MergedData1)
table(MergedData1$RATEHAPrec2)
table(MergedData1$WithoutMWithR)
summary(model6)
#Marginal effects and estimated marginal means from regression models
# ggeffects package computes predicted marginal means for the response
plot(ggemmeans(model6, "WithoutMWithR"))
ggemmeans(model6, "WithoutMWithR")


(ci <- confint(model6))
exp(cbind(coef(model6),t(ci)))



model7 <- lm(MergedData1$wealth_index~MergedData1$WithoutMWithR, 
             w=MergedData1$mweights, data=MergedData1)
table(MergedData1$wealth_index)
table(MergedData1$WithoutMWithR)
summary(model7)
#Marginal effects and estimated marginal means from regression models
# ggeffects package computes predicted marginal means for the response
plot(ggemmeans(model7, "WithoutMWithR"))
ggemmeans(model7, "WithoutMWithR")




######################################################
#Matching migrant HHs with remittances VS non-migrant HHs
######################################################
title: "Migration Matching analysis code"

library(Matching)
library(MASS)
library(RTextTools)
library(haven)
library(descr)
library(Matching)
library(rgenoud)
library(psych)
library(sandwich)
library(lmtest)
library(multiwayvcov)
library(survey)
library(margins)
library(ggeffects)
library(sjmisc)
library(ggplot2)
library(effects)
library(MASS)
library(ggeffects)
library(emmeans)
library(erer)
library(sjPlot)
library(dplyr)
library(janitor)
library(dplyr)
library(ggplot2)


rm(list = ls())
getwd()
setwd("D:\\anano\\Policy Briefs\\migration\\migration moneytot DKRA-it")
CB2019<-read_dta("CB2019_Georgia_response_30Jan2020 - UPDATED.dta")
designgeoun <- svydesign(id=~PSU,weights=~HHWT, strat=~SUBSTRATUM, data=CB2019)
names(designgeoun)


sum(is.na(CB2019$MRVSNM))
sum(is.na(CB2019$SUBSTRATUM))
sum(is.na(CB2019$ETHNIC))
sum(is.na(CB2019$HHSIZE))
sum(is.na(CB2019$RELCONDrec))
sum(is.na(CB2019$RESPEDU))
sum(is.na(CB2019$shareFemale))
sum(is.na(CB2019$ECONSTNrec))
sum(is.na(CB2019$RATEHAPrec))
cols2=c('SUBSTRATUM', 'ETHNIC', 'HHSIZE','RESPEDU', 'shareFemale')
sum(is.na(CB2019[,cols2]))

#MRVSNM - 1 is with remittances, 0 is without remittances
CB2019_cln <- CB2019[!is.na(CB2019$MRVSNM),]
dim(CB2019)
dim(CB2019_cln)
CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$RATEHAPrec),]

CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$ETHNIC),]
dim(CB2019_cln)

CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$RELCONDrec),]
dim(CB2019_cln)

CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$RESPEDU),]
dim(CB2019_cln)

#economic situation
CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$ECONSTNrec),]
dim(CB2019_cln)


CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$wealth_index),]
dim(CB2019_cln)


#Estimand --	A character string for the estimand. 
##The default estimand is "ATT", the sample average treatment effect for the treated. 
#"ATE" is the sample average treatment effect (for all), and 
#"ATC" is the sample average treatment effect for the controls.
#M is a scalar for the number of matches which should be found. The default is one-to-one matching. Also see the ties option.

#Tr=CB2019_cln$MRVSNM - 1 if "With remittances from emigrant", 0 if "without"
CB2019_cln$SUBSTRATUM<-as_factor(CB2019_cln$SUBSTRATUM)
CB2019_cln$ETHNIC<-as_factor(CB2019_cln$ETHNIC)
CB2019_cln$RESPEDU<-as_factor(CB2019_cln$RESPEDU)
CB2019_cln$ECONSTNrec<-as_factor(CB2019_cln$ECONSTNrec)
CB2019_cln$MONYTOTcontextrec<-as_factor(CB2019_cln$MONYTOTcontextrec)
CB2019_cln$SPENDMOrec<-as_factor(CB2019_cln$SPENDMOrec)
CB2019_cln$RATEHAPrec<-as_factor(CB2019_cln$RATEHAPrec1)
CB2019_cln$RATEHAPrec<-as_factor(CB2019_cln$wealth_index)


psmodel<-glm(CB2019_cln$MRVSNM~ CB2019_cln$SUBSTRATUM
             + CB2019_cln$ETHNIC + CB2019_cln$HHSIZE 
             + CB2019_cln$RESPEDU + CB2019_cln$shareFemale, family = "binomial")
hist(psmodel$fitted.values)
CB2019_cln$prop.score<-psmodel$fitted.values


X <- cbind(CB2019_cln$SUBSTRATUM, CB2019_cln$ETHNIC,
           CB2019_cln$RESPEDU, CB2019_cln$HHSIZE,
           CB2019_cln$shareFemale,  
           CB2019_cln$prop.score) 
BalanceMat <- cbind(CB2019_cln$SUBSTRATUM, CB2019_cln$ETHNIC,
                    CB2019_cln$RESPEDU, CB2019_cln$HHSIZE,
                    CB2019_cln$shareFemale,  
                    CB2019_cln$prop.score)
dim(BalanceMat)

#Now we are looking for optimal weights
genout <- GenMatch(Tr=CB2019_cln$MRVSNM, X=X,
                   BalanceMatrix=BalanceMat, estimand = "ATE", 
                   pop.size = 1000, wait.generations = 1)

#Matching
mout<- Match(Tr=CB2019_cln$MRVSNM, X=X, 
             Weight.matrix = genout)
summary(mout)


#Determining if balance has actually obtained or not on th variables of interst
mb<- MatchBalance(CB2019_cln$MRVSNM ~ CB2019_cln$SUBSTRATUM
                  + CB2019_cln$ETHNIC + CB2019_cln$HHSIZE 
                  + CB2019_cln$RESPEDU + CB2019_cln$shareFemale
                  +CB2019_cln$prop.score, 
                  match.out = mout, nboots = 1000)
summary(mb)
treated<-CB2019_cln[mout$index.treated,]
treated$mweights<-mout$weights
treated$mweights

control<-CB2019_cln[mout$index.control,]
control$mweights<-mout$weights
control$mweights

MergedData <- merge(treated, control, all=TRUE)

glm(MergedData$RELCONDrec~MergedData$MRVSNM,w=MergedData$mweights)
str(MergedData$MRVSNM)

#######################################
#######################################

##saving merged data as dta file 
setwd("D:\\anano\\Policy Briefs\\migration\\migration moneytot DKRA-it\\Matching Migrant HHs with Remittances VS non-migrant HHs")
merged_data_2 = data.frame(MergedData)
names(merged_data_2)
getwd()
names(merged_data_2) = gsub("[.]", "_", names(merged_data_2))
new_path_2 = file.path(getwd(), "mergedMhhrVSNMhh_withoutWI.dta")

write_dta(merged_data_2, new_path_2)


#######################################
#######################################
#######################################

MergedData2<-read_dta("mergedMhhrVSNMhh_withoutWI.dta")

table(MergedData2$MRVSNM)
MergedData2$MRVSNM<-as.factor(MergedData2$MRVSNM)

MergedData2$RELCONDrec<-as.factor(MergedData2$RELCONDrec)
MergedData2$SUBSTRATUM<-as_factor(MergedData2$SUBSTRATUM)
MergedData2$ETHNIC<-as_factor(MergedData2$ETHNIC)
MergedData2$RESPEDU<-as_factor(MergedData2$RESPEDU)
MergedData2$ECONSTNrec<-as_factor(MergedData2$ECONSTNrec)
MergedData2$MONYTOTcontext<-as_factor(MergedData2$MONYTOTcontext)
MergedData2$SPENDMOrec<-as_factor(MergedData2$SPENDMOrec)
MergedData2$RATEHAPrec<-as_factor(MergedData2$RATEHAPrec1)

MergedData2$RELCONDrec<-factor(MergedData2$RELCONDrec, levels = c(1,2,3), labels = c("Relatively poor", "Fair", "Relatively good"))
str(MergedData2$WithoutMWithR)
MergedData2$WithoutMWithR<-as_factor(MergedData2$WithoutMWithR)
MergedData2$RATEHAPrec1<-as_factor(MergedData2$RATEHAPrec1)
MergedData2 <- MergedData2[!is.na(MergedData2$wealth_index),]




###################################
#checkbalance
###################################
#MergedData2$MRVSNM<-as.numeric(MergedData2$MRVSNM)

#checkbalance <- MatchBalance(MergedData2$MRVSNM ~ MergedData2$SUBSTRATUM
 #                            + MergedData2$ETHNIC + MergedData2$HHSIZE 
  #                           + MergedData2$RESPEDU + MergedData2$shareFemale 
   #                          + MergedData2$prop_score, 
    #                         nboots = 1000)

###################################
#checkbalance
###################################

#model 1 - C14 real economic condition
model1 <- polr(MergedData2$RELCONDrec~MergedData2$MRVSNM, w=MergedData2$mweights, data=MergedData2)
table(MergedData2$RELCONDrec)
table(MergedData2$MRVSNM)

margs<-ggemmeans(model1, terms = c("MRVSNM"))
plot(margs)
margs


#RELCONDrec 1 "Relatively poor", 2 "Fair", 3 "Relatively good"
summary(model1)
#Marginal effects and estimated marginal means from regression models
# ggeffects package computes predicted marginal means for the response
plot(ggemmeans(model1, "MRVSNM"))
ggemmeans(model1, "MRVSNM")


(ci <- confint(model1))
exp(cbind(coef(model1),t(ci)))
control$mweights
mout$weights




#model 2 - hh income
##1 "Up to USD 100", 2 "USD 101-400",3 "More than USD 400"
model2 <- polr(MergedData2$MONYTOTcontext~MergedData2$MRVSNM, 
               w=MergedData2$mweights, data=MergedData2)
table(MergedData2$MONYTOTcontext)
table(MergedData2$MRVSNM)
summary(model2)
#Marginal effects and estimated marginal means from regression models
# ggeffects package computes predicted marginal means for the response
plot(ggemmeans(model2, "MRVSNM"))
ggemmeans(model2, "MRVSNM")


(ci <- confint(model2))
exp(cbind(coef(model2),t(ci)))


#model 3 - hh spending
##1 "Spending-Up to USD 100", 2 "Spending-USD 101-400",3 "Spending-More than USD 400"
model3 <- polr(MergedData2$SPENDMOrec~MergedData2$MRVSNM, 
               w=MergedData2$mweights, data=MergedData2)
table(MergedData2$SPENDMOrec)
table(MergedData2$MRVSNM)
summary(model3)
#Marginal effects and estimated marginal means from regression models
# ggeffects package computes predicted marginal means for the response
plot(ggemmeans(model3, "MRVSNM"))
ggemmeans(model3, "MRVSNM")


(ci <- confint(model3))
exp(cbind(coef(model3),t(ci)))
control$mweights
mout$weights



#model 4 - C1 economic situation
##4 "Can afford to buy expensive durables",3 "Enough for food/clothes",
##2 "Enough for food only", 1 "Not enough for food"

model4 <- polr(MergedData2$ECONSTNrec~MergedData2$MRVSNM, 
               w=MergedData2$mweights, data=MergedData2)
table(MergedData2$ECONSTNrec)
table(MergedData2$MRVSNM)
summary(model4)
#Marginal effects and estimated marginal means from regression models
# ggeffects package computes predicted marginal means for the response
plot(ggemmeans(model4, "MRVSNM"))
ggemmeans(model4, "MRVSNM")



(ci <- confint(model4))
exp(cbind(coef(model4),t(ci)))
control$mweights
mout$weights
ctable <- coef(summary(model4))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE)*2
(ctable <- cbind(ctable, "p value"=p))



model5 <- polr(MergedData2$RATEHAPrec1~MergedData2$MRVSNM, 
               w=MergedData2$mweights, data=MergedData2)
table(MergedData2$RATEHAPrec1)
table(MergedData2$MRVSNM)
summary(model5)
#Marginal effects and estimated marginal means from regression models
# ggeffects package computes predicted marginal means for the response
plot(ggemmeans(model5, "MRVSNM"))
ggemmeans(model5, "MRVSNM")


(ci <- confint(model5))
exp(cbind(coef(model5),t(ci)))




model6 <- lm(MergedData2$wealth_index~MergedData2$MRVSNM, 
             w=MergedData2$mweights, data=MergedData2)
table(MergedData2$wealth_index)
table(MergedData2$MRVSNM)
summary(model6)
#Marginal effects and estimated marginal means from regression models
# ggeffects package computes predicted marginal means for the response
plot(ggemmeans(model6, "MRVSNM"))
ggemmeans(model6, "MRVSNM")





######################################################
#Matching migrant HHs with Remmitances VS migrant HHs without Remittances
######################################################
title: "Migration Matching analysis code"

library(Matching)
library(MASS)
library(RTextTools)
library(haven)
library(descr)
library(Matching)
library(rgenoud)
library(psych)
library(sandwich)
library(lmtest)
library(multiwayvcov)
library(survey)
library(margins)
library(ggeffects)
library(sjmisc)
library(ggplot2)
library(effects)
library(MASS)
library(ggeffects)
library(emmeans)

rm(list = ls())
getwd()
setwd("D:\\anano\\Policy Briefs\\migration\\migration moneytot DKRA-it")
CB2019<-read_dta("CB2019_Georgia_response_30Jan2020 - UPDATED.dta")
designgeoun <- svydesign(id=~PSU,weights=~HHWT, strat=~SUBSTRATUM, data=CB2019)
names(designgeoun)


sum(is.na(CB2019$MRVSMWR))
sum(is.na(CB2019$SUBSTRATUM))
sum(is.na(CB2019$ETHNIC))
sum(is.na(CB2019$HHSIZE))
sum(is.na(CB2019$RELCONDrec))
sum(is.na(CB2019$RESPEDU))
sum(is.na(CB2019$shareFemale))
sum(is.na(CB2019$ECONSTNrec))
sum(is.na(CB2019$RATEHAPrec))
cols2=c('SUBSTRATUM', 'ETHNIC', 'HHSIZE','RESPEDU', 'shareFemale')
sum(is.na(CB2019[,cols2]))

#MRVSMWR - 1 is with remittances, 0 is without remittances
CB2019_cln <- CB2019[!is.na(CB2019$MRVSMWR),]
dim(CB2019)
dim(CB2019_cln)
CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$RATEHAPrec),]

CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$ETHNIC),]
dim(CB2019_cln)

CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$RELCONDrec),]
dim(CB2019_cln)

CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$RESPEDU),]
dim(CB2019_cln)

#economic situation
CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$ECONSTNrec),]
dim(CB2019_cln)

CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$wealth_index),]
dim(CB2019_cln)



#Estimand --	A character string for the estimand. 
##The default estimand is "ATT", the sample average treatment effect for the treated. 
#"ATE" is the sample average treatment effect (for all), and 
#"ATC" is the sample average treatment effect for the controls.
#M is a scalar for the number of matches which should be found. The default is one-to-one matching. Also see the ties option.

#Tr=CB2019_cln$MRVSMWR - 1 if "With remittances from emigrant", 0 if "without"
CB2019_cln$SUBSTRATUM<-as_factor(CB2019_cln$SUBSTRATUM)
CB2019_cln$ETHNIC<-as_factor(CB2019_cln$ETHNIC)
CB2019_cln$RESPEDU<-as_factor(CB2019_cln$RESPEDU)
CB2019_cln$ECONSTNrec<-as_factor(CB2019_cln$ECONSTNrec)
CB2019_cln$MONYTOTcontext<-as_factor(CB2019_cln$MONYTOTcontext)
CB2019_cln$SPENDMOrec<-as_factor(CB2019_cln$SPENDMOrec)
CB2019_cln$RATEHAPrec<-as_factor(CB2019_cln$RATEHAPrec1)
CB2019_cln$RATEHAPrec<-as_factor(CB2019_cln$wealth_index)

psmodel<-glm(CB2019_cln$MRVSMWR~ CB2019_cln$SUBSTRATUM
             + CB2019_cln$ETHNIC + CB2019_cln$HHSIZE 
             + CB2019_cln$RESPEDU + CB2019_cln$shareFemale, family = "binomial")
hist(psmodel$fitted.values)
CB2019_cln$prop.score<-psmodel$fitted.values


X <- cbind(CB2019_cln$SUBSTRATUM, CB2019_cln$ETHNIC,
           CB2019_cln$RESPEDU, CB2019_cln$HHSIZE,
           CB2019_cln$shareFemale, CB2019_cln$prop.score) 
BalanceMat <- cbind(CB2019_cln$SUBSTRATUM, CB2019_cln$ETHNIC,
                    CB2019_cln$RESPEDU, CB2019_cln$HHSIZE,
                    CB2019_cln$shareFemale, CB2019_cln$prop.score)
dim(BalanceMat)

#Now we are looking for optimal weights
genout <- GenMatch(Tr=CB2019_cln$MRVSMWR, X=X,
                   BalanceMatrix=BalanceMat, estimand = "ATE", 
                   pop.size = 1000, wait.generations = 1)

#Matching
mout<- Match(Tr=CB2019_cln$MRVSMWR, X=X, 
             Weight.matrix = genout)
summary(mout)


#Determining if balance has actually obtained or not on th variables of interst
mb<- MatchBalance(CB2019_cln$MRVSMWR ~ CB2019_cln$SUBSTRATUM
                  + CB2019_cln$ETHNIC + CB2019_cln$HHSIZE 
                  + CB2019_cln$RESPEDU + CB2019_cln$shareFemale
                  +CB2019_cln$prop.score, 
                  match.out = mout, nboots = 1000)

summary(mb)
treated<-CB2019_cln[mout$index.treated,]
treated$mweights<-mout$weights
treated$mweights

control<-CB2019_cln[mout$index.control,]
control$mweights<-mout$weights
control$mweights

MergedData <- merge(treated, control, all=TRUE)
glm(MergedData$RELCONDrec~MergedData$MRVSMWR,w=MergedData$mweights)



#######################################
#######################################

##saving merged data as dta file 
setwd("D:\\anano\\Policy Briefs\\migration\\migration moneytot DKRA-it\\Matching HHs with Remittances VS migrant HHs without Remittances")
merged_data = data.frame(MergedData)
names(merged_data)
getwd()
names(merged_data) = gsub("[.]", "_", names(merged_data))
new_path = file.path(getwd(), "mergedhhrVSmhhwr_withoutWI.dta")

write_dta(merged_data, new_path)


#######################################
#######################################
#######################################

MergedData3<-read_dta("mergedhhrVSmhhwr_withoutWI.dta")

MergedData3$MRVSMWR<-as.factor(MergedData3$MRVSMWR)
MergedData3$RELCONDrec<-as.factor(MergedData3$RELCONDrec)
MergedData3$SUBSTRATUM<-as_factor(MergedData3$SUBSTRATUM)
MergedData3$ETHNIC<-as_factor(MergedData3$ETHNIC)
MergedData3$RESPEDU<-as_factor(MergedData3$RESPEDU)
MergedData3$ECONSTNrec<-as_factor(MergedData3$ECONSTNrec)
MergedData3$MONYTOTcontext<-as_factor(MergedData3$MONYTOTcontext)
MergedData3$SPENDMOrec<-as_factor(MergedData3$SPENDMOrec)
MergedData3$RATEHAPrec<-as_factor(MergedData3$RATEHAPrec1)

MergedData3$RELCONDrec<-factor(MergedData3$RELCONDrec, levels = c(1,2,3), labels = c("Relatively poor", "Fair", "Relatively good"))
str(MergedData3$WithoutMWithR)
MergedData3$WithoutMWithR<-as_factor(MergedData3$WithoutMWithR)
MergedData3$RATEHAPrec1<-as_factor(MergedData3$RATEHAPrec1)




###################################
#checkbalance
###################################
#MergedData3$MRVSMWR<-as.numeric(MergedData3$MRVSMWR)

#checkbalance <- MatchBalance(MergedData3$MRVSMWR ~ MergedData3$SUBSTRATUM
 #                            + MergedData3$ETHNIC + MergedData3$HHSIZE 
  #                           + MergedData3$RESPEDU + MergedData3$shareFemale 
   #                          + MergedData3$prop_score, 
    #                         nboots = 1000)

###################################
#checkbalance
###################################


#model 1 - C14 real economic condition
sum(is.na(MergedData3$RELCONDrec))
model1 <- polr(MergedData3$RELCONDrec~MergedData3$MRVSMWR, w=MergedData3$mweights, data=MergedData3)
table(MergedData3$RELCONDrec)
table(MergedData3$MRVSMWR)

#RELCONDrec 1 "Relatively poor", 2 "Fair", 3 "Relatively good"
summary(model1)
#Marginal effects and estimated marginal means from regression models
# ggeffects package computes predicted marginal means for the response
plot(ggemmeans(model1, "MRVSMWR"))
ggemmeans(model1, "MRVSMWR")

(ci <- confint(model1))
exp(cbind(coef(model1),t(ci)))
control$mweights
mout$weights

ctable1 <- coef(summary(model1))
p <- pnorm(abs(ctable1[, "t value"]), lower.tail = FALSE)*2
(ctable1 <- cbind(ctable1, "p value"=p))


#model 2 - hh income
##1 "Up to USD 100", 2 "USD 101-400",3 "More than USD 400"
model2 <- polr(MergedData3$MONYTOTcontext~MergedData3$MRVSMWR, 
               w=MergedData3$mweights, data=MergedData3)
table(MergedData3$MONYTOTcontext)
table(MergedData3$MRVSMWR)
summary(model2)
#Marginal effects and estimated marginal means from regression models
# ggeffects package computes predicted marginal means for the response
plot(ggemmeans(model2, "MRVSMWR"))
ggemmeans(model2, "MRVSMWR")

(ci <- confint(model2))
exp(cbind(coef(model2),t(ci)))

ctable <- coef(summary(model2))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE)*2
(ctable <- cbind(ctable, "p value"=p))
ctable

#model 3 - hh spending
##1 "Spending-Up to USD 100", 2 "Spending-USD 101-400",3 "Spending-More than USD 400"
model3 <- polr(MergedData3$SPENDMOrec~MergedData3$MRVSMWR, 
               w=MergedData3$mweights, data=MergedData3)
table(MergedData3$SPENDMOrec)
table(MergedData3$MRVSMWR)
summary(model3)
#Marginal effects and estimated marginal means from regression models
# ggeffects package computes predicted marginal means for the response
plot(ggemmeans(model3, "MRVSMWR"))
ggemmeans(model3, "MRVSMWR")

(ci <- confint(model3))
exp(cbind(coef(model3),t(ci)))
control$mweights
mout$weights

#model 4 - C1 economic situation
##4 "Can afford to buy expensive durables",3 "Enough for food/clothes",
##2 "Enough for food only", 1 "Not enough for food"


model4 <- polr(MergedData3$ECONSTNrec~MergedData3$MRVSMWR, 
               w=MergedData3$mweights, data=MergedData3)
table(MergedData3$ECONSTNrec)
table(MergedData3$MRVSMWR)
summary(model4)
#Marginal effects and estimated marginal means from regression models
# ggeffects package computes predicted marginal means for the response
plot(ggemmeans(model4, "MRVSMWR"))
ggemmeans(model4, "MRVSMWR")



(ci <- confint(model4))
exp(cbind(coef(model4),t(ci)))


model5 <- polr(MergedData3$RATEHAPrec1~MergedData3$MRVSMWR, 
               w=MergedData3$mweights, data=MergedData3)
table(MergedData3$RATEHAPrec1)
table(MergedData3$MRVSMWR)
summary(model5)
#Marginal effects and estimated marginal means from regression models
# ggeffects package computes predicted marginal means for the response
plot(ggemmeans(model5, "MRVSMWR"))
ggemmeans(model5, "MRVSMWR")


#wealth index
################################
#wealth index3
MergedData3$wealth_index2 <- 0

MergedData3[MergedData3$wealth_index %in% c(1,2,3), "wealth_index2"] <- 1
MergedData3[MergedData3$wealth_index %in% c(4,5,6,7), "wealth_index2"] <- 2
MergedData3[MergedData3$wealth_index %in% c(8,9,10), "wealth_index2"] <- 3

MergedData3$wealth_index2
table(MergedData3$wealth_index2)
MergedData3$wealth_index2<-as_factor(MergedData3$wealth_index2)
###################################


#wealth index-OLS
model6 <- lm(MergedData3$wealth_index~MergedData3$MRVSMWR, 
             w=MergedData3$mweights, data=MergedData3)
table(MergedData3$wealth_index)
table(MergedData3$MRVSMWR)
summary(model6)
#Marginal effects and estimated marginal means from regression models
# ggeffects package computes predicted marginal means for the response
plot(ggemmeans(model6, "MRVSMWR"))
ggemmeans(model6, "MRVSMWR")





######################################################
#Migrant HHs without remittances to non-migrant HHs
######################################################
title: "Migration Matching analysis code"

library(Matching)
library(MASS)
library(RTextTools)
library(haven)
library(descr)
library(Matching)
library(rgenoud)
library(psych)
library(sandwich)
library(lmtest)
library(multiwayvcov)
library(survey)
library(margins)
library(ggeffects)
library(sjmisc)
library(ggplot2)
library(effects)
library(MASS)
library(sjPlot)
library(emmeans)

rm(list = ls())
getwd()
setwd("D:\\anano\\Policy Briefs\\migration\\migration moneytot DKRA-it")
CB2019<-read_dta("CB2019_Georgia_response_30Jan2020 - UPDATED.dta")
designgeoun <- svydesign(id=~PSU,weights=~HHWT, strat=~SUBSTRATUM, data=CB2019)
names(designgeoun)


sum(is.na(CB2019$MWRVSNM))
sum(is.na(CB2019$SUBSTRATUM))
sum(is.na(CB2019$ETHNIC))
sum(is.na(CB2019$HHSIZE))
sum(is.na(CB2019$RELCONDrec))
sum(is.na(CB2019$RESPEDU))
sum(is.na(CB2019$shareFemale))
sum(is.na(CB2019$ECONSTNrec))
sum(is.na(CB2019$RATEHAPrec))
cols2=c('SUBSTRATUM', 'ETHNIC', 'HHSIZE','RESPEDU', 'shareFemale')
sum(is.na(CB2019[,cols2]))

#MWRVSNM - 1 is with remittances, 0 is without remittances
CB2019_cln <- CB2019[!is.na(CB2019$MWRVSNM),]
dim(CB2019)
dim(CB2019_cln)
CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$RATEHAPrec),]

CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$ETHNIC),]
dim(CB2019_cln)

CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$RELCONDrec),]
dim(CB2019_cln)

CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$RESPEDU),]
dim(CB2019_cln)

#economic situation
CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$ECONSTNrec),]
dim(CB2019_cln)

CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$wealth_index),]
dim(CB2019_cln)


#Estimand --	A character string for the estimand. 
##The default estimand is "ATT", the sample average treatment effect for the treated. 
#"ATE" is the sample average treatment effect (for all), and 
#"ATC" is the sample average treatment effect for the controls.
#M is a scalar for the number of matches which should be found. The default is one-to-one matching. Also see the ties option.

#Tr=CB2019_cln$MWRVSNM - 1 if "With remittances from emigrant", 0 if "without"
CB2019_cln$SUBSTRATUM<-as_factor(CB2019_cln$SUBSTRATUM)
CB2019_cln$ETHNIC<-as_factor(CB2019_cln$ETHNIC)
CB2019_cln$RESPEDU<-as_factor(CB2019_cln$RESPEDU)
CB2019_cln$ECONSTNrec<-as_factor(CB2019_cln$ECONSTNrec)
CB2019_cln$MONYTOTcontext<-as_factor(CB2019_cln$MONYTOTcontext)
CB2019_cln$SPENDMOrec<-as_factor(CB2019_cln$SPENDMOrec)
CB2019_cln$RATEHAPrec<-as_factor(CB2019_cln$RATEHAPrec1)
CB2019_cln$RATEHAPrec<-as_factor(CB2019_cln$wealth_index)


psmodel<-glm(CB2019_cln$MWRVSNM~ CB2019_cln$SUBSTRATUM
             + CB2019_cln$ETHNIC + CB2019_cln$HHSIZE 
             + CB2019_cln$RESPEDU + CB2019_cln$shareFemale, family = "binomial")
hist(psmodel$fitted.values)
CB2019_cln$prop.score<-psmodel$fitted.values


X <- cbind(CB2019_cln$SUBSTRATUM, CB2019_cln$ETHNIC,
           CB2019_cln$RESPEDU, CB2019_cln$HHSIZE,
           CB2019_cln$shareFemale,  CB2019_cln$prop.score) 
BalanceMat <- cbind(CB2019_cln$SUBSTRATUM, CB2019_cln$ETHNIC,
                    CB2019_cln$RESPEDU, CB2019_cln$HHSIZE,
                    CB2019_cln$shareFemale, CB2019_cln$prop.score)
dim(BalanceMat)

#Now we are looking for optimal weights
genout <- GenMatch(Tr=CB2019_cln$MWRVSNM, X=X,
                   BalanceMatrix=BalanceMat, estimand = "ATE", 
                   pop.size = 1000, wait.generations = 1)

#Matching
mout<- Match(Tr=CB2019_cln$MWRVSNM, X=X, 
             Weight.matrix = genout)
summary(mout)


#Determining if balance has actually obtained or not on th variables of interst
mb<- MatchBalance(CB2019_cln$MWRVSNM ~ CB2019_cln$SUBSTRATUM
                  + CB2019_cln$ETHNIC + CB2019_cln$HHSIZE 
                  + CB2019_cln$RESPEDU + CB2019_cln$shareFemale
                  +CB2019_cln$prop.score, 
                  match.out = mout, nboots = 1000)
summary(mb)
treated<-CB2019_cln[mout$index.treated,]
treated$mweights<-mout$weights
treated$mweights

control<-CB2019_cln[mout$index.control,]
control$mweights<-mout$weights
control$mweights

MergedData <- merge(treated, control, all=TRUE)


glm(MergedData$RELCONDrec~MergedData$MWRVSNM,w=MergedData$mweights)
#######################################
#######################################

##saving merged data as dta file 
setwd("D:\\anano\\Policy Briefs\\migration\\migration moneytot DKRA-it\\Matching Migrant HHs without Remittances VS non-migrant HHs")
merged_data = data.frame(MergedData)
names(merged_data)
getwd()
names(merged_data) = gsub("[.]", "_", names(merged_data))
new_path = file.path(getwd(), "mergedmhhwrVSnmhhwr_withoutWI.dta")

write_dta(merged_data, new_path)


#######################################
#######################################
#######################################

MergedData4<-read_dta("mergedmhhwrVSnmhhwr_withoutWI.dta")

MergedData4$RELCONDrec<-as.factor(MergedData4$RELCONDrec)
MergedData4$SUBSTRATUM<-as_factor(MergedData4$SUBSTRATUM)
MergedData4$ETHNIC<-as_factor(MergedData4$ETHNIC)
MergedData4$RESPEDU<-as_factor(MergedData4$RESPEDU)
MergedData4$ECONSTNrec<-as_factor(MergedData4$ECONSTNrec)
MergedData4$MONYTOTcontext<-as_factor(MergedData4$MONYTOTcontext)
MergedData4$SPENDMOrec<-as_factor(MergedData4$SPENDMOrec)
MergedData4$RATEHAPrec<-as_factor(MergedData4$RATEHAPrec1)

MergedData4$RELCONDrec<-factor(MergedData4$RELCONDrec, levels = c(1,2,3), labels = c("Relatively poor", "Fair", "Relatively good"))
str(MergedData4$WithoutMWithR)
MergedData4$WithoutMWithR<-as_factor(MergedData4$WithoutMWithR)
MergedData4$RATEHAPrec1<-as_factor(MergedData4$RATEHAPrec1)
MergedData4$MWRVSNM<-as_factor(MergedData4$MWRVSNM)


###################################
#checkbalance
###################################
#MergedData4$MWRVSNM<-as.numeric(MergedData4$MWRVSNM)
#checkbalance <- MatchBalance(MergedData4$MWRVSNM ~ MergedData4$SUBSTRATUM
 #                            + MergedData4$ETHNIC + MergedData4$HHSIZE 
  #                           + MergedData4$RESPEDU + MergedData4$shareFemale 
   #                          + MergedData4$prop_score, 
    #                         nboots = 1000)

###################################
#checkbalance
###################################



#model 1 - C14 real economic condition

model1 <- polr(MergedData4$RELCONDrec~MergedData4$MWRVSNM, w=MergedData4$mweights, data=MergedData4)

table(MergedData4$RELCONDrec)
table(MergedData4$MWRVSNM)

#RELCONDrec 1 "Relatively poor", 2 "Fair", 3 "Relatively good"
summary(model1)
#Marginal effects and estimated marginal means from regression models
# ggeffects package computes predicted marginal means for the response
plot(ggemmeans(model1, "MWRVSNM"))

ggemmeans(model1, "MWRVSNM")
(ci <- confint(model1))
exp(cbind(coef(model1),t(ci)))
control$mweights
mout$weights

ctable1 <- coef(summary(model1))
p <- pnorm(abs(ctable1[, "t value"]), lower.tail = FALSE)*2
(ctable1 <- cbind(ctable1, "p value"=p))

plot(ggemmeans(model1, "MWRVSNM"))
##for multinomial logistic
sjPlot::plot_model(model1,
                   # type="eff",
                   title = "kitxvis teqsti",
                   colors = c("seagreen3", "deepskyblue2", "orchid"),
                   show.values = TRUE,
                   value.offset = .4,
                   value.size = 6,
                   dot.size = 4,
                   line.size = 1.5,
                   vline.color = "red",
                   width = 0.5)
+
  font_size(title = 18, labels.y = 14)

#margplot
margs<-ggemmeans(model1, terms = c("MWRVSNM"))
plot(margs)
margs
margplot<-ggplot(margs, aes(x, predicted))+
  geom_point(aes(color=x))+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, color=x), width=0.1, size=1, position=position_dodge(width=0.4))+
  labs(title ="satauri",
       x="x cvladi",
       y="Predicted probability")+
  theme_ggeffects()+theme(legend.position = "none")+
  theme(legend.position="top")
margplot

#model 2 - hh income
##1 "Up to USD 100", 2 "USD 101-400",3 "More than USD 400"
model2 <- polr(MergedData4$MONYTOTcontext~MergedData4$MWRVSNM, 
               w=MergedData4$mweights, data=MergedData4)
table(MergedData4$MONYTOTcontext)
table(MergedData4$MWRVSNM)
summary(model2)
#Marginal effects and estimated marginal means from regression models
# ggeffects package computes predicted marginal means for the response
plot(ggemmeans(model2, "MWRVSNM"))
ggemmeans(model2, "MWRVSNM")
(ci <- confint(model2))
exp(cbind(coef(model2),t(ci)))

ctable <- coef(summary(model2))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE)*2
(ctable <- cbind(ctable, "p value"=p))
ctable

#model 3 - hh spending
##1 "Spending-Up to USD 100", 2 "Spending-USD 101-400",3 "Spending-More than USD 400"
model3 <- polr(MergedData4$SPENDMOrec~MergedData4$MWRVSNM, 
               w=MergedData4$mweights, data=MergedData4)
table(MergedData4$SPENDMOrec)
table(MergedData4$MWRVSNM)
summary(model3)
#Marginal effects and estimated marginal means from regression models
# ggeffects package computes predicted marginal means for the response
plot(ggemmeans(model3, "MWRVSNM"))
ggemmeans(model3, "MWRVSNM")
(ci <- confint(model3))
exp(cbind(coef(model3),t(ci)))
control$mweights
mout$weights

#model 4 - C1 economic situation
##4 "Can afford to buy expensive durables",3 "Enough for food/clothes",
##2 "Enough for food only", 1 "Not enough for food"

model4 <- polr(MergedData4$ECONSTNrec~MergedData4$MWRVSNM, 
               w=MergedData4$mweights, data=MergedData4)
table(MergedData4$ECONSTNrec)
table(MergedData4$MWRVSNM)
summary(model4)
#Marginal effects and estimated marginal means from regression models
# ggeffects package computes predicted marginal means for the response
plot(ggemmeans(model4, "MWRVSNM"))
ggemmeans(model4, "MWRVSNM")
(ci <- confint(model4))
exp(cbind(coef(model4),t(ci)))




#wealth index-OLS
model6 <- lm(MergedData4$wealth_index~MergedData4$MWRVSNM, 
             w=MergedData4$mweights, data=MergedData4)
table(MergedData4$wealth_index)
table(MergedData4$MWRVSNM)
summary(model6)
#Marginal effects and estimated marginal means from regression models
# ggeffects package computes predicted marginal means for the response
plot(ggemmeans(model6, "MWRVSNM"))
ggemmeans(model6, "MWRVSNM")




######################################################
#draft who are more likely to have close relatives abroad and remittances
######################################################
title: "WHO RECEIVES REMITTANCES?"

library(Matching)
library(MASS)
library(RTextTools)
library(haven)
library(descr)
library(Matching)
library(rgenoud)
library(psych)
library(sandwich)
library(lmtest)
library(multiwayvcov)
library(survey)
library(margins)
library(ggeffects)
library(sjmisc)
library(ggplot2)
library(effects)
library(MASS)
library(ggeffects)
library(emmeans)
library(erer)
library(sjPlot)

rm(list = ls())
getwd()
setwd("D:\\anano\\Policy Briefs\\migration\\CB 2019")
CB2019<-read_dta("CB2019_Georgia_response_30Jan2020.dta")
designgeoun <- svydesign(id=~PSU,weights=~HHWT, strat=~SUBSTRATUM, data=CB2019)
names(designgeoun)


sum(is.na(CB2019$WithoutMWithR))
sum(is.na(CB2019$SUBSTRATUM))
sum(is.na(CB2019$ETHNIC))
sum(is.na(CB2019$HHSIZE))
sum(is.na(CB2019$RELCONDrec))
sum(is.na(CB2019$RESPEDU))
sum(is.na(CB2019$shareFemale))
sum(is.na(CB2019$ECONSTNrec))
sum(is.na(CB2019$RATEHAPrec))
cols2=c('SUBSTRATUM', 'ETHNIC', 'HHSIZE','RESPEDU', 'shareFemale')
sum(is.na(CB2019[,cols2]))


#WithoutMWithR - 1 is with remittances, 0 is without remittances
CB2019_cln <- CB2019[!is.na(CB2019$WithoutMWithR),]
dim(CB2019)
dim(CB2019_cln)
CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$RATEHAPrec),]

CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$ETHNIC),]
dim(CB2019_cln)

CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$RELCONDrec),]
dim(CB2019_cln)

CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$RESPEDU),]
dim(CB2019_cln)

#economic situation
CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$ECONSTNrec),]
dim(CB2019_cln)

#close relative living abroad
CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$CLRLABRrec),]
dim(CB2019_cln)


#remittances from abroad
CB2019_cln <- CB2019_cln[!is.na(CB2019_cln$INCSOUABrec),]
dim(CB2019_cln)

#Estimand --	A character string for the estimand. 
##The default estimand is "ATT", the sample average treatment effect for the treated. 
#"ATE" is the sample average treatment effect (for all), and 
#"ATC" is the sample average treatment effect for the controls.
#M is a scalar for the number of matches which should be found. The default is one-to-one matching. Also see the ties option.

#Tr=CB2019_cln$WithoutMWithR - 1 if "With remittances from emigrant", 0 if "without"
CB2019_cln$SUBSTRATUM<-as_factor(CB2019_cln$SUBSTRATUM)
CB2019_cln$ETHNIC<-as_factor(CB2019_cln$ETHNIC)
CB2019_cln$RESPEDU<-as_factor(CB2019_cln$RESPEDU)
CB2019_cln$ECONSTNrec<-as_factor(CB2019_cln$ECONSTNrec)
CB2019_cln$MONYTOTrec<-as_factor(CB2019_cln$MONYTOTrec)
CB2019_cln$SPENDMOrec<-as_factor(CB2019_cln$SPENDMOrec)
CB2019_cln$RATEHAPrec<-as_factor(CB2019_cln$RATEHAPrec1)
CB2019_cln$RATEHAPrec<-as_factor(CB2019_cln$CLRLABRrec)
CB2019_cln$RATEHAPrec<-as_factor(CB2019_cln$INCSOUABrec)


#who are more likely to have close relative currently living abroad


table(CB2019_cln$RATEHAPrec)
table(CB2019_cln$settypeREC)
relatmodelmigrants<-glm(CLRLABRrec~ settypeREC
                        + ETHNICREC + HHASIZE_r 
                        + children_dummy,  data=CB2019_cln)
ggemmeans(relatmodelmigrants, "settypeREC")
ggemmeans(relatmodelmigrants, "ETHNICREC")
ggemmeans(relatmodelmigrants, "HHASIZE_r")
ggemmeans(relatmodelmigrants, "children_dummy")
table(CB2019_cln$ETHNICREC)


relatmodelremitt<-glm(INCSOUABrec~ settypeREC
                      + ETHNICREC + HHASIZE_r 
                      + children_dummy,  data=CB2019_cln)
ggemmeans(relatmodelremitt, "settypeREC")
ggemmeans(relatmodelremitt, "ETHNICREC")
ggemmeans(relatmodelremitt, "HHASIZE_r")
ggemmeans(relatmodelremitt, "children_dummy")
table(CB2019_cln$ETHNICREC)
table(CB2019_cln$HHASIZE_r)

