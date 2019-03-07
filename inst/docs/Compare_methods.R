# Summarize Parameter Effect on Model Accuracy
# Gabriel Odom
# 2018-07-05

# Given the DMR method comparison results in the DMRcompare package, assess the
#   strengths and directions between the model parameters and performance
#   metrics for each of the four methods.
library(DMRcompare)
library(tidyverse)


######  DMRcate  ##############################################################
data("dmrcateRes_df")
colnames(dmrcateRes_df)
# Model components:
#   Predictors: lambda and C
#   Responses : power, AuPR, precision, mcc, and F1
#   Covariates: delta, nCPG_med, nCPG_q3
resultsDMRcate_df <-
  dmrcateRes_df %>%
  select(power, AuPR, precision, mcc, F1,
         lambda, C,
         delta, seed, nCPG_med, nCPG_q3) %>%
  mutate(nCPG_med = as.numeric(nCPG_med)) %>%
  mutate(nCPG_q3 = as.numeric(nCPG_q3))

###  Replace Missing Values with 0  ###
# For each of the responses, 0 is as bad as it gets. Therefore, we can replace
#   all NAs and NaNs with 0. However, precision start at 1 and never reach 0,
#   so we replace all of the 0s for precision with 1.
resultsDMRcate_df <- resultsDMRcate_df[complete.cases(resultsDMRcate_df), ]


###  Correlation Matrix  ###
# First we build a correlation matrix to inspect possible multicollinearity.
round(cor(resultsDMRcate_df), 2)
# We see that all of the responses (except for precision) are nearly perfectly
#   correlated with each other. Precision is rather strongly negatively
#   correlated with the other output variables, but not perfectly (-0.6). This
#   means we can simply model the AuPR instead. Further, it appears that
#   neither C nor lambda are correlated with the outcome at all. These two nCPG
#   values are also highly correlated. Thankfully, the seed has no effect.
resultsDMRcate_df <-
  resultsDMRcate_df %>%
  select(AuPR, lambda, C, delta, nCPG_q3)
cor(resultsDMRcate_df)
# We've removed any multicollinearity issues.
cor(resultsDMRcate_df, method = "spearman")
# Now, the correlation test:
cor.test(resultsDMRcate_df$lambda,
         resultsDMRcate_df$AuPR,
         method = "spearman")$p.value

cor.test(resultsDMRcate_df$C,
         resultsDMRcate_df$AuPR,
         method = "spearman")$p.value
# Both are related.


###  Statistical Prediction for AuPR  ###
summary(resultsDMRcate_df$AuPR)
plot(density(resultsDMRcate_df$AuPR))
# The AuPR values can only range from 0-1, so I think we should be using
#   beta regression if the residuals of OLS do not appear normal.
dmrcateLin_mod <- lm(AuPR ~ ., data = resultsDMRcate_df)
par(mfrow = c(2, 2))
plot(dmrcateLin_mod)
# Well that's a resounding "NO". However, we see that observation 180 is a
#   severe outlier. Let's check it.
dmrcateRes_df[178:180, ]
# I don't see anything wrong with it, but  it's apparently very influential.

# Beta regression it is. The default link is logit.
library(betareg)
resultsDMRcate_df <-
  resultsDMRcate_df %>%
  mutate(sigma = lambda / C)
dmrcateBeta_mod <- betareg(AuPR ~ ., data = resultsDMRcate_df)
plot(dmrcateBeta_mod)
summary(dmrcateBeta_mod)
# The directions of the relationships jive with the Spearman correlations. What
#   we can then say is this: given delta and the number of CPGs, for each unit
#   increase in lambda, logit(AuPR) decreases by 0.18892654; for each unit
#   increase of C, logit(AuPR) increases by 0.25862122. Also, the interaction
#   between C and lambda (sigma) is significant: as lambda / C increases, AuPR
#   decreases.

# Let's have a sanity check:
par(mfrow = c(1,1))
boxplot(resultsDMRcate_df$AuPR ~ resultsDMRcate_df$C)
# I'll buy the results that C is positively related with AuPC, but not that the
#   relationship is logit-linear
boxplot(resultsDMRcate_df$AuPR ~ resultsDMRcate_df$lambda)
# Same as before, but with the negative relationship.
boxplot(resultsDMRcate_df$AuPR ~ resultsDMRcate_df$sigma)
# Smaller lambda / C is better


######  ProbeLasso  #####################################################
data("probeLassoRes_df")
colnames(probeLassoRes_df)
# Model components:
#   Predictors: adjPval, mLassoRad, and minDmrSep
#   Responses : power, AuPR, precision, mcc, and F1
#   Covariates: delta, nCPG_med, nCPG_q3
resultsPL_df <-
  probeLassoRes_df %>%
  select(power, AuPR, precision, mcc, F1,
         adjPval, mLassoRad, minDmrSep,
         delta, seed, nCPG_med, nCPG_q3) %>%
  mutate(nCPG_med = as.numeric(nCPG_med)) %>%
  mutate(nCPG_q3 = as.numeric(nCPG_q3))

###  Replace Missing Values with 0  ###
# For each of the responses, 0 is as bad as it gets. Therefore, we can replace
#   all NAs and NaNs with 0. However, precision start at 1 and never reach 0,
#   so we replace all of the 0s for precision with 1.
resultsPL_df <- resultsPL_df[complete.cases(resultsPL_df), ]

###  Correlation Matrix  ###
# First we build a correlation matrix to inspect possible multicollinearity.
round(cor(resultsPL_df), 2)
# Basically the same story as DMRcate. Because we have the mean Lasso radius as
#   one of our predictors, we will take the median number of CPGs instead of
#   upper quartile as our block effect.
resultsPL_df <-
  resultsPL_df %>%
  select(AuPR, adjPval, mLassoRad, minDmrSep, delta, nCPG_q3)
# We've removed any multicollinearity issues, except for the nCPG with mLasRad
cor(resultsPL_df, method = "spearman")
# Now, the correlation test:
cor.test(resultsPL_df$adjPval,
         resultsPL_df$AuPR,
         method = "spearman")$p.value # related

cor.test(resultsPL_df$mLassoRad,
         resultsPL_df$AuPR,
         method = "spearman")$p.value # related

cor.test(resultsPL_df$minDmrSep,
         resultsPL_df$AuPR,
         method = "spearman")$p.value # not related


###  Statistical Prediction for AuPR  ###
plLin_mod <- lm(AuPR ~ ., data = resultsPL_df)
par(mfrow = c(2, 2))
plot(plLin_mod)
# Well that's a resounding "NO". However, we see that observation 420 is a
#   potential outlier.
dmrcateRes_df[418:420, ]
# I don't see any problems.

resultsPL_df <-
  resultsPL_df %>%
  mutate(pVal_X_mLasRd = adjPval * mLassoRad)
plBeta_mod <- betareg(AuPR ~ ., data = resultsPL_df)
plot(plBeta_mod)
summary(plBeta_mod)
coef(plBeta_mod)
# The direction of the adjPval and mLassoRd relationships jive with their
#   Spearman correlations, as does the "not a relationship" result for
#   minDmrSep. What we can then say is this: given delta and the number of
#   CPGs, for each unit increase in adjPval, logit(AuPR) increases by 0.96451899;
#   for each unit increase of mLassoRd, logit(AuPR) increases by 0.00230829
#   (but remember that this effect may be occluded by the nCPG block).

# Let's have a sanity check:
par(mfrow = c(1,1))
unique(resultsPL_df$adjPval)
boxplot(resultsPL_df$AuPR ~ resultsPL_df$adjPval)
# This shows no effect.
boxplot(resultsPL_df$AuPR ~ resultsPL_df$mLassoRad)
# This effect is clear: increasing the mean lasso radius increases the AuPR. We
#   take this result with a grain of salt, however: obviously if we test a
#   larger area, we will find more stuff. Interpeting the mean lasso radius
#   fairly (conditional on the number of CPGs), may prove to be quite difficult
boxplot(resultsPL_df$AuPR ~ resultsPL_df$pVal_X_mLasRd)
# What I see here: the interaction is simply showing how the p-value make the
#   lasso radius effect more pronounced the smaller it is



######  Bumphunter  ###########################################################
data("bumphunterRes_df")
colnames(bumphunterRes_df)
# Model components:
#   Predictors: cutoffQ and maxGap
#   Responses : power, AuPR, precision, mcc, and F1
#   Covariates: delta, nCPG_med, nCPG_q3
resultsBump_df <-
  bumphunterRes_df %>%
  select(power, AuPR, precision, mcc, F1,
         cutoffQ, maxGap,
         delta, seed, nCPG_med, nCPG_q3) %>%
  mutate(nCPG_med = as.numeric(nCPG_med)) %>%
  mutate(nCPG_q3 = as.numeric(nCPG_q3))

###  Replace Missing Values with 0  ###
# For each of the responses, 0 is as bad as it gets. Therefore, we can replace
#   all NAs and NaNs with 0. However, precision start at 1 and never reach 0,
#   so we replace all of the 0s for precision with 1.
resultsBump_df    <- resultsBump_df[complete.cases(resultsBump_df), ]

###  Correlation Matrix  ###
# First we build a correlation matrix to inspect possible multicollinearity.
round(cor(resultsBump_df), 2)
resultsBump_df <-
  resultsBump_df %>%
  select(AuPR, cutoffQ, maxGap, delta, nCPG_med)
# We've removed any multicollinearity issues.
cor(resultsBump_df, method = "spearman")
# Now, the correlation test:
cor.test(resultsBump_df$cutoffQ,
         resultsBump_df$AuPR,
         method = "spearman")$p.value

cor.test(resultsBump_df$maxGap,
         resultsBump_df$AuPR,
         method = "spearman")$p.value
# Both are related


###  Statistical Prediction for AuPR  ###
bumpLin_mod <- lm(AuPR ~ ., data = resultsBump_df)
par(mfrow = c(2, 2))
plot(bumpLin_mod)
# Well that's a resounding "NO".

bumpBeta_mod <- betareg(AuPR ~ ., data = resultsBump_df)
plot(bumpBeta_mod)
summary(bumpBeta_mod)
coef(bumpBeta_mod)
# The cutoffQ parameter results agree with the Spearman correlation. The maxGap
#   parameter is not significant in the model, but the sign changes.

# Let's have a sanity check:
par(mfrow = c(1,1))
boxplot(resultsBump_df$AuPR ~ resultsBump_df$cutoffQ)
# The smaller quantiles yield better performance
boxplot(resultsBump_df$AuPR ~ resultsBump_df$maxGap)
# There is no effect.



######  Comb-p  ###############################################################
data("combpRes_df")
colnames(combpRes_df)
# Model components:
#   Predictors: combSeed and combDist
#   Responses : power, AuPR, precision, mcc, and F1
#   Covariates: delta, nCPG_med, nCPG_q3
resultsComb_df <-
  combpRes_df %>%
  select(power, AuPR, precision, mcc, F1,
         combSeed, combDist,
         delta, seed, nCPG_med, nCPG_q3) %>%
  mutate(nCPG_med = as.numeric(nCPG_med)) %>%
  mutate(nCPG_q3 = as.numeric(nCPG_q3))

###  Replace Missing Values with 0  ###
# For each of the responses, 0 is as bad as it gets. Therefore, we can replace
#   all NAs and NaNs with 0. However, precision start at 1 and never reach 0,
#   so we replace all of the 0s for precision with 1.
resultsComb_df    <- resultsComb_df[complete.cases(resultsComb_df), ]


###  Correlation Matrix  ###
# First we build a correlation matrix to inspect possible multicollinearity. We
#   will keep nCPG_q3 because nCPG_med has near-constant variance. Half of the
#   values are 7, and another third of the values are 8.
round(cor(resultsComb_df), 2)
resultsComb_df <-
  resultsComb_df %>%
  select(AuPR, combSeed, combDist, delta, nCPG_med)
# We've removed any multicollinearity issues, except for with nCPGs
cor(resultsComb_df, method = "spearman")
# Now, the correlation test:
cor.test(resultsComb_df$combSeed,
         resultsComb_df$AuPR,
         method = "spearman")$p.value # not related

cor.test(resultsComb_df$combDist,
         resultsComb_df$AuPR,
         method = "spearman")$p.value # related


###  Statistical Prediction for AuPR  ###
combLin_mod <- lm(AuPR ~ ., data = resultsComb_df)
par(mfrow = c(2, 2))
plot(combLin_mod)
# Well that's a resounding "NO".

combBeta_mod <- betareg(AuPR ~ ., data = resultsComb_df)
plot(combBeta_mod)
summary(combBeta_mod)
coef(combBeta_mod)
# The directions and significances of the parameters agree with the Spearman
#   correlations.

# Let's have a sanity check:
par(mfrow = c(1,1))
boxplot(resultsComb_df$AuPR ~ resultsComb_df$combSeed)
# There is no effect.
boxplot(resultsComb_df$AuPR ~ resultsComb_df$combDist)
# Increasing combDist yields better performance.
