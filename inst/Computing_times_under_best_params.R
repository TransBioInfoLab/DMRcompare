# Computing Times under Optimal Parameter Settings for each Design Point
# Gabriel Odom
# 2018-07-06

# Summary: now that I've done all this work, my computer is too slow, so Lily
#   re-ran this whole script on her machine.

######  Best-Performing DMRcate  ##############################################
library(DMRcompare)
data("betaVals_mat")
data("startEndCPG_df")
data("cpgLocation_df")

# Given the table of performance by parameter setting (report in inst/docs/),
#   we construct the following design matrix:
# Actually, we are going with only one parameter set: 500, 5

# Calculate the execution times for these designs. The overall wrapper will be
#   serial, but the internal calls will use 8 cores.
# JUST KIDDING. The freaking dmrcate function doesn't support parallelization
#   on Windows.
delta_num <- c(0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4)
a1 <- Sys.time()
times1Num_ls <- lapply(delta_num, function(delta){

  ###  Set Parameters  ###
  lambda <- 500
  C_int  <- 5

  ###  Replicate over Seeds  ###
  seeds_int <- c(100, 210, 330, 450, 680)
  innerTimes_num <- numeric(length = length(seeds_int))

  for(j in seq_along(seeds_int)){

    treatment_ls <- SimulateData(beta_mat = betaVals_mat,
                                 AclustCPG_df = startEndCPG_df,
                                 delta_num = delta,
                                 seed_int = seeds_int[j])
    betas_df <- treatment_ls$simBetaVals_df

    suppressMessages(
      res_ls <- RunDMRcate(betaVals_mat = betas_df,
                           cpgLocation_df = cpgLocation_df,
                           lambda_int = lambda,
                           C_int = C_int)
    )
    innerTimes_num[j] <- res_ls[[2]]

  }

  innerTimes_num

})
Sys.time() - a1 # 0.7034515 hrs

###  Fix First Run  ###
# The first run has require() calls:
treatment_ls <- SimulateData(beta_mat = betaVals_mat,
                             AclustCPG_df = startEndCPG_df,
                             delta_num = 0.025,
                             seed_int = 100)
betas_df <- treatment_ls$simBetaVals_df
res_ls <- RunDMRcate(betaVals_mat = betas_df,
                     cpgLocation_df = cpgLocation_df,
                     lambda_int = 500,
                     C_int = 5)

times1Num_ls[[1]][1] <- res_ls[[2]]


###  Summary  ###
names(times1Num_ls) <- as.character(delta_num)
out_ls <- list()
out_ls[[1]] <- dplyr::bind_rows(
  lapply(times1Num_ls, function(x){
    data.frame(Mean = mean(x), StdDev = sd(x))
  }),
  .id = "Delta"
)

# Delta    Mean    StdDev
# 0.025  62.846  1.539896
# 0.050  62.326  1.696019
# 0.100  62.828  1.920201
# 0.150  61.846  1.766332
# 0.200  62.566  1.666652
# 0.300  65.908  2.119733
# 0.400  66.488  2.067878



######  Best-Performing ProbeLasso  ###########################################
# We are going with only one parameter set: 0.1, 1000, 500
a2 <- Sys.time()
times2Num_ls <- lapply(delta_num, function(delta){

  ###  Set Parameters  ###
  adjPval   <- 0.05
  mLassoRad <- 1000
  minDmrSep <- 1000

  ###  Replicate over Seeds  ###
  seeds_int <- c(100, 210, 330, 450, 680)
  innerTimes_num <- numeric(length = length(seeds_int))

  for(j in seq_along(seeds_int)){

    treatment_ls <- SimulateData(beta_mat = betaVals_mat,
                                 AclustCPG_df = startEndCPG_df,
                                 delta_num = delta,
                                 seed_int = seeds_int[j])
    betas_df <- treatment_ls$simBetaVals_df

    res_ls <- RunProbeLasso(betaVals_mat = betas_df,
                            cpgLocation_df = cpgLocation_df,
                            adjPvalProbe_num = adjPval,
                            meanLassoRadius_int = mLassoRad,
                            minDmrSep_int = minDmrSep)
    innerTimes_num[j] <- res_ls[[2]]

  }

  innerTimes_num

})
Sys.time() - a2 # 0.4875482 hrs


###  Summary  ###
names(times2Num_ls) <- as.character(delta_num)
out_ls[[2]] <- dplyr::bind_rows(
  lapply(times2Num_ls, function(x){
    data.frame(Mean = mean(x), StdDev = sd(x))
  }),
  .id = "Delta"
)

# Delta    Mean     StdDev
# 0.025  46.204  4.3422321
# 0.050  48.436  1.0597075
# 0.100  48.604  2.0908204
# 0.150  48.828  1.2204794
# 0.200  48.920  1.2730475
# 0.300  46.366  0.9054999
# 0.400  46.264  1.3178885



######  Best-Performing Bumphunter  ###########################################

# We are going with only two parameter sets: 0.9, 200 for delta < 0.2; 0.9, 500
#   for delta >= 0.2
a3 <- Sys.time()
times3Num_ls <- lapply(delta_num, function(delta){

  ###  Set Parameters  ###
  cutoffQ <- 0.95
  maxGap  <- 250

  ###  Replicate over Seeds  ###
  seeds_int <- c(100, 210, 330, 450, 680)
  innerTimes_num <- numeric(length = length(seeds_int))

  for(j in seq_along(seeds_int)){

    treatment_ls <- SimulateData(beta_mat = betaVals_mat,
                                 AclustCPG_df = startEndCPG_df,
                                 delta_num = delta,
                                 seed_int = seeds_int[j])
    betas_df <- treatment_ls$simBetaVals_df

    ###  Data Wrangling  ###
    mergedBetas_df <- merge(betas_df, cpgLocation_df,
                            by.x = "row.names",
                            by.y = "ILMNID")

    cpgInfo_df <- subset(mergedBetas_df, select = c("chr", "MAPINFO"))
    cpgInfo_df$chr <- substr(cpgInfo_df$chr, 4, 6)

    betaSorted_df <- mergedBetas_df
    row.names(betaSorted_df) <- betaSorted_df$Row.names
    betaSorted_df$Row.names <-
      betaSorted_df$chr <-
      betaSorted_df$MAPINFO <-
      NULL
    betaSorted_mat <- as.matrix(betaSorted_df)

    # Clean memory
    rm(treatment_ls, betas_df, mergedBetas_df, betaSorted_df)


    ###  Run  ###
    suppressMessages(
      res_ls <- RunBumphunter(betaVals_mat = betaSorted_mat,
                              chromos_char = cpgInfo_df$chr,
                              chromPosit_num = cpgInfo_df$MAPINFO,
                              cpgLocation_df = cpgLocation_df,
                              pickCutoffQ_num = cutoffQ,
                              maxGap_int = maxGap,
                              numCores = 1)
    )
    innerTimes_num[j] <- res_ls[[2]]

  }

  innerTimes_num

})
Sys.time() - a3 # 5.037842 hrs


###  Summary  ###
names(times3Num_ls) <- as.character(delta_num)
out_ls[[3]] <- dplyr::bind_rows(
  lapply(times3Num_ls, function(x){
    data.frame(Mean = mean(x), StdDev = sd(x))
  }),
  .id = "Delta"
)

# Delta     Mean    StdDev
# 0.025  526.688  20.78909
# 0.050  526.500  16.71671
# 0.100  517.552  18.29805
# 0.150  511.940  16.64845
# 0.200  506.506  19.60679
# 0.300  499.946  16.37973
# 0.400  494.084  13.30777


###  Save Times  ###
names(out_ls) <- c("DMRcate", "ProbeLasso", "Bumphunter")
outTimes_df <- dplyr::bind_rows(out_ls, .id = "Method")
write.csv(outTimes_df,
          file = "./inst/resultsData/DMRMethodTimes_Lily.csv",
          row.names = FALSE)
