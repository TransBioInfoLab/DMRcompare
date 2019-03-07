
######  Test  #################################################################
library(DMRcompare)
library(doParallel)
# Load Data
data("betaVals_mat")
data("cpgLocation_df")
data("startEndCPG_df")

###  Single-Design Test  ###
# Because the bumphunter algorithm itself is written in parallel, we will not
#   have a "parallel" version of this function.
WriteBumphunterResults(beta_mat = betaVals_mat,
                       CPGs_df = cpgLocation_df,
                       Aclusters_df = startEndCPG_df,
                       deltas_num = 0.40,
                       seeds_int = 100,
                       cutoffQ_num = 0.9,
                       maxGap_int = 1000,
                       resultsDir = "inst/comparison_results/")


###  Multi-Param Test  ###
# This also tests garbage collection
getwd()
a <- Sys.time()
WriteBumphunterResults(beta_mat = betaVals_mat,
                       CPGs_df = cpgLocation_df,
                       Aclusters_df = startEndCPG_df,
                       deltas_num = 0.40,
                       seeds_int = 100,
                       cutoffQ_num = 0.95,
                       resultsDir = "inst/comparison_results/")
Sys.time() - a # 5.49612 min
# This will write 5 files. At it's most expensive, this required 20Gb of RAM

###  All  ###
a <- Sys.time()
WriteBumphunterResults(beta_mat = betaVals_mat,
                       CPGs_df = cpgLocation_df,
                       Aclusters_df = startEndCPG_df,
                       numCores = 20,
                       resultsDir = "../DMRcomparison_results/")
Sys.time() - a
# 600 files in 12.03613 hrs (data.table::fwrite); 12.24827 hrs (base::saveRDS)
