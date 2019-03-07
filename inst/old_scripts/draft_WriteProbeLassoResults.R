
######  Test  #################################################################
library(DMRcompare)
# Load Data
data("betaVals_mat")
data("cpgLocation_df")
data("startEndCPG_df")

# Single-design Test
WriteProbeLassoResults(beta_mat = betaVals_mat,
                       CPGs_df = cpgLocation_df,
                       Aclusters_df = startEndCPG_df,
                       parallel = FALSE,
                       deltas_num = 0.40,
                       seeds_int = 100,
                       pVals_num = 0.1,
                       aveLassoRad_int = 375,
                       minDmrSep_int = 1000,
                       resultsDir = "inst/comparison_results/")



######  Parallel Test  ########################################################
library(DMRcompare)
library(doParallel)
# Load Data
data("betaVals_mat")
data("cpgLocation_df")
data("startEndCPG_df")

getwd()
a <- Sys.time()
WriteProbeLassoResults(beta_mat = betaVals_mat,
                       CPGs_df = cpgLocation_df,
                       Aclusters_df = startEndCPG_df,
                       deltas_num = 0.40,
                       seeds_int = 100,
                       aveLassoRad_int = 1000,
                       resultsDir = "inst/comparison_results/")
Sys.time() - a # 3.946578 min
# This will write 20 files. At it's most expensive, this required 52Gb of RAM

# Test garbage collection
a <- Sys.time()
WriteProbeLassoResults(beta_mat = betaVals_mat,
                       CPGs_df = cpgLocation_df,
                       Aclusters_df = startEndCPG_df,
                       deltas_num = c(0.3, 0.4),
                       seeds_int = 100,
                       pVals_num = 0.01,
                       resultsDir = "inst/comparison_results/")
Sys.time() - a # 5.589961 min for 15 parameter points over 2 deltas
# Garbage collection via stopCluster() is working fine. I think we're ready to
#   leave it running over night.

###  All  ###
a <- Sys.time()
WriteProbeLassoResults(beta_mat = betaVals_mat,
                       CPGs_df = cpgLocation_df,
                       Aclusters_df = startEndCPG_df,
                       resultsDir = "../DMRcomparison_results/")
Sys.time() - a
# 2400 files in 4.302536 hrs (data.table::fwrite); 4.295948 hrs (base::saveRDS)
