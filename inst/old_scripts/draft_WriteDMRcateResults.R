
######  Test  #################################################################
library(DMRcompare)
# Load Data
data("betaVals_mat")
data("cpgLocation_df")
data("startEndCPG_df")

# Run the print version (computing commented out, just to see the file names)
WriteDMRcateResults(beta_mat = betaVals_mat,
                    CPGs_df = cpgLocation_df,
                    Aclusters_df = startEndCPG_df)
# This would have written 8 * 5 * 5 * 5 = 1000 files

WriteDMRcateResults(beta_mat = betaVals_mat,
                    CPGs_df = cpgLocation_df,
                    Aclusters_df = startEndCPG_df,
                    deltas_num = 0.4,
                    seeds_int = 100,
                    lambdas_num = 1000,
                    Cs_int = 2,
                    resultsDir = "inst/comparison_results/")
# We had to add DMRcatedata to the Depends field (and use the @import tag)

getwd()
a <- Sys.time()
WriteDMRcateResults(beta_mat = betaVals_mat,
                    CPGs_df = cpgLocation_df,
                    Aclusters_df = startEndCPG_df,
                    deltas_num = c(0, 0.4),
                    seeds_int = c(100, 210),
                    lambdas_num = c(200, 1000),
                    Cs_int = c(1, 5),
                    resultsDir = "inst/comparison_results/")
Sys.time() - a # 17.3511 min
# This will write 16 files


######  Parallel Test  ########################################################
library(DMRcompare)
library(doParallel)
# Load Data
data("betaVals_mat")
data("cpgLocation_df")
data("startEndCPG_df")

getwd()
a <- Sys.time()
WriteDMRcateResults(beta_mat = betaVals_mat,
                    CPGs_df = cpgLocation_df,
                    Aclusters_df = startEndCPG_df,
                    numCores = 18,
                    deltas_num = 0.3,
                    seeds_int = 100,
                    lambdas_num = c(200, 250, 500, 750),
                    Cs_int = 1:4,
                    resultsDir = "inst/comparison_results/")
Sys.time() - a # 3.813431, 3.420486 min
# This will write 16 files. At it's most expensive, this required 56Gb of RAM

# Test garbage collection
a <- Sys.time()
WriteDMRcateResults(beta_mat = betaVals_mat,
                    CPGs_df = cpgLocation_df,
                    Aclusters_df = startEndCPG_df,
                    numCores = 18,
                    deltas_num = c(0.1, 0.15, 0.2),
                    seeds_int = 100,
                    lambdas_num = c(200, 250, 500, 750),
                    Cs_int = 1:4,
                    resultsDir = "../DMRcomparison_results/")
Sys.time() - a # 10.01768 min for 16 parameter points over the 3 small deltas
# Garbage collection via stopCluster() is working fine. I think we're ready to
#   leave it running over night.

###  All  ###
a <- Sys.time()
WriteDMRcateResults(beta_mat = betaVals_mat,
                    CPGs_df = cpgLocation_df,
                    Aclusters_df = startEndCPG_df,
                    numCores = detectCores() - 2,
                    seeds_int = c(330, 450, 680),
                    resultsDir = "../DMRcomparison_results/")
Sys.time() - a
# 1000 files in 3.245189 hrs (data.table::fwrite)
