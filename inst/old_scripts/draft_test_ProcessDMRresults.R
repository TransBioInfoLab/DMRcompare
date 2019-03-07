ProcessDMR_Results <- function(methodFile, aclust_df) {

  # browser()

  # # This was the original function header and first line. We will take in an
  # #   Aclust results data frame instead of a file. (There are 40 unique such
  # #   Aclust data frames, so we have no need of importing them 4k times.)
  # ProcessDMR_Results <- function(methodFile, aclustFile) {
  #
  #   aclustFile <- read.csv (paste0("DATA_simulated/", aclustFile), row.names = 1)
  #
  # result.one <- readRDS (paste0 ("results/", methodFile))

  result.one <- readRDS(paste0("../DMRcomparison/results/", methodFile))

  ranges  <- result.one [[1]]
  elapsedtime <- result.one [[2]]

  if (!is.null(ranges)) {

    ############ query = significant DMRs ----------------------------------------------------------------------

    # need to limit to min.cpgs > 4 and pval < 0.05
    sig.ranges <- ranges[ranges$dmr.n.cpgs > 4 & ranges$dmr.pval < 0.05, ]

    require(IRanges)
    require(GenomicRanges)
    sig.Iranges <- IRanges(sig.ranges$dmr.start, sig.ranges$dmr.end)
    query <- GRanges(seqnames = sig.ranges$dmr.chr, ranges = sig.Iranges)

    ############ subject = Aclusters ----------------------------------------------------------------
    clusters.df <- aclust_df[, c("Clusternumber", "chromosome", "start_position", "end_position", "actual")]
    clusters.unique <- unique(clusters.df)

    aclust.Iranges <- IRanges(clusters.unique$start_position, clusters.unique$end_position)
    subject <- GRanges(seqname = clusters.unique$chromosome, ranges = aclust.Iranges)

    ############## overlap
    overlapinfo <- as.data.frame(findOverlaps(query, subject, type = "any", select = "all"))
    overlapinfo$dmr.order <- as.numeric(overlapinfo$queryHits)
    overlapinfo$aclust.order <- as.numeric(overlapinfo$subjectHits)

    ############### merge with dmr info
    results <- sig.ranges

    results$dmr.row <- as.numeric(seq.int(nrow(results)))
    results$predicted <- "positive"

    dmrs.overlap <- merge(x = overlapinfo, y = results, by.x = "dmr.order", by.y = "dmr.row", all = TRUE)

    ############### merge with aclust info
    aclust <- clusters.unique

    aclust$aclust.row <- as.numeric(seq.int(nrow(aclust)))

    dmrs.aclust.overlap <- merge(x = dmrs.overlap, y = aclust,
                                 by.x = "aclust.order", by.y = "aclust.row", all = TRUE)

    all <- dmrs.aclust.overlap
  }

  if (is.null(ranges)) {

    clusters.df <- aclust_df[, c("Clusternumber", "chromosome", "start_position", "end_position", "actual")]
    clusters.unique <- unique(clusters.df)


    all <- clusters.unique
    all$predicted <- "negative"


  }

  ################## classify into different status

  all$predicted[is.na(all$predicted)] = "negative"
  all$actual[is.na(all$actual)] = "negative"

  all$status[all$actual == "positive" & all$predicted == "positive"] = "TP"
  all$status[all$actual == "positive" & all$predicted == "negative"] = "FN"

  all$status[all$actual == "negative" & all$predicted == "positive"] = "FP"
  all$status[all$actual == "negative" & all$predicted == "negative"] = "TN"

  # write.csv(all, file = paste0(result.dir, methods[d], "_mergedresult_with_aclust_for_miu", mval_each, "_rep", r, ".csv"), row.names = F)

  ####################### Frequency count of each status - based on unique aclusters, for power calculation

  all.status <- all[, c("Clusternumber", "status")]
  nrow(all.status)  ##3065

  all.status.unique <- unique(all.status)
  dim(all.status.unique)  ##3063

  status.count <- as.data.frame(table(all.status.unique$status))

  require(tidyr)
  temp.one <- spread(data = status.count, key = Var1, value = Freq)
  if (length(grep("FP", colnames(temp.one))) == 0) {
    temp.one$FP = 0
  }
  if (length(grep("TP", colnames(temp.one))) == 0) {
    temp.one$TP = 0
  }
  if (length(grep("TN", colnames(temp.one))) == 0) {
    temp.one$TN = 0
  }
  if (length(grep("FN", colnames(temp.one))) == 0) {
    temp.one$FN = 0
  }

  temp.one$stat <- "power"
  temp.one$value <- temp.one$TP/(temp.one$TP + temp.one$FN)
  temp.one$based.on <- "aclusters"
  temp.one$time <- elapsedtime
  temp.one$total <- temp.one$TP + temp.one$FN

  temp.one <- temp.one[, c("time", "FN", "FP", "TN", "TP", "stat", "value", "based.on", "total")]

  ####################### Frequency count of each status - based on unique DMRs, for precision calculation
  if (!is.null(ranges)) {

    all.dmr.status <- all[, c("dmr.order", "status")]
    nrow(all.dmr.status)

    all.dmr.status.unique <- unique(all.dmr.status)
    dim(all.dmr.status.unique)

    status.count <- as.data.frame(table(all.dmr.status.unique$status))

    require(tidyr)
    temp.two <- spread(data = status.count, key = Var1, value = Freq)
    if (length(grep("FP", colnames(temp.two))) == 0) {
      temp.two$FP = 0
    }
    if (length(grep("TP", colnames(temp.two))) == 0) {
      temp.two$TP = 0
    }
    if (length(grep("TN", colnames(temp.two))) == 0) {
      temp.two$TN = 0
    }
    if (length(grep("FN", colnames(temp.two))) == 0) {
      temp.two$FN = 0
    }

    temp.two$stat <- "precision"
    temp.two$value <- temp.two$TP/(temp.two$TP + temp.two$FP)

    temp.two$TN <- "NA"
    temp.two$FN <- "NA"

    temp.two$based.on <- "dmrs"
    temp.two$time <- elapsedtime
    temp.two$total <- temp.two$TP + temp.two$FP

    temp.two <- temp.two[, c("time", "FN", "FP", "TN", "TP", "stat", "value", "based.on", "total")]
    temp.both <- rbind(temp.one, temp.two)
  }

  if (is.null(ranges)) {
    temp.both <- temp.one
  }

  # temp.both$method <- method
  # temp.both$mu <- mval_each
  # temp.both$rep <- r

  return(list(temp.both, all))
}


######  Test  #################################################################
# Load data
library(DMRcompare)
data("betaVals_mat")
data("startEndCPG_df")

# This first step is exploratory only. I have no idea what this function does,
#   so I need to figure that out first. The inputs appear to be file names, so
#   I will have to change those to be appropriate (as my results files are
#   named differently, and I don't have any saved Aclust files).
# The results file we will inspect is:
#   "ProbeLassoResults_delta0.4_seed210_adjPvalProbe0.01_meanLassoRd1000_minDmrSep1000.csv"
# The corresponding Aclust results will be:
treatment_ls <- SimulateData(beta_mat = betaVals_mat,
                             AclustCPG_df = startEndCPG_df,
                             delta_num = 0.4,
                             seed_int = 210)
trueClusters_df <- treatment_ls$simAclusters_df

# # Test
# testPath <- "ProbeLassoResults_delta0.4_seed210_adjPvalProbe0.01_meanLassoRd1000_minDmrSep1000.csv"
# ProcessDMR_Results(methodFile = testPath,
#                    aclust_df = trueClusters_df)

# Well, the results I need require the computing time as well, and my .csv files
#   don't save this. Further, for files less than 0.5Mb, the saveRDS() and
#   readRDS() can be faster than even fwrite() / fread(), all while using less
#   hard disk space (the compression ratio for raw numeric data is very high).

# Test Again
testPath <- "ProbeLasso_m0.4_rep2_adjPvalProbe0.01_meanLassoRd1000_minDmrSep1000.RDS"
ProcessDMR_Results(methodFile = testPath,
                   aclust_df = trueClusters_df)

# intermediary results:
all2 <- CleanDMRResults(results_ls = result.one, Aclusters_df = aclust_df)
temp.both2 <- SummarizeDMRResults(cleanDMR_df = all2, time_num = result.one[[2]])
# We pass. We can't exactly write a test for this because the format of the
#   results will be different. However, this should be close to a proper test:


###  Original  ###
testPath <- "ProbeLasso_m0.4_rep2_adjPvalProbe0.01_meanLassoRd1000_minDmrSep1000.RDS"
oldOut_ls <- ProcessDMR_Results(methodFile = testPath,
                                aclust_df = trueClusters_df)

###  New  ###
res_ls <- readRDS(paste0("../DMRcomparison/results/", testPath))
newOut1_df <- CleanResults(dmrResults_ls = res_ls,
                           Aclusters_df = trueClusters_df)
newOut2_df <- SummarizeResults(cleanDMR_df = newOut1_df,
                               time_num = res_ls[[2]])

###  Compare  ###
all.equal(oldOut_ls[[2]], newOut1_df)
all.equal(as.numeric(oldOut_ls[[1]][1, 2:5]), as.numeric(newOut2_df[, 2:5]))
# We can't make this test test official because it will fail on any other
#   machine (the paths will be different).
