context("RunBumphunter")

combine.dmrs.cpgs <- function(dmrs, cpgs, sig.threshold) {

  require(GenomicRanges)

  # 1. make ranges of dmr info
  dmrs <- dmrs[dmrs$dmr.pval < sig.threshold ,]

  sig.ranges <- IRanges(dmrs$dmr.start, dmrs$dmr.end)
  dmr.ranges <- GRanges (seqname=dmrs$dmr.chr, ranges=sig.ranges)

  # 2. make ranges of cpgs info
  temp.cpg.ranges <- IRanges(cpgs$MAPINFO, cpgs$MAPINFO)
  cpg.ranges <-GRanges (seqname=cpgs$chr, ranges=temp.cpg.ranges)

  #3. find overlaps
  dmr.cpgs.overlaps <- as.data.frame(findOverlaps(query=dmr.ranges, subject=cpg.ranges, type="any"))

  dmr.cpgs.overlaps$dmr.order <- as.numeric(dmr.cpgs.overlaps$queryHits)
  dmr.cpgs.overlaps$cpg.order <- as.numeric(dmr.cpgs.overlaps$subjectHits)

  # 4. merge with dmr info

  dmrs$dmr.row <- as.numeric(seq.int(nrow(dmrs)))

  dmrs.info <- merge (x=dmr.cpgs.overlaps, y=dmrs, by.x="dmr.order", by.y="dmr.row")

  # 5. merge with cpgs info

  cpgs$row <- as.numeric(seq.int(nrow(cpgs)))

  dmr.cpg.overlap.info <- merge (dmrs.info, cpgs, by.x="cpg.order", by.y="row")

  dmr.cpg.overlap.info <- dmr.cpg.overlap.info[order(dmr.cpg.overlap.info$dmr.order) ,]

  # 6. add number of cpgs
  ncpgs <- as.data.frame(table (dmr.cpg.overlap.info$dmr.order))

  dmr.cpg.overlap.ncpgs <-merge (dmr.cpg.overlap.info, ncpgs, by.x="dmr.order", by.y="Var1")


  return (dmr.cpg.overlap.ncpgs)
}


runBumphunter <- function(beta.values, mval_each, r,
                          pickCutoffQ_num,
                          maxGap_int,
                          dmr.sig.threshold = 0.05,
                          min.cpgs = 5,
                          numCores = 2) {

  # source ("scripts_6-22-2018/3_combineDMRcpgs_ncpgs.R")

  ptm <- proc.time()
  require(bumphunter)
  require(minfi)
  require(doParallel)

  myBetas <- as.matrix(beta.values)
  myMs <- logit2(myBetas)

  type <- factor(sub(".*-", "", colnames(myMs)))
  design1 <- model.matrix(~type)

  # Gabriel O.: my machine does not have 22 cores, only 20.
  # registerDoParallel(cores = 22)
  registerDoParallel(cores = numCores)
  dmrcoutput <- tryCatch(bumphunter(myMs, design = design1, chr = temp.location$CHR,
                                    pos = temp.location$MAPINFO, cluster = NULL, coef = 2,
                                    cutoff = NULL, pickCutoff = TRUE,
                                    pickCutoffQ = pickCutoffQ_num,
                                    maxGap = maxGap_int,
                                    nullMethod = "permutation",
                                    smooth = FALSE, smoothFunction = locfitByCluster,
                                    useWeights = FALSE, permutations = NULL, B = 10,
                                    verbose = TRUE, type = "M"), error = function(e1) {
                                      return(NULL)
                                    })

  bumphunter:::foreachCleanup()

  elapsedtime <- proc.time() - ptm

  # extract results
  if (!is.null(dmrcoutput)) {
    results.df <- dmrcoutput$table

    # add variables to be processed later - this part is unique to each method
    results.df$dmr.pval <- results.df$p.valueArea

    results.df$dmr.chr <- paste0("chr", as.character(results.df$chr))
    results.df$dmr.start <- results.df$start
    results.df$dmr.end <- results.df$end

    # get ncpgs for each dmr
    temp <- combine.dmrs.cpgs(dmrs = results.df, cpgs = cpg.location, sig.threshold = dmr.sig.threshold)
    temp.ncpgs <- unique(subset(temp, select = c(dmr.chr, dmr.start, dmr.end, Freq)))

    results.df2 <- merge(results.df, temp.ncpgs, by = c("dmr.chr", "dmr.start", "dmr.end"))
    results.df2$dmr.n.cpgs <- results.df2$Freq

    results.df2 <- results.df2[results.df2$dmr.n.cpgs >= min.cpgs, ]
    results.df2 <- subset(results.df2, select = -Freq)

    if (nrow(results.df2) == 0) {
      results.df2 <- NULL
    }

  }

  # if no predicted cluster
  if (is.null(dmrcoutput)) {
    results.df2 <- NULL
  }

  res <- list(results.df2, elapsedtime[3])
  return(res)
}


######  Test  ######
# Load Data
data("betaVals_mat")
data("cpgLocation_df")
data("startEndCPG_df")

# Simulate a Data Set
treatment_ls <- SimulateData(beta_mat = betaVals_mat,
                             Aclusters_df = startEndCPG_df,
                             delta_num = 0.4,
                             seed_int = 100)


###  Setup from Main Script  ###
temp <- treatment_ls$simBetaVals_df
cpg.location <- cpgLocation_df

# get location info for the cpgs
temp$cpg <- rownames(temp)

temp.location <- merge(temp, cpg.location,
                       by.x = "cpg", by.y = "ILMNID")
rownames(temp.location) <- temp.location$cpg
temp.location$CHR <- substr(temp.location$chr, 4, 6)

location.index <- c(grep("CHR", colnames(temp.location)),
                    grep("MAPINFO", colnames(temp.location)),
                    grep("cpg", colnames(temp.location)),
                    grep("chr", colnames(temp.location)))

#take only the annotated beta values
beta.matrix <- as.matrix (temp.location[, -location.index])

# Gabriel O.: literally, the beta.matrix is just the temp matrix ordered by the
#   CPGs (rownames). The chromosomes and MAP values do not match, however. I'll
#   try to make this cleaner:
mergedBetas_df <- merge(treatment_ls$simBetaVals_df,
                        cpgLocation_df,
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


###  Compare Functions  ###
# Old functions
set.seed(100)
a <- Sys.time()
oldRes <- runBumphunter(beta.values = beta.matrix,
                        mval_each, r,
                        pickCutoffQ_num = 0.9,
                        maxGap_int = 200)
Sys.time() - a # 1.691564 min


# New function
set.seed(100)
b <- Sys.time()
newRes <- RunBumphunter(betaVals_mat = betaSorted_mat,
                        chromos_char = cpgInfo_df$chr,
                        chromPosit_num = cpgInfo_df$MAPINFO,
                        cpgLocation_df = cpgLocation_df,
                        pickCutoffQ_num = 0.9,
                        maxGap_int = 200,
                        numCores = 2)
Sys.time() - b # 1.669844 min

# test
test_that("New RunBumphunter output equal to legacy output", {

  expect_equal(oldRes[[1]], newRes[[1]])

})
