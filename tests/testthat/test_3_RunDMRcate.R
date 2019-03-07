context("RunDMRcate")

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


runDMRcate <- function (beta.values, mval_each, r,
                        lambda_int, C_int,
                        dmr.sig.threshold = 0.05, min.cpgs = 5)
{
  # source ("scripts_6-22-2018/3_combineDMRcpgs_ncpgs.R")

  ptm <- proc.time()
  require(DMRcate)
  myBetas<-as.matrix(beta.values)
  myMs <- logit2(myBetas)

  type <- factor(sub(".*-", "", colnames(myMs)))
  design1 <- model.matrix(~ type)

  myannotation <- cpg.annotate("array", myMs, what="M", arraytype = "450K",
                               analysis.type="differential", design=design1, coef=2, fdr=0.05)

  dmrcoutput <- tryCatch(dmrcate(myannotation, lambda = lambda_int, C = C_int), error = function(e1) {return(NULL)})

  elapsedtime <- proc.time() - ptm

  # extract results
  #if any predicted cluster is identified from DMRcate
  if(!is.null(dmrcoutput))
  {
    # browser()

    results.df <- data.frame (extractRanges(dmrcoutput, genome = "hg19"))

    # add variables to be processed later - this part is unique to each method
    results.df$dmr.pval <- results.df$Stouffer

    results.df$dmr.chr <- results.df$seqnames
    results.df$dmr.start <- results.df$start
    results.df$dmr.end   <- results.df$end

    # get ncpgs for each dmr
    temp <- combine.dmrs.cpgs (dmrs = results.df, cpgs = cpg.location, sig.threshold = dmr.sig.threshold)
    temp.ncpgs <- unique(subset (temp, select = c(dmr.chr, dmr.start, dmr.end, Freq) ) )

    results.df2 <- merge (results.df, temp.ncpgs, by = c("dmr.chr", "dmr.start", "dmr.end"))
    results.df2$dmr.n.cpgs <- results.df2$Freq

    results.df2 <- results.df2[results.df2$dmr.n.cpgs >= min.cpgs ,]
    results.df2 <- subset (results.df2, select = -Freq)

    # output
    # if (nrow(results.df2) > 0 ){
    #   write.csv (results.df2, file = paste0 ("results/",
    #                                          methods[d], "_results", mval_each,
    #                                          "_rep",r,
    #                                          "_lambda",params_mat[i, 1],
    #                                          "_C", params_mat[i, 2],
    #                                          ".csv"),
    #              row.names = FALSE)
    # }

    if (nrow(results.df2) == 0 ) {
      results.df2 <- NULL
    }

  }
  # if no predicted cluster
  if(is.null(dmrcoutput) )
  {
    results.df2 <- NULL
  }

  res <- list(results.df2, elapsedtime[3])
  return(res)
}

######  Test  ######
library(DMRcompare)
# Load Data
data("betaVals_mat")
data("cpgLocation_df")
data("startEndCPG_df")

# Simulate a Data Set
treatment_ls <- SimulateData(beta_mat = betaVals_mat,
                             Aclusters_df = startEndCPG_df,
                             delta_num = 0.4,
                             seed_int = 100)

# Old functions
cpg.location <- cpgLocation_df
oldRes <- runDMRcate(beta.values = treatment_ls$simBetaVals_df,
                     lambda_int = 1000, C_int = 2)

# New function
newRes <- RunDMRcate(betaVals_mat = treatment_ls$simBetaVals_df,
                     cpgLocation_df = cpgLocation_df,
                     lambda_int = 1000, C_int = 2)

# test
test_that("New RunDMRcate output equal to legacy output", {

  expect_equal(oldRes[[1]], newRes[[1]])

})
