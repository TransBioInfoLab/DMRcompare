#' Build a List of DMR-Overlap Lists
#'
#' @description Given a directory of best-performing results files from one of
#'    the simulation functions (\code{\link{WriteDMRcateResults}},
#'    \code{\link{WriteProbeLassoResults}},
#'    \code{\link{WriteBumphunterResults}}, or results from the \code{Comb-p}
#'    method in \code{Python}), import the raw data files and DMR overlap lists
#'    via the \code{\link[GenomicRanges]{GRanges}} and
#'    \code{\link[IRanges]{IRanges}} functions.
#'
#' @param bestResultsDir The name of the directory where the method results from
#'    the best-performing parameter settings are stored. For the full design we
#'    have included (\code{delta = c(0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4)} and
#'    \code{seed = c(100, 210, 330, 450, 680)}), this directory should contain
#'    35 \code{.RDS} files per method.
#' @param delta A treatment size corresponding to one of the simulations with
#'    completed results files in the \code{bestResultsDir} directory.
#' @param seed A seed value corresponding to one of the simulations with
#'    completed results files in the \code{bestResultsDir} directory.
#' @param CPGs_df An annotation table that indicates locations of CpGs.
#'    This data frame has CPG IDs as the rows with matching chromosome and
#'    location info in the columns. Specifically, the columns are: \code{ILMNID}
#'     - the CPG ID; \code{chr} - the chromosome label; and \code{MAPINFO} -
#'    the chromosome location. An example is given in the \code{cpgLocation_df}
#'    data set. This data set is only necessary if the results directory
#'    contains Comb-p results with the specified \code{delta} and \code{seed}
#'    values.
#' @param min.cpgs The minimum number of CPGs before we consider a result
#'    significant. Defaults to 5. This argument is only required if the results
#'    directory contains Comb-p results with the specified \code{delta} and
#'    \code{seed} values.
#'
#' @details This function is called internally by the \code{\link{PlotOverlaps}}
#'    function.
#'
#' @return A list of DMR overlaps to pass to the
#'    \code{\link[ChIPpeakAnno]{makeVennDiagram}} function.
#'
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#'
#' @export
#'
#' @keywords internal
#'
#' @examples
#'   # Called internally by the PlotOverlaps() function.
#' \dontrun{
#' 
#'   data("cpgLocation_df")
#' 
#'    BuildOverlaps(
#'      bestResultsDir = "best_cases_results/",
#'      delta = 0.4,
#'      seed = 100,
#'      CPGs_df = cpgLocation_df
#'    )
#' }
#'   
BuildOverlaps <- function(bestResultsDir,
                          delta = c(0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4),
                          seed = c(100, 210, 330, 450, 680),
                          CPGs_df,
                          min.cpgs = 5){

  ### List Results Files  ###
  fileNames_char <- list.files(bestResultsDir)
  targetNames_char <- paste0("delta", delta, "_seed", seed)
  correctFiles_idx <- grep(targetNames_char, fileNames_char)
  correctNames_char <- fileNames_char[correctFiles_idx]
  names(correctNames_char) <- sapply(
    strsplit(correctNames_char, "_"),
    function(x){
      gsub(pattern = "Results", replacement = "", x[1])
    }
  )

  ### Load and Clean Appropriate Files  ###
  allRes_ls <- lapply(correctNames_char, function(x){

    res_ls <- readRDS(paste0(bestResultsDir, x))

    if(grepl("Combp", x)){

      results_df <- StandardizeOutput(
        methodOut_df = res_ls[[1]],
        method = "Comb_p",
        cpgLocation_df = CPGs_df
      )

      # The raw Comb-p results were filtered to > 1, not > 4.
      results_df <- results_df[results_df$dmr.n.cpgs >= min.cpgs, ]
      res_ls[[1]] <- results_df

    }

    results_df <- res_ls[[1]]
    keepRows <- rowSums(is.na(results_df)) != ncol(results_df)
    results_df <- results_df[keepRows, ]

    if(nrow(results_df) > 0){

      results_IRanges <- IRanges(
        start = results_df$dmr.start,
        end = results_df$dmr.end
      )
      GRanges(
        seqnames = results_df$dmr.chr,
        ranges = results_IRanges
      )

    }

  })

  ###  Return  ###
  attr(allRes_ls, "delta") <- delta
  attr(allRes_ls, "repl")  <- which(c(100, 210, 330, 450, 680) %in% seed)
  allRes_ls

}
