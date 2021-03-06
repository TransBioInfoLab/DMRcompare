% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4_simulate_and_save_ProbeLasso_results.R
\name{WriteProbeLassoResults}
\alias{WriteProbeLassoResults}
\title{Calculate and Save ProbeLasso Method Results for Specified Design Points}
\usage{
WriteProbeLassoResults(beta_mat, CPGs_df, Aclusters_df, parallel = TRUE,
  numCores = detectCores() - 2, deltas_num = c(0, 0.025, 0.05, 0.1,
  0.15, 0.2, 0.3, 0.4), seeds_int = c(100, 210, 330, 450, 680),
  pVals_num = c(0.001, 0.01, 0.05, 0.1), aveLassoRad_int = c(375, 700,
  1000), minDmrSep_int = c(200, 250, 500, 750, 1000),
  resultsDir = "DMRcate_compare/", verbose = !parallel)
}
\arguments{
\item{beta_mat}{A beta value matrix for methylation samples from a
   450k methylation array with CpG IDs as the row names and sample IDs as
   the column names. An example is given in the \code{betaVals_mat} data set.
#'}

\item{CPGs_df}{An annotation table that indicates locations of CpGs.
This data frame has CpG IDs as the rows with matching chromosome and
location info in the columns. Specifically, the columns are: \code{ILMNID}
 - the CPG ID; \code{chr} - the chromosome label; and \code{MAPINFO} -
the chromosome location. An example is given in the \code{cpgLocation_df}
data set.}

\item{Aclusters_df}{A data frame of beta values and CpG information for
clusters of CpGs over a 450k methylation array. The rows correspond to the
CPGs. The columns have information on the cluster number, chromosome,
cluster start and end locations, and the beta values for each subject
grouped by some clinical indicator (e.g. case v. control). An example is
given in the \code{startEndCPG_df} data set. This data set can be
generated by the file \code{/inst/1_Aclust_data_import.R}}

\item{parallel}{Should computing be completed over multiple computing cores?
Defaults to \code{TRUE}.}

\item{numCores}{If \code{parallel}, how many cores should be used? Defaults
to two less than the number of available cores (as calculated by the
\code{\link[parallel]{detectCores}} function).}

\item{deltas_num}{A vector of treatment sizes: non-negative real numbers to
add to the beta values within randomly-selected clusters for a single
class of subjects. This artifically creates differentially-methylated
regions (DMRs).}

\item{seeds_int}{A vector of seed values passed to the
\code{\link[base]{Random}} function to enable reproducible results}

\item{pVals_num}{A vector of the minimum thresholds of significance for
probes to be includede in DMRs, passed through the
\code{\link{RunProbeLasso}} function to the \code{\link[ChAMP]{champ.DMR}}
function.}

\item{aveLassoRad_int}{A vector of radii around each differential
methylation position to detect DMR, passed to the
\code{\link[ChAMP]{champ.DMR}} function.}

\item{minDmrSep_int}{A vector of the minimum seperation (bp) values between
neighbouring DMRs, passed to the \code{\link[ChAMP]{champ.DMR}} function.}

\item{resultsDir}{Where should the results be saved? Defaults to
\code{DMRcate_compare/}.}

\item{verbose}{Should the function print progress messages? Defaults to
\code{TRUE} only if \code{parallel = FALSE}.}
}
\value{
Saves output files in the specified results directory.
}
\description{
Given a set of design points (treatment effect size to be added
   and number of repetitions), simulate simulate methylation data with DMRs
   and apply the ProbeLasso method (via the \code{\link[ChAMP]{champ.DMR}}
   function) to them. Write the results to a file.
}
\details{
This function creates matrices of all combinations of design points
   and all combinations of parameters. For each combination, this function
   executes the internal \code{\link{RunProbeLasso}} function and saves the
   results as a compressed \code{.RDS} file.
}
\examples{
\dontrun{
   data("betaVals_mat")
   data("cpgLocation_df")
   data("startEndCPG_df")

   WriteProbeLassoResults(
     beta_mat = betaVals_mat,
     CPGs_df = cpgLocation_df,
     Aclusters_df = startEndCPG_df
   )
}
}
