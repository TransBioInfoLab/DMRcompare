% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3_standard_method_output.R
\name{StandardizeOutput}
\alias{StandardizeOutput}
\title{Standardize DMR Method Comparison Output}
\usage{
StandardizeOutput(methodOut_df, method = c("DMRcate", "ProbeLasso",
  "Bumphunter", "Comb_p"), cpgLocation_df, dmr.sig.threshold = 0.05,
  min.cpgs = 5)
}
\arguments{
\item{methodOut_df}{An object returned by \code{\link[DMRcate]{dmrcate}},
\code{\link[ChAMP]{champ.DMR}} (with \code{method = "ProbeLasso"}), or
\code{\link[bumphunter]{bumphunter}}.}

\item{method}{A character string matching \code{"DMRcate"},
\code{"ProbeLasso"}, or \code{"Bumphunter"}. Partial matching is handled
via the \code{\link[base]{match.arg}} function.}

\item{cpgLocation_df}{A data frame matching chromosomes to CPG names and
locations. This is given in the \code{cpgLocation_df} data set.}

\item{dmr.sig.threshold}{Significance level to select regions (with
\code{dmr.pval} less than the specified value) passed to the internal
\code{\link{MergeDMRsWithCPGs}} function.}

\item{min.cpgs}{The minimum number of CPGs before we consider a result
significant. Defaults to 5.}
}
\value{
A data frame of ...
}
\description{
Given the output from either \code{\link[DMRcate]{dmrcate}},
   \code{\link[ChAMP]{champ.DMR}} (for the ProbeLasso method), or
   \code{\link[bumphunter]{bumphunter}}, standardize the $p$-value, start
   and end of each DMR, and chromosome columns.
}
\examples{
  # Called internally by the RunDMRcate(), RunProbeLasso(), RunBumphunter(),
  #   ProcessCombpResults(), BuildOverlaps(), and BuildPRcurve() functions.
}
\keyword{internal}
