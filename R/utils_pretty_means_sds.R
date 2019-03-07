#' Calculate Mean and Standard Deviation
#'
#' @description Given a vector, calculate its mean and standard deviation, then
#'    display the results in a "pretty" manner.
#'
#' @param vec A numeric vector to summarize
#' @param sigFigsMean How many significant figures to display the mean? Defaults
#'    to 2.
#' @param sigFigsSD How many significant figures to display the standard
#'    deviations? Defaults to 2.
#' @param charact Should the results be returned as a character string in
#'    "\code{xx.xx (x.xx)}" format? Defaults to \code{TRUE}. If \code{FALSE},
#'    this function will return a data frame with column names \code{mean} and
#'    \code{std_dev}.
#' @param na.rm Should missing values be removed? Defaults to \code{TRUE}.
#'
#' @details This function is for help in building large tables of output wherein
#'    the user may need many entries in the form "mean (std dev)".
#'
#' @return A character string in the form of "mean (sd)", if \code{charact} is
#'    set to \code{TRUE}. Otherwise, a data frame of the mean and standard
#'    deviation of \code{vec}.
#'
#' @importFrom stats sd
#'
#' @export
#'
#' @keywords internal
#'
#' @examples
#'    # Called by the script in inst/docs/Method_compare_report_20180705.Rmd
#'    CalcMeanSD(1:10)
CalcMeanSD <- function(vec,
                       sigFigsMean = 2,
                       sigFigsSD = 2,
                       charact = TRUE,
                       na.rm = TRUE){

  xBar_num <- round(mean(vec, na.rm = na.rm), digits = sigFigsMean)
  s_num    <- round(sd(vec, na.rm = na.rm), digits = sigFigsSD)


  if(charact){
    paste0(xBar_num, " (", s_num, ")")
  } else {
    data.frame(mean = xBar_num, std_dev = s_num)
  }

}
