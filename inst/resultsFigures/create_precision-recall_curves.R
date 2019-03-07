
######  Test PR Curve Function  ###############################################
# migrated to R/6_Build_Precision-Recall_Curve_List.R



######  Build Plotting Function  ##############################################
# migrated to R/6_Plot_Precision-Recall_Curves.R



######  Calculate PR-Curve Lists  #############################################
seeds_int <- c(100, 210, 330, 450, 680)

###  Delta = 0.025  ###
prCurves0.025_ls <- lapply(seeds_int, function(seed){
  BuildPRcurve(
    bestResultsDir = "inst/extdata/best_cases_results/",
    delta = 0.025, seed = seed
  )
})

pdf("inst/resultsFigures/PR_curve_mu0.025.pdf", width = 8.5, height = 11)
par(mfrow = c(3, 2))
lapply(prCurves0.025_ls, PlotPRCurve)
dev.off()


###  Delta = 0.05  ###
prCurves0.05_ls <- lapply(seeds_int, function(seed){
  BuildPRcurve(
    bestResultsDir = "inst/extdata/best_cases_results/",
    delta = 0.05, seed = seed
  )
})

pdf("inst/resultsFigures/PR_curve_mu0.05.pdf", width = 8.5, height = 11)
par(mfrow = c(3, 2))
lapply(prCurves0.05_ls, PlotPRCurve)
dev.off()


###  Delta = 0.1  ###
prCurves0.1_ls <- lapply(seeds_int, function(seed){
  BuildPRcurve(
    bestResultsDir = "inst/extdata/best_cases_results/",
    delta = 0.1, seed = seed
  )
})

pdf("inst/resultsFigures/PR_curve_mu0.1.pdf", width = 8.5, height = 11)
par(mfrow = c(3, 2))
lapply(prCurves0.1_ls, PlotPRCurve)
dev.off()


###  Delta = 0.15  ###
prCurves0.15_ls <- lapply(seeds_int, function(seed){
  BuildPRcurve(
    bestResultsDir = "inst/extdata/best_cases_results/",
    delta = 0.15, seed = seed
  )
})

pdf("inst/resultsFigures/PR_curve_mu0.15.pdf", width = 8.5, height = 11)
par(mfrow = c(3, 2))
lapply(prCurves0.15_ls, PlotPRCurve)
dev.off()


###  Delta = 0.2  ###
prCurves0.2_ls <- lapply(seeds_int, function(seed){
  BuildPRcurve(
    bestResultsDir = "inst/extdata/best_cases_results/",
    delta = 0.2, seed = seed
  )
})

pdf("inst/resultsFigures/PR_curve_mu0.2.pdf", width = 8.5, height = 11)
par(mfrow = c(3, 2))
lapply(prCurves0.2_ls, PlotPRCurve)
dev.off()


###  Delta = 0.3  ###
prCurves0.3_ls <- lapply(seeds_int, function(seed){
  BuildPRcurve(
    bestResultsDir = "inst/extdata/best_cases_results/",
    delta = 0.3, seed = seed
  )
})

pdf("inst/resultsFigures/PR_curve_mu0.3.pdf", width = 8.5, height = 11)
par(mfrow = c(3, 2))
lapply(prCurves0.3_ls, PlotPRCurve)
dev.off()


###  Delta = 0.4  ###
prCurves0.4_ls <- lapply(seeds_int, function(seed){
  BuildPRcurve(
    bestResultsDir = "inst/extdata/best_cases_results/",
    delta = 0.4, seed = seed
  )
})

pdf("inst/resultsFigures/PR_curve_mu0.4.pdf", width = 8.5, height = 11)
par(mfrow = c(3, 2))
lapply(prCurves0.4_ls, PlotPRCurve)
dev.off()


######  PR Curves for One Seed  ###############################################
library(DMRcompare)
data("betaVals_mat")
data("startEndCPG_df")
data("cpgLocation_df")

# delta_num <- c(0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4)
delta_num <- c(0.05, 0.4)

prCurves_ls <- lapply(delta_num, function(x){
  BuildPRcurve(
    bestResultsDir = "inst/extdata/best_cases_results/",
    delta = x, seed = 100
  )
})

PlotPRCurve(prCurves_ls[[2]])


# PDF
pdf("inst/resultsFigures/PR_curve_seed100.pdf", width = 11, height = 8.5)
par(mfrow = c(2, 2))
lapply(prCurves_ls, PlotPRCurve)
dev.off()


# JPEG
jpeg("inst/resultsFigures/PR_curve_mu0.05_seed100.jpeg", width = 600, height = 400)
PlotPRCurve(prCurves_ls[[1]])
dev.off()

jpeg("inst/resultsFigures/PR_curve_mu0.4_seed100.jpeg", width = 600, height = 400)
PlotPRCurve(prCurves_ls[[2]])
dev.off()


# TIFF
tiff(
  'PR_curve_mu0.05_seed100.tiff',
  units = "in", width = 5.5, height = 4.25, res = 300
)
PlotPRCurve(prCurves_ls[[1]], plotTitle = "Fig 4(A)     mu = 0.05")
dev.off()

tiff(
  'PR_curve_mu0.4_seed100.tiff',
  units = "in", width = 5.5, height = 4.25, res = 300
)
PlotPRCurve(prCurves_ls[[2]], plotTitle = "Fig 4(B)     mu = 0.4")
dev.off()
