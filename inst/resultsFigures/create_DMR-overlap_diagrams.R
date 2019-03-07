######  Overlap List Function  ################################################
# migrated to R/6_Build_DMR-Overlap_List.R


######  Create a Venn Diagram Plotting Function  ##############################
# migrated to R/6_Plot_DMR-Overlaps_Venn.R


######  Test: Build Venn Diagrams  ############################################

# PDF
PlotOverlaps(
  bestResultsDir = "inst/extdata/best_cases_results/",
  figFileName = "inst/resultsFigures/testVenn_allDesigns2"
)

PlotOverlaps(
  bestResultsDir = "inst/extdata/best_cases_results/",
  figFileName = "inst/resultsFigures/Venn_Diagrams_for_seed_100.pdf",
  seeds_int = 100
)


# JPEG
PlotOverlaps(
  bestResultsDir = "inst/extdata/best_cases_results/",
  figFileName = "inst/resultsFigures/Venn_Diagram_delta0.05_seed_100.jpeg",
  device = jpeg,
  delta_num = 0.05,
  seeds_int = 100
)

PlotOverlaps(
  bestResultsDir = "inst/extdata/best_cases_results/",
  figFileName = "inst/resultsFigures/Venn_Diagram_delta0.4_seed_100.jpeg",
  device = jpeg,
  delta_num = 0.4,
  seeds_int = 100
)


# TIFF
PlotOverlaps(
  bestResultsDir = "inst/extdata/best_cases_results/",
  figFileName = "Venn_Diagram_delta0.05_seed_100.tiff",
  device = tiff,
  plotTitle = "Fig 6(A)     mu = 0.05",
  delta_num = 0.05,
  seeds_int = 100,
  units = "in", width = 5.5, height = 4.25, res = 300
)

PlotOverlaps(
  bestResultsDir = "inst/extdata/best_cases_results/",
  figFileName = "Venn_Diagram_delta0.4_seed_100.tiff",
  device = tiff,
  plotTitle = "Fig 6(B)     mu = 0.4",
  delta_num = 0.4,
  seeds_int = 100,
  units = "in", width = 5.5, height = 4.25, res = 300
)
