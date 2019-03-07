# Save all Analysis and Time Elapsed files as .RDS files in a single directory
# Gabriel Odom
# 2018-07-04

# yep, I'm working on a holiday...

trueSeeds_int <- c(100, 210, 330, 450, 680)

######  File Structure  #######################################################
# Zhen has saved the comb-p results in nested files: 
#   1. Top Level: seed-0.001, seed-0.01, seed-0.05, seed-0.1
#   2. Second Level: dist-200, dist-250, dist-500, dist-750, dist-1000
#   3. Third Level: linear_model_result_for_mu<DELTA>_repetition_<SEED>, for
#      every combination of delta = (0, 0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4)
#      and seed = (100, 210, 330, 450, 680).
#   4. Results Level: the results file is stored in "Train.anno.hg19.BED" and
#      the time elapsed is stored in "run_time.txt". Once inside these results
#      directories, we need to read in both files, create a list object which
#      contains them, then save this list as a compressed RData file into a top
#      level directory called "Results_Lists".

###  Top-Level Folder Names  ###
level1Files_char <- list.dirs()
resDirs_char <- level1Files_char[grep("linear", level1Files_char)]


###  Second-Level Comb-p Parameter Names  ###
level2Params_ls <- strsplit(resDirs_char, "/")

ExtractParamLevels <- function(splitNames_ls, param){
  
  posIdx <- grep(param, splitNames_ls[[1]])
  point_char <- sapply(splitNames_ls, `[`, posIdx)
  point_num <- as.numeric(
    gsub(pattern = param, replacement = "", point_char)
  )
  
  sort(unique(point_num))
  
}

combSeeds_num <- ExtractParamLevels(level2Params_ls, "seed-")
combDists_int <- ExtractParamLevels(level2Params_ls, "dist-")


###  Third-Level Design Point Names  ###
level3designs_char <- sapply(level2Params_ls, `[`, 4)
level3designs_ls <- strsplit(level3designs_char, "_")

mu_num <- ExtractParamLevels(level3designs_ls, "mu")
reps_int <- as.numeric(unique(sapply(level3designs_ls, `[`, 7)))

######  for-loop Setup  #######################################################
designPts_mat <- expand.grid(reps_int, mu_num)
paramsGrid_mat <- expand.grid(combSeeds_num, combDists_int)

a <- Sys.time()
for(i in 1:nrow(designPts_mat)){
  
  repetition  <- designPts_mat[i, 1]
  mu <- designPts_mat[i, 2]
  
  ###  Inner for-loop  ###
  for(j in 1:nrow(paramsGrid_mat)){
    
    seed <- paramsGrid_mat[j, 1]
    dist  <- paramsGrid_mat[j, 2]
    
    
    ###  Create Directory Name  ###
    dir_char <- paste0("seed-", seed,
                       "/dist-", dist,
                       "/linear_model_result_for",
                       "_mu", mu,
                       "_repetition_", repetition)
    
    
    ###  Read Files  ###
    inFileName_char <- paste0(dir_char, "/Train.anno.hg19.BED")
    if(file.exists(inFileName_char)){
      
      combp_in <- read.delim(inFileName_char, header = TRUE)
      
      inNames_char <- colnames(combp_in)
      inNames_char[1] <- "chrom"
      colnames(combp_in) <- inNames_char
      
    } else {
      
      combp_in <- data.frame(
        chrom     = NA_character_,
        start     = NA_integer_,
        end       = NA_integer_,
        min_p     = NA_real_,
        n_probes  = NA_integer_,
        z_p       = NA_real_,
        z_sidak_p = NA_real_,
        refGene_name     = NA_character_,
        refGene_distance = NA_character_,
        refGene_feature  = NA_character_,
        cpgIslandExt_name     = NA_character_,
        cpgIslandExt_distance = NA_integer_,
        cpgIslandExt_feature  = NA_character_
      )
      
    }
    
    timeFileName_char <- paste0(dir_char, "/run_time.txt")
    elapsedTime_in <- read.table(timeFileName_char) / 1000
    
    
    ###  Write Data File  ###
    outFileName <- paste0("CombpResults",
                          "_delta", mu,
                          "_seed", trueSeeds_int[repetition],
                          "_combSeed", seed,
                          "_combDist", dist,
                          ".RDS")
    out_ls <- list(combp_in, unlist(elapsedTime_in))
    saveRDS(out_ls, file = paste0("Results_Lists/", outFileName))
    
  } # END for(j)
  print(i)
  
} # END for(i)
Sys.time() - a # 23.20102 sec