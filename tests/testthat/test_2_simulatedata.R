context("Simulating Data")

simulateData <- function(beta.value.all, res_aclcpg, beta_start_columnid, exsampleno, ctsampleno, miuv, seed_value, totalori, trimNames = FALSE) {

  colm <- ncol(res_aclcpg)  ##newly added the line

  totalclsno <- max(res_aclcpg$Clusternumber)  ###find the total number of clusters initially from Aclust results

  set.seed(seed_value)  ##set the seed value from the 'rp'-th number of element from the vector 'seed_values'

  randclss <- sample(totalclsno, totalori)  ###randomly pick totalori=500 clusters

  res_aclcpg_repeach <- res_aclcpg[which(res_aclcpg$Clusternumber %in% randclss), ]  ##full datainfo of randomly picked 500 clusters

  bval_res <- res_aclcpg_repeach[, beta_start_columnid:colm]  ##extract only betavalues for selected clusters


  ##### check average b-val/m-val group-wise in each cpg, and add mu value###########

  if (miuv > 0) {
    ### for all cpgs of totalori=500 random picked clusters
    for (i in 1:nrow(bval_res))
    {
      avgbetagr1 <- mean(as.numeric(bval_res[i, 1:exsampleno]))
      avgbetagr2 <- mean(as.numeric(bval_res[i, (exsampleno + 1):ncol(bval_res)]))

      if (avgbetagr1 > avgbetagr2) {
        bval_res[i, 1:exsampleno] <- as.numeric(bval_res[i, 1:exsampleno]) + miuv
      } else {
        bval_res[i, (exsampleno + 1):ncol(bval_res)] <- as.numeric(bval_res[i, (exsampleno + 1):ncol(bval_res)]) + miuv
      }
    }
    tt125 = as.matrix(bval_res)
    tt125[tt125 == 1] <- 0.999  ###if resultant bval_res after adding miuv is greater than 1, replace it by 1 (max beta score)


    tt125[tt125 > 1] <- 0.999

    bval_res <- as.data.frame(tt125)

  }

  res_aclcpg_repeach$actual <- "positive"  ##newly added

  ## newly added two lines
  res_aclcpg_repeach_remain <- res_aclcpg[!(rownames(res_aclcpg) %in% res_aclcpg_repeach$cpg), ]  ##full datainfo of remaining cpgs, i.e. except randomly picked 500 clusters
  res_aclcpg_repeach_remain$actual <- "negative"  ##newly added

  bval_res_remain <- beta.value.all[!(rownames(beta.value.all) %in% res_aclcpg_repeach$cpg), ]  ##full datainfo (betavalues) of remaining cpgs, i.e. except randomly picked 500 clusters

  ## newly added the line
  res_aclcpg.res.with_act_inter_stat <- rbind(res_aclcpg_repeach, res_aclcpg_repeach_remain)  ###combine the cpgs (belonging to selected clusters) and remianing cpgs (belonging to non-selected clusters)

  bval.res <- rbind(bval_res, bval_res_remain)  ###combine the mu-value addeded cpgs (belonging to selected clusters) and mu-value not added cpgs (belonging to non-selected clusters)

  ### reorder the 'actual' column to previous of betavalues ##newly added
  res_aclcpg.res.with_act_inter_stat <- as.data.frame(res_aclcpg.res.with_act_inter_stat)

  data_impvars_aclust <- res_aclcpg.res.with_act_inter_stat[c("Clusternumber", "cpg", "CHR", "MAPINFO", "start_position", "end_position",
                                                              "coordinate_37", "chromosome", "actual")]

  data_rest_aclust <- res_aclcpg.res.with_act_inter_stat[setdiff(names(res_aclcpg.res.with_act_inter_stat), c("Clusternumber", "cpg",
                                                                                                              "CHR", "MAPINFO", "start_position", "end_position", "coordinate_37", "chromosome", "actual"))]

  res_aclcpg.res.with_act_inter_stat_ordered <- cbind(data_impvars_aclust, data_rest_aclust)

  ### reorder the rows (cpgs) according to clusternumber
  res_aclcpg.res.with_act_inter_stat_ordered2 <- res_aclcpg.res.with_act_inter_stat_ordered[order(res_aclcpg.res.with_act_inter_stat_ordered$Clusternumber),
                                                                                            ]

  # Added by Gabriel O. on 2018-06-26 in order to match output for cleaned data
  if(trimNames){

    tt12 <- as.character(colnames(bval_res))
    v <- substr(colnames(bval_res), 7, nchar(tt12))  ######deleting 'GSE100' from the column names
    colnames(bval.res) <- v

  }

  if (miuv == 0) {
    res_aclcpg.res.with_act_inter_stat_ordered2$actual = "negative"
  }

  return(list(bval.res, res_aclcpg.res.with_act_inter_stat_ordered2))  ##updated line
}

# Test
data("betaVals_mat")
data("startEndCPG_df")

a <- Sys.time()
oldRes <- simulateData(beta.value.all = betaVals_mat,
                       res_aclcpg = startEndCPG_df,
                       beta_start_columnid = 9,
                       exsampleno = 7,
                       miuv = 0.4,
                       seed_value = 100,
                       totalori = 500)
Sys.time() - a # 2.631118 sec

b <- Sys.time()
newRes <- SimulateData(beta_mat = betaVals_mat,
                       Aclusters_df = startEndCPG_df,
                       delta_num = 0.4,
                       seed_int = 100)
Sys.time() - b # 1.441503 sec: over 1.8 times faster

test_that("New simulation data equal to legacy simulation data", {

  expect_equal(oldRes[[1]], newRes[[1]])
  expect_equal(oldRes[[2]], newRes[[2]])

})

