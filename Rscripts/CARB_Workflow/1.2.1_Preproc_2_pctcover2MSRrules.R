# pctcover2MSRrules
# 
# Uses an input basin boundary map and set of lifeform percent cover maps (tree, shrub, herb, other)
# using these inputs it creates a set of MSR rules and a corresponding map
library(rhutils)
library(terra)
library(data.table)
source("r/4.0_MSR_thinning_functions.R")

# ------------------------------ INPUTS ------------------------------
lpc_data_path = "../data/Landcover/2019/" 
lpc_pattern = "h003v008"
mask_map_path = "preprocessing/spatial90m/basin.tif"
round_value = 10
output_rules_map = "preprocessing/spatial90m/rules_90m.tif"
output_rules_file = "preprocessing/rules/LPC_90m.rules"
# --------------------------------------------------------------------

# If your basin covers multiple lpc tiles you may need to edit this function
lpc_map_list = lpc_extract(lpc_data_path = lpc_data_path, lpc_pattern = lpc_pattern, mask_map_path = mask_map_path)

lpc_df = lpc_as_df(lpc_map_list)

lpc_bin_df = lpc_bin(lpc_df = lpc_df, round_value = round_value)

outlist = make_std_rules(lpc_bin_df)
lpc_rules_df = outlist[[1]] # cover data and rule ID for each patch
unique_rules = outlist[[2]] # unique cover combinations/rules for the basin

rules_map_std = make_rules_map(lpc_rules_df, mask_map_path, lpc_map_list)
writeRaster(rules_map_std, output_rules_map, overwrite = T)

write_rules_file(
  unique_rules,
  output_rules_file,
  vegid_tree = 1,
  vegid_shrub = 5,
  vegid_herb = 3,
  vegid_other = 4,
  strata_ct = NULL
)







# ---------------------------OLD CODE EVERYTHING IS IN FUNCTIONS NOW ------------------------------
old = F
if (old) {
  # ------------------------------ IMPORT DATA - EDIT HERE ------------------------------
  write_files = F
  # OUTPUT MAP NAMES
  # veg cover map to be used with standard (non-MSR) RHESSys
  output_std_vegcover_map = "preprocessing/spatial90m/LPC_veg_cover.tif"
  # map of (patch family MSR) rules based on LPC data, distinct rules defining differences in cover combinations (tree, shrub, herb, other)
  output_rules_map = "preprocessing/spatial90m/rules_90m.tif"
  # rules file defining the patch family rules used in the above map, input when running RHESSysPreprocessing for MSR world setup
  output_rules_file = "preprocessing/rules/LPC_90m.rules"
  
  # mask map eg basin map
  mask_map = rast("preprocessing/spatial90m/basin.tif")
  mask_vect = as.polygons(mask_map)
  
  # tc tree cover, sc shrub cover, hc herb cover, oc other cover
  tc_in = rast("../data/Landcover/2019/LPC_h003v008_TC_2019_V01.tif")
  sc_in = rast("../data/Landcover/2019/LPC_h003v008_SC_2019_V01.tif")
  hc_in = rast("../data/Landcover/2019/LPC_h003v008_HC_2019_V01.tif")
  oc_in = rast("../data/Landcover/2019/LPC_h003v008_OC_2019_V01.tif")
  
  # IF BASIN COVERS TWO GRIDS
  # tc_in1 = rast("../data/Landcover/2019/LPC_h003v009_TC_2019_V01.tif")
  # sc_in1 = rast("../data/Landcover/2019/LPC_h003v009_SC_2019_V01.tif")
  # hc_in1 = rast("../data/Landcover/2019/LPC_h003v009_HC_2019_V01.tif")
  # oc_in1 = rast("../data/Landcover/2019/LPC_h003v009_OC_2019_V01.tif")
  # 
  # tc_in2 = rast("../data/Landcover/2019/LPC_h003v010_TC_2019_V01.tif")
  # sc_in2 = rast("../data/Landcover/2019/LPC_h003v010_SC_2019_V01.tif")
  # hc_in2 = rast("../data/Landcover/2019/LPC_h003v010_HC_2019_V01.tif")
  # oc_in2 = rast("../data/Landcover/2019/LPC_h003v010_OC_2019_V01.tif")
  # 
  # tc_in = terra::merge(tc_in1, tc_in2)
  # sc_in = terra::merge(sc_in1, sc_in2)
  # hc_in = terra::merge(hc_in1, hc_in2)
  # oc_in = terra::merge(oc_in1, oc_in2)
  
  # ------------------------------ STOP EDIT ------------------------------
  mask_proj = project(mask_vect, tc_in)
  
  plot(tc_in$LPC_h003v008_TC_2019_V01_1)
  plot(mask_proj, add=T)
  
  # project
  tc_proj = project(tc_in, mask_map)
  sc_proj = project(sc_in, mask_map)
  hc_proj = project(hc_in, mask_map)
  oc_proj = project(oc_in, mask_map)
  # crop 
  tc_crop = crop(tc_proj, mask_map)
  sc_crop = crop(sc_proj, mask_map)
  hc_crop = crop(hc_proj, mask_map)
  oc_crop = crop(oc_proj, mask_map)
  # mask everything
  tc = mask(tc_crop, mask_map)
  sc = mask(sc_crop, mask_map)
  hc = mask(hc_crop, mask_map)
  oc = mask(oc_crop, mask_map)
  
  
  # write maps - REPEAT FOR DIFFERENT RES IF NEEDED
  if (write_files) {
    writeRaster(tc, "preprocessing/spatial_source/LPC_basin_TC_2019_V01.tif", overwrite = T)
    writeRaster(sc, "preprocessing/spatial_source/LPC_basin_SC_2019_V01.tif", overwrite = T)
    writeRaster(hc, "preprocessing/spatial_source/LPC_basin_HC_2019_V01.tif", overwrite = T)
    writeRaster(oc, "preprocessing/spatial_source/LPC_basin_OC_2019_V01.tif", overwrite = T)
  }
  
  tc = rast("preprocessing/spatial_source/LPC_basin_TC_2019_V01.tif")
  sc = rast("preprocessing/spatial_source/LPC_basin_SC_2019_V01.tif")
  hc = rast("preprocessing/spatial_source/LPC_basin_HC_2019_V01.tif")
  oc = rast("preprocessing/spatial_source/LPC_basin_OC_2019_V01.tif")
  
  mean_dt = data.table(values(tc$LPC_1)[!is.na(values(tc$LPC_1))],
                       values(sc$LPC_1)[!is.na(values(sc$LPC_1))],
                       values(hc$LPC_1)[!is.na(values(hc$LPC_1))],
                       values(oc$LPC_1)[!is.na(values(oc$LPC_1))])
  names(mean_dt) = c("Tree","Shrub","Herb","Other")
  
  mean_dt$id = 1:nrow(mean_dt)
  mean_dt$sum = mean_dt$Tree + mean_dt$Shrub + mean_dt$Herb + mean_dt$Other
  
  par(mfrow = c(2,2))
  hist(mean_dt$Tree)
  hist(mean_dt$Shrub)
  hist(mean_dt$Herb)
  hist(mean_dt$Other)
  
  # ------------------------------ BINNING FOR STANDARD VEG COVER ------------------------------
  mean_dt$maxcol = max.col((mean_dt[,c(1:4)]))
  #mean_dt$max_veg_par = dplyr::recode(mean_dt$maxcol, `1` = 7, `2` = 50, `3` = 71, `4` = 31)
  # new veg IDs - evergreen 1, shrub 5, grass 3, noveg 4 c(1, 5, 3, 4)
  mean_dt$max_veg_par = dplyr::recode(mean_dt$maxcol, `1` = 1, `2` = 5, `3` = 3, `4` = 4)
  
  cover_pct_mean = mean(rowSums(mean_dt[,1:3]))
  cover_pct_median = median(rowSums(mean_dt[,1:3]))
  
  if (write_files) {
    # writeLines(text = paste0("Mean pct cover:",cover_pct_mean, "\nMedian pct cover:",cover_pct_median), "preprocessing/spatial90m/LPC_coverfrac.txt")
  }
  maxrule = names(which.max(summary(as.factor(mean_dt$max_veg_par))))
  
  # if these have different lengths/missing values, choose the one missing the most
  LPC_veg_cover_map = tc$LPC_1
  length(values(LPC_veg_cover_map)[!is.na(values(LPC_veg_cover_map))]) == length(mean_dt$max_veg_par)
  
  values(LPC_veg_cover_map)[!is.na(values(LPC_veg_cover_map))] <- mean_dt$max_veg_par
  values(LPC_veg_cover_map)[is.na(values(LPC_veg_cover_map)) & !is.na(values(mask_map))] = as.numeric(maxrule)
  names(LPC_veg_cover_map) = "LPC_veg_cover"
  
  # ------------------------------ OUTPUT STANDARD VEG COVER MAP ------------------------------
  if (write_files) {
    writeRaster(LPC_veg_cover_map, output_std_vegcover_map, overwrite = T)
  }
  
  # ------------------------------ BINNING FOR MSR RULES ------------------------------
  # binning - tree, binned to nearest 10
  mean_dt$Tree_rnd = round(mean_dt$Tree, digits = -1)
  mean_dt$Shrub_rnd = round(mean_dt$Shrub, digits = -1)
  mean_dt$Herb_rnd = round(mean_dt$Herb, digits = -1)
  mean_dt$Other_rnd = round(mean_dt$Other, digits = -1)
  
  rnd2x = F
  if (rnd2x) {
    rnd = 5
    mean_dt$Tree_rnd = round(mean_dt$Tree/rnd)*rnd
    mean_dt$Shrub_rnd = round(mean_dt$Shrub/rnd)*rnd
    mean_dt$Herb_rnd = round(mean_dt$Herb/rnd)*rnd
    mean_dt$Other_rnd = round(mean_dt$Other/rnd)*rnd
  }
  
  # IF just tree + shrub + herb > 100%, remove from herb
  mean_dt$Herb_rnd_crt =  ifelse(mean_dt$Tree_rnd + mean_dt$Shrub_rnd + mean_dt$Herb_rnd > 100, 100 - (mean_dt$Tree_rnd + mean_dt$Shrub_rnd),mean_dt$Herb_rnd)
  summary(mean_dt$Tree_rnd + mean_dt$Shrub_rnd + mean_dt$Herb_rnd)
  summary(mean_dt$Tree_rnd + mean_dt$Shrub_rnd + mean_dt$Herb_rnd_crt)
  
  # IF tree+shrub+herb+other > 100% remove from other
  mean_dt$Other_rnd_crt =  ifelse(mean_dt$Tree_rnd + mean_dt$Shrub_rnd + mean_dt$Herb_rnd_crt + mean_dt$Other_rnd > 100, 
                                  100 - (mean_dt$Tree_rnd + mean_dt$Shrub_rnd + mean_dt$Herb_rnd_crt), mean_dt$Other_rnd)
  # IF < 100%, ADD to other
  mean_dt$Other_rnd_crt =  ifelse(mean_dt$Tree_rnd + mean_dt$Shrub_rnd + mean_dt$Herb_rnd_crt + mean_dt$Other_rnd < 100, 
                                  100 - (mean_dt$Tree_rnd + mean_dt$Shrub_rnd + mean_dt$Herb_rnd_crt), mean_dt$Other_rnd_crt)
  
  summary(mean_dt$Tree_rnd + mean_dt$Shrub_rnd + mean_dt$Herb_rnd_crt + mean_dt$Other_rnd_crt)
  
  binned_dt = mean_dt[,c("id", "Tree_rnd","Shrub_rnd","Herb_rnd_crt","Other_rnd_crt")]
  
  # ------------------------------ MAKE MSR RULES ------------------------------
  binned_dt_unique = unique(binned_dt[,-1])
  binned_dt_unique$rule = 1:nrow(binned_dt_unique)
  binned_dt = binned_dt[binned_dt_unique, on = .(Tree_rnd, Shrub_rnd, Herb_rnd_crt, Other_rnd_crt),][order(id)]
  
  # IF THERES MISSING DATA
  # FILL MISSING DATA FROM PCT COVER - use most common rule, lazy choice
  maxrule = which.max(summary(as.factor(binned_dt$rule)))
  
  # if these have different lengths/missing values, choose the one missing the most
  length(unique(length(values(tc$LPC_1)), length(values(sc$LPC_1)), length(values(hc$LPC_1)),length(values(oc$LPC_1)))) == 1
  
  rules_map = tc$LPC_1
  length(values(rules_map))
  
  values(rules_map)[!is.na(values(rules_map))] <- binned_dt$rule
  values(rules_map)[is.na(values(rules_map)) & !is.na(values(mask_map))] = as.numeric(maxrule)
  
  scengen = F
  if (scengen) {
    # generate scenarios for creepgrass
    scens = list()
    cg = c(5, 10, 15, 20, 25, 30, 35, 40)
    for (i in seq_along(cg)) {
      scens[[i]] = binned_dt_unique
      scens[[i]]$BErmd = 0
      scens[[i]]$Herb_rnd_crt = scens[[i]]$Herb_rnd_crt + cg[i]
      scens[[i]]$BErmd = ifelse(scens[[i]]$Other_rnd_crt > cg[i], 0, cg[i] - scens[[i]]$Other_rnd_crt)
      scens[[i]]$Other_rnd_crt = scens[[i]]$Other_rnd_crt - cg[i]
      scens[[i]]$Other_rnd_crt[scens[[i]]$Other_rnd_crt < 0] = 0
      scens[[i]]$Shrub_rnd = scens[[i]]$Shrub_rnd - scens[[i]]$BErmd
      scens[[i]]$SHrmd = ifelse(scens[[i]]$Shrub_rnd < 0, scens[[i]]$Shrub_rnd, 0)
      scens[[i]]$Shrub_rnd[scens[[i]]$Shrub_rnd<0] = 0
      scens[[i]]$Herb_rnd_crt = scens[[i]]$Herb_rnd_crt + scens[[i]]$SHrmd
    }
  }
  
  # ------------------------------ WRITE FILES ------------------------------
  if (write_files) {
    # ------------------------------ OUTPUT MSR RULES MAP ------------------------------
    writeRaster(rules_map, output_rules_map, overwrite = T)
    #writeRaster(rules_map, "preprocessing/spatial180m/rules_180m.tif", overwrite = T)
    
    # ------------------------------ OUTPUT MSR RULES FILE ------------------------------
    file_out = output_rules_file
    fcon = file(file_out,open = "wt")
    #veg = c(7, 50, 71, 31)
    veg = c(1, 5, 3, 4)
    strata = c(1,1,1,1)
    
    for (i in 1:nrow(binned_dt_unique)) {
      include = which(binned_dt_unique[i,-5] != 0)
      text_out = c(paste0("ID\t",binned_dt_unique[i,"rule"]),
                   paste0("subpatch_count\t", length(include)),
                   paste0("_patch"),
                   paste0("\tpct_family_area\tvalue\t", paste(binned_dt_unique[i,..include]/100,collapse = " | ")),
                   paste0("\t_canopy_strata\t",paste(strata[include],collapse = " | ")),
                   paste0("\t\tveg_parm_ID\t",paste(veg[include],collapse = " | ")),
                   ""
      )
      writeLines(text = text_out, con = fcon)
    }
    close(fcon)
  }
  
  # ------------------------------ PLOTS + EVALUATION ------------------------------
  plots = F
  if (plots) {
    
    # DIFFERENCES FROM BASELINE
    mean_dt$tree_diff = mean_dt$Tree_rnd - mean_dt$Tree
    mean_dt$shrub_diff = mean_dt$Shrub_rnd - mean_dt$Shrub
    mean_dt$herb_diff = mean_dt$Herb_rnd_crt - mean_dt$Herb
    mean_dt$other_diff = mean_dt$Other_rnd_crt - mean_dt$Other
    
    par(mfrow = c(2,2))
    hist(mean_dt$tree_diff)
    hist(mean_dt$shrub_diff)
    hist(mean_dt$herb_diff)
    hist(mean_dt$other_diff)
    
    # PLOTS for report
    mean_cover = c(tc[[1]], sc[[1]], hc[[1]], oc[[1]])
    mean_cover = trim(mean_cover)
    names(mean_cover) = c("Tree Cover", "Shrub Cover", "Herb Cover", "Other Cover")
    plot(mean_cover)
    plot(trim(rules_map))
    #plot(trim(rules_map2))
    
    diff_cover = mean_cover
    names(diff_cover) = c("Tree Cover Change", "Shrub Cover Change", "Herb Cover Change", "Other Cover Change")
    values(diff_cover[[1]])[!is.na(values(diff_cover[[1]]))] = mean_dt$tree_diff
    values(diff_cover[[2]])[!is.na(values(diff_cover[[2]]))] = mean_dt$shrub_diff
    values(diff_cover[[3]])[!is.na(values(diff_cover[[3]]))] = mean_dt$herb_diff
    values(diff_cover[[4]])[!is.na(values(diff_cover[[4]]))] = mean_dt$other_diff
    
    plot(diff_cover)
    
    summary(mean_dt)
    # investigate some of the edge cases
    mean_dt[mean_dt$other_diff >= 13 | mean_dt$other_diff <= -13,]
    
    mean_dt[mean_dt$other_diff >= 6.5 | mean_dt$other_diff <= -6.5,]
    
    # When rounding to 10's I get some reaching up to ~ +/- 15 for the other cover, but for tree, shrub, and herb eveerything is < +/- 5
    # the cause of the larger differences looks to be cases where multiple initial values are just on the edge in terms
    # of rounding, and all go the same direction, 
    # EXAMPLE:
    # tree, shrub, herb, and other have starting values: 26.16355 15.63966 26.22934 33.457813
    # rounding to the nearest 10, results in: 30        20       30        30, totalling 110
    # Since other is the first category to modify, it gets reduced to 20 to make the total cover correct
    # combined with the initial reduction from roudning, the total change to other is -13.4...
    # All the largest changes follow this trend of multiple values being rounded notably in the same direction
    # Also, initial values may not originally sum to 100 which can add to these errors
    # 
    # When rounding to 5s, the errors get smaller, max change again is for other cover, and isfrom -8.7 to + 9
    # These are smaller than the errors when rounding to 10, but not scaling (eg the errors don't nevessarily get smaller at the same rate as the rounding)
    
  }
  
  
  # ------------------------------ Write Rules function ------------------------------
  write_rules_file = function(rules, file_out, veg_ids, strata_ct) {
    fcon = file(file_out,open = "wt")
    vals = as.data.frame(rules)[!names(rules) %in% c("id", "rule")]
    
    for (i in 1:nrow(rules)) {
      inc_cols = which(vals[i, ] != 0)
      text_out = c(paste0("ID\t",rules[i,"rule"]),
                   paste0("subpatch_count\t", length(inc_cols)),
                   paste0("_patch"),
                   paste0("\tpct_family_area\tvalue\t", paste(vals[i,inc_cols]/100,collapse = " | ")),
                   paste0("\t_canopy_strata\t",paste(strata_ct[inc_cols],collapse = " | ")),
                   paste0("\t\tveg_parm_ID\t",paste(veg_ids[inc_cols],collapse = " | ")),
                   ""
      )
      writeLines(text = text_out, con = fcon)
    }
    close(fcon)
  }
  
  # ================================================================================
  # ------------------------------ MSR THINNING RULES ------------------------------
  # ================================================================================
  thin = F
  if (thin) {
    thin_dt = binned_dt_unique[,c("rule","Tree_rnd","Shrub_rnd","Herb_rnd_crt","Other_rnd_crt")]
    names(thin_dt) = c("rule","Tree","Shrub","Herb","Other")
    
    # GOAL: Remove 40% overstory biomass, 90% understory
    # from each rule type with overstory, add a second overstory patch with RELATIVE coverage split of 40/60
    # same for understory, but 90/10
    
    # these are modifications of the PERCENT COVER
    thin_dt$Tree_treat_thin = thin_dt$Tree * 0.4
    thin_dt$Tree_treat_remain = thin_dt$Tree * 0.6
    thin_dt$Tree = NULL
    
    thin_dt$Shrub_treat_thin = thin_dt$Shrub * 0.9
    thin_dt$Shrub_treat_remain = thin_dt$Shrub * 0.1
    thin_dt$Shrub = NULL
    
    thin_dt = thin_dt[,c("rule","Tree_treat_thin","Tree_treat_remain","Shrub_treat_thin","Shrub_treat_remain",
                         "Herb", "Other")]
    veg = c(1, 1, 5, 5, 3, 4)
    strata = rep(1,length(veg))
    file_out = "preprocessing/rules/LPC_90m_thin1.rules"
    
    write_rules_file(thin_dt, file_out, veg, strata)
    
  }
  
  mean_dt$HO_comb_rnd = round(mean_dt$Other + mean_dt$Herb, digits = -1)
  # VERSION 2 with COMBINED HERB AND OTHER
  mean_dt$HO_comb_rnd_crt =  ifelse(mean_dt$Tree_rnd + mean_dt$Shrub_rnd + mean_dt$HO_comb_rnd > 100, 100 - (mean_dt$Tree_rnd + mean_dt$Shrub_rnd), mean_dt$HO_comb_rnd)
  summary(mean_dt$Tree_rnd + mean_dt$Shrub_rnd +  mean_dt$HO_comb_rnd)
  summary(mean_dt$Tree_rnd + mean_dt$Shrub_rnd + mean_dt$HO_comb_rnd_crt)
  
  mean_dt$HO_comb_rnd_crt =  ifelse(mean_dt$Tree_rnd + mean_dt$Shrub_rnd + mean_dt$HO_comb_rnd_crt < 100, 
                                    100 - (mean_dt$Tree_rnd + mean_dt$Shrub_rnd), mean_dt$HO_comb_rnd_crt)
  summary(mean_dt$Tree_rnd + mean_dt$Shrub_rnd + mean_dt$HO_comb_rnd_crt)
  
  # 36 bins instead of 85
  binned_dt2 = mean_dt[,c("id", "Tree_rnd","Shrub_rnd","HO_comb_rnd_crt")]
  binned_dt_unique2 = unique(binned_dt2[,-1])
  binned_dt_unique2$rule = 1:nrow(binned_dt_unique2)
  binned_dt2 = binned_dt2[binned_dt_unique2, on = .(Tree_rnd, Shrub_rnd, HO_comb_rnd_crt),][order(id)]
  
  # IF THERES MISSING DATA
  # FILL MISSING DATA FROM PCT COVER - use most common rule, lazy choice
  maxrule2 = which.max(summary(as.factor(binned_dt2$rule)))
  
  rules_map2 = tc$LPC_1
  length(values(rules_map2))
  values(rules_map2)[!is.na(values(rules_map2))] <- binned_dt2$rule
  values(rules_map2)[is.na(values(rules_map2)) & !is.na(values(mask_map))] = as.numeric(maxrule2)
  
  file_out = "preprocessing/rules/ward30m_v2.rules"
  fcon = file(file_out,open = "wt")
  veg = c(7, 50, 71)
  strata = c(1,1,1)
  
  for (i in 1:nrow(binned_dt_unique2)) {
    include = which(binned_dt_unique2[i,-4] != 0)
    text_out = c(paste0("ID\t",binned_dt_unique[i,"rule"]),
                 paste0("subpatch_count\t", length(include)),
                 paste0("_patch"),
                 paste0("\tpct_family_area\tvalue\t", paste(binned_dt_unique2[i,..include]/100,collapse = " | ")),
                 paste0("\t_canopy_strata\t",paste(strata[include],collapse = " | ")),
                 paste0("\t\tveg_parm_ID\t",paste(veg[include],collapse = " | ")),
                 ""
    )
    writeLines(text = text_out, con = fcon)
  }
}

