# fun_spinup


# add level index to worldfile dataframe
world_add_level_i = function(world) {
  # add numeric level col
  world$level_i = world$level
  l = c("world","basin","hillslope","zone","patch","canopy_strata")
  li = c(1,2,3,4,5,6)
  world[.(from = l, to = li), on = paste0("level_i", "==from"), ("level_i") := i.to]
  world$level_i = as.numeric(world$level_i)
  world$i = 1:length(world$vars)
  return(world)
}

# add vegparmID col to worldfile dataframe
world_add_patch_vegparmIDcol = function(world) {
  # ---------- add patch and vegparm columns ----------
  world$pID = world$ID
  world$pID[!(world$level == "patch" | world$level == "canopy_strata")] = NA
  world$pID[!is.na(world$pID) & world$level == "canopy_strata"] = 
    substr(world$pID[!is.na(world$pID) & world$level == "canopy_strata"], 0, nchar(world$pID[!is.na(world$pID) & world$level == "canopy_strata"]) - 1 )
  
  in_vegparm = world[world$vars == "veg_parm_ID", "values"]
  in_patches = unique(world[world$level == "patch", "ID"])
  in_veg_patches = data.table(in_vegparm, in_patches)
  names(in_veg_patches) = c("vegparm", "pID")
  world = merge.data.table(world, in_veg_patches, by = "pID", all.x = TRUE, sort = F)
  return(world)
}

# extract a world based on a target unique index ID
extract_world = function(world, target_unique_ID) {
  # use to look up parent/child levels
  IDworld = world[
    which(world$vars == "world_ID" | world$vars == "basin_ID" | world$vars == "hillslope_ID" |
            world$vars == "zone_ID" | world$vars == "patch_ID" | world$vars == "canopy_strata_ID"), 
    c("pID", "ID", "unique_ID", "level_i", "vegparm", "i")]
  
  target = IDworld[unique_ID == target_unique_ID,]
  # Parents should only be single levels, ie if target is a patch, you only need the single encompassing zone
  l_parent = c(1:6)[c(1:6) < target$level_i]
  unid_parent = rep(NA, length(l_parent))
  for (l in l_parent) {
    unid_parent[l] = IDworld[level_i == l, ][i < target$i,][which.min(abs(target$i - i)), "unique_ID"]
  }
  # child could be multiple units, need everything until the next element at the same level as target
  l_child = c(1:6)[c(1:6) > target$level_i]
  unid_child = list()
  max_i = IDworld[IDworld$level_i == (target$level_i) & i > target$i, min(i)]
  for (l in  seq_along(l_child)) {
    unid_child[l] = IDworld[level_i == l_child[l], ][i > target$i & i < max_i, "unique_ID"]
  }
  unid_extract = sort(c(unlist(unid_parent), target_unique_ID, unlist(unid_child)))
  world_extract = world[unique_ID %in% unid_extract,]
  IDworld_extract = IDworld[IDworld$unique_ID %in% unid_extract, ]
  # correct number of sub units and zone areas
  IDworld_extract$ct = NA
  zoneareas = data.table(old = world_extract[world_extract$level == "zone" & world_extract$vars == "area", values])
  zoneareas$new = NA
  zi = 0
  for (ind in which(IDworld_extract$level_i < 6)) {
    max_i = IDworld_extract[IDworld_extract$level_i == IDworld_extract$level_i[ind] & i > IDworld_extract$i[ind], suppressWarnings(min(i))]
    IDworld_extract$ct[ind] = length(IDworld_extract[IDworld_extract$level_i == (IDworld_extract$level_i[ind] + 1),][i > IDworld_extract$i[ind] & i < max_i, i])
    if (IDworld_extract$level_i[ind] == 4) {
      zi = zi + 1
      zoneareas$new[zi] = sum(as.numeric(world_extract[world_extract$vars == "area" & 
                                                         world_extract$level_i == (IDworld_extract$level_i[ind] + 1),][i > IDworld_extract$i[ind] & i < max_i, values]))
    }
  }
  world_extract$values[grepl("^num_\\w", world_extract$vars)] = as.character(IDworld_extract$ct[which(IDworld_extract$level_i < 6)])
  world_extract[world_extract$level == "zone" & world_extract$vars == "area", "values"] = as.character(zoneareas$new)
  
  return(world_extract)
}


make_1pflow = function(world) {
  hill = unique(world$ID[world$level == "hillslope"])
  patch = unique(world$pID[!is.na(world$pID)])
  zone = unique(world$ID[world$level == "zone"])
  xcen = world$values[world$level == "patch" & world$vars == "x"]
  ycen = world$values[world$level == "patch" & world$vars == "y"]
  zcen = world$values[world$level == "patch" & world$vars == "z"]
  area = world$values[world$level == "patch" & world$vars == "area"]
  flow_out = paste0(length(hill),"\n",
                    hill, "\t", length(patch),"\n",
                    patch, "\t", zone, "\t", hill, "\t", xcen, "\t", ycen, "\t", zcen, "\t", area, "\t", area, "\t", 1, "\t", area, "\t", 0)
  return(flow_out)
}

# # ASSUMES 1 HILLSLOPE
# make_1pflow = function(world) {
#   hill = unique(world$ID[world$level == "hillslope"])
#   patch = unique(world$pID[!is.na(world$pID)])
#   zone = unique(world$ID[world$level == "zone"])
#   xcen = world$values[world$level == "patch" & world$vars == "x"]
#   ycen = world$values[world$level == "patch" & world$vars == "y"]
#   zcen = world$values[world$level == "patch" & world$vars == "z"]
#   area = world$values[world$level == "patch" & world$vars == "area"]
#   flow_out = paste0(length(hill),"\n",
#                     hill, "\t", length(patch),"\n",
#                     patch, "\t", zone, "\t", hill, "\t", xcen, "\t", ycen, "\t", zcen, "\t", area, "\t", area, "\t", 1, "\t", area, "\t", 0)
#   return(flow_out)
# }