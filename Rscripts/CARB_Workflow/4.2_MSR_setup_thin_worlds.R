# MSR_4_init_thin_worlds
library(RHESSysIOinR)
library(rhutils)
library(tidyverse)
library(data.table)


thinrules = c(2,3,4,5,6)

template = "preprocessing/template/carb_msr.template"

type = "raster"
typepars = "preprocessing/spatial90m/"
overwrite = TRUE
streams = "streams.tif"
unique_strata_ID = TRUE
seq_patch_IDs = F

for (i in seq_along(thinrules)) {
  name = paste0("preprocessing/preprocess_out/Ward_msr90m_thin",thinrules[i])
  asprules = paste0("preprocessing/rules/LPC_90m_thin",thinrules[i],".rules")
  
  RHESSysPreprocess(template = template,
                    name = name,
                    type = type,
                    typepars = typepars,
                    streams = streams,
                    overwrite = overwrite,
                    asprules = asprules,
                    unique_strata_ID = unique_strata_ID,
                    seq_patch_IDs = seq_patch_IDs)
  
  file.copy(from = paste0(name,".world"), to = file.path("worldfiles",basename(paste0(name,".world"))), overwrite = T)
}



# This script transposes values between two MSR worldfiles, with the input worldfile containing 
# spun-up values you want to move to (a copy of) the second. This transposition is based on veg paramater ID,
# and assumes the input world only has 1 of each veg parameter, and then transposes each set of patch and strata
# variables based on veg parm, with the output/destination worldfile having 1 or more of each veg parm. Can 
# account for variable patch family compositions within the world.
# 
# Ex: 
# input world patches: ID 1, 2, 3          | veg parms 11, 22, 33
# output world patches: ID 1, 2, 3, 4, 5   | veg parms 11, 11, 22, 22, 33
# input world patch and strata vars from ID 1 get transposed to output world patches 1 & 2
# 

# ==================== START I/O ====================
# world_in_path = "worldfiles/Ward_msr90m_soilvegspun.world"
# world_out_path = "worldfiles/Ward_msr90m_thin2.world"
# world_output_filename = "worldfiles/Ward_msr90m_thin2_init.world"
# ==================== END I/O ====================

# --------- shift vars from soil spun up world to world w 4 msr patches -----------
shift_msr_world_vals = function(world_in_path, world_out_path, world_output_filename) {
  
  # ---------- read worldfiles ----------
  world_in = as.data.table(read_world(worldfile = world_in_path))
  world_out = as.data.table(read_world(worldfile = world_out_path))

  # ---------- vars to change ----------
  vars_zonehill = c("gw.storage", "gw.NO3")
  #vars_zone = c("gw.storage", "gw.NO3", "e_horizon", "w_horizon")
  vars_list = c("Ksat_vertical", "mpar", "rz_storage", "unsat_storage", "sat_deficit", "snowpack.water_equivalent_depth", "snowpack.water_depth", "snowpack.T", 
                "snowpack.surface_age", "snowpack.energy_deficit", "litter.cover_fraction","litter.rain_stored", "litter_cs.litr1c", "litter_ns.litr1n", 
                "litter_cs.litr2c", "litter_cs.litr3c", "litter_cs.litr4c", "soil_cs.soil1c", "soil_ns.sminn", "soil_ns.nitrate", "soil_cs.soil2c", "soil_cs.soil3c", 
                "soil_cs.soil4c", "rootzone.depth", "snow_stored", "rain_stored", "cs.cpool", "cs.leafc", "cs.dead_leafc", 
                "cs.leafc_store", "cs.leafc_transfer", "cs.live_stemc", "cs.livestemc_store", "cs.livestemc_transfer", "cs.dead_stemc", "cs.deadstemc_store", 
                "cs.deadstemc_transfer", "cs.live_crootc", "cs.livecrootc_store", "cs.livecrootc_transfer", "cs.dead_crootc", "cs.deadcrootc_store", 
                "cs.deadcrootc_transfer", "cs.frootc", "cs.frootc_store", "cs.frootc_transfer", "cs.cwdc", "epv.prev_leafcalloc", "ns.npool", "ns.leafn", "ns.dead_leafn", 
                "ns.leafn_store", "ns.leafn_transfer", "ns.live_stemn", "ns.livestemn_store", "ns.livestemn_transfer", "ns.dead_stemn", "ns.deadstemn_store", 
                "ns.deadstemn_transfer", "ns.live_crootn", "ns.livecrootn_store", "ns.livecrootn_transfer", "ns.dead_crootn", "ns.deadcrootn_store", 
                "ns.deadcrootn_transfer", "ns.frootn", "ns.frootn_store", "ns.frootn_transfer", "ns.cwdn", "ns.retransn", "epv.wstress_days", "epv.max_fparabs", 
                "epv.min_vwc")
  
  # ---------- check patch and family counts ----------
  in_pID = unique(world_in[world_in$level == "patch", "ID"])
  #in_famID = unique(world_in[world_in$vars == "family_ID", "values"])
  # ^^ DOESNT WORK BECAUSE FAM ID WASNT OUTPUT W WORLD STATE
  in_famID = unique(substr(in_pID$ID, 0, nchar(in_pID$ID) - 2))
  out_pID = unique(world_out[world_out$level == "patch", "ID"])
  out_famID = unique(world_out[world_out$level == "patch" & world_out$vars == "family_ID", "values"])
  cat("          Input | Output\nPatches:", length(in_pID$ID), " | ", length(out_pID$ID), "\nFamilies:", length(in_famID), " | ", length(out_famID$values))

  # ---------- add family ID column ----------
  world_in$pID = world_in$ID
  world_in$pID[!(world_in$level == "patch" | world_in$level == "canopy_strata")] = NA
  world_in$pID[!is.na(world_in$pID) & world_in$level == "canopy_strata"] = 
    substr(world_in$pID[!is.na(world_in$pID) & world_in$level == "canopy_strata"], 0, nchar(world_in$pID[!is.na(world_in$pID) & world_in$level == "canopy_strata"]) - 1 )
  world_in$famID = world_in$pID
  world_in$famID[!is.na(world_in$famID)] = substr(world_in$famID[!is.na(world_in$famID)], 0, nchar(world_in$famID[!is.na(world_in$famID)]) - 2 )
  
  world_out$pID = world_out$ID
  world_out$pID[!(world_out$level == "patch" | world_out$level == "canopy_strata")] = NA
  world_out$pID[!is.na(world_out$pID) & world_out$level == "canopy_strata"] = 
    substr(world_out$pID[!is.na(world_out$pID) & world_out$level == "canopy_strata"], 0, nchar(world_out$pID[!is.na(world_out$pID) & world_out$level == "canopy_strata"]) - 1 )
  world_out$famID = world_out$pID
  world_out$famID[!is.na(world_out$famID)] = substr(world_out$famID[!is.na(world_out$famID)], 0, nchar(world_out$famID[!is.na(world_out$famID)]) - 2 )
  
  # ---------- add vegparm column  ----------
  in_vegparm = world_in[world_in$vars == "veg_parm_ID", "values"]
  in_vegparm$values = as.numeric(in_vegparm$values)
  in_vegparm$values[in_vegparm$values > 20] = in_vegparm$values[in_vegparm$values > 20] - 20
  out_patches = unique(world_out[world_out$level == "patch", "ID"])
  in_patches = unique(world_in[world_in$level == "patch", "ID"])
  in_veg_patches = data.table(in_vegparm, in_patches)
  names(in_veg_patches) = c("vegparm", "pID")
  world_in = merge.data.table(world_in, in_veg_patches, by = "pID", all.x = TRUE, sort = F)
  
  out_vegparm = world_out[world_out$vars == "veg_parm_ID", "values"]
  out_vegparm$values = as.numeric(out_vegparm$values)
  out_vegparm$values[out_vegparm$values > 20] = out_vegparm$values[out_vegparm$values > 20] - 20
  out_patches = unique(world_out[world_out$level == "patch", "ID"])
  out_veg_patches = data.table(out_vegparm, out_patches)
  names(out_veg_patches) = c("vegparm", "pID")
  world_out = merge.data.table(world_out, out_veg_patches, by = "pID", all.x = TRUE, sort = F)
  
  # ---------- init world for merge ----------
  world_merge = world_out
  
  # ---------- merge zone vars ----------
  world_merge[world_merge$level %in% c("zone", "hillslope")  & world_merge$vars %in% vars_zonehill, c("values","vars","level","ID")] = 
    merge.data.table(world_out[world_out$level %in% c("zone", "hillslope")  & world_out$vars %in% vars_zonehill, c("vars","level","ID")], 
                     world_in[world_in$level %in% c("zone", "hillslope")  & world_in$vars %in% vars_zonehill, c("values","vars","level","ID")], 
                     by = c("level", "ID", "vars"), sort = F, all.x = T)
  
  # ---------- merge patch+strata vars ----------
  world_merge[!is.na(world_merge$famID) & world_merge$vars %in% vars_list, c("famID", "vegparm", "vars", "level","ID","unique_ID","pID","values")] = 
    merge.data.table(world_out[!is.na(world_out$famID) & world_out$vars %in% vars_list,c("vars","level","ID","unique_ID","famID","pID","vegparm")], 
                     world_in[!is.na(world_in$famID) & world_in$vars %in% vars_list,c("values","vars","famID","vegparm")], 
                     by = c("famID", "vegparm", "vars"), sort = F, all.x = T)
  
  # ---------- Write output world ----------
  write_world(world_merge, world_output_filename)
  
  return("Success!")
}

shift_msr_world_vals("worldfiles/Ward_msr90m_soilvegspun.world", "worldfiles/Ward_msr90m_thin2.world", "worldfiles/Ward_msr90m_thin2_init.world")
shift_msr_world_vals("worldfiles/Ward_msr90m_soilvegspun.world", "worldfiles/Ward_msr90m_thin3.world", "worldfiles/Ward_msr90m_thin3_init.world")
shift_msr_world_vals("worldfiles/Ward_msr90m_soilvegspun.world", "worldfiles/Ward_msr90m_thin4.world", "worldfiles/Ward_msr90m_thin4_init.world")
shift_msr_world_vals("worldfiles/Ward_msr90m_soilvegspun.world", "worldfiles/Ward_msr90m_thin5.world", "worldfiles/Ward_msr90m_thin5_init.world")
shift_msr_world_vals("worldfiles/Ward_msr90m_soilvegspun.world", "worldfiles/Ward_msr90m_thin6.world", "worldfiles/Ward_msr90m_thin6_init.world")


