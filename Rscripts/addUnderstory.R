#title: "AddUnderstory"
#authors: "ErinJHanan and RyanBart"

# Read in worldfile
setwd("...")
jc=read.table("world.x", head=F) # enter the name of your worldfile here

# Add Understory for Pine and Deciduous
# Veg types            Understory
# Pine = 42            shrub_understory = 50
# Grass = 71           noveg_understory = 51
# Shrub = 52           noveg_understory = 51
# Deciduous = 41       shrub_understory = 49
# Water = 11           noveg_understory = 51
# Rock = 31            noveg_understory = 51


# Determine highest canopy_strata_ID so understory IDs will start higher
# In worldfile:
# :v/canopy_strata_ID/d   #deletes all lines that don't contain canopy_strata_ID
# w canopy_strata_ID_world.x  #use your worldfile name as suffix
# Read into R:
canopy_ID=read.table("canopy_strata_ID_world.x", head=F) #use the same name as line 21 above
summary(canopy_ID)
# Max = 235421


# ---------
# Modified from Ryan Bart:
# Find every patch that has pine, add understory, then repeat routine (below) for each veg type
# need to change in/out file and veg IDs for each cycle
# Add an understory canopy to a worldfile

world_name_in <- "input_worldfile"
world_name_out <- "output_worldfile"

canopy_default_ID <- 42
understory_canopy_default_ID <- 50
understory_canopy_strata_ID <- 1000000    # should be larger than highest canopy_strata_ID

add_understory <- function(world_name_in,
                           world_name_out,
                           canopy_default_ID,
                           understory_canopy_default_ID,
                           understory_canopy_strata_ID){
  
  # Read in worldfile
  worldfile <- read.table(world_name_in, header = FALSE, stringsAsFactors = FALSE)
  
  change_flag <- 0
  aa <- 1
  while(aa < nrow(worldfile)){
    #for (aa in seq_len(nrow(worldfile))){
    if (aa%%1000 == 0 ){print(paste("Worldfile line", aa, "of", nrow(worldfile)))} # Counter
    
    if (worldfile[aa,2] == "canopy_strata_ID" && (worldfile[aa+1,1] == canopy_default_ID)){
      change_flag <- 1
      worldfile[aa-1,1] <- 2      # Changes num_canopy_strata from 1 to 2
    }
    
    if (change_flag == 1 && worldfile[aa,2] == "n_basestations"){
      
      # Make space for chunk
      worldfile[seq(aa+55,nrow(worldfile)+55),] <- worldfile[seq(aa,nrow(worldfile)),]
      
      worldfile[aa+1,] <- data.frame(V1=understory_canopy_strata_ID,V2="canopy_strata_ID", stringsAsFactors=FALSE)[1,]
      worldfile[aa+2,] <- data.frame(V1=understory_canopy_default_ID,V2="veg_parm_ID", stringsAsFactors=FALSE)[1,]
      worldfile[aa+3,] <- data.frame(V1=0.1,V2="cover_fraction", stringsAsFactors=FALSE)[1,]
      worldfile[aa+4,] <- data.frame(V1=0,V2="gap_fraction", stringsAsFactors=FALSE)[1,]
      worldfile[aa+5,] <- data.frame(V1=0,V2="rootzone_depth", stringsAsFactors=FALSE)[1,]
      worldfile[aa+6,] <- data.frame(V1=0,V2="snow_stored", stringsAsFactors=FALSE)[1,]
      worldfile[aa+7,] <- data.frame(V1=0,V2="rain_stored", stringsAsFactors=FALSE)[1,]
      worldfile[aa+8,] <- data.frame(V1=0,V2="cs.cpool", stringsAsFactors=FALSE)[1,]
      worldfile[aa+9,] <- data.frame(V1=0,V2="cs.leafc", stringsAsFactors=FALSE)[1,]
      worldfile[aa+10,] <- data.frame(V1=0,V2="cs.dead_leafc", stringsAsFactors=FALSE)[1,]
      worldfile[aa+11,] <- data.frame(V1=0,V2="cs.leafc_store", stringsAsFactors=FALSE)[1,]
      worldfile[aa+12,] <- data.frame(V1=0,V2="cs.leafc_transfer", stringsAsFactors=FALSE)[1,]
      worldfile[aa+13,] <- data.frame(V1=0,V2="cs.live_stemc", stringsAsFactors=FALSE)[1,]
      worldfile[aa+14,] <- data.frame(V1=0,V2="cs.livestemc_store", stringsAsFactors=FALSE)[1,]
      worldfile[aa+15,] <- data.frame(V1=0,V2="cs.livestemc_transfer", stringsAsFactors=FALSE)[1,]
      worldfile[aa+16,] <- data.frame(V1=0,V2="cs.dead_stemc", stringsAsFactors=FALSE)[1,]
      worldfile[aa+17,] <- data.frame(V1=0,V2="cs.deadstemc_store", stringsAsFactors=FALSE)[1,]
      worldfile[aa+18,] <- data.frame(V1=0,V2="cs.deadstemc_transfer", stringsAsFactors=FALSE)[1,]
      worldfile[aa+19,] <- data.frame(V1=0,V2="cs.live_crootc", stringsAsFactors=FALSE)[1,]
      worldfile[aa+20,] <- data.frame(V1=0,V2="cs.livecrootc_store", stringsAsFactors=FALSE)[1,]
      worldfile[aa+21,] <- data.frame(V1=0,V2="cs.livecrootc_transfer", stringsAsFactors=FALSE)[1,]
      worldfile[aa+22,] <- data.frame(V1=0,V2="cs.dead_crootc", stringsAsFactors=FALSE)[1,]
      worldfile[aa+23,] <- data.frame(V1=0,V2="cs.deadcrootc_store", stringsAsFactors=FALSE)[1,]
      worldfile[aa+24,] <- data.frame(V1=0,V2="cs.deadcrootc_transfer", stringsAsFactors=FALSE)[1,]
      worldfile[aa+25,] <- data.frame(V1=0,V2="cs.frootc", stringsAsFactors=FALSE)[1,]
      worldfile[aa+26,] <- data.frame(V1=0,V2="cs.frootc_store", stringsAsFactors=FALSE)[1,]
      worldfile[aa+27,] <- data.frame(V1=0,V2="cs.frootc_transfer", stringsAsFactors=FALSE)[1,]
      worldfile[aa+28,] <- data.frame(V1=0,V2="cs.cwdc", stringsAsFactors=FALSE)[1,]
      worldfile[aa+29,] <- data.frame(V1=0,V2="epv.prev_leafcalloc", stringsAsFactors=FALSE)[1,]
      worldfile[aa+30,] <- data.frame(V1=0,V2="ns.npool", stringsAsFactors=FALSE)[1,]
      worldfile[aa+31,] <- data.frame(V1=0,V2="ns.leafn", stringsAsFactors=FALSE)[1,]
      worldfile[aa+32,] <- data.frame(V1=0,V2="ns.dead_leafn", stringsAsFactors=FALSE)[1,]
      worldfile[aa+33,] <- data.frame(V1=0,V2="ns.leafn_store", stringsAsFactors=FALSE)[1,]
      worldfile[aa+34,] <- data.frame(V1=0,V2="ns.leafn_transfer", stringsAsFactors=FALSE)[1,]
      worldfile[aa+35,] <- data.frame(V1=0,V2="ns.live_stemn", stringsAsFactors=FALSE)[1,]
      worldfile[aa+36,] <- data.frame(V1=0,V2="ns.livestemn_store", stringsAsFactors=FALSE)[1,]
      worldfile[aa+37,] <- data.frame(V1=0,V2="ns.livestemn_transfer", stringsAsFactors=FALSE)[1,]
      worldfile[aa+38,] <- data.frame(V1=0,V2="ns.dead_stemn", stringsAsFactors=FALSE)[1,]
      worldfile[aa+39,] <- data.frame(V1=0,V2="ns.deadstemn_store", stringsAsFactors=FALSE)[1,]
      worldfile[aa+40,] <- data.frame(V1=0,V2="ns.deadstemn_transfer", stringsAsFactors=FALSE)[1,]
      worldfile[aa+41,] <- data.frame(V1=0,V2="ns.live_crootn", stringsAsFactors=FALSE)[1,]
      worldfile[aa+42,] <- data.frame(V1=0,V2="ns.livecrootn_store", stringsAsFactors=FALSE)[1,]
      worldfile[aa+43,] <- data.frame(V1=0,V2="ns.livecrootn_transfer", stringsAsFactors=FALSE)[1,]
      worldfile[aa+44,] <- data.frame(V1=0,V2="ns.dead_crootn", stringsAsFactors=FALSE)[1,]
      worldfile[aa+45,] <- data.frame(V1=0,V2="ns.deadcrootn_store", stringsAsFactors=FALSE)[1,]
      worldfile[aa+46,] <- data.frame(V1=0,V2="ns.deadcrootn_transfer", stringsAsFactors=FALSE)[1,]
      worldfile[aa+47,] <- data.frame(V1=0,V2="ns.frootn", stringsAsFactors=FALSE)[1,]
      worldfile[aa+48,] <- data.frame(V1=0,V2="ns.frootn_store", stringsAsFactors=FALSE)[1,]
      worldfile[aa+49,] <- data.frame(V1=0,V2="ns.frootn_transfer", stringsAsFactors=FALSE)[1,]
      worldfile[aa+50,] <- data.frame(V1=0,V2="ns.cwdn", stringsAsFactors=FALSE)[1,]
      worldfile[aa+51,] <- data.frame(V1=0,V2="ns.retransn", stringsAsFactors=FALSE)[1,]
      worldfile[aa+52,] <- data.frame(V1=0,V2="epv.wstress_days", stringsAsFactors=FALSE)[1,]
      worldfile[aa+53,] <- data.frame(V1=0,V2="epv.max_fparabs", stringsAsFactors=FALSE)[1,]  
      worldfile[aa+54,] <- data.frame(V1=0,V2="epv.min_vwc", stringsAsFactors=FALSE)[1,]
      worldfile[aa+55,] <- data.frame(V1=0,V2="n_basestations", stringsAsFactors=FALSE)[1,]
      
      change_flag <- 0
      understory_canopy_strata_ID <- understory_canopy_strata_ID + 1
      aa <- aa + 55
    }
    aa <- aa + 1
  }
  # Write new file
  worldfile$V1 <- format(worldfile$V1, scientific = FALSE)
  write.table(worldfile, file = world_name_out, row.names = FALSE, col.names = FALSE, quote=FALSE, sep="  ")
}


add_understory(world_name_in,
               world_name_out,
               canopy_default_ID,
               understory_canopy_default_ID,
               understory_canopy_strata_ID)

