#title: "ChangeCoverFraction"
#authors: "ErinJHanan and Ryan Bart"

worldIn <- "world.jc.h77.s100.v200.suFix"
worldOut <- "world.jc.h77.s100.v200.suFix.60"

# changes cover fraction for all veg types

changeFraction <- function(worldIn,
                           worldOut){
  
  # Read in worldfile
  worldfile <- read.table(worldIn, header = FALSE, stringsAsFactors = FALSE)
  aa <- 1
  for (aa in seq_len(nrow(worldfile))){
    if (aa%%1000 == 0 ){print(paste("Worldfile line", aa, "of", nrow(worldfile)))} # Counter
    
    if (worldfile[aa,2] == "cover_fraction") { 
      worldfile[aa,1] <- 0.6      # Changes cover fraction
    }
    aa <- aa + 1
  }
  # Write new file
  worldfile$V1 <- format(worldfile$V1, scientific = FALSE)
  write.table(worldfile, file = worldOut, row.names = FALSE, col.names = FALSE, quote=FALSE, sep="  ")
}

changeFraction(worldIn,
               worldOut)  

###################################################################################################
# ----------
#If want to change only for specific veg types

worldIn <- "world.jc.h77.s100.v200.suFix"
worldOut <- "world.jc.h77.s100.v200.suFix.60"

changeFractionPartial <- function(worldIn,
                           worldOut){
  
  # Read in worldfile
  worldfile <- read.table(worldIn, header = FALSE, stringsAsFactors = FALSE)
  aa <- 1
  for (aa in seq_len(nrow(worldfile))){
    if (aa%%1000 == 0 ){print(paste("Worldfile line", aa, "of", nrow(worldfile)))} # Counter
    
    if (worldfile[aa,2] == "veg_parm_ID" && worldfile[aa,1] == 41) { # put in the ID for the canopy type you want to change
      worldfile[aa+1,1] <- 0.6      # Changes cover fraction
    }
    aa <- aa + 1
  }
  # Write new file
  worldfile$V1 <- format(worldfile$V1, scientific = FALSE)
  write.table(worldfile, file = worldOut, row.names = FALSE, col.names = FALSE, quote=FALSE, sep="  ")
}

changeFractionPartial(worldIn,
               worldOut)  








  
  
