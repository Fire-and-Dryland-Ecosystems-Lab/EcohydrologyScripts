# Change the cover fraction for all patches in a worldfile

setwd("...")

worldIn <- "world.jc.h77.s100.v200.suFix"
worldOut <- "world.jc.h77.s100.v200.suFix.60"

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

  
  
