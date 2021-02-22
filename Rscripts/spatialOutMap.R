#title: "spatialOutMap"
#author: "ErinJHanan"

# read in patch daily results, uses readin_rhessys_output function
setwd("...")
output = readin_rhessys_output("yourOutputPrefix", p=1)

# subset individual days ... do for each day you output
tmp = subset(output$pdg,output$pdg$year == 2001)
out.day.pdg = subset(tmp, tmp$month == 10) #name based on date for desired spatial output
# ...  

# Export for GRASS
setwd("...")
tmp2 = sprintf("%d:%d:%f", out.day.pdg$patchID, out.day.pdg$stratumID, out.day.pdg$proj_lai)
write.table(tmp2, file="lai.txt", quote=F, row.names=FALSE, col.names=F) #name based on date for desired spatial output


# To read it back into GRASS:
# Note: p.rip30.up270.cl = name of my patch map
# cd ../auxdata
# r.recode input=p.rip30.up270.cl output=lai.new.2014 rules=lai.new.2001.txt
