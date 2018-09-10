# From N. Tague
# PASTE THE COMMANDS BELOW IN R TO CREATE THE WATER BALANCE FUNCTION

watbal = function(rb) {
  rb$watbal.tmp=with(rb,precip-streamflow-trans-evap)
  rb$sd=with(rb,sat_def-rz_storage-unsat_stor)
  tmp=diff(rb$sd)
  tmp=c(0,tmp)
  rb$sddiff=tmp
  tmp=diff(rb$snowpack)
  tmp=c(0,tmp)
  rb$snodiff=tmp
  tmp=diff(rb$detention_store)
  tmp=c(0,tmp)
  rb$detdiff=tmp
  tmp=diff(rb$litter_store)
  tmp=c(0,tmp)
  rb$litdiff=tmp
  tmp=diff(rb$canopy_store)
  tmp=c(0,tmp)
  rb$candiff=tmp
  tmp=diff(rb$gw.storage)
  tmp=c(0,tmp)
  rb$gwdiff=tmp
  rb$watbal=with(rb,watbal.tmp+sddiff-snodiff-detdiff-litdiff-candiff-gwdiff)
  rb
}

#To use this function in an R data set:
#If your basin daily output is called (for example) "bd" type "bd=watbal(bd)".
#This will add a number of fields to your basin daily output file, the last of which will be called "watbal". This is your water balance error. You should see a minor numerical error here even if your water balances (on the order of 10^-6). If your watbal values are negative then water inputs are less than water outputs, and vice versa.

