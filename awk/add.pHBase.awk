# This will add a pH base station record to the worldfile (use when running with pH timeseries)

{a = 0;}
($2 == "Ksat_vertical") {printf("%f  %s\n",$1=0.5,$2); a=1;}
(a == 0) {printf("%s   %s\n",$1,$2);}
