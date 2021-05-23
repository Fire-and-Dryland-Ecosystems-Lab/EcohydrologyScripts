# Changes vertical drainage, e.g. to generate hydrophobicity. This may be applied to worldfiles for multiple years to enable hydrophobicity to return to prefire conditions over time.

{a = 0;}
($2 == "Ksat_vertical") {printf("%f  %s\n",$1=0.10,$2); a=1;}
(a == 0) {printf("%s   %s\n",$1,$2);}
