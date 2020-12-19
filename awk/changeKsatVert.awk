# Changes vertical drainage, e.g. to generate hydrophobicity

{a = 0;}
($2 == "Ksat_vertical") {printf("%f  %s\n",$1=0.10,$2); a=1;}
(a == 0) {printf("%s   %s\n",$1,$2);}
