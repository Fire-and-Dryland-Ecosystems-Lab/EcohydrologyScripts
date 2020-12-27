# Changes aspect for all patches in the worldfile, useful when running single patch

{a = 0;}
($2 == "aspect") {printf("%f  %s\n",$1=0.00,$2); a=1;}
(a == 0) {printf("%s   %s\n",$1,$2);}
