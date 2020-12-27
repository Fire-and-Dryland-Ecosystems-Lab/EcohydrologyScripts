# Changes the veg def file in worldfile header. Useful if need to change several at a time

{a = 0;}
($1 == "../defs/veg_chaparral.def7.mod") {printf("%s  %s\n",$1="../defs/veg_chaparral_seed",$2); a=1;}
(a == 0) {printf("%s   %s\n",$1,$2);}
