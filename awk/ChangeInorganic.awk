# zero out ammonium and nitrate pools after spinup, give a small amount of seed C

{a = 0;}
($2 == "soil_ns.sminn")  {printf("%f	%s\n",0.0001,$2); a=1;}
($2 == "soil_ns.nitrate")  {printf("%f	%s\n",0.0002,$2); a=1;}
($2 == "cs.cpool")  {printf("%f	%s\n",0.004,$2); a=1;}
($2 == "ns.npool")  {printf("%f	%s\n",0.002,$2); a=1;}

(a == 0) {printf("%s	%s\n",$1,$2);}
