# Use to add fire ID to a worldfile

BEGIN {h=0;}
{a = 0; }
($2 == "landuse_default_ID") {printf("%f	%s\n1.0 fire_default_ID\n",$1,$2); a=1;}
(a == 0) {printf("%s	%s\n",$1,$2);}
