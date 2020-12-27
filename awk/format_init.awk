# Creates awk scripts for initializing vegetation in a worldfile (using previous worldfile). This is for a single vegetation type.
# Prior to running this, must extract appropriate lines from prev. worldfile :75,118 w spinup_vals.txt
# Then run awk -f format_init.awk < spinup_vals.txt > initialize.awk

{printf("\n($2 == \""$2"\") {printf(\"%%f %%s\\n\",$1="$1",$2); a=1;}"); 
}
