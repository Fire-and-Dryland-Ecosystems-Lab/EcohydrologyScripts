# Creates awk scripts for initializing vegetation in a worldfile (using previous worldfile). This is for multiple vegetation types.
# Prior to running this, must extract appropriate lines from prev. worldfile for EACH vegetation type 
# e.g., :75,118 w conifer.txt ... etc.
# Then run awk -f format_init2.awk < conifer.txt > initialize.conifer.awk

{printf("\n($2 == \""$2"\") && (h==1) {printf(\"%%f %%s\\n\",$1="$1",$2); a=1;}"); 
}
