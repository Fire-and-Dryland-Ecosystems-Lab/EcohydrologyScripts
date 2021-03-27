#use format_loc.awk to generate another awk that will search for line number and replace the x or y with the new value. 
#To run this:
#awk -f format_loc <jc.xloc> replace_xloc.awk
#awk -f format_loc <jc.yloc> replace_yloc.awk
#additional instructions in Bear notes
{printf("\n($1 == \""$1"\") {printf(\"%%d %%f %%s\\n\", $1, $2="$2",$3); a=1;}"); 
}
