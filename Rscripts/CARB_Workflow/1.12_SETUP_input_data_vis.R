# create a pdf with all input maps - maybe some extras for soil and landcover

library(RHESSysPreprocessing)
library(rhutils)
library(terra)


map_dir = "preprocessing/whitebox/"
dest = "preprocessing"
out_name = "map_plots"
map_exts = c(".tif|.tiff")
pdfwidth = 7 
pdfheight = 7

template = "preprocessing/template/walla_std.template"
# template_read(template = "")

input_map_plotpdf(map_dir = map_dir, template = template)


# check template
template = "preprocessing/template/walla_std.template"


check_template(template)
