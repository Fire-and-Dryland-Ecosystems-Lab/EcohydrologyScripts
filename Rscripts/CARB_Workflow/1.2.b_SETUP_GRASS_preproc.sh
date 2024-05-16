#!/bin/bash
# input one $dem_file
# input two threhold to define subbains
# input three and four basin outlet lat and lon

echo "import DEM raster"
# remove mask if needed
r.mask -r

set inputThreshold=500

echo "CHANGE DEM"

r.in.gdal -e input="preprocessing/spatial90m/dem.tif" output="dem30m"

g.region rast=dem30m -p

echo "fill depression of dem"
r.fill.dir input=dem30m output=dem30f dir=dir30m

echo "calculate the flow -s(higherthan6.4) means D8 single direction"
r.watershed --overwrite elevation=dem30f accumulation=acc30m drainage=drain30m 
r.topidx input=dem30f output=wetness_index30m --overwrite

echo "calculate the slope and aspect"
r.slope.aspect -e --overwrite elevation=dem30f slope=slope30m aspect=aspect30m

echo "calcualte the east and west horizon, distance is sampling distance step coef"
echo "check this website: https://github.com/RHESSys/RHESSys/wiki/East-and-West-Horizon-Maps"

r.horizon --overwrite -d elevation=dem30f direction=180 output="west"
r.horizon --overwrite -d elevation=dem30f direction=0 output="east"
r.mapcalc expression=e_horizon=sin(east_000)
r.mapcalc expression=w_horizon=sin(west_180)

echo "calculate subbasins and hillslopes, the theshold is important"
echo "v.info map=gauge_loc@PERMANENT"
echo "find the basin shape from usgs, input into grass, then overlap them with stream_1000, pay attention to elevation and stream order and null"
r.water.outlet --overwrite input=drain30m@PERMANENT output=basin coordinates=
r.mask -r
r.mask basin

echo "threshold is %inputThreshold%"
r.watershed --o elevation=dem30f threshold=%inputThreshold% basin="subbasins" stream="streams" half_basin="hillslope"

echo "calculate other maps and variables"

r.mapcalc "slope30m_fill = if(isnull(slope30m),0.143,slope30m)" --overwrite
r.mapcalc "aspect30m_fill = if(isnull(aspect30m),abs(drain30m)*45,aspect30m)" --overwrite


echo "import soils layers, originally from R-polaris workflow"

r.import input=C:\Users\burke\Documents\CARB\Pitman\preprocessing\spatial_source\POLARISOut\mean\clay\0_5\lat3738_lon-120-119.tif output=clay_0_5
r.import input=C:\Users\burke\Documents\CARB\Pitman\preprocessing\spatial_source\POLARISOut\mean\sand\0_5\lat3738_lon-120-119.tif output=sand_0_5

r.soils.texture sand=sand_0_5@PERMANENT clay=clay_0_5@PERMANENT scheme=C:\Users\burke\Documents\CARB\data\USDA.dat output=soil_texture

echo "import the soil texture"

r.proj input=soil_texture location=soil mapset=rhessys output=soil_texture method=nearest --v --o


echo "CHANGE BASE NAME HERE"

set base="C:/Users/burke/Documents/Carb/Ward/preprocessing/spatial"
set base="C:/Users/burke/Documents/Carb/Ward/preprocessing/spatial90m"
set base="C:/Users/burke/Documents/Carb/Ward/preprocessing/spatial180m"


r.out.gdal in=dem30f output="%base%/dem.tif" format=GTiff --o
r.out.gdal in=basin output="%base%/basin.tif" format=GTiff --o
r.out.gdal in=subbasins output="%base%/subbasins.tif" format=GTiff --o
r.out.gdal in=hillslope output="%base%/hillslope.tif" format=GTiff --o
r.out.gdal in=streams output="%base%/streams.tif" format=GTiff --o
r.out.gdal in=slope30m output="%base%/slope.tif" format=GTiff --o
r.out.gdal in=aspect30m output="%base%/aspect.tif" format=GTiff --o
r.out.gdal in=e_horizon output="%base%/e_horizon.tif" format=GTiff --o
r.out.gdal in=w_horizon output="%base%/w_horizon.tif" format=GTiff --o
r.out.gdal in=soil_texture output="%base%/soil_texture.tif" format=GTiff --o

echo "done"