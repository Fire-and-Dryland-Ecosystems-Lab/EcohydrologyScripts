# NETCDF Header Formating

Author: Ashley Cale
Created: Apr 12, 2021 6:59 PM
Tags: #Climate Change #NETCDF #RHESSys

This workflow is for formatting the headers of netcdf files (.nc) from [MACA Statistical Downscaling Method](https://climate.northwestknowledge.net/MACA/data_portal.php) for reading into rhessys. Netcdf’s downloaded from MACA are downscaled  GCM’s and are in gridmet resolution. You shouldn’t have to change your basestation file in order to run RHESSys with these files.

I recommend making a series of .nc files as you go _test1.nc, _test2.nc, etc. You can delete extra's later, this makes it easy to keep track of your changes and see if anything didn't work. If you're courageous you can overwrite your files as you go (leave out "outfile.nc.)  What are you trying to prove though? 

This workflow uses the [NCO package](http://nco.sourceforge.net/). You'll also need the netcdf libraries in order to use ncdump command. 

To see the header of your netcdf:

```bash
ncdump -h infile.nc
```

To see lon, lat or other variables:

```bash
ncdump -vlat [infile.nc](http://infile.nc) 
```

*If you ask me for questions I will always help! Will I also be googling in another tab? Yes...but we'll figure it out together :)*

## Get Files

Download .nc files from [MACA Statistical Downscaling Method](https://climate.northwestknowledge.net/MACA/data_portal.php).
At the time of writing this, MACA doesn’t have all variables to support expanded climate variables in RHESSys.

## Convert longitude coordinates

MACA files come raw with longitude coordinates with range (0,360). RHESSys needs (180,-180)

```bash
ncap2 -O -s 'where(lon>180) lon=lon-360' [infile.nc](http://infile.nc/) [outfile.nc](http://outfile.nc/)
```

## Time Dimension

MACA outputs dimension associated with time as "time" we need it to be "day" 

```bash
ncrename -d time,day infile.nc outfile.nc
```

## Time Variable

MACA outputs time dimension variable as "time" we need it to be "day" 

```bash
ncrename -v time,day [infile.nc](http://infile.nc) outfile.nc
```

## Coordinates for Climate Variable

If you're looking at a precipitation netcdf, your climate variable will be "precipitation" or "precipitation amount" (something like that). In your header there will be a precipitation:coordinates = "time lat lon" line we need to change it to "day lat lon" 

```bash
ncatted -a coordinates,precipitation,o,c,"day lat lon" [in](http://inline.nc)file.nc [outfile.nc](http://outfile.nc) 
```

You'll edit that for temperature

## Change Time Dimension to Double Variable

You'll notice that your time variable is a float variable. We need this to be a Double variable to match the lat and lon. 

```bash
ncap2 -s 'day=double(day)' [in](http://inline.nc)file.nc outfile.nc
```