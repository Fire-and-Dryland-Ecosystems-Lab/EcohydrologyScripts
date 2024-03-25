#define MAXFILENAME 200
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <limits.h>
#define MAXID 50000

int searchid(int *array, int data, int dim) {
	int i = 0;
	int b = 0;
	while(!b && i < dim) {
		if(array[i] == data) {
			b = 1;
			break;
		}
		i++;
	}
	if(b == 1) return i;
	else return -1;
}

int main(int argc, char *argv[]) {
	float **data; //[cellid][var]: 0:x 1:y 2:lai 3:z 4:lat 5:lont
	int *cellidnum;
	FILE *fin[3]; //0:lai 1:dem 2:coordinate
	FILE *fcellid, *fcor, *fout;
	int rows, cols, cellmin, cellmax;
	int i, j, k, itmp, id, cell;
	int **cellid, cells;
	char buff[1000];
	float ftmp, ftmp2;
	float lat, lont, x, y;
	char gridl[200];
	setvbuf(stdout, NULL, _IONBF, 0);

	if(argc < 2) {
		fprintf(stderr, "usage:%s <cellid_asc> <lai_asc> <dem_asc> <coordinate_info> <basin_name_or_path_to_met> <out_base> <year_start>\n", argv[0]);
		return 1;
	}

	if((fcellid = fopen(argv[1], "r")) == NULL) {
		fprintf(stderr, "Can't open cellid file:%s\n", argv[1]);
		return 1;
	}

	if((fout = fopen(argv[6],"w")) == NULL) {
		fprintf(stderr, "Cannot create base file:%s\n", argv[6]);
		return 1;
	}

	/****read grid cellid****/
	buff[0] = '\0';
	fgets(gridl, 100, fcellid);
	sscanf(gridl, "%s %i", buff, &cols);
	fgets(gridl, 100, fcellid);
	sscanf(gridl, "%s %i", buff, &rows);
	fgets(gridl, 100, fcellid);
	fgets(gridl, 100, fcellid);
	fgets(gridl, 100, fcellid);
	fgets(gridl, 100, fcellid);

	cellid = (int **)malloc(rows * sizeof(int *));
	if(cellid == NULL) {
		fprintf(stderr,  "out of memory\n");
		return 1;
	}
	for(i = 0; i < rows; i++) {
		cellid[i] = (int *)malloc(cols * sizeof(int));
		if(cellid[i] == NULL) {
			fprintf(stderr,  "out of memory\n");
		return 1;
		}
	}
	printf("malloc sucess\n");
	//cellmin = INT_MAX;
	//cellmax = INT_MIN;
	for (i = 0; i < rows; i++) {
		for(j = 0; j < cols; j++) {
			cellid[i][j] = -9999;
			fscanf(fcellid, "%i", &cellid[i][j]);
			//if(cellid[i][j]>0 && cellid[i][j]<cellmin) cellmin = cellid[i][j];
			//else if(cellid[i][j]>0 && cellid[i][j]>cellmax) cellmax = cellid[i][j];
		}
		fgets(buff, 100, fcellid);
	}
	fclose(fcellid);
	//cells = cellmax-cellmin+1;
	cells = 0;
	for (i=0; i < rows; i++) {
		for(j=0; j < cols; j++) {
			if(cellid[i][j] != -9999) {
				cells++;
			}
		}
    }

	printf("cells:%i \n", cells);

	cellidnum = (int *) malloc(cells*sizeof(int));
	data = (float **) malloc(cells * sizeof(float *));
	if(data == NULL || cellidnum == NULL) {
		fprintf(stderr,  "out of memory\n");
        return 1;
    }
	
	for(i = 0; i < cells; i++) {
		data[i] = (float *) malloc(6 * sizeof(float));
		if(data[i] == NULL) {
			fprintf(stderr, "out of memory\n");
		return 1;
		}
	}
	
	//initialize data
	for(i = 0; i < cells; i++) {
		for(j = 0; j < 6; j++) data[i][j] = -9999.0;
	}
	cell = 0;
	
	for (i = 0; i < rows; i++) {
        for(j = 0; j < cols; j++) {
			if(cellid[i][j] != -9999) {
				cellidnum[cell] = cellid[i][j];
				cell++;
			}
        }
    }

	printf("reading head of lai, dem, and cor\n");
	//read lai, dem, cor's head
	for(i = 0; i < 3; i++) {
		if((fin[i] = fopen(argv[i + 2], "r")) == NULL) {
			fprintf(stderr, "cannot open file:%s\n", argv[i+2]);	
			return 1;
		}
		if(i < 2) {
			for(j = 0; j < 6; j++) fgets(buff, 200, fin[i]); //ascii grid head
		}
		//else fgets(buff, 200, fin[i]); //cor has no head
	}
	
	//read lai, dem
	for(k = 0; k < 2; k++) {
		cell = 0;
		for (i = 0; i < rows; i++) {
			for(j = 0; j < cols; j++) {
				fscanf(fin[k], "%f", &ftmp);
				itmp = cellid[i][j];
				if(itmp != -9999){
					data[cell][k+2] = ftmp;
					cell++;
				}
			}
			fgets(buff, 100, fin[k]);
		}
	fclose(fin[k]);
	}
	
	//read cor
	while(!feof(fin[2])) {
		id = -9999;
		fscanf(fin[2], "%f %f %f %f %f %f", &ftmp, &ftmp2, &lat, &lont, &x, &y);
		fgets(buff, 200, fin[2]);
		id = (int)(ftmp);
		cell = searchid(cellidnum, id, cells);
		printf("id:%i cell:%i\n", id, cell);
		if(cell != -1) {
			data[cell][0] = x;
			data[cell][1] = y;
			data[cell][4] = lat;
			data[cell][5] = lont;
		}
	}
	fclose(fin[2]);
	
	//outfile
	fprintf(fout, "%d grid_cells\n", cells);
    fprintf(fout, "0.041667 location_searching_distance\n");
    fprintf(fout, "%s year_start_index\n", argv[7]);
    fprintf(fout, "0.0 day_offset\n");
    fprintf(fout, "1.0 leap_year_include\n");
    fprintf(fout, "0.001 precip_multiplier\n");
	fprintf(fout, "0.01 rhum_multiplier\n");
    fprintf(fout, "K temperature_unit\n");
    fprintf(fout, "lon netcdf_var_x\n");
    fprintf(fout, "lat netcdf_var_y\n");
    fprintf(fout, "%sCL_tasmax_1979_2020.nc netcdf_tmax_filename\n", argv[5]);
    fprintf(fout, "air_temperature netcdf_var_tmax\n");
    fprintf(fout, "%sCL_tasmin_1979_2020.nc netcdf_tmin_filename\n", argv[5]);
    fprintf(fout, "air_temperature netcdf_var_tmin\n");
    fprintf(fout, "%sCL_pr_1979_2020.nc netcdf_rain_filename\n", argv[5]);
    fprintf(fout, "precipitation_flux netcdf_var_rain\n");
    fprintf(fout, "%sCL_huss_1979_2020.nc netcdf_huss_filename\n", argv[5]);
    fprintf(fout, "specific_humidity netcdf_var_huss\n");	
    fprintf(fout, "%sCL_rmax_1979_2020.nc netcdf_rmax_filename\n", argv[5]);
    fprintf(fout, "relative_humidity netcdf_var_rmax\n");	
    fprintf(fout, "%sCL_rmin_1979_2020.nc netcdf_rmin_filename\n", argv[5]);
    fprintf(fout, "relative_humidity netcdf_var_rmin\n");		
    fprintf(fout, "%sCL_rsds_1979_2020.nc netcdf_rsds_filename\n", argv[5]);
    fprintf(fout, "surface_downwelling_shortwave_flux_in_air netcdf_var_rsds\n");	
    fprintf(fout, "%sCL_was_1979_2020.nc netcdf_was_filename\n", argv[5]);
    fprintf(fout, "wind_speed netcdf_var_was\n");	
    
	for(i = 0; i < cells; i++) {
		if(cellidnum[i] != -9999.0) {
			fprintf(fout, "%i base_station_id\n", cellidnum[i]);
			fprintf(fout, "%f lon\n", data[i][0]);
			fprintf(fout, "%f lat\n", data[i][1]);
			fprintf(fout, "%f z_coordinate\n", data[i][3]);
			//fprintf(fout, "lon netcdf_var_x\n");
			//fprintf(fout, "lat netcdf_var_y\n");
			//fprintf(fout, "%f lon_coordinate\n", data[i][5]);
			//fprintf(fout, "%f lat_coordinate\n", data[i][4]);
			fprintf(fout, "%f effective_lai\n", data[i][2]);
			fprintf(fout, "10.0 screen_height\n");
			//for 1/24-degree met data
			
			//fprintf(fout, "0.045 location_searching_distance\n");
			
			
			//fprintf(fout, "1.0 day_offset\n");
			
			//fprintf(fout, "%sYakima_UofI_7914_tmax_conv.nc netcdf_tmax_filename\n", argv[5]);
			
			
			//fprintf(fout, "%sYakima_UofI_7914_tmin_conv.nc netcdf_tmin_filename\n", argv[5]);
			
			
			//fprintf(fout, "%sYakima_UofI_7914_prec_conv.nc netcdf_rain_filename\n", argv[5]);
			
			
		}
	}
	fclose(fout);
	
	for(i = 0; i < cells; i++) free(data[i]);
	free(cellidnum);
	free(data);
}		
	

