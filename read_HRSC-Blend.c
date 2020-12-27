#include<stdio.h>
#include<stdlib.h>
#include<string.h>

// grid size of source data (Mars_HRSC_MOLA_BlendDEM_Global_200mp_v2_gdal.csv)
#define NX 106694
#define NY 53347
// latitude and longitude ranges to read
#define LAT_MIN 0.0
#define LAT_MAX 1.0
#define LON_MIN -170.0
#define LON_MAX -169.0

// function to read BlendDEM csv data
// --------------------------------------------------------------------------------
void read_HRSC_Blend(short *dat, float *lat, float *lon, int nx, int ny){

  FILE *fp;
  int i, j, k;
  int ix, iy;
  float x, y, z;
  fp=fopen("../HRSC-Blend/Mars_HRSC_MOLA_BlendDEM_Global_200mp_v2_gdal.csv","r");

  // grid number of LON_MIN
  ix = (int)((float)NX/360.0*(LON_MIN+180.0));
  // grid number of LAT_MAX
  iy = (int)((float)NY/180.0*(90.0-LAT_MAX));

  for (j=0; j<iy+ny; j++){
    for (i=0; i<NX; i++){
      // read data from the csv file
      fscanf(fp, "%f %f %f", &x, &y, &z);
      // store data if ix<=i<ix+nx & iy<=j<iy+ny
      if (i >= ix && i < ix+nx && j >= iy){
        k = (i-ix) + (j-iy)*nx;
        dat[k] = (short)z;
        lon[k] = x;
        lat[k] = y;
      }
    }
  }
  fclose(fp);
}
// --------------------------------------------------------------------------------

void main(void)
{
  short *datr;
  float *lat,*lon;
  int nx, ny;
  int i,j,k;
  FILE *fp;

  // number of grid for longitude (from LON_MIN to LON_MAX)
  nx = (int)((float)NX/360.0*(LON_MAX-LON_MIN));
  // number of grid for latitude (from LAT_MAX to LAT_MIN)
  ny = (int)((float)NY/180.0*(LAT_MAX-LAT_MIN));
  
  // read data
  datr = (short *)malloc(sizeof(short)*nx*ny);
  lat = (float *)malloc(sizeof(float)*nx*ny);
  lon = (float *)malloc(sizeof(float)*nx*ny);
  read_HRSC_Blend(datr,lat,lon,nx,ny);

  // output data
  fp = fopen("output_HRSC.dat", "w");
  for (j=0;j<ny;j++){
    for (i=0;i<nx;i++){
      k=i+j*nx;
      fprintf(fp,"%d %d %f %f %d\n", i,j,lon[k],lat[k],datr[k]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  

}

