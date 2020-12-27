#include <stdio.h>

void main(void){

  FILE *fpi, *fpo;
  char *fi = "../HRSC-Blend/Mars_HRSC_MOLA_BlendDEM_Global_200mp_v2_gdal.csv";
  char *fo = "output.log";
  double x, y, z, xp, yp, zp;
  unsigned int m=0; 
  unsigned int n=0; 
  unsigned int i=0;
  double x1=1e10, x2=-1e10; 

  // open files
  fpi = fopen( fi, "r" );
  fpo = fopen( fo, "w" );
  fprintf(fpo,"#   IX    X[IX]     NY   min(Y)   max(Y)\n");

  // read until the end of file is detected
  while(!feof(fpi)){

     // read one line
     fscanf(fpi, "%lf %lf %lf", &x, &y, &z); 

     if (i){ // skip the first line
             // check changes in y value
	     if(y != yp || feof(fpi)){
                     // output y value (n and yp), number of x data, min and max value of x
		     printf("%6d %8.3f %6d %8.3f %8.3f\n", n, yp, m, x1, x2);
		     fprintf(fpo, "%6d %8.3f %6d %8.3f %8.3f\n", n, yp, m, x1, x2);
		     m = 0;
		     n ++;
		     x1 = 1e10;
		     x2 = -1e10;
	     }
     }
     i = 1;
     // store the prevous values
     xp = x;
     yp = y;
     zp = z;
     if (x > x2) x2 = x;
     if (x < x1) x1 = x;
     m ++;


  }				  
  // close files
  fclose( fpi );
  fclose( fpo );
}

