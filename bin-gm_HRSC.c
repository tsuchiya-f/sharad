// **********************************************************************
// *** bin-gm: ?????
// ***   2019/12/23  modified by YK --- Addition of 'Latitude selection'
// **********************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tiffio.h"
      
#define CV (299792458.0)
#define DR (0.0375e-6*CV/2)
#define RM (3396000)

#define NX (106694) //11520*4
#define NY ( 53347) // 5632*4

#define DL (3)
#define DEG (M_PI/180.0)
#define M (10)

// **** Only "LAT = LAT_MIN - LAT-MAX" can pass. ****
#define LAT_MAX (45.0)
#define LAT_MIN (15.0)
#define LON_MAX (-20.0)
#define LON_MIN (-60.0)

// Vector normaliztion
double norm(double *a,double *b) {
  double r;
  int i;
  r=sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
  for (i=0;i<3;i++) b[i]=a[i]/r;
  return r;
}

// Inner product
double dot(double *a,double *b) {
  return(a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
}

// Outer product w/o normalization
void cross(double *a,double *b,double *c) {
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[2]*b[0]-a[0]*b[2];
  c[2]=a[0]*b[1]-a[1]*b[0];
}

// Outer product with normalization
void cross_n(double *a,double *b,double *c) {
  double cc;
  int i;
  cross(a,b,c);
  cc=sqrt(dot(c,c));
  for (i=0;i<3;i++) c[i]/=cc;
}

//-----------------------------------------------------

int main(int argc,char **argv) {
  static float dat[3604];
  static char dat2[24];
  FILE *fp;
  short *datr;
  static double datr0[M+1][M+1];
  static double datr1[M+1][M+1];
  //short *datt;
  int n;
  int ix;
  int iy;
  int ix0;
  int ix1;
  int iy0;
  int iy1;
  int jx;
  int jx0;
  int jx1;
  double rr_ar;
  static double e_sc[3];
  //static double r_sc[3];
  double rr_sc;
  static double e_su[3];
  static double r_su[3];
  double rr_su;
  static double r_su1[3];
  static double r_su2[3];
  static double r_su3[3];
  static double r_su4[3];
  double vv;
  static double n1[3];
  static double n2[3];
  static double n3[3];
  static double n4[3];
  static double nv[3];
  static char s[64];
    
  double lo;
  double la;
  double lo0;
  double la0;
  int i,j,k;
  int ns;
  int ne;
  FILE *fp1;
  int kx;
  int ky;
  double p1;
  double p2;
  double p3;
  double p4;
  double p5;
  double p6;
  static double v1[3];
  double rr;
  int ii;
  
  int nx, ny;
  double lon1, lon2, lat1, lat2;
  double lon_nbin, lat_nbin;
  float  f_val;

  sscanf(argv[1],"%d",&ns);
  sscanf(argv[2],"%d",&ne);

  // read HRSC Blend data ------------------------------------------
  fp=fopen("output_HRSC.bin","rb");
  fread(&nx, sizeof(int), 1, fp);      // number of bin of longitude
  fread(&ny, sizeof(int), 1, fp);      // number of bin of latitude
  fread(&f_val, sizeof(float), 1, fp); lon2 = (double)f_val; // minimum longitude
  fread(&f_val, sizeof(float), 1, fp); lon1 = (double)f_val; // maximum longitude
  fread(&f_val, sizeof(float), 1, fp); lat1 = (double)f_val; // maximum latitude
  fread(&f_val, sizeof(float), 1, fp); lat2 = (double)f_val; // minimum latitude
  if (lon1 < 0.0) lon1 += 360.0;
  if (lon2 < 0.0) lon2 += 360.0;
  fprintf(stderr, " longitude: %d %lf %lf\n",nx,lon1,lon2);
  fprintf(stderr, " latitude : %d %lf %lf\n",ny,lat1,lat2);
  
  datr = (short *)malloc(sizeof(short)*nx*ny);
  fread(datr, sizeof(short)*nx*ny, 1, fp);
  
  //printf("%d\n",datr[0]);
  fclose(fp);

  //----------------------------------------------------------------

  // read radergram data from stdin --------------------------------
  // read radergram file name (*.img)
  fread(dat2,sizeof(char),24,stdin);

  // read radergram data
  for (n=0;n<ns;n++) {
    fread(dat,sizeof(float),3604,stdin);if (feof(stdin)) exit(0);
  }
  for (n=ns;n<ne;n++) {
    fread(dat,sizeof(float),3604,stdin);if (feof(stdin)) break;
    rr_sc=dat[3603]*1e3;   // R of S/C
    rr_ar=dat[3602]*1e3;   // R of Areoid
    lo0=dat[3601];         // Longitude of S/C
    la0=dat[3600];         // Latitude of S/C

    // **** Only "LAT = LAT_MIN to LAT-MAXi & LON = LON_MIN to LON_MAX" can pass.  ****
    if (la0 > lat1 || la0 < lat2 || lo0 > lon1 || lo0 < lon2){
      fprintf(stderr,"# %s %lf %lf %lf %lf  ==> SKIP !!!\n",
		dat2,dat[3600],dat[3601],dat[3602],dat[3603]);
      continue;
    }
  
    fprintf(stderr,"# %s %lf %lf %lf %lf\n",
	    dat2,dat[3600],dat[3601],dat[3602],dat[3603]);

    lon_nbin = (float)nx/(lon1-lon2);   // number of grid (longitude, HRSC data)
    ix0=floor((lo0 - lon2 - DL)*lon_nbin); if (ix0<1   ) ix0=1;
    ix1=floor((lo0 - lon2 + DL)*lon_nbin); if (ix1>nx-2) ix1=nx-1;
    lat_nbin = (float)ny/(lat1-lat2);   // nuber of grid (latitude, HRSC data)
    iy0=floor((lat1 - la0 - DL)*lat_nbin); if (iy0<1   ) iy0=1;
    iy1=floor((lat1 - la0 + DL)*lat_nbin); if (iy1>ny-2) iy1=ny-2; //NY NY

    e_sc[0]=cos(la0*DEG)*cos(lo0*DEG);
    e_sc[1]=cos(la0*DEG)*sin(lo0*DEG);
    e_sc[2]=sin(la0*DEG);
    //for (i=0;i<3;i++) r_sc[i]=rr_sc*e_sc[i];

    sprintf(s,"./data/%04d.gm",n);
    fp=fopen(s,"w");
    fwrite(e_sc  ,sizeof(double),3,fp);
    fwrite(&rr_sc,sizeof(double),1,fp);
    fwrite(&rr_ar,sizeof(double),1,fp);

    //fprintf(stderr,"%d %d %d %d\n",ix0,ix1,iy0,iy1);
    for (ix=ix0;ix<ix1;ix++) {
      if ((ix%100)==0) fprintf(stderr,"(ix0,ix,ix1) = (%d %d %d)\n",ix0,ix,ix1);
      for (iy=iy0;iy<iy1;iy++) {
	//fprintf(stderr,"%d %d %d %d %d %d\n",ix0,ix,ix1,iy0,iy,iy1);
	jx=ix;
	while (jx<0)   jx+=ny; //NY
	while (jx>=nx) jx-=nx; //NX
	//jx0=ix-1;
	//while (jx0<0)   jx0+=NX;
	//while (jx0>=NX) jx0-=NX;
	jx1=ix+1;
	while (jx1<0)   jx1+=nx; //NX
	while (jx1>=nx) jx1-=nx; //NX

	p1=datr[jx*ny+iy  ];p2=datr[jx1*ny+iy  ]; //NY NY
	p3=datr[jx*ny+iy+1];p4=datr[jx1*ny+iy+1]; //NY NY
	//printf("%lf %lf %lf %lf\n", p1,p2,p3,p4);
	for (kx=0;kx<=M;kx++) {
	  p5=p1+(p2-p1)*kx/(double)M;
	  p6=p3+(p4-p3)*kx/(double)M;
	  for (ky=0;ky<=M;ky++) {
	    datr0[kx][ky]=p5+(p6-p5)*ky/(double)M;
	  }
	}
	for (ky=0;ky<=M;ky++) {
	  p5=p1+(p3-p1)*ky/(double)M;
	  p6=p2+(p4-p2)*ky/(double)M;
	  for (kx=0;kx<=M;kx++) {
	    datr1[kx][ky]=p5+(p6-p5)*kx/(double)M;
	  }
	}

	for (kx=0;kx<=M;kx++) {
	  for (ky=0;ky<=M;ky++) {
	    datr0[kx][ky]=0.5*datr0[kx][ky]+0.5*datr1[kx][ky];
	  }
	}
	for (kx=0;kx<M;kx++) {
	  for (ky=0;ky<M;ky++) {
	    datr1[kx][ky]
	      =0.25*datr0[kx  ][ky]+0.25*datr0[kx  ][ky+1]
	      +0.25*datr0[kx+1][ky]+0.25*datr0[kx+1][ky+1];
	  }
	}
	
	/*
	if ((ix==ix0+10)&&(iy==iy0+10)) {
	  fp1=fopen("ll","w");
	  fprintf(fp1,"# %lf %lf %lf %lf\n",p1,p2,p3,p4);
	  for (kx=0;kx<=100;kx++) {
	    for (ky=0;ky<=100;ky++) {
	      fprintf(fp1,"%d %d %lf\n",kx,ky,datr0[kx][ky]);
	    }
	    fprintf(fp1,"\n");
	  }
	  fclose(fp1);
	}
	*/
	
	for (kx=0;kx<M;kx++) {
	  for (ky=0;ky<M;ky++) {
	    rr_su=RM+datr1[kx][ky];
	    lo = lon2 + (ix+(kx+0.5)/(double)M) / lon_nbin;
	    la = lat1 - (iy+(ky+0.5)/(double)M) / lat_nbin;
	    e_su[0]=cos(la*DEG)*cos(lo*DEG);
	    e_su[1]=cos(la*DEG)*sin(lo*DEG);
	    e_su[2]=sin(la*DEG);
	    for (i=0;i<3;i++) r_su[i]=rr_su*e_su[i];
	    //rr=(r_su[0]-r_sc[0])*(r_su[0]-r_sc[0])
	    //+(r_su[1]-r_sc[1])*(r_su[1]-r_sc[1])
	    //+(r_su[2]-r_sc[2])*(r_su[2]-r_sc[2]);
	    //rr=sqrt(rr);
	    //ii=1800-(rr_sc-rr-rr_ar)/DR;  // r0=rr+i*D+ra
	    //if (ii<0)    ii=0;
	    //if (ii>3599) ii=3599;

	    for (i=0;i<3;i++) v1[i]=r_su[i]-rr_sc*e_sc[i];
	    rr=sqrt(dot(v1,v1));
	    ii=1800+((rr_sc-rr_ar)-rr)/DR;
	    if ((0<=ii)&&(ii<3600)) {
	    
	      vv=RM+datr0[kx][ky];
	      lo = lon2 + (ix+kx/(double)M) / lon_nbin;
	      la = lat1 - (iy+ky/(double)M) / lat_nbin;
	      r_su1[0]=vv*cos(la*DEG)*cos(lo*DEG);
	      r_su1[1]=vv*cos(la*DEG)*sin(lo*DEG);
	      r_su1[2]=vv*sin(la*DEG);
	      
	      vv=RM+datr0[kx][ky+1];
	      lo = lon2 + (ix+kx/(double)M)     / lon_nbin;
	      la = lat1 - (iy+(ky+1)/(double)M) / lat_nbin;
	      r_su2[0]=vv*cos(la*DEG)*cos(lo*DEG);
	      r_su2[1]=vv*cos(la*DEG)*sin(lo*DEG);
	      r_su2[2]=vv*sin(la*DEG);
	      
	      vv=RM+datr0[kx+1][ky+1];
	      lo = lon2 + (ix+(kx+1)/(double)M) / lon_nbin;
	      la = lat1 - (iy+(ky+1)/(double)M) / lat_nbin;
	      r_su3[0]=vv*cos(la*DEG)*cos(lo*DEG);
	      r_su3[1]=vv*cos(la*DEG)*sin(lo*DEG);
	      r_su3[2]=vv*sin(la*DEG);
	      
	      vv=RM+datr0[kx+1][ky];
	      lo = lon2 + (ix+(kx+1)/(double)M) / lon_nbin;
	      la = lat1 - (iy+ky/(double)M)     / lat_nbin;
	      r_su4[0]=vv*cos(la*DEG)*cos(lo*DEG);
	      r_su4[1]=vv*cos(la*DEG)*sin(lo*DEG);
	      r_su4[2]=vv*sin(la*DEG);
	      
	      for (i=0;i<3;i++) r_su1[i]-=r_su[i];
	      for (i=0;i<3;i++) r_su2[i]-=r_su[i];
	      for (i=0;i<3;i++) r_su3[i]-=r_su[i];
	      for (i=0;i<3;i++) r_su4[i]-=r_su[i];
	      cross_n(r_su1,r_su2,n1);
	      cross_n(r_su2,r_su3,n2);
	      cross_n(r_su3,r_su4,n3);
	      cross_n(r_su4,r_su1,n4);
	      for (i=0;i<3;i++) n1[i]=n1[i]+n2[i]+n3[i]+n4[i];
	      norm(n1,nv);

	      fwrite(e_su  ,sizeof(double),3,fp);
	      fwrite(&rr_su,sizeof(double),1,fp);
	      fwrite(nv    ,sizeof(double),3,fp);
	    }
	  }
	}
	//printf("%lf ",lo);
	//printf("%lf ",la);
	//printf("%d ",datt[jx*NY+iy]);
	//printf("%d ",datr[jx*NY+iy]);
	//printf("%lf ",sqrt(rr));
	//printf("%d ",ii);
	//printf("%lf ",10*log10(dat[ii]));
	//printf("\n");
      }
      //printf("\n");
    }
    fclose(fp);
  }
  exit(0);
}
