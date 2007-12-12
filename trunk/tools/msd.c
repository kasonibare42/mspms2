/* 
   Mean-squared displacement calculation, straightforward algorithm.
   Reads in an arbitrary number of configuration snapshots 
   in XYZ format, and computes msd(t).

   Cameron F. Abrams

   Written for the course CHE 800-002, Molecular Simulation
   Spring 0304

   compile using "gcc -o msd.x msd.c -lm"

   Drexel University, Department of Chemical Engineering
   Philadelphia
   (c) 2004
*/
/* 
   revised by Yang Wang
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* Prints usage information */
void usage ( void ) {
  fprintf(stdout,"msd usage:\n");
  fprintf(stdout,"msd [options]\n\n");
  fprintf(stdout,"Options:\n");
  fprintf(stdout,"\t -fn [string]\t\tFile name\n");
  fprintf(stdout,"\t -nmole \t\tnumber of molecules\n");
  fprintf(stdout,"\t -natom \t\tnumber of atoms in one molecule, default is 1\n");
  fprintf(stdout,"\t -for total,span,interval\t\tloop control\n");
  fprintf(stdout,"\t -mddt \t\tMD time step, default is 1.0\n");
  fprintf(stdout,"\t -h           \t\tPrint this info.\n");
}

/* Reads in a configuration snapshot in XYZ format (such
   as that created by mdlj.c) */
int xyz_in (FILE * fp, double * rx, double * ry, double * rz, 
	    int * N) {
  int has_vel = 0;

  fread(rx,sizeof(double),*N,fp);
  fread(ry,sizeof(double),*N,fp);
  fread(rz,sizeof(double),*N,fp);

  return has_vel;
}

void com(double *rx, double *ry, double *rz, 
	int N, int nmole, int natom, double *atm_wgt, double weight,
	double *rmx, double *rmy, double *rmz, 
	double *cx, double *cy, double *cz)
{
    int i, j, k;
    *cx = 0.0; 
    *cy = 0.0;
    *cz = 0.0;
    for (i=0;i<nmole;i++)
    {
	rmx[i] = 0.0;
	rmy[i] = 0.0;
	rmz[i] = 0.0;
	for (j=0;j<natom;j++)
	{
	    rmx[i] += rx[i*natom+j]*atm_wgt[j];
	    rmy[i] += ry[i*natom+j]*atm_wgt[j];
	    rmz[i] += rz[i*natom+j]*atm_wgt[j];
	}

	rmx[i] = rmx[i]/weight;
	rmy[i] = rmy[i]/weight;
	rmz[i] = rmz[i]/weight;
    }
    for (i=0;i<nmole;i++)
    {
	*(cx) += rmx[i];
	*(cy) += rmy[i];
	*(cz) += rmz[i];
    }
}

int main ( int argc, char * argv[] ) {

  double * rx0, * ry0, * rz0;
  double * rx, * ry, * rz;
  double *rmx0, *rmy0, *rmz0;
  double *rmx, *rmy, *rmz;
  double cx0, cy0, cz0;
  double cx, cy, cz;
  int N;
  double * sdx, * sdy, * sdz;
  double * sdcx, * sdcy, * sdcz;
  int t, dt;
  double md_time_step = 1.0;
  double L=0.0, V=0.0;
  double dr=0.1, r, vb, nid;
  double rc2 = 10.0, rho=0.85;
  int start=0,stop=0,step=1,ngr=0;
  int total, span, interval;
  int last_origin;
  int number_of_origins;
  int nmole;
  int natom = 1;
  double atm_wgt[100];
  double weight;
  char useless[100];
  int cnt;
  int i, ii;

  FILE * fp;

  char *fn;

  /* Here we parse the command line arguments;  If
   you add an option, document it in the usage() function! */
  for (i=1;i<argc;i++) 
  {
    if (!strcmp(argv[i],"-for"))
      sscanf(argv[++i],"%i,%i,%i",&total,&span,&interval);
    else if (!strcmp(argv[i],"-fn")) fn=argv[++i];
    else if (!strcmp(argv[i],"-nmole")) nmole=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-natom")) natom = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-mddt")) md_time_step=atof(argv[++i]);
    else if (!strcmp(argv[i],"-h")) 
    {
      usage(); exit(0);
    }
    else 
    {
      fprintf(stderr,"Error: Command-line argument '%s' not recognized.\n", argv[i]);
      exit(-1);
    }
  }

  N = nmole*natom;
  if (natom>1)
  {
      for (i=0;i<natom;i++)
      {
	  fprintf(stderr,"Atom %i weight = ",i);
	  atm_wgt[i] = atof(gets(useless));
	  fprintf(stderr,"echo %f\n",atm_wgt[i]);
      }
  }
  else
  {
      atm_wgt[0] = 1.0;
  }
  weight = 0.0;
  for (i=0;i<natom;i++)
      weight += atm_wgt[i];
  fprintf(stderr,"molecule weight = %f\n",weight);

  last_origin = total-span;
  number_of_origins = (total-span)/interval + 1;
  fprintf(stderr,"# Expecting %i frames...\n",total);
  fprintf(stderr,"# Expecting %i time origins...\n",number_of_origins);

  /* Allocate the trajectory arrays */
  rx0 = (double*)malloc(N*sizeof(double));
  ry0 = (double*)malloc(N*sizeof(double));
  rz0 = (double*)malloc(N*sizeof(double));
  rx = (double*)malloc(N*sizeof(double));
  ry = (double*)malloc(N*sizeof(double));
  rz = (double*)malloc(N*sizeof(double));

  /* Allocate the molecule arrays */
  rmx0 = (double*)malloc(nmole*sizeof(double));
  rmy0 = (double*)malloc(nmole*sizeof(double));
  rmz0 = (double*)malloc(nmole*sizeof(double));
  rmx = (double*)malloc(nmole*sizeof(double));
  rmy = (double*)malloc(nmole*sizeof(double));
  rmz = (double*)malloc(nmole*sizeof(double));

  /* Allocate and initialize the squared-displacement teams */
  sdx=(double*)calloc(span+1,sizeof(double));
  sdy=(double*)calloc(span+1,sizeof(double));
  sdz=(double*)calloc(span+1,sizeof(double));
  sdcx=(double*)calloc(span+1,sizeof(double));
  sdcy=(double*)calloc(span+1,sizeof(double));
  sdcz=(double*)calloc(span+1,sizeof(double));


  /* Compute the mean-squared displacement */
  fprintf(stderr,"# computing...\n");fflush(stderr);
  for (t=0;t<=last_origin;t+=interval) 
  { 
      fp = fopen(fn,"rb");
      /* ignore t frames */
      for (ii=0;ii<t;ii++)
	  xyz_in(fp,rx0,ry0,rz0,&N);
      /* read the current frame */
      fprintf(stderr,"# Reading in trajectory %i ...",t);fflush(stderr);
      xyz_in(fp,rx0,ry0,rz0,&N);
      com(rx0,ry0,rz0,N,nmole,natom,atm_wgt,weight,rmx0,rmy0,rmz0,&cx0,&cy0,&cz0);
      /* loop through possible intervals */
      cnt = 1;
      for (dt=1;dt<=span;dt++) 
      {
	  xyz_in(fp,rx,ry,rz,&N);
	  com(rx,ry,rz,N,nmole,natom,atm_wgt,weight,rmx,rmy,rmz,&cx,&cy,&cz);
	  for (i=0;i<nmole;i++) 
	  {
	      sdx[dt] += (rmx[i] - rmx0[i])*(rmx[i] - rmx0[i]);
	      sdy[dt] += (rmy[i] - rmy0[i])*(rmy[i] - rmy0[i]);
	      sdz[dt] += (rmz[i] - rmz0[i])*(rmz[i] - rmz0[i]);
	  }
	  sdcx[dt] +=(cx-cx0)*(cx-cx0);
	  sdcy[dt] +=(cy-cy0)*(cy-cy0);
	  sdcz[dt] +=(cz-cz0)*(cz-cz0);
	  cnt++;
      }
      fclose(fp);
      fprintf(stderr,"done: %i frames.\n",cnt);fflush(stderr);
  }
  for (t=0;t<=span;t++) 
  {
      sdx[t] /= nmole;
      sdy[t] /= nmole;
      sdz[t] /= nmole;
      sdx[t] /= number_of_origins;
      sdy[t] /= number_of_origins;
      sdz[t] /= number_of_origins;
      sdcx[t] /= nmole;
      sdcy[t] /= nmole;
      sdcz[t] /= nmole;
      sdcx[t] /= number_of_origins;
      sdcy[t] /= number_of_origins;
      sdcz[t] /= number_of_origins;
      fprintf(stdout,"%.1lf %.8lf %.8lf %.8lf %.8lf %.8lf %.8lf\n",
	      t*step*md_time_step,sdx[t],sdy[t],sdz[t],
	      sdcx[t],sdcy[t],sdcz[t]);
  }
}
