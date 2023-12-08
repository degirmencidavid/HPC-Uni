/*==========================================================//
//  ANSI-C code (un optimised) for PX425 assignment 2 2021  //
//  Evolves a function u in 2D via finite differences       //
//  of the 2D Burgers Equation (simplified)                 //
//                                                          //
//  d u          du   du           d^2       d^2            //
//  ----- = -u ( -- + -- ) + nu ( ---- u  +  ---- u )       //
//  d t          dx   dy          dx^2       dy^2           //
//                                                          //
//  Based originally on code created by D. Quigley          //
//  Adapted by N. Hine                                      //
//==========================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "makePNG.h"     /* For visualisation */
#include "mt19937ar.h"   /* Random number generator */
#include <string.h>

/* Function prototypes for memory management routines */
void allocate2d(double ***a,int num_rows,int num_cols);
void free2d(double ***a,int num_rows);

int main () {

  /* Function u on new and current grid */
  double **u_new, **u;

  /* Approximate Laplacian */
  //double Lapl, grad;

  /* Number of grid points */
  int Nx = 256;
  int Ny = 256;

  /* Loop counters */
  int ix,iy,istep;

  /* Filename to which the grid is drawn */
  int  isnap=0;
  char filename[25];

  /*--------------*/
  /* Initial time */
  /*--------------*/
  clock_t t1 = clock();

  /*------------------------------------*/
  /* Initialise random number generator */
  /*------------------------------------*/
  unsigned long seed = 120549784972;
  init_genrand(seed);

  /*--------------------------*/
  /* Set grid spacing         */
  /*--------------------------*/
  double dx = 1.0;
  double dy = 1.0;

  /*-----------------------------*/
  /* Set timestep and K          */
  /*-----------------------------*/
  double dt = 0.001;
  double nu = 5.0;

  /*--------------------------*/
  /* Number of steps to run   */
  /*--------------------------*/
  int nstep = 10000;

  /*--------------------------------------*/
  /* Allocate memory for a bunch of stuff */
  /*--------------------------------------*/
  allocate2d(&u,Nx,Ny);
  allocate2d(&u_new,Nx,Ny);

  /*--------------------------------*/
  /* Initialise with random numbers */
  /*--------------------------------*/
  for(ix=0;ix<Nx;ix++) {
    for(iy=0;iy<Ny;iy++) {
      u[ix][iy] = 2.0*genrand() - 1.0;
    }
  }
      
  /*------------------------------------*/
  /* Write an image of the initial grid */
  /*------------------------------------*/
  int stepincr = 10; 
  sprintf(filename,"snapshot%08d.png",isnap); 
  writePNG(filename,u,Nx,Ny); 
  isnap++;

  /* setup time */
  clock_t t2 = clock();
  printf("Setup time                    : %15.6f seconds\n",(double)(t2-t1)/(double)CLOCKS_PER_SEC);
  t1 = t2;

  /*===============================*/
  /* BEGIN SECTION TO BE OPTIMISED */
  /*===============================*/
 
  /*------------------------------------------*/
  /* Loop over the number of output timesteps */
  /*------------------------------------------*/
  
//to avoid repeated calculation and declaration
double dtnu = dt*nu;
//this int i is actually not necessary, but I've used it and it's not worth it to go back and rename every instance of it for no improvement
int i;

  
  for (istep=1;istep<nstep;istep++) {
	  

	for (i=1;i<Nx-1;i++) {
		
		///NB///
		//u_new[ix][iy] = u[ix][iy]*(1-dt*grad) + dtnu*Lapl; but I have eliminated the need of Lap1 and grad being calculated on separate lines
		
		///points away from boundary///
		for (iy=1;iy<Ny-1;iy++) {
			u_new[i][iy] = u[i][iy]*(1-dt*(-u[i-1][iy] - u[i][iy-1] + u[i][iy+1] + u[i+1][iy])*0.5) + dtnu*(u[i-1][iy] + u[i][iy-1] - 4.0*u[i][iy] + u[i][iy+1] + u[i+1][iy]);
		}
		///END///
		
		///ix==0 boundary///
		//away from iy==Ny-1 boundary, iterating through iy=1 to iy=Ny-2
		u_new[0][i] = u[0][i]*(1-dt*(-u[0][i-1] + u[0][i+1] + u[1][i] - u[Nx-1][i])*0.5) + dtnu*(u[0][i-1] -4.0*u[0][i] + u[0][i+1] + u[1][i] + u[Nx-1][i]);
		///END///	
		
		
		///iy==0 boundary///
		//away from ix==Nx-1 boundary, iterating through ix=1 to ix=Nx-2
		u_new[i][0] = u[i][0]*(1-dt*(-u[i-1][0] + u[i][1] - u[i][Ny-1] + u[i+1][0])*0.5) + dtnu*(u[i-1][0] -4.0*u[i][0] + u[i][1] + u[i][Ny-1] + u[i+1][0]);
		///END///		
		
		///ix==Nx-1 boundary///
		//iy!=Ny-1,0 , iterating through iy=1 to iy=Ny-2
		u_new[Nx-1][i] = u[Nx-1][i]*(1-dt*(u[0][i] - u[Nx-2][i] - u[Nx-1][i-1] + u[Nx-1][i+1])*0.5) + dtnu*(u[0][i] + u[Nx-2][i] + u[Nx-1][i-1] -4.0*u[Nx-1][i] + u[Nx-1][i+1]);
		///END///
		
		///iy==Ny-1 boundary///
		//ix!=Nx-1,0, iterating through ix=1 to ix=Nx-2
		u_new[i][Ny-1] = u[i][Ny-1]*(1-dt*(-u[i-1][Ny-1] + u[i][0] - u[i][Ny-2] + u[i+1][Ny-1])*0.5) + dtnu*(u[i-1][Ny-1] + u[i][0] + u[i][Ny-2] -4.0*u[i][Ny-1] + u[i+1][Ny-1]);
		///END///
	}

	
	
	///non loop///
	
		
	///bottom left (ix,iy==0)///
	u_new[0][0] = u[0][0]*(1-dt*(u[0][1] - u[0][Ny-1] + u[1][0] - u[Nx-1][0])*0.5) + dtnu*(-4.0*u[0][0] + u[0][1] + u[0][Ny-1] + u[1][0] + u[Nx-1][0]);
	///END///
	
	
	///ix==0 boundary///
	//iy==Ny-1 boundary
	//iy = Ny-1;
	u_new[0][Ny-1] = u[0][Ny-1]*(1-dt*(u[0][0] - u[0][Ny-2] + u[1][Ny-1] - u[Nx-1][Ny-1])*0.5) + dtnu*(u[0][0] + u[0][Ny-2] -4.0*u[0][Ny-1] + u[1][Ny-1] + u[Nx-1][Ny-1]);
	///END///


	///iy==0 boundary///
	//ix==Nx-1 boundary
	u_new[Nx-1][0] = u[Nx-1][0]*(1-dt*(u[0][0] - u[Nx-2][0] + u[Nx-1][1] - u[Nx-1][Ny-1])*0.5) + dtnu*(u[0][0] + u[Nx-2][0] -4.0*u[Nx-1][0] + u[Nx-1][1] + u[Nx-1][Ny-1]);
	///END///
	
	
	
	///ix==Nx-1,iy==Ny-1 (top right)
	u_new[Nx-1][Ny-1] = u[Nx-1][Ny-1]*(1-dt*(u[0][Ny-1] - u[Nx-2][Ny-1] + u[Nx-1][0] - u[Nx-1][Ny-2])*0.5) + dtnu*(u[0][Ny-1] + u[Nx-2][Ny-1] + u[Nx-1][0] + u[Nx-1][Ny-2] -4.0*u[Nx-1][Ny-1]);
	///END///
	
	///memcpy u_new into u///
	//iterate through slices to fill u with u_new, memcpy(u, u_new, Nx*Ny*sizeof(double)) would be nice but not possible with the way memory has been allocated
	for(iy=0;iy<Ny;iy++) {	
		memcpy(u[iy], u_new[iy], Nx*sizeof(double*));
    }
	


    /*-----------------------------*/
    /* Snapshots of grid to file   */
    /*-----------------------------*/
    if ( istep==isnap)  { 
        sprintf(filename,"snapshot%08d.png",isnap); 
        writePNG(filename,u,Nx,Ny); 
        isnap *= stepincr; 
    }  

  }

  /*=============================*/
  /* END SECTION TO BE OPTIMISED */
  /*=============================*/
  
  /* calculation time */
  t2 = clock();
  printf("Time taken for %8d steps : %15.6f seconds\n",nstep,(double)(t2-t1)/(double)CLOCKS_PER_SEC);
    
  /*----------------------------------*/
  /* Write an image of the final grid */
  /*----------------------------------*/
  sprintf(filename,"snapshot%08d.png",istep); 
  writePNG(filename,u,Nx,Ny); 

  /*--------------------------------------------*/
  /* Write final time-evolved solution to file. */
  /*--------------------------------------------*/
  FILE *fp = fopen("final_grid.dat","w");
  if (fp==NULL) printf("Error opening final_grid.dat for output\n");

  for(ix=0;ix<Nx-1;ix++) {
    for(iy=0;iy<Ny-1;iy++) {
      /* x and y at the current grid points */
      double x = dx*(double)ix;
      double y = dy*(double)iy;
      fprintf(fp,"%8.4f %8.4f %8.4e\n",x,y,u[ix][iy]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);

  /* Release memory */
  free2d(&u,Nx);
  free2d(&u_new,Nx);
  
  return 0;
  
}
  


/*===========================================*/
/* Auxilliary routines for memory management */ 
/*===========================================*/
void allocate2d(double ***a,int Nx,int Ny) {

  double **b_loc; 

  b_loc = (double **)calloc(Nx,sizeof(double *));
  if (b_loc==NULL) printf("malloc error in allocate2d\n"); 

  int iy;
  for (iy=0;iy<Nx;iy++) {
    
    b_loc[iy] = (double *)calloc(Ny,sizeof(double));
    if (b_loc[iy]==NULL) printf("malloc error for row %d of %d in allocate2d\n",iy,Nx);

  }

  *a = b_loc;

}

void free2d(double ***a,int Nx) {

  int iy;

  double **b_loc = *a;

  /* Release memory */
  for (iy=0;iy<Nx;iy++) { 
    free(b_loc[iy]);
  }
  free(b_loc);
  *a = b_loc;

}