/*=========================================================/
/  Comms routines for PX425 assignment 4. Contains all     /
/  routines which interact with MPI libraries. Many of     /
/  these are currently incomplete and will work only in    /
/  serial. You will need to correct this.                  /
/                                                          /
/  Original code created by N. Hine  - November 2021       /
/  (based on previous code by D. Quigley)                  /
/=========================================================*/
#include "mpi.h"
#include "grid.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "comms.h"


int p;               /* Number of processor       */
int my_rank;         /* Rank of current processor */
int my_rank_cart;    // cartesian grid rank ""

MPI_Comm cart_comm;  /* Cartesian communicator    */

/* Coordinates of current rank in the processor grid */
int my_rank_coords[2];

/* Ranks of neighbours to the current processor (left, right, down, up) */
int my_rank_neighbours[4];

/* Time and initialisation and shutdown */
double t1,t2;

void comms_initialise(int argc, char **argv) {
  /*==================================================================/
  / Function to initialise MPI, get the communicator size p and the   /
  / rank my_rank of the current process within that communicator.     /
  /-------------------------------------------------------------------/
  / N. Hine (based on code by D. Quigley) - Univ. Warwick             /
  /==================================================================*/
  int proot;                   /* square root of p */

  //initialise mpi
  MPI_Init(&argc, &argv);
  //store rank of current processor in my_rank
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  //store sizie of communicator in p
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  //printf("%d: sss (p=%d)\n", my_rank, p);


  /* Start the timer. Set t1 using MPI_Wtime() which returns a double */
  //as stated above
  t1 = MPI_Wtime();


  /* Check that we have a square number of processors */
  proot = (int)sqrt((double)p+0.5);
  if (proot*proot!=p) {
    if (my_rank==0) {
      printf("Number of processors must be an exact square!\n");
      exit(EXIT_FAILURE); 
    }
  }

  return;

}

void comms_processor_map() {
  /*====================================================================/
  / Function to map our p processors into a 2D Cartesian grid of      /
  / dimension proot by proot where proot = sqrt(p).                   /
  /                                                                   /
  / Should populate the arrays my_rank_cooords, which contains the    /
  / location of the current MPI task within the processor grid, and   /
  / my_rank_neighbours, which contains (in the order left, right,     /
  / down and up) ranks of neighbouring MPI tasks on the grid with     /
  / which the current task will need to communicate.                  /
  /-------------------------------------------------------------------/
  / N. Hine (based on code by D. Quigley) - University of Warwick     /
  /==================================================================*/
  
  /* Information for setting up a Cartesian communicator */
  int ndims = 2;
  int reorder = 0;
  int pbc [2] = {1,1};
  int dims[2];
    
  /* Local variables */
  int proot;            /* square root of p */

  /* Square root of number of processors */
  proot = (int)sqrt((double)p+0.5);

  /* Dimensions of Cartesian communicator */
  dims[x] = proot; dims[y] = proot;

  //////

  //create cartesian communicator
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, pbc, reorder, &cart_comm);

  //set communicator grid rank
  MPI_Comm_rank(cart_comm, &my_rank_cart);

  //current rank coordinates in frid
  MPI_Cart_coords(cart_comm, my_rank_cart, 2, my_rank_coords);

  //store neighbours in halo ring
  MPI_Cart_shift(cart_comm, x, 1, &my_rank_neighbours[left], &my_rank_neighbours[right]);
  MPI_Cart_shift(cart_comm, y, 1, &my_rank_neighbours[down], &my_rank_neighbours[up]);

  //print
  /*
  printf("my_rank: %d; my_rank_coords: (%d,%d)\n", my_rank, my_rank_coords[0], my_rank_coords[1]);
  printf("my_rank: %d; my_rank_neighbours: left: %d, right=%d, down=%d, up=%d\n", my_rank, my_rank_neighbours[left], my_rank_neighbours[right],  my_rank_neighbours[down], my_rank_neighbours[up]);
  */

  return;

  
} 

void comms_get_global_mag(double local_mag,double *global_mag) {
  /*==================================================================/
  / Function to compute the glocal magnetisation of the grid by       /
  / averaging over all values of local_mag, and storing the result    /
  / in global_mag.                                                    /
  /-------------------------------------------------------------------/
  / N. Hine (based on code by D. Quigley) - University of Warwick     /
  /==================================================================*/

  /* This is only correct on one processor. You will need        */
  /* to use a collective communication routine to correct this.  */
  *global_mag = local_mag;


  /* Insert collective communication operation here */
  MPI_Reduce(&local_mag, global_mag, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


  *global_mag = *global_mag/(double)p;

} 

void comms_halo_swaps() {
  /*==================================================================/
  / Function to send boundary spins on each side of the local grid    /
  / to neighbour processors, and to receive from those processors the /
  / halo information needed to perform computations involving spins   /
  / on the boundary processors grid.                                  /
  /-------------------------------------------------------------------/
  / N. Hine (based on code by D. Quigley) - University of Warwick     /
  /==================================================================*/

  /* Send and receive buffers */
  int *sendbuf, *recvbuf;
    
  /* MPI Status */
  MPI_Status status;

  int ix,iy;      /* Loop counters */

  /* If running on 1 processor copy boundary elements into opposite halo */
  if (p==1) {

    for (iy = 0 ; iy < grid_domain_size ; iy++) {
      grid_halo[right][iy] = grid_spin[iy][0]; 
      grid_halo[left][iy]  = grid_spin[iy][grid_domain_size-1];
    }

    for (ix = 0 ; ix < grid_domain_size ; ix++) {
      grid_halo[up][ix]   = grid_spin[0][ix];
      grid_halo[down][ix] = grid_spin[grid_domain_size-1][ix];
    }

    return; /* Do not do any comms */
  }

  /* Allocate buffers */
  sendbuf = (int *)malloc(grid_domain_size*sizeof(int));
  if (sendbuf==NULL) { 
    printf("Error allocating sendbuf in comms_halo_swaps\n");
    exit(EXIT_FAILURE);
  }
  recvbuf = (int *)malloc(grid_domain_size*sizeof(int));
  if (recvbuf==NULL) { 
    printf("Error allocating recvbuf in comms_halo_swaps\n");
    exit(EXIT_FAILURE);
  }

  /* Send left hand boundary elements of grid_spin to my_rank_neighbours[left] 
     and receive from my_rank_neighbours[right] into the appropriate part
     of grid_halo. Remember to use the appropriate communicator. */


  //put left boundary of grid_spin into sendbuf
  for(iy=0;iy<grid_domain_size;iy++) {
    sendbuf[iy] = grid_spin[iy][0];
  }

  /* Insert MPI calls here to implement this swap. Use sendbuf and recvbuf */
  /* REDUNDANT, MPI_Sendrecv does all of this
  //initialise mpi status,request for nonblocking Isend,Irecv
  MPI_Status send_stat_0;
  MPI_Status recv_stat_0;

  MPI_Request send_req_0;
  MPI_Request recv_req_0;
  //actual sending and receiving (nonblocking)
  MPI_Isend(sendbuf, grid_domain_size, MPI_INT, my_rank_neighbours[left], my_rank_neighbours[left]+600, MPI_COMM_WORLD, &send_req_0);
  MPI_Irecv(recvbuf, grid_domain_size, MPI_INT, my_rank_neighbours[right], my_rank+600, MPI_COMM_WORLD, &recv_req_0);

  //wait for these to be completed
  MPI_Wait(&send_req_0, &send_stat_0);
  MPI_Wait(&recv_req_0, &recv_stat_0);
  */

  //MPI_Sendrecv isn't like MPI_Send followed by MPI_Recv
  //it's like MPI_Isend followed by MPI_Irecv followed by the two MPI_Wait s
  //so nonblocking sending and receiving is achieved, which is essential to
  //avoid deadlocking

  //does the same as commented out section above, status can be reused
  MPI_Sendrecv(sendbuf,grid_domain_size, MPI_INT, my_rank_neighbours[left], my_rank_neighbours[left]+600,recvbuf, grid_domain_size, MPI_INT, my_rank_neighbours[right], my_rank+600, MPI_COMM_WORLD, &status);

  //put elements of recvbuf into grid_halo as right hand boundary
  for(iy=0;iy<grid_domain_size;iy++) {
    grid_halo[right][iy] = recvbuf[iy];
  }


  /* Send right hand boundary elements of grid_spin to my_rank_neighbours[right]
     and receive from my_rank_neighbours[left] into the appropriate part 
     of grid_halo. Remember to use the appropriate communicator. */

  //put right boundary of grid_spin into sendbuf
  for(iy=0;iy<grid_domain_size;iy++) {
    sendbuf[iy] = grid_spin[iy][grid_domain_size-1];
  }

  /* Insert MPI calls here to implement this swap. Use sendbuf and recvbuf */

  /* REDUNDANT, MPI_Sendrecv does all of this
  //MPI_Status can be reused, it is a struct containing only ints
  MPI_Status send_stat_1;
  MPI_Status recv_stat_1;
  //MPI_requests can be reused 
  //looking at mpi.h, MPI_Request is just an integer, so after you
  //finish using it (with MPI_Wait), it is safe to reuse, not that there's any benefit
  MPI_Request send_req_1;
  MPI_Request recv_req_1;
  
  
  //send recv
  MPI_Isend(sendbuf, grid_domain_size, MPI_INT, my_rank_neighbours[right], my_rank_neighbours[right]+700, MPI_COMM_WORLD, &send_req_0);
  MPI_Irecv(recvbuf, grid_domain_size, MPI_INT, my_rank_neighbours[left], my_rank+700, MPI_COMM_WORLD, &recv_req_0);

  //wait for these to be completed
  MPI_Wait(&send_req_0, &send_stat_1);
  MPI_Wait(&recv_req_0, &recv_stat_1);
  */

  //does the same as commented out section above, status can be reused
  MPI_Sendrecv(sendbuf, grid_domain_size, MPI_INT, my_rank_neighbours[right], my_rank_neighbours[right]+700, recvbuf, grid_domain_size, MPI_INT, my_rank_neighbours[left], my_rank+700, MPI_COMM_WORLD, &status);

  //put elements of recvbuf into left boundary of grid_halo
  for(iy=0;iy<grid_domain_size;iy++) {
    grid_halo[left][iy] = recvbuf[iy];
  }


  /* Send bottom boundary elements of grid_spin to my_rank_neighbours[down]
     and receive from my_rank_neighbours[up] into the appropriate part
     of grid halo. Remember to use the appropriate communicator.  */

  

  //put lower boundary of grid_spin into sendbuf
  for(ix=0;ix<grid_domain_size;ix++) {
    sendbuf[ix] = grid_spin[0][ix];
  }
  
  /* Insert MPI calls here to implement this swap. Use sendbuf and recvbuf */

  /*
  MPI_Status send_stat_2;
  MPI_Status recv_stat_2;
  MPI_Request send_req_2;
  MPI_Request recv_req_2;

  MPI_Isend(sendbuf, grid_domain_size, MPI_INT, my_rank_neighbours[down], my_rank_neighbours[down]+800, MPI_COMM_WORLD, &send_req_0);
  MPI_Irecv(recvbuf, grid_domain_size, MPI_INT, my_rank_neighbours[up], my_rank+800, MPI_COMM_WORLD, &recv_req_0);

  MPI_Wait(&send_req_0, &send_stat_2);
  MPI_Wait(&recv_req_0, &recv_stat_2);
  */

  //does the same as commented out section above, status can be reused
  MPI_Sendrecv(sendbuf, grid_domain_size, MPI_INT, my_rank_neighbours[down], my_rank_neighbours[down]+800, recvbuf, grid_domain_size, MPI_INT, my_rank_neighbours[up], my_rank+800, MPI_COMM_WORLD, &status);
  
  //put elements of recvbuf into up boundary of grid_halo
  for(ix=0;ix<grid_domain_size;ix++) {
    grid_halo[up][ix] = recvbuf[ix];
  }

  /* Send top boundary elements of grid_spin to my_rank_neighbours[up]
     and receive from my_rank_neighbours[down] into the appropriate part
     of grid halo. Remember to use the appropriate communicator. */

  //put upper boundary of grid_spin into sendbuf
  for(ix=0;ix<grid_domain_size;ix++) {
    sendbuf[ix] = grid_spin[grid_domain_size-1][ix];
  }

  /* Insert MPI call or calls here to implement this swap. Use sendbuf and recvbuf */

  /*
  MPI_Status send_stat_3;
  MPI_Status recv_stat_3;
  MPI_Request send_req_3;
  MPI_Request recv_req_3;

  MPI_Isend(sendbuf, grid_domain_size, MPI_INT, my_rank_neighbours[up], my_rank_neighbours[up]+900, MPI_COMM_WORLD, &send_req_0);
  MPI_Irecv(recvbuf, grid_domain_size, MPI_INT, my_rank_neighbours[down], my_rank+900, MPI_COMM_WORLD, &recv_req_0);

  MPI_Wait(&send_req_0, &send_stat_3);
  MPI_Wait(&recv_req_0, &recv_stat_3);
  */

  //does the same as commented out section above, status can be reused
  MPI_Sendrecv(sendbuf, grid_domain_size, MPI_INT, my_rank_neighbours[up], my_rank_neighbours[up]+900, recvbuf, grid_domain_size, MPI_INT, my_rank_neighbours[down], my_rank+900, MPI_COMM_WORLD, &status);

  //put elements of recvbuf into down boundary of grid_halo
  for(ix=0;ix<grid_domain_size;ix++) {
    grid_halo[down][ix] = recvbuf[ix];
  }

  /* Release memory */
  free(sendbuf);
  free(recvbuf);

  return;

}


void comms_get_global_grid() {
  /*==================================================================/
  / Function to collect all contributions to the global grid onto     /
  / rank zero for visualisation.                                      /
  /-------------------------------------------------------------------/
  / N. Hine (based on code by D. Quigley) - University of Warwick     /
  /==================================================================*/

  /* comms buffer */
  int *combuff;

  /* MPI Status */
  MPI_Status status;

  /* Information on the remote domain */
  int remote_domain_start[2] = {0,0};

  /* Loop counters and error flag */
  int ix,iy,ixg,iyg,ip;

  /* Just point at local grid if running on one processor */
  if (p==1) {
    global_grid_spin=grid_spin;
    return;
  }

  if (my_rank==0) {

    /* Rank 0 first fills out its part of the global grid */
    for (iy = 0 ; iy < grid_domain_size ; iy++) {
      for (ix = 0 ; ix < grid_domain_size ; ix++) {
	
	/* Global indices */
	ixg = ix+grid_domain_start[x];
	iyg = iy+grid_domain_start[y];
	
	global_grid_spin[iyg][ixg] = grid_spin[iy][ix];
      
      }
    }

  }


  /* Remove the following line when you have inserted appropriate MPI calls below */
  //return;

  /* Allocate buffer */
  combuff = (int *)malloc(grid_domain_size*sizeof(int));
  if (combuff==NULL) { 
    printf("Error allocating combuff in comms_get_global_grid\n");
    exit(EXIT_FAILURE);
  }

  if (my_rank==0) {
    
    /* Now loops over all other ranks receiving their data */
    for (ip = 1; ip < p ; ip++) {

      /* First receive remote_domain_start from rank ip */
      /* Insert an appropriate MPI call here */
      MPI_Recv(&remote_domain_start, 2, MPI_INT, ip, 400, MPI_COMM_WORLD, &status);


      /* Loop over rows within a domain */
      for (iy = 0 ; iy < grid_domain_size ; iy++) {
	
        /* Receive this row from rank ip */ 
        /* Insert appropriate MPI call here */
        MPI_Recv(combuff, grid_domain_size, MPI_INT, ip, 500, MPI_COMM_WORLD, &status);

        for (ix = 0 ; ix < grid_domain_size ; ix++) {

          /* Global indices */
          ixg = ix+remote_domain_start[x];
          iyg = iy+remote_domain_start[y];
        
          /* Store in global_grid_spin */
          global_grid_spin[iyg][ixg] = combuff[ix];


        } /* elements in row */

      } /* rows */

    } /* processors */

  } else {

    /* All other processors must send the data rank 0 needs */

    /* Send grid_domain_start to rank 0 */
    /* Insert appropriate MPI call here */
    MPI_Send(&grid_domain_start, 2, MPI_INT, 0, 400, MPI_COMM_WORLD);

    /* Loop over rows in the domain, sending them to rank 0*/
    for (iy = 0 ; iy < grid_domain_size ; iy++) {

      /* Insert appropriate MPI call here */
      MPI_Send(&grid_spin[iy][0], grid_domain_size, MPI_INT, 0, 500, MPI_COMM_WORLD);

    }
    
  }
  
  /* Free memory */
  free(combuff);

  return;
  
}
 

void comms_finalise() {
  /*==================================================================/
  / Function to finalise MPI functionality and exit cleanly           /
  /-------------------------------------------------------------------/
  / N. Hine (based on code by D. Quigley) - University of Warwick     /
  /==================================================================*/
  
  /* Measure the time t2 using MPI_Wtime() which returns a double */
  t2 = MPI_Wtime();

  //set condition to p>=1 so time can be displayed for running 1 MPI task
  if (my_rank==0 && p>=1 ) {
    printf("Total time elapsed since MPI initialised :  %12.6f s\n",t2-t1);
  }

  /* Shutdown MPI - insert appropriate call here */
  MPI_Finalize();

  return;

}

			     
  



