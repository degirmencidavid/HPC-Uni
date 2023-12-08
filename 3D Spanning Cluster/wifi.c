#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <unistd.h>
#include <time.h>
#include "mt19937ar.h"

/* Typedefs */

/* variables for a "router", containing position & distance from center,
   router radius, and cluster number */
struct router {
  double x,y,z;
  double m,r;
  int cluster;
};

/* variables describing a domain decomposition of a volume filled with routers,
   containing sizes, number of cells, pointers to the routers in each cell,
   total number of routers, total number of clusters, and whether the
   clusters span the space */
struct cell_domain {
  double S, Lx, Ly, Lz;
  int nx, ny, nz, nc;
  struct router*** cell_rtr;
  int *cell_nrtr;
  int nrtr;
  int ncluster;
  int spanning_cluster;
};

/* Function Prototypes */
int read_args(int argc, char **argv, double* R, double* S, double* D, int* N,
              double* P, double* Q, int* M, int* seedtime);
int isprint(int opt);
void find_cluster(int Nrtr, struct router** Rtr, int i);;
int connected(struct router* ra, struct router* rb);
int find_spanning_cluster(struct cell_domain* dom);
int merge_clusters(int nra, struct router** ra, int nrb, struct router** rb);
int count_clusters(struct cell_domain* dom);
void generate_routers(int* Nrtr, struct router** Rtr, double S,double R,double P);
void create_domain_decomp(struct cell_domain* dom, int Nrtr, struct router* Rtr,
                          int nx, int ny, int nz);
void find_all_clusters(struct cell_domain* dom);
void destroy_domain_decomp(struct cell_domain* dom);

/* Main Routine */
int main (int argc, char **argv)  {

  double R, S, Sinit, P, Pinit;
  int seedtime=0;

  /* Default Values */
  Sinit = 100; Pinit = 0.7; R = 1.0;
  int nP = 1; double deltaP = 0.015;
  int nS = 1; double deltaS = 5.0;
  int cellmin=20, cellmax=20;
  cellmin = (int)round(Sinit/16)*8;
  cellmax = cellmin;
  //printf("%d %d", cellmin, cellmax);

  /* Read command line arguments */
  int ra = read_args(argc,argv,&R,&Sinit,&deltaS,&nS,&Pinit,&deltaP,&nP,
                     &seedtime);

  /* Echo command line arguments to stdout */
  if (!ra) {
    printf("# Command line arguments successfully processed\n");
  }
  if (ra) {
    printf("Error Processing command line arguments\n"); return EXIT_FAILURE;
  }

  /* Seed random number generator */
  unsigned long seed=20350292;
  if (seedtime) {
    seed = time(NULL);
    printf("# Using time-based random seed %ld\n",seed);
  }
  init_genrand(seed);

  /* Loop variables */
  int i,j,iP,iS,irun;
  /* Router-related variables */
  int Nrtr=0;
  struct router* Rtr = NULL;
  /* Cell decomposition information */
  struct cell_domain dom;
  /* Number of runs for this invocation of the program */
  int nruns = nP*nS;

  /* Loop over values of filling fraction P and system size S */
  for (irun=0;irun<nruns;irun++) {
     /* Find values of iP and iS for this run */
     iP = irun%nP;
     iS = (irun-iP)/nP;
     /* Find values of P and S for this run */
     P = Pinit + iP*deltaP; 
     S = Sinit + iS*deltaS;
     dom.S = S; dom.Lx = 2.0*S; dom.Ly = 2.0*S; dom.Lz = 2.0*S;
     /* Generate randomly-placed routers in the domain */
     generate_routers(&Nrtr,&Rtr,S,R,P);
     /* Output sizes and volume fraction (+newline if multiple cell sizes) */
     printf("S = %6.2f P = %8.6f ",S,P);
     if (cellmin!=cellmax) printf("\n");
     /* Loop over domain decomposition grid sizes */
     for (i=cellmin;i<=cellmax;i++){
        if (cellmin!=cellmax) printf("ncells = %3d ",i);
        /* Initialise the domain decomposition structure */
        create_domain_decomp(&dom,Nrtr,Rtr,i,i,i);
        /* Find clusters in cells, merge between cells, count clusters
           and find spanning cluster if it exists */
        find_all_clusters(&dom);
        /* Write results to stdout */
        printf(": %6d clusters, ",dom.ncluster);
        if (dom.spanning_cluster==0) {
           printf("none spanning\n");
        }
        else {
           printf("%7d spans\n",dom.spanning_cluster);
        }
        /* remove storage associated with domain decomposition and
           reset router cluster values */
        destroy_domain_decomp(&dom);
        for (j=0;j<Nrtr;j++)
        {
           Rtr[j].cluster = 0;
        }
     }
     free(Rtr);
  }
 
  return(EXIT_SUCCESS);

}

/* Set up number of routers for a given volume fraction, and
   generate random coordinates for each one, setting cluster
   value to 0 to indicate cluster not found yet */
void generate_routers(int* Nrtr, struct router** Rtr,
                      double S,double R,double P)
{
  *Nrtr = (int)(S*S*S*P/(R*R*R));
  *Rtr = malloc((*Nrtr)*sizeof(struct router));
 
  int i;
  double mag,x,y,z;
  for (i=0;i<(*Nrtr);i++)
  {
    /* Find random coordinates inside the sphere of radius S 
       whose centre is at (S,S,S) */
    mag = S*4.0;
    while (mag>S) {
      x = genrand()*S*2.0; y=genrand()*S*2.0; z=genrand()*S*2.0;
      mag = sqrt((x-S)*(x-S)+(y-S)*(y-S)+(z-S)*(z-S));
    }
    (*Rtr)[i].x = x; (*Rtr)[i].y = y;  (*Rtr)[i].z = z; (*Rtr)[i].m = mag;
    (*Rtr)[i].r = R; (*Rtr)[i].cluster = 0;
  }
}

/* Set up a cartesian grid of cells as a "domain decomposition" for helping
   find clusters of routers which span the space. Assign each router to
   a cell, and set up lists of pointers to the routers in each cell */
void create_domain_decomp(struct cell_domain* dom, int Nrtr, struct router* Rtr,
                          int nx, int ny, int nz)
{
  
  dom->nx = nx; dom->ny = ny; dom->nz = nz;
  dom->nc = nx*ny*nz;
  dom->cell_nrtr = malloc(dom->nx * dom->ny * dom->nz * sizeof(int));
  dom->cell_rtr = malloc(dom->nx * dom->ny * dom->nz * sizeof(struct router*));
  dom->nrtr = Nrtr;
  dom->spanning_cluster = 0;
  
  /* set counts of each cell's routers to zero */
  int ix, iy, iz, ic, i;
  for (ix=0;ix<dom->nx;ix++) {
    for (iy=0;iy<dom->ny;iy++) { 
      for (iz=0;iz<dom->nz;iz++) {
        ic = ix+dom->ny*iy+dom->ny*dom->nz*iz;
        dom->cell_nrtr[ic] = 0;
      }
    }
  }

  /* count the routers associated with each cell */
  for (i=0;i<Nrtr;i++)
  {
    ix = floor(dom->nx * Rtr[i].x / dom->Lx);
    iy = floor(dom->ny * Rtr[i].y / dom->Ly);
    iz = floor(dom->nz * Rtr[i].z / dom->Lz);
    ic = ix+dom->ny*iy+dom->ny*dom->nz*iz;
    dom->cell_nrtr[ic]++;
  }

  /* allocate memory for storage of each cell's set of pointers to routers,
     then re-set the counts of numbers of routers to zero */
  for (ic=0;ic<dom->nc;ic++) {
    dom->cell_rtr[ic] = malloc(dom->cell_nrtr[ic]*sizeof(struct router*));
    dom->cell_nrtr[ic] = 0;
  }

  /* fill up the list of routers associated with each cell, counting them
     again as we go for indexing purposes */
  for (i=0;i<Nrtr;i++)
  {
    ix = floor(dom->nx * Rtr[i].x / dom->Lx);
    iy = floor(dom->ny * Rtr[i].y / dom->Ly);
    iz = floor(dom->nz * Rtr[i].z / dom->Lz);
    ic = ix+dom->ny*iy+dom->ny*dom->nz*iz;
    dom->cell_rtr[ic][dom->cell_nrtr[ic]] = &Rtr[i];
    dom->cell_nrtr[ic]++;
  }
}

/* deallocate storage associated with a domain decomposition struct */
void destroy_domain_decomp(struct cell_domain* dom)
{
  int ic;
  for (ic=0;ic<dom->nc;ic++) {
    free(dom->cell_rtr[ic]);
  }
  free(dom->cell_rtr);
  free(dom->cell_nrtr);
}

typedef struct cube_region {
  int startX, startY, startZ;
  int size;
} cube_t;

int to1d(int x, int y, int z, int domSize) {
  return x + y * domSize + z * domSize * domSize;
}

void merge_subcube(struct router*** routers, int* nrouters, cube_t* cube, int domSize) {
  int changed = 1;
  int neighb, neighbourIndex;
  int dx,dy,dz;
  int cx,cy,cz;
  int index;
  while (changed) {
    changed = 0;
    for (cx = cube->startX; cx < cube->startX + cube->size; cx++) {
      for (cy = cube->startY; cy < cube->startY + cube->size; cy++) {
        for (cz = cube->startZ; cz < cube->startZ + cube->size; cz++) {
          index = to1d(cx, cy, cz, domSize);
          for (neighb = 14; neighb < 27; neighb++) {
            dx = neighb % 3 - 1;
            dy = ((neighb - dx) / 3) % 3 - 1;
            dz = (neighb - dx - 3 * dy) / 9 - 1;

            if ((cx + dx >= cube->startX + cube->size) || (cy + dy >= cube->startY + cube->size) || (cz + dz >= cube->startZ + cube->size)) continue;
            if ((cx + dx < 0) || (cy + dy < 0) || (cz + dz < 0)) continue;
            // find index of neighbour cell
            neighbourIndex = index + to1d(dx, dy, dz, domSize);
            changed = changed + merge_clusters(nrouters[index],
              routers[index], nrouters[neighbourIndex], routers[neighbourIndex]);
          }
        }
      }
    }
  }
}

void join_cubes_x(struct router*** routers, int* nrouters, cube_t* first, int direction, int domSize) {
  int cx = first->startX + (direction == 1 ? first->size - 1 : 0);
  int cy, cz;
  int i1, i2;
  int ox = cx + direction;
  int oy, oz;
  int dchange, dy, dz;
  int changed = 1;
  while (changed) {
    changed = 0;
    for (cy = first->startY; cy < first->startY + first->size; cy++) {
      for (cz = first->startZ; cz < first->startZ + first->size; cz++) {
        i1 = to1d(cx, cy, cz, domSize);
        for (dchange = 0; dchange < 9; dchange++) {
          dy = dchange % 3;
          dz = dchange / 3;
          oy = cy - 1 + dy;
          oz = cz - 1 + dz;
          if (oy < 0 || oz < 0 || oy >= domSize || oz >= domSize) continue;
          i2 = to1d(ox, oy, oz, domSize);
          changed += merge_clusters(nrouters[i1], routers[i1], nrouters[i2], routers[i2]);
        }
      }
    }
  }
}

void join_cubes_y(struct router*** routers, int* nrouters, cube_t* first, int direction, int domSize) {
  int cy = first->startY + (direction == 1 ? first->size - 1 : 0);
  int cx, cz;
  int i1, i2;
  int oy = cy + direction;
  int ox, oz;
  int dchange, dx, dz;
  int changed = 1;
  while (changed) {
    changed = 0;
    for (cx = first->startX; cx < first->startX + first->size; cx++) {
      for (cz = first->startZ; cz < first->startZ + first->size; cz++) {
        i1 = to1d(cx, cy, cz, domSize);
        for (dchange = 0; dchange < 9; dchange++) {
          dx = dchange % 3;
          dz = dchange / 3;
          ox = cx - 1 + dx;
          oz = cz - 1 + dz;
          if (ox < 0 || oz < 0 || ox >= domSize || oz >= domSize) continue;
          i2 = to1d(ox, oy, oz, domSize);
          changed += merge_clusters(nrouters[i1], routers[i1], nrouters[i2], routers[i2]);
        }
      }
    }
  }
}

void join_cubes_z(struct router*** routers, int* nrouters, cube_t* first, int direction, int domSize) {
  int cz = first->startZ + (direction == 1 ? first->size - 1 : 0);
  int cx, cy;
  int i1, i2;
  int oz = cz + direction;
  int ox, oy;
  int dchange, dx, dy;
  int changed = 1;
  while (changed) {
    changed = 0;
    for (cx = first->startX; cx < first->startX + first->size; cx++) {
      for (cy = first->startY; cy < first->startY + first->size; cy++) {
        i1 = to1d(cx, cy, cz, domSize);
        for (dchange = 0; dchange < 9; dchange++) {
          dx = dchange % 3;
          dy = dchange / 3;
          ox = cx - 1 + dx;
          oy = cy - 1 + dy;
          if (ox < 0 || oy < 0 || ox >= domSize || oy >= domSize) continue;
          i2 = to1d(ox, oy, oz, domSize);
          changed += merge_clusters(nrouters[i1], routers[i1], nrouters[i2], routers[i2]);
        }
      }
    }
  }
}

/*
static int cube_adj[8][3] = {
  { 1, 2, 4 },
  { 0, 3, 5 },
  { 3, 0, 6 },
  { 2, 1, 7 },
  { 5, 6, 0 },
  { 4, 7, 1 },
  { 7, 4, 2 },
  { 6, 5, 3 }
};
*/

static int cube_directions[8][3] = {
  { 1, 1, 1 },
  { -1, 1, 1 },
  { 1, -1, 1 },
  { -1, -1, 1 },
  { 1, 1, -1 },
  { -1, 1, -1 },
  { 1, -1, -1 },
  { -1, -1, -1 }
};

void join_cubes(struct router*** routers, int* nrouters, cube_t** cubes, int start, int domSize, cube_t* bounding) {
  int i, xd, yd, zd;
  for (i = start; i < start + 8; i++) {
    xd = cube_directions[i - start][0];
    yd = cube_directions[i - start][1];
    zd = cube_directions[i - start][2];
    if (xd != 0) {
      join_cubes_x(routers, nrouters, cubes[i], xd, domSize);
    }
    if (yd != 0) {
      join_cubes_y(routers, nrouters, cubes[i], yd, domSize);
    }
    if (zd != 0) {
      join_cubes_z(routers, nrouters, cubes[i], zd, domSize);
    }
  }

  bounding->startX = cubes[start]->startX;
  bounding->startY = cubes[start]->startY;
  bounding->startZ = cubes[start]->startZ;
  bounding->size = 2 * cubes[start]->size;
}

cube_t* gen_cube(int startX, int startY, int startZ, int size) {
  cube_t* cube = (cube_t*)malloc(sizeof(cube_t));
  cube->startX = startX;
  cube->startY = startY;
  cube->startZ = startZ;
  cube->size = size;
  return cube;
}

void split(cube_t** regions, cube_t* region, int offset) {
  int c;
  int splitSize = region->size / 2;
  for (c = 0; c < 8; c++) {
    int startX = region->startX + (c & 1) * splitSize;
    int startY = region->startY + (c >> 1 & 1) * splitSize;
    int startZ = region->startZ + (c >> 2 & 1) * splitSize;
    regions[c + offset] = gen_cube(startX, startY, startZ, splitSize);
  }
}

void print_cube(cube_t* cube) {
  printf("cube { %d, %d, %d, %d }\n", cube->startX, cube->startY, cube->startZ, cube->size);
}

/* subroutine to identify all clusters in a domain decomposition structure */
void find_all_clusters(struct cell_domain* dom)
{

  double t1 = omp_get_wtime();
  int cl = 1;
  int ic,i;
  
  /* loop over all cells of the domain decomposition, then loop over the routers
     in each cell. If cluster has not yet been identified, find all the 
     connected routers within this cell's list of routers (fast!) */
  for (ic=0;ic<dom->nc;ic++) {
    for (i=0;i<dom->cell_nrtr[ic];i++) {
      if (dom->cell_rtr[ic][i]->cluster==0) {
        dom->cell_rtr[ic][i]->cluster = cl;
        find_cluster(dom->cell_nrtr[ic],dom->cell_rtr[ic],i);
        cl++;
      }
    }
  }
  double t2 = omp_get_wtime();

  /* merge clusters between cells if they are connected. Start from first
     cell and move outwards, checking the "outward" half of the set of nearest
     neighbour cells. Always retain lower numbered cluster to prevent circular
     "flows" of cluster value */
  /* keep repeating loop until nothing changes any more
  while(changed) {
    changed = 0;
    for (ic=0;ic<dom->nx*dom->ny*dom->nz;ic++) {
      ix = ic%dom->nx;
      iy = ((ic-ix)/dom->nx)%dom->ny;
      iz = (ic-ix-dom->ny*iy)/(dom->nx*dom->ny);
      for (neighb=14;neighb<27;neighb++) {
        dx = neighb%3 - 1;
        dy = ((neighb-dx)/3)%3 - 1;
        dz = (neighb-dx-3*dy)/9 - 1;
        if ((ix+dx>=dom->nx)||(iy+dy>=dom->ny)||(iz+dz>=dom->nz)) continue;
        if ((ix+dx<0)||(iy+dy<0)||(iz+dz<0)) continue;
        icp = ic + dx + dom->ny*dy + dom->ny*dom->nz*dz;
        changed = changed + merge_clusters(dom->cell_nrtr[ic],
            dom->cell_rtr[ic], dom->cell_nrtr[icp], dom->cell_rtr[icp]);
      }
    }
  } */

  // even cells
  if (dom->nx % 2 == 0) {
    int numthreads = omp_get_max_threads();
    int maxCubes = dom->nx * dom->ny * dom->nz;
    int splits = 0;
    if (numthreads > 1 && numthreads <= 8 && maxCubes >= 8) {
      splits = 1;
    } else if (numthreads > 8 && maxCubes >= 64) {
      splits = 2;
    }
    int ncubes = (int)pow(8, splits);
    cube_t domCube = { .startX = 0, .startY = 0, .startZ = 0, .size = dom->nx };
    cube_t** regions = (cube_t**)malloc(ncubes * sizeof(cube_t*));
    regions[0] = &domCube;
    if (splits >= 1) {
      split(regions, &domCube, 0);

      if (splits == 2) {
        int c = 0;
        for (c = 0; c < 8; c++) {
          int inv = 8 - c - 1;
          split(regions, regions[c], inv * 8);
        }
      }
    }

    for (int a = 0; a < ncubes; a++) {
      //print_cube(regions[a]);
    }

    int ncube = 0;
    #pragma omp parallel for
    for (ncube = 0; ncube < ncubes; ncube++) {
      //printf("thread num %d\n", omp_get_thread_num());
      cube_t* cube = regions[ncube];
      merge_subcube(dom->cell_rtr, dom->cell_nrtr, cube, dom->nx);
    }

    if (splits == 2) {
      int nc;
      cube_t** supercubes = (cube_t**)malloc(sizeof(cube_t*) * 8);
      for (nc = 0; nc < 8; nc++) {
        supercubes[nc] = (cube_t*)malloc(sizeof(cube_t));
        join_cubes(dom->cell_rtr, dom->cell_nrtr, regions, nc * 8, dom->nx, supercubes[nc]);
        //print_cube(supercubes[nc]);
      }
      //printf("joined supercube super\n");
      //fflush(stdout);
      //dealloc regions
      regions = supercubes;

      int flip;
      for (flip = 0; flip < 4; flip++) {
        cube_t* swap = regions[flip];
        regions[flip] = regions[7 - flip];
        regions[7 - flip] = swap;
      }
    }
    
    cube_t bounding;
    
    //printf("%p\n", &bounding);
    join_cubes(dom->cell_rtr, dom->cell_nrtr, regions, 0, dom->nx, &bounding);
    // remember to dealloc
  }

  //printf("done\n");
  //fflush(stdout);


  double t3 = omp_get_wtime();
  dom->ncluster = count_clusters(dom);

  //printf("counted\n");
  //fflush(stdout);
  dom->spanning_cluster = find_spanning_cluster(dom);
  //printf("spanning\n");
  //fflush(stdout);
  double t4 = omp_get_wtime();
  /* print timings - you may need to disable this in parallel if
     it is causing problems */
  printf("(%8.4f %8.4f %8.4f sec) ",
      (double)(t2-t1), 
      (double)(t3-t2), 
      (double)(t4-t3)/*,
      (double)(t4-t1)*/);
}

/* recursive subroutine that finds all the "connected" routers in a list,
   starting from a specific router i, and sets their cluster value to
   match that of the starting router */ 
void find_cluster(int Nrtr, struct router** Rtr, int i)
{
  int j;
  for (j=0;j<Nrtr;j++)
  {
    if (connected(Rtr[i],Rtr[j])&&(Rtr[j]->cluster!=Rtr[i]->cluster))
    {
      Rtr[j]->cluster = Rtr[i]->cluster;
      find_cluster(Nrtr,Rtr,j);
    }
  }
}

/* function to check if two routers are "connected", ie their separation is
   less than the sum of their radii */
int connected(struct router* ra, struct router* rb)
{
   if ((ra->x-rb->x)*(ra->x-rb->x) + 
       (ra->y-rb->y)*(ra->y-rb->y) +
       (ra->z-rb->z)*(ra->z-rb->z) <= (ra->r+rb->r)*(ra->r+rb->r))
   { return 1;}
   else
   { return 0;}
}

/* function to merge the clusters associated with two lists of routers if
   pairs of them are closer than the sum of their radii. Retains lower-
   numbered cluster and converts higher-numbered cluster to match */
int merge_clusters(int nra, struct router** ra, int nrb, struct router** rb)
{
  int i,j,k,cl,changed;
  changed = 0;
  /* Loop over routers i, j in the two cells */
  for (i=0;i<nra;i++) {
    for (j=0;j<nrb;j++) {
      /* search for pairs of routers whose cluster values are not equal */
      if (ra[i]->cluster!=rb[j]->cluster) {
        /* if they are closer than the sum of their radii ...*/
        if (connected(ra[i],rb[j])) {
          /* convert whichever cluster has the higher index to match the
             lower index */
          if (rb[j]->cluster > ra[i]->cluster) {
             /* cluster in cell A has lower index */
            cl = rb[j]->cluster;
            for (k=0;k<nrb;k++) {
              if (rb[k]->cluster==cl) {
                rb[k]->cluster = ra[i]->cluster;
              }
            }
          }  /* else cluster in cell B has lower index */
          else {
            cl = ra[i]->cluster;
            for (k=0;k<nra;k++) {
              if (ra[k]->cluster==cl) {
                ra[k]->cluster = rb[j]->cluster;
              }
            }
          }
          /* set flag to remember that something changed within this call so
             that we can halt merging once nothing changes any more */
          changed = 1;
        }
      }
    }
  }
  return changed;
}

/* Count the number of unique clusters in the list of routers */
int count_clusters(struct cell_domain* dom)
{
  int count, found;
  int* counted = malloc(dom->nrtr*sizeof(int));
  int ic,i,j;
  count = 0;
  for (ic=0;ic<dom->nc;ic++) {
    for (i=0;i<dom->cell_nrtr[ic];i++) {
      found = 0;
      /* check if we have already counted this cluster */
      for (j=0;j<count;j++) {
        if (dom->cell_rtr[ic][i]->cluster==counted[j]) {
          found = 1;
          break;
        }
      }
      /* if not, add it to the list and increment the count */
      if (!found) {
        counted[count] = dom->cell_rtr[ic][i]->cluster;
        count++;
      }
    }
  }
  free(counted);
  return count;
}

/* check if there are clusters extending to the surface in each octant */
int find_spanning_cluster(struct cell_domain* dom)
{
  int* octant_check_list[8];
 
  int ioct,ic,cl,i,icl,jcl;
  int spanning_cluster = 0;
  double x,y,z,r,m;

  /* Set up and initialise storage for finding spanning cluster */
  for (ioct=0;ioct<8;ioct++) {
    octant_check_list[ioct] = malloc(dom->ncluster*sizeof(int));
    for (i=0;i<dom->ncluster;i++) {
       octant_check_list[ioct][i] = -1;
    }
  }

  /* Loop over cells */
  for (ic=0;ic<dom->nc;ic++) {
     /* Loop over routers in each cell */
     for (i=0;i<dom->cell_nrtr[ic];i++) {
        x = dom->cell_rtr[ic][i]->x;
        y = dom->cell_rtr[ic][i]->y;
        z = dom->cell_rtr[ic][i]->z;
        m = dom->cell_rtr[ic][i]->m;
        r = dom->cell_rtr[ic][i]->r;
        cl = dom->cell_rtr[ic][i]->cluster;
        /* If router touches the sphere edge... */
        if (m+r>dom->S) {
           ioct = 0;
           /* Calculate which octant router is in */
           if (x>dom->S) ioct++;
           if (y>dom->S) ioct+=2;
           if (z>dom->S) ioct+=4;
           jcl = 0;
           /* Check if this cluster has been found for this octant */
           for(icl=0;icl<dom->ncluster;icl++) {
              /* If we reach the end of the list of clusters, add this one at 
                 the end and stop looking */
              if (octant_check_list[ioct][icl]==-1) {
                 octant_check_list[ioct][icl] = cl;
                 jcl = 1;
                 break;
              }
              /* If this cluster has already been listed, stop looking */
              if (octant_check_list[ioct][icl] == cl) {break;}
           }
        }
     }
  }
  int sum;
  /* Check if the same cluster appears in all 8 lists */
  for (icl=0;icl<dom->ncluster;icl++) {
     sum = 0;
     /* Get next cluster from octant 0, quit if blank */
     cl = octant_check_list[0][icl];
     if (cl==-1) break;
     /* Check all 8 octants to check this cluster appears in each */
     for (ioct=0;ioct<8;ioct++) {
        for (jcl=0;jcl<dom->ncluster;jcl++) {
           if (octant_check_list[ioct][jcl]==cl) {
              sum++;
              break;
           }
        }
     } 
     if (sum==8) { spanning_cluster++;} 
  }
  /* Free storage associated with spanning cluster check */
  for (ioct=0;ioct<8;ioct++) {
    free(octant_check_list[ioct]);
  }
  return spanning_cluster;
}

/* Parse Command line arguments */
int read_args(int argc, char **argv, double* R, double* S, double* D, int* N,
              double* P, double* Q, int* M, int* seedtime)
{
  int index;
  int c;
  opterr = 0;

  /* Process all flags found in argc */
  while ((c = getopt(argc, argv, "R:S:D:N:P:Q:M:ts")) != -1)
    switch (c)
      {
      case 't':
        *seedtime = 1;
        break;
      case 'R':
        *R = atof(optarg);
        if (*R<=0.0) {
           opterr = 1; 
           fprintf (stderr, "Radius argument could not be read: %s\n", optarg);
        }
        break;
      case 'S':
        *S = atof(optarg);
        if (*S<=0.0) {
           opterr = 1; 
           fprintf (stderr, "Station size argument could not be read: %s\n", optarg);
        }
        break;
      case 'P':
        *P = atof(optarg);
        if (*P<=0.0) {
           opterr = 1; 
           fprintf (stderr, "Volume fraction argument could not be read: %s\n", optarg);
        }
        break;
      case 'D':
        *D = atof(optarg);
        if (*D<0.0) {
           opterr = 1; 
           fprintf (stderr, "Station size step argument could not be read: %s\n", optarg);
        }
        break;
      case 'Q':
        *Q = atof(optarg);
        if (*Q<0.0) {
           opterr = 1; 
           fprintf (stderr, "Volume fraction step argument could not be read: %s\n", optarg);
        }
        break;
      case 'N':
        *N = atoi(optarg);
        if (*N<=0.0) {
           opterr = 1; 
           fprintf (stderr, "Number of station size steps argument could not be read: %s\n", optarg);
        }
        break;
      case 'M':
        *M = atof(optarg);
        if (*M<=0.0) {
           opterr = 1; 
           fprintf (stderr, "Number of volume fraction steps argument could not be read: %s\n", optarg);
        }
        break;
      case '?':
        if ((optopt == 'R')||(optopt=='P')||(optopt=='S'))
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint(optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,"Unknown option character `\\x%x'.\n",optopt);
        return 1;
      default:
        abort ();
      }

  /* List unrecognised arguments */
  for (index = optind; index < argc; index++)
    printf ("Non-option argument %s\n", argv[index]);

  return opterr;
}
