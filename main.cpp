#include "main.h"


int main(int argc, char *argv[])
{
  double dt = 0.01; // timestep size
  double T = 0.0; // physical time

  // Global interior grid
  int ggrid[NDIMS] = {32,32,32}; // global interior grid pts on entire domain

  // Initialize MPI
  int numprocs, numprocs_dir[NDIMS], rank;
  int periodic[NDIMS], reorder = 0;
  int nbrs[2*NDIMS], coords[NDIMS];
  int igrid[NDIMS]; 

  MPI_Comm cartcomm;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  // Processor grid
  for (int i =0; i<NDIMS; ++i)
  {
    if (NDIMS==3)
      numprocs_dir[i] = cbrt(numprocs);
    else if (NDIMS==2)
      numprocs_dir[i] = sqrt(numprocs);

    periodic[i] = 1; // bdry flag for MPI cartesian grid 
    igrid[i] = ggrid[i]/numprocs_dir[i]; // interior grid pts on each proc 
  }

  MPI_Cart_create(MPI_COMM_WORLD, NDIMS, numprocs_dir, periodic,
                  reorder, &cartcomm);
  MPI_Comm_rank(cartcomm, &rank);
  MPI_Cart_coords(cartcomm, rank, NDIMS, coords);

  // Find the neighbors
  for (int i=0; i<NDIMS; ++i)
  {
    int ierr = MPI_Cart_shift(cartcomm, i, 1, &nbrs[2*i], &nbrs[2*i+1]);
  }

  // Initialize grid
  Grid coarse_grid(ggrid[0], ggrid[1], ggrid[2]);

  // Initialize velocity
  Field velocity(3, &coarse_grid);
  velocity.initialize();
  // Fill in the boundary cells
  velocity.update_bdry('PERIODIC');

  // Initialize pressure
  Field pressure(1, &coarse_grid);
  pressure.initialize();
  // Fill in the boundary cells
  pressure.update_bdry('PERIODIC');

  // Initialize time integration object
  //PRAO: Which of the two, MPI_COMM_WORLD or cartcomm, be passed into this function?
  TimeIntegration ts(&velocity, &pressure, dt, MPI_COMM_WORLD, coords, igrid, ggrid);

  // Propagate in time
  for (int itr = 0; itr < nsteps; itr++)
  {
   // Pick a time-stepping scheme out of EULER, RK2 and RK4 (later)
   ts.time_stepping('EULER'); 

   // Increment physical time  
   T += dt;  
  }

  // Finalize MPI
   MPI_Finalize();
  
  return 0;
}
