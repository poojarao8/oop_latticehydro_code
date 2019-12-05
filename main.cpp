#include "main.h"

using namespace std;

int main(int argc, char *argv[])
{
  // Initialize the grid
  Grid coarse_grid(32, 32, 32);

  // Initialize velocity
  Field velocity(3, &coarse_grid);
  velocity.initialize();
  velocity.update_bdry('PERIODIC');

  // Initialize pressure
  Field pressure(1, &coarse_grid);
  pressure.initialize();
  pressure.update_bdry('PERIODIC');

  // Initialize a time integration object
  double dt = 0.01;
  TimeIntegration ts(velocity, pressure, dt, 'EULER');

  // Propagate in time
  for (int itr = 0; itr < nsteps; itr++)
  {
   // Pick a time-stepping scheme out of EULER, RK2 and RK4 (later)
   //ts.forward_euler(); 

   // Increment the time  
  
  }

  // Initialize MPI
  int numprocs, numprocs_dir[NDIMS], rank;
  int periodic[NDIMS], reorder = 0;
  int nbrs[2*NDIMS], coords[NDIMS];

  MPI_Comm cartcomm;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  for (int i =0; i<NDIMS; ++i)
  {
    if (NDIMS==3)
      numprocs_dir[i] = cbrt(numprocs);
    else if (NDIMS==2)
      numprocs_dir[i] = sqrt(numprocs);

    periodic[i] = 1; 
  }

  MPI_Cart_create(MPI_COMM_WORLD, NDIMS, numprocs_dir, periodic,
                  reorder, &cartcomm);
  MPI_Comm_rank(cartcomm, &rank);
  MPI_Cart_coords(cartcomm, rank, NDIMS, coords);

  // find the neighbors
  for (int i=0; i<NDIMS; ++i)
  {
    int ierr = MPI_Cart_shift(cartcomm, i, 1, &nbrs[2*i], &nbrs[2*i+1]);
  }

  // communicate using MPI
  // call the solver
  
  // Finalize MPI
   MPI_Finalize();

  return 0;
}
