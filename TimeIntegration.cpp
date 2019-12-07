#include "main.h"

TimeIntegration::TimeIntegration(Field* obj_vel, Field* obj_pres, double delta_t, MPI_Comm mpi_communicator, int* mpi_coords, int* int_grid, int* glo_grid)
{
  dt = delta_t;
  vel = obj_vel; 
  pres = obj_pres;
 
  L = vel->obj->L;
  W = vel->obj->W;
  H = vel->obj->H;

  mpi_comm = mpi_communicator;
  coords = mpi_coords;
  igrid = int_grid;
  ggrid = glo_grid;
  cout << "TimeIntegration object is being created" << endl;
}

TimeIntegration::~TimeIntegration(void)
{
  cout << "TimeIntegration object is being deleted" << endl;
}

//PRAO: this is temporary
// w refers to x or y or z components
// w is 0 for scalars
int TimeIntegration::I(int w, int i, int j, int k, int NSIZE)
{
  return (i*W*H*NSIZE + j*H*NSIZE + k*NSIZE + w);
}

//PRAO: this is temporary
// Function for calculating divergence (which is the boundary map of the 1-chain.) Outputs to a 0 chain of dimension L,W,H
void TimeIntegration::bd10(double* in_arr, double* out)
{
  for(int i=NGUARD; i<L+NGUARD; i++) {
    for(int j=NGUARD; j<W+NGUARD; j++) {
      for(int k=NGUARD; k<H+NGUARD; k++) {

        out[I(X,i,j,k,1)] = in_arr[I(X,i-1,j,k,3)] - in_arr[I(X,i+1,j,k,3)]
                        + in_arr[I(Y,i,j-1,k,3)] - in_arr[I(Y,i,j+1,k,3)]
                        + in_arr[I(Z,i,j,k-1,3)] - in_arr[I(Z,i,j,k+1,3)];
      }
    }
  }
}

// Computation of the non-linear term d(Vf dot vf).
void TimeIntegration::dVfvf(double* out)
{
  for(int i=NGUARD; i<L+NGUARD; i++) {
    for(int j=NGUARD; j<W+NGUARD; j++) {
      for(int k=NGUARD; k<H+NGUARD; k++) {

        double recipr = 1.0/(2.0*vel->obj->dx);
        int ind = vel->I(X,i,j,k);
        out[ind] = recipr*(  vel->arr[vel->I(X,i+1,j,k)] * vel->arr[vel->I(X,i+1,j,k)]
                   - vel->arr[vel->I(X,i-1,j,k)] * vel->arr[vel->I(X,i-1,j,k)] 
                   + vel->arr[vel->I(X,i,j+1,k)] * vel->arr[vel->I(Y,i,j+1,k)]
                   - vel->arr[vel->I(X,i,j-1,k)] * vel->arr[vel->I(Y,i,j-1,k)]
                   + vel->arr[vel->I(X,i,j,k+1)] * vel->arr[vel->I(Z,i,j,k+1)] 
                   - vel->arr[vel->I(X,i,j,k-1)] * vel->arr[vel->I(Z,i,j,k-1)]  );

        out[ind+1] = recipr*(  vel->arr[vel->I(Y,i+1,j,k)] * vel->arr[vel->I(X,i+1,j,k)]
                     - vel->arr[vel->I(Y,i-1,j,k)] * vel->arr[vel->I(X,i-1,j,k)]
                     + vel->arr[vel->I(Y,i,j+1,k)] * vel->arr[vel->I(Y,i,j+1,k)]
                     - vel->arr[vel->I(Y,i,j-1,k)] * vel->arr[vel->I(Y,i,j-1,k)]
                     + vel->arr[vel->I(Y,i,j,k+1)] * vel->arr[vel->I(Z,i,j,k+1)]
                     - vel->arr[vel->I(Y,i,j,k-1)] * vel->arr[vel->I(Z,i,j,k-1)]  );

        out[ind+2] = recipr*(  vel->arr[vel->I(Z,i+1,j,k)] * vel->arr[vel->I(X,i+1,j,k)]
                     - vel->arr[vel->I(Z,i-1,j,k)] * vel->arr[vel->I(X,i-1,j,k)]
                     + vel->arr[vel->I(Z,i,j+1,k)] * vel->arr[vel->I(Y,i,j+1,k)]
                     - vel->arr[vel->I(Z,i,j-1,k)] * vel->arr[vel->I(Y,i,j-1,k)]
                     + vel->arr[vel->I(Z,i,j,k+1)] * vel->arr[vel->I(Z,i,j,k+1)]
                     - vel->arr[vel->I(Z,i,j,k-1)] * vel->arr[vel->I(Z,i,j,k-1)]  );
      }
    }
  }
 
}

//Leray projection using multigrid solver (Boomer AMG) in HYPRE
void TimeIntegration::pressure_solve(double* div_vstar, double* pressure)
{
  // setup hypre data structures here
  HYPRE_StructGrid     grid;
  HYPRE_StructStencil  stencil;
  HYPRE_StructMatrix   A;
  HYPRE_StructVector   b;
  HYPRE_StructVector   x;
  HYPRE_StructSolver   solver;

  /* Create an empty 3D grid object */
  HYPRE_StructGridCreate(mpi_comm, NDIMS, &grid);

  /* Add boxes to the grid */
  // Find x-index of 3d proc_id
  int ilower[NDIMS] = {};
  int iupper[NDIMS] = {};
  for (int i=0; i<NDIMS; ++i)
  {
    ilower[i] = coords[i]*igrid[i];
    iupper[i] = coords[i]*igrid[i] + (igrid[i]-1);
  }

  HYPRE_StructGridSetExtents(grid, ilower, iupper);
 
 /* This is a collective call finalizing the grid assembly.
     The grid is now "ready to be used" */
  HYPRE_StructGridSetPeriodic(grid, ggrid);
  //HYPRE_StructMatrixSetNumGhost(A, NGUARD); // Not sure what it does
  HYPRE_StructGridAssemble(grid);

  /* 2. Define the discretization stencil */
   {
      /* Create an empty 3D, 7-pt stencil object */
      int nentries = 7; // number of stencil points
      HYPRE_StructStencilCreate(NDIMS, nentries, &stencil);

      /* Define the geometry of the stencil. Each represents a
         relative offset (in the index space). */
      {
         //FIXME: C++ doesn't like variable arguments in static arrays...
         int offsets[7][3] = {{0,0,0}, {-1,0,0}, {1,0,0}, {0,-1,0}, {0,1,0}, {0,0,-1}, {0,0,1}};

         /* Assign each of the 5 stencil entries */
         for (int entry = 0; entry < nentries; entry++)
            HYPRE_StructStencilSetElement(stencil, entry, offsets[entry]);
      }
   }

   /* 3. Set up a Struct Matrix */
  {
    /* Create an empty matrix object */
    HYPRE_StructMatrixCreate(mpi_comm, grid, stencil, &A);

    /* Indicate that the matrix coefficients are ready to be set */
    HYPRE_StructMatrixInitialize(A);

    /* Set the matrix coefficients.  Each processor assigns coefficients
         for the boxes in the grid that it owns. Note that the coefficients
         associated with each stencil entry may vary from grid point to grid
         point if desired.  Here, we first set the same stencil entries for
         each grid point.  Then we make modifications to grid points near
         the boundary. */

    {
      int nentries = 7;
      int stencil_indices[7] = {0,1,2,3,4,5,6}; /* labels for the stencil entries -
                                               these correspond to the offsets
                                               defined above */
      //grid_points is interior only
      int grid_points = igrid[0]*igrid[1]*igrid[2];
      double values[grid_points*nentries];

      /* each grid point has 7 stencil entries */
      for (int i = 0; i < grid_points; i++)
      {
         int index = i*nentries;
         values[index] = 6.0; // the diagonal entry is 1st
         for (int j = 1; j < nentries; j++)
            values[index+j] = -1.0;
      }

      HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, nentries,
                                     stencil_indices, values);
    }

    /* This is a collective call finalizing the matrix assembly.
       The matrix is now ``ready to be used'' */
      HYPRE_StructMatrixAssemble(A);
  }


/* 4. Set up Struct Vectors for b and x.  Each processor sets the vectors
      corresponding to its boxes. */
   {
      /* Create an empty vector object */
      HYPRE_StructVectorCreate(mpi_comm, grid, &b);
      HYPRE_StructVectorCreate(mpi_comm, grid, &x);

      /* Indicate that the vector coefficients are ready to be set */
      HYPRE_StructVectorInitialize(b);
      HYPRE_StructVectorInitialize(x);

      /* Set the vector coefficients */
      {
        // grid_points is interior only
         // FIXME: So here we are copying the interior only points into the values
         // array which will be used to set b, so we just need to reverse
         // engineer this below
         int grid_points = igrid[0]*igrid[1]*igrid[2];
         double values[grid_points]; /* grid points */
         int cnt = 0;
         for (int i = NGUARD; i < igrid[0]+NGUARD; i++)
         for (int j = NGUARD; j < igrid[1]+NGUARD; j++)
         for (int k = NGUARD; k < igrid[2]+NGUARD; k++)
         {
           int scalarInd = I(0, i, j, k, 1);
           // Should be negative since the matrix is set up as -Lap(phi) and we are solving
           // Lap(phi) = f(x)
           // the dx^2 is here, this assumes single h spacing, for 2h spacing, this needs
           // to be -4dx^2
           values[cnt] = -(this->vel->obj->dx*this->vel->obj->dx)*div_vstar[scalarInd];
           cnt = cnt + 1;
         }

         HYPRE_StructVectorSetBoxValues(b, ilower, iupper, values);

         for (int i = 0; i < grid_points; i ++)
            values[i] = 0.0;
         HYPRE_StructVectorSetBoxValues(x, ilower, iupper, values);

      }

      /* This is a collective call finalizing the vector assembly.
         The vectors are now "ready to be used" */
      HYPRE_StructVectorAssemble(b);
      HYPRE_StructVectorAssemble(x);
   }

    /* 5. Set up and use a solver (See the Reference Manual for descriptions
      of all of the options.) */
   {
      /* Create an empty PCG Struct solver */
      HYPRE_StructPCGCreate(mpi_comm, &solver);

      /* Set some parameters */
      HYPRE_StructPCGSetTol(solver, 1.0e-06); /* convergence tolerance */
      HYPRE_StructPCGSetPrintLevel(solver, 2); /* amount of info. printed */

      /* Setup and solve */
      HYPRE_StructPCGSetup(solver, A, b, x);
      HYPRE_StructPCGSolve(solver, A, b, x);
   }

  // put solution x into a press array, which has buffer
  // So at this point we have our solution in x, which is the size of 
  // the interior. 
   int grid_points = igrid[0]*igrid[1]*igrid[2]; 
   double temp[grid_points];
   HYPRE_StructVectorGetBoxValues(x, ilower, iupper, temp);

   int cnt = 0;
   for (int i = NGUARD; i < igrid[0]+NGUARD; i++)
   for (int j = NGUARD; j < igrid[1]+NGUARD; j++)
   for (int k = NGUARD; k < igrid[2]+NGUARD; k++)
   {
     // PRAO: not sure how to do the index
     int ind = I(0, i, j, k, 1); //scalar index
     /* get the local solution */
     //printf("BEFORE ind=%d, cnt=%d\n", ind, cnt);
     pres->arr[ind] = temp[cnt];
     //printf("AFTER temp[cnt]=%f\n", temp[cnt]);
     cnt = cnt + 1;
   }

  // print the hypre matrix
  int ierr = HYPRE_StructMatrixPrint("A.out", A, 0);
  ierr = HYPRE_StructVectorPrint("x.out", x, 0);
  ierr = HYPRE_StructVectorPrint("b.out", b, 0);

  /* Free memory */
  HYPRE_StructGridDestroy(grid);
  HYPRE_StructStencilDestroy(stencil);
  HYPRE_StructMatrixDestroy(A);
  HYPRE_StructVectorDestroy(b);
  HYPRE_StructVectorDestroy(x);
  HYPRE_StructPCGDestroy(solver);

}

// Computes (Vf wedge vf) - nu*Laplacian(V)
void TimeIntegration::nav_stoke(double out[]) 
{
  double *Vfvf_res = new double[vel->ARR_SIZE];
  dVfvf(Vfvf_res); // calculate the non-linear term
  vel->laplacian(Vfvf_res); // add the laplacian term
}

void TimeIntegration::time_stepping(char method)
{
  switch(method)
  {
    default:
    case 'EULER':
       cout << "Time integration using Forward Euler method" << endl;
       forward_euler();
       break;
    case  'RK2' :
    case  'RK4' :
       cout << "Time integration using 4th order Runge-Kutta" << endl;
       runge_kutta4();
       break;
    case 'AB4':
       cout << "Time integration using 4th order Adams-Bashford" << endl;
       break;
  }
}

void TimeIntegration::forward_euler() 
{ 
  double *R = new double[vel->ARR_SIZE]; //vector field
  double *vstar = new double[vel->ARR_SIZE]; //vector field
  double *div_vstar = new double[pres->ARR_SIZE]; //scalar field

  nav_stoke(R); // sum of non-linear and laplacian terms
  for (int i = 0; i < vel->ARR_SIZE; i++) vstar[i] = vel->arr[i] - R[i] * dt;

  // calculate divergence of the intermediate velocity
  //PRAO: How do I make it so bd10 can calculate the divergence of the intermediare velocity?
  bd10(vstar, div_vstar); 

  // pass the divergence of vstar in the pressure solve
  pressure_solve(div_vstar, pres->arr);    
}


void TimeIntegration::runge_kutta4()
{
  double *R = new double[vel->ARR_SIZE];
  double *temp = new double[vel->ARR_SIZE]; 
  double *Total = new double[vel->ARR_SIZE];
  nav_stoke(R);
  for (int i = 0; i < vel->ARR_SIZE; i++) temp[i] = vel->arr[i] + R[i] * dt * 0.5;
  for (int i = 0; i < vel->ARR_SIZE; i++)
    Total[i] = R[i];  // separate loop to make compiler to use SIMD
  nav_stoke(R);
  for (int i = 0; i < vel->ARR_SIZE; i++) temp[i] = vel->arr[i] + R[i] * dt * 0.5;
  for (int i = 0; i < vel->ARR_SIZE; i++) Total[i] += 2 * R[i];
  nav_stoke(R);
  for (int i = 0; i < vel->ARR_SIZE; i++) temp[i] = vel->arr[i] + R[i] * dt;
  for (int i = 0; i < vel->ARR_SIZE; i++) Total[i] += 2 * R[i];
  nav_stoke(R);
  for (int i = 0; i < vel->ARR_SIZE; i++) Total[i] += R[i];
  for (int i = 0; i < vel->ARR_SIZE; i++) vel->arr[i] += Total[i] * (dt / 6.0);
  T += dt;
}

