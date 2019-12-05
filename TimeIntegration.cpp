#include "main.h"
#include "Field.h"

using namespace std;

TimeIntegration::TimeIntegration(Field* obj_vel, Field* obj_pres, double delta_t, char method_name):
L(vel->obj->L),W(vel->obj->W),H(vel->obj->H)
{
  T = 0.0;
  dt = delta_t;
  vel = obj_vel; 
  pres = obj_pres;
 
  L = vel->obj->L;
  W = vel->obj->W;
  H = vel->obj->H;
 
  int LL = L/2;
  int WW = W/2;
  int HH = H/2;
  

  double *Vfvf_res = new double[vel->ARR_SIZE];
  cout << "TimeIntegration object is being created" << endl;
}

TimeIntegration::~TimeIntegration(void)
{
  cout << "Field object is being deleted" << endl;
}


// Computation of the non-linear term d(Vf dot vf).
void TimeIntegration::dVfvf(double out[])
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

//Leray projection using multigrid solver.
void TimeIntegration::pressure_solve(MPI_Comm mpi_comm, double* div_vstar, double* pressure)
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
  /*for (int i=0; i<NDIMS; ++i)
  {
    ilower[i] = coords[i]*nblock[i];
    iupper[i] = coords[i]*nblock[i] + (nblock[i]-1);
  }*/

  //HYPRE_StructGridSetExtents(grid, ilower, iupper);
 
 /* This is a collective call finalizing the grid assembly.
     The grid is now "ready to be used" */
  //HYPRE_StructGridSetPeriodic(grid, phyGrid);
  //HYPRE_StructMatrixSetNumGhost(A, NGUARD); // Not sure what it does
  //HYPRE_StructGridAssemble(grid);

}

// Computes (Vf wedge vf) - nu*Laplacian(V)
void TimeIntegration::nav_stoke(double out[]) 
{
  dVfvf(Vfvf_res); // calculate the non-linear term
  vel->laplacian(Vfvf_res); // add the laplacian term
}

void TimeIntegration::time_stepping(char method)
{
  switch(method)
  {
    default:
    case 'EULER':
       cout << "Time integration using Forward Euler emthod" << endl;
       forward_euler();
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
  double *R = new double[vel->ARR_SIZE];
  double *vstar = new double[vel->ARR_SIZE]; 

  nav_stoke(R); // sum of non-linear and laplacian terms
  for (int i = 0; i < vel->ARR_SIZE; i++) vstar[i] = vel->arr[i] - R[i] * dt;
  // calculate divergence of the intermediate velocity
  //PRAO: How do I make it so bd10 can calculate the divergence of the intermediare velocity?
  //bd10(vstar, div_vstar); 
  // pass the divergence of vstar in the pressure solve
  T += dt;
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

