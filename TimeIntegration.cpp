#include "main.h"
#include "Field.h"

using namespace std;

TimeIntegration::TimeIntegration(Field* obj_vel, Field* obj_pres, double delta_t, char method):
L(vel->obj->L),W(vel->obj->W),H(vel->obj->H)
{
  T = 0.0;
  dt = delta_t;
  vel = obj_vel; 
  pres = obj_pres;
 
  L = vel->obj->L;
  W = vel->obj->W;
  H = vel->obj->H;
 
  constexpr int LL = L/2;
  constexpr int WW = W/2;
  constexpr int HH = H/2;
  

  double *Vfvf_res = new double[vel->ARR_SIZE];
  double *pressure_old = new double[pres->ARR_SIZE];
  double *pressure_lap = new double[pres->ARR_SIZE];
  double *pressure_lap_old = new double[pres->ARR_SIZE];

  //using level_t = level<Real, LL, WW, HH, _LEVELS>;

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
/*void TimeIntegration::proj(double out[])
{
  level_t& b = *(new level_t);
  level_t& x = *(new level_t);

  double *temp= new double[pres->ARR_SIZE];
  double tot;
  double avg;
  double thresh = 1.0e-15;

  // Vfvf is the non-linear term
  // bd10 is the div function for 
  vel->bd10_NL(Vfvf_res, pressure_lap);

  for(int i=0; i<(pres->ARR_SIZE); i++) temp[i] = pressure_lap[i];
  for(int i=0; i<(pres->ARR_SIZE); i++) pressure_lap[i] -= pressure_lap_old[i];

  int octal;
  int K;
  for(int b1=0; b1<2; b1++) {for(int b2=0; b2<2; b2++) {for(int b3=0; b3<2; b3++) {
  octal=b1*4+b2*2+b3; //Loop over 8 octal
  tot=0.0;
    for(int i=0; i<pres->obj->L/2; i++) {
      for(int j=0; j<pres->obj->W/2; j++) {
        for(int k=0; k<pres->obj->H/2; k++) {

          b.dat.at(i,j,k) = pressure_lap[(i*2+b1)*W*H+(j*2+b2)*H+k*2+b3];
          tot += b.dat.at(i,j,k);
        }
      }
    }
    if (abs(tot)>1.0e-17) {thresh= abs(tot)*10;}
    else  {thresh= 1.0e-15;}

    //cout << "octal: " << octal <<endl;
    //cout << "total: " <<tot <<endl;
    K = Slash(b, x, 15, 500, thresh);
    //cout<<"Multigrid solver applied "<<K<<" times"<<endl;
    for(int i=0; i<pres->obj->L/2; i++) {
      for(int j=0; j<pres->obj->W/2; j++) {
        for(int k=0; k<pres->obj->H/2; k++) {
          pres->arr[(i*2+b1)*pres->obj->W*pres->obj->H+(j*2+b2)*pres->obj->H+k*2+b3] = x.dat.at(i,j,k);
        }
      }
    }
  }}}

  for(int i=0; i<pres->ARR_SIZE; i++) pressure_lap[i] += pressure_lap_old[i];
  for(int i=0; i<pres->ARR_SIZE; i++) pressure_lap_old[i] = temp[i];
  for(int i=0; i<pres->ARR_SIZE; i++) pressure_old[i] = pres->arr[i];

  pres->d01(out);

  for(int i=0; i< vel->ARR_SIZE; i++) out[i] = Vfvf_res[i] - out[i];

  delete &x;
  delete &b;
  delete[] temp;

}
*/
// Computes proj(Vf wedge vf) - nu*Laplacian(V), where proj is the Leray projection.
void TimeIntegration::nav_stoke(double out[]) 
{
  dVfvf(Vfvf_res); //Vfvf_proj is used as a temporary array before applying diffuse()
  vel->laplacian(Vfvf_res); //add the Laplace term. This can be done either before or after the Leray projection.
  //proj(out); // pressure projection
}

void TimeIntegration::time_stepping(char method)
{
  switch(method)
  {
    default:
    case 'RK4' :
       cout << "Time integration using 4th order Runge-Kutta" << endl;
       runge_kutta4();
       break;
    case 'AB4' :
       cout << "Time integration using 4th order Adams-Bashford" << endl;
       break;
  }
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

