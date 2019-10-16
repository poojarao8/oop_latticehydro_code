#include "main.h"
#include "Field.h"

using namespace std;

TimeIntegration::TimeIntegration(Field* obj_vec, Field* obj_scl, double delta_t, char method)
{
  dt = delta_t;
  vobj = obj_vec; 
  sobj = obj_scl;
  L = vobj->obj->L;
  W = vobj->obj->W;
  H = vobj->obj->H;

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

        double recipr = 1.0/(2.0*vobj->obj->dx);
        int ind = vobj->I(X,i,j,k);
        out[ind] = recipr*(  vobj->arr[vobj->I(X,i+1,j,k)] * vobj->arr[vobj->I(X,i+1,j,k)]
                           - vobj->arr[vobj->I(X,i-1,j,k)] * vobj->arr[vobj->I(X,i-1,j,k)]
                           + vobj->arr[vobj->I(X,i,j+1,k)] * vobj->arr[vobj->I(Y,i,j+1,k)]
                           - vobj->arr[vobj->I(X,i,j-1,k)] * vobj->arr[vobj->I(Y,i,j-1,k)]
                           + vobj->arr[vobj->I(X,i,j,k+1)] * vobj->arr[vobj->I(Z,i,j,k+1)]
                           - vobj->arr[vobj->I(X,i,j,k-1)] * vobj->arr[vobj->I(Z,i,j,k-1)]  );

        out[ind+1] = recipr*(  vobj->arr[vobj->I(Y,i+1,j,k)] * vobj->arr[vobj->I(X,i+1,j,k)]
                             - vobj->arr[vobj->I(Y,i-1,j,k)] * vobj->arr[vobj->I(X,i-1,j,k)]
                             + vobj->arr[vobj->I(Y,i,j+1,k)] * vobj->arr[vobj->I(Y,i,j+1,k)]
                             - vobj->arr[vobj->I(Y,i,j-1,k)] * vobj->arr[vobj->I(Y,i,j-1,k)]
                             + vobj->arr[vobj->I(Y,i,j,k+1)] * vobj->arr[vobj->I(Z,i,j,k+1)]
                             - vobj->arr[vobj->I(Y,i,j,k-1)] * vobj->arr[vobj->I(Z,i,j,k-1)]  );

        out[ind+2] = recipr*(  vobj->arr[vobj->I(Z,i+1,j,k)] * vobj->arr[vobj->I(X,i+1,j,k)]
                             - vobj->arr[vobj->I(Z,i-1,j,k)] * vobj->arr[vobj->I(X,i-1,j,k)]
                             + vobj->arr[vobj->I(Z,i,j+1,k)] * vobj->arr[vobj->I(Y,i,j+1,k)]
                             - vobj->arr[vobj->I(Z,i,j-1,k)] * vobj->arr[vobj->I(Y,i,j-1,k)]
                             + vobj->arr[vobj->I(Z,i,j,k+1)] * vobj->arr[vobj->I(Z,i,j,k+1)]
                             - vobj->arr[vobj->I(Z,i,j,k-1)] * vobj->arr[vobj->I(Z,i,j,k-1)]  );
      }
    }
  }
 
}

/*

//Leray projection using multigrid solver.
void proj(double deg1ch[], double out[], double pressure[], double pressure_old[], double solution_old[]){
  level_t& b = *(new level_t);
  level_t& x = *(new level_t);
  double *temp= new double[L*W*H];
  double tot;
  double avg;
  double thresh = 1.0e-15;
  bd10(deg1ch,pressure);

  for(int i=0; i<(L*W*H); i++) temp[i]=pressure[i];
  for(int i=0; i<(L*W*H); i++) pressure[i] -= pressure_old[i];

int octal;
int K;
for(int b1=0; b1<2; b1++) {for(int b2=0; b2<2; b2++) {for(int b3=0; b3<2; b3++) {
octal=b1*4+b2*2+b3; //Loop over 8 octal
tot=0.0;
  for(int i=0; i<LL; i++) {
    for(int j=0; j<WW; j++) {
      for(int k=0; k<HH; k++) {

        b.dat.at(i,j,k)=pressure[(i*2+b1)*W*H+(j*2+b2)*H+k*2+b3];
        tot += b.dat.at(i,j,k);
      }
    }
  }
if (abs(tot)>1.0e-17) {thresh= abs(tot)*10;}
else  {thresh= 1.0e-15;}

//cout << "octal: " << octal <<endl;
//cout << "total: " <<tot <<endl;
  K=Slash(b, x, 15, 500, thresh);
  //cout<<"Multigrid solver applied "<<K<<" times"<<endl;
    for(int i=0; i<LL; i++) {
      for(int j=0; j<WW; j++) {
        for(int k=0; k<HH; k++) {
          pressure[(i*2+b1)*W*H+(j*2+b2)*H+k*2+b3]=  x.dat.at(i,j,k);
        }
      }
    }
  }}}

  for(int i=0; i<(L*W*H); i++) pressure[i]+=solution_old[i];
  for(int i=0; i<(L*W*H); i++) pressure_old[i]=temp[i];
  for(int i=0; i<(L*W*H); i++) solution_old[i]=pressure[i];

    d01(pressure,out);

    for(int i=0; i<ARR_SIZE; i++) out[i] = deg1ch[i] - out[i];

    delete &x;
    delete &b;
    delete[] temp;

  }

// Computes proj(Vf wedge vf) - nu*Laplacian(V), where proj is the Leray projection.
void nav_stoke(double h,double nu,double V[], double Vfvf_res[], double out[], double pressure[], double pressure_old[], double solution_old[]) {
      double recipr = 1.0/(2.0*h);

      dVfvf(recipr,V, Vfvf_res); //Vfvf_proj is used as a temporary array before applying diffuse()
      addLaplacian(V,Vfvf_res,nu,h); //add the Laplace term. This can be done either before or after the Leray projection.
      proj(Vfvf_res, out, pressure, pressure_old, solution_old);
    }


*/
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
 /* 
  nav_stoke(h, nu, V, Vfvf_res, R, scal, pressure_old,
      solution_old);  // scal is pressure.
  for (int i = 0; i < ARR_SIZE; i++) temp[i] = V[i] + R[i] * dt * 0.5;
  for (int i = 0; i < ARR_SIZE; i++)
    Total[i] = R[i];  // separate loop to make compiler to use SIMD
  nav_stoke(h, nu, temp, Vfvf_res, R, scal, pressure_old, solution_old);
  for (int i = 0; i < ARR_SIZE; i++) temp[i] = V[i] + R[i] * dt * 0.5;
  for (int i = 0; i < ARR_SIZE; i++) Total[i] += 2 * R[i];
  nav_stoke(h, nu, temp, Vfvf_res, R, scal, pressure_old, solution_old);
  for (int i = 0; i < ARR_SIZE; i++) temp[i] = V[i] + R[i] * dt;
  for (int i = 0; i < ARR_SIZE; i++) Total[i] += 2 * R[i];
  nav_stoke(h, nu, temp, Vfvf_res, R, scal, pressure_old, solution_old);
  for (int i = 0; i < ARR_SIZE; i++) Total[i] += R[i];
  for (int i = 0; i < ARR_SIZE; i++) V[i] += Total[i] * (dt / 6);
  t += dt;
  */
}

