#include <iostream>
#include "main.h"

using namespace std;
using level_t = level<Real, LL, WW, HH, _LEVELS>;

// Computation of the non-linear term d(Vf dot vf).
void dVfvf(double recipr,double arr[], double out[])
{
  for(int i=0; i<L; i++) {
    for(int j=0; j<W; j++) {
      for(int k=0; k<H; k++) {
        out[I(X,i,j,k)] = recipr*(arr[I(X,(i!=(L-1)) ? (i+1) : 0,j,k)]*arr[I(X,(i!=(L-1)) ? (i+1) : 0,j,k)] - arr[I(X,(i!=0) ? (i-1) : (L-1),j,k)]*arr[I(X,(i!=0) ? (i-1) : (L-1),j,k)] + arr[I(X,i,(j!=(W-1)) ? (j+1) : 0,k)]*arr[I(Y,i,(j!=(W-1)) ? (j+1) : 0,k)] - arr[I(X,i,(j!=0) ? (j-1) : (W-1),k)]*arr[I(Y,i,(j!=0) ? (j-1) : (W-1),k)] + arr[I(X,i,j,(k!=(H-1)) ? (k+1) : 0)]*arr[I(Z,i,j,(k!=(H-1)) ? (k+1) : 0)] - arr[I(X,i,j,(k!=0) ? (k-1) : (H-1))]*arr[I(Z,i,j,(k!=0) ? (k-1) : (H-1))]);
        out[I(Y,i,j,k)] = recipr*(arr[I(Y,(i!=(L-1)) ? (i+1) : 0,j,k)]*arr[I(X,(i!=(L-1)) ? (i+1) : 0,j,k)] - arr[I(Y,(i!=0) ? (i-1) : (L-1),j,k)]*arr[I(X,(i!=0) ? (i-1) : (L-1),j,k)] + arr[I(Y,i,(j!=(W-1)) ? (j+1) : 0,k)]*arr[I(Y,i,(j!=(W-1)) ? (j+1) : 0,k)] - arr[I(Y,i,(j!=0) ? (j-1) : (W-1),k)]*arr[I(Y,i,(j!=0) ? (j-1) : (W-1),k)] + arr[I(Y,i,j,(k!=(H-1)) ? (k+1) : 0)]*arr[I(Z,i,j,(k!=(H-1)) ? (k+1) : 0)] - arr[I(Y,i,j,(k!=0) ? (k-1) : (H-1))]*arr[I(Z,i,j,(k!=0) ? (k-1) : (H-1))]);
        out[I(Z,i,j,k)] = recipr*(arr[I(Z,(i!=(L-1)) ? (i+1) : 0,j,k)]*arr[I(X,(i!=(L-1)) ? (i+1) : 0,j,k)] - arr[I(Z,(i!=0) ? (i-1) : (L-1),j,k)]*arr[I(X,(i!=0) ? (i-1) : (L-1),j,k)] + arr[I(Z,i,(j!=(W-1)) ? (j+1) : 0,k)]*arr[I(Y,i,(j!=(W-1)) ? (j+1) : 0,k)] - arr[I(Z,i,(j!=0) ? (j-1) : (W-1),k)]*arr[I(Y,i,(j!=0) ? (j-1) : (W-1),k)] + arr[I(Z,i,j,(k!=(H-1)) ? (k+1) : 0)]*arr[I(Z,i,j,(k!=(H-1)) ? (k+1) : 0)] - arr[I(Z,i,j,(k!=0) ? (k-1) : (H-1))]*arr[I(Z,i,j,(k!=0) ? (k-1) : (H-1))]);
      }
    }
  }
}

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


int checkCFL(double dt, double h, double nu, double umax)
{
   int CF_FLAG = 0;
   // double umax =  maxElement(arr);
   double dt_cfl = std::min(h/umax, h*h/nu);
   //PRAO: Dennis wanted to put in the extra factor of 4 in checking cfl
   if (dt_cfl < dt/4.0 ) CF_FLAG = 1;

   return CF_FLAG;
}
