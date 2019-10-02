#include "main.h"
#include "initialize_var.h"

void initialize_var(double* vel, double h, double eta, double coef)
{
  for(int i=0; i<LL; i++) {
    for(int j=0; j<WW; j++) {
      for(int k=0; k<HH; k++) {

        double id = static_cast<double> (i-LL/2);
        double jd = static_cast<double> (j-WW/2);
        double kd = static_cast<double> (k-HH/2);
 
       vel[I(X,i,j,k,LL,WW,HH)] = coef*exp(-eta*(h/2.0)*(h/2.0)*((id*id)+(jd*jd)+(kd*kd)));
      }
    }
  }
}
