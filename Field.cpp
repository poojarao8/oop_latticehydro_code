#include <iostream>
#include "Field.h"

using namespace std;

Field::Field(int nsize, Grid &obj)
{ 
  ARR_SIZE = nsize*obj.grid_pts;
  cout << "Grid object is being created" << endl;  
}

Field::~Field(void)
{
  cout << "Field object is being deleted" << endl;
}

void Field::initialize(Grid &obj)
{
  for(int i=0; i<obj.L; i++) {
    for(int j=0; j<obj.W; j++) {
      for(int k=0; k<obj.H; k++) {

        double id = static_cast<double> (i-obj.L/2);
        double jd = static_cast<double> (j-obj.W/2);
        double kd = static_cast<double> (k-obj.H/2);

       //vel[I(X,i,j,k,LL,WW,HH)] = coef*exp(-eta*(h/2.0)*(h/2.0)*((id*id)+(jd*jd)+(kd*kd)));
      }
    }
  }
 
}

