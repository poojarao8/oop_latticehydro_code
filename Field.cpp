#include "main.h"
#include "Field.h"

using namespace std;

Field::Field(int nsize, Grid &obj)
{
  NSIZE = nsize; 
  ARR_SIZE = nsize*obj.grid_pts;
  cout << "Grid object is being created" << endl;  
}

Field::~Field(void)
{
  cout << "Field object is being deleted" << endl;
}

// w refers to x or y or z components
// 0 for scalars
//PRAO: Is there a sensible way to avoid passing in grid obj every time 
// indexing function is evoked (it is evoked a lot of times)
// Should I make Field a derived class of Grid??
int Field::I(int w, int i, int j, int k, Grid &obj)
{
  return (i*obj.W*obj.H*NSIZE + j*obj.H*NSIZE + k*NSIZE + w);
}

void Field::initialize(Grid &obj)
{
  for(int i=0; i<obj.L; i++) {
    for(int j=0; j<obj.W; j++) {
      for(int k=0; k<obj.H; k++) {

        double id = static_cast<double> (i-obj.L/2);
        double jd = static_cast<double> (j-obj.W/2);
        double kd = static_cast<double> (k-obj.H/2);

      
        I(X,i,j,k) = 2;
        //vel[I(X,i,j,k)] = coef*exp(-eta*(obj.dx/2.0)*(obj.dx/2.0)
          //                   *((id*id)+(jd*jd)+(kd*kd)));
      }
    }
  }
 
}

