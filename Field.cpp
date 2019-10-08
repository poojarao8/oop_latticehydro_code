#include "main.h"
#include "Field.h"

using namespace std;

Field::Field(int nsize, Grid *test, int btype)
{
  obj = test;
  BTYPE = btype;
  NSIZE = nsize;
  int BUFF_SIZE = (obj->L+2*NGUARD)*(obj->W+2*NGUARD)*(obj->H+2*NGUARD);
  ARR_SIZE = nsize*BUFF_SIZE; 
  arr = new double[ARR_SIZE];
  cout << "Grid object is being created" << endl;  
}

Field::~Field(void)
{
  cout << "Field object is being deleted" << endl;
}

// w refers to x or y or z components
// w is 0 for scalars
int Field::I(int w, int i, int j, int k)
{
  return (i*obj->W*obj->H*NSIZE + j*obj->H*NSIZE + k*NSIZE + w);
}

void Field::initialize()
{
  for(int i=0; i<obj->L; i++) {
    for(int j=0; j<obj->W; j++) {
      for(int k=0; k<obj->H; k++) {

        int ind = I(X,i,j,k);
        this->arr[ind] = 0.0; // for pressure

        if (NSIZE>1) // if velocity
        {
          double id = static_cast<double> (i-obj->L/2.0);
          double jd = static_cast<double> (j-obj->W/2.0);
          double kd = static_cast<double> (k-obj->H/2.0);

          this->arr[ind] = coef*exp(-eta*(obj->dx)*(obj->dx)
                             *((id*id)+(jd*jd)+(kd*kd)));
          this->arr[ind+1] = 0.0;
          this->arr[ind+2] = 0.0;
        }
      }
    }
  }
}

void update_guard()
{

}


/*Boundary operator of a two chain. This is the curl operator. Since the * operator between 2-chain and 1-chain don't do anything, we are really using this as *del.i NOTE: This method is for updating the interior cells only. Having a separate boundary update makes it easier to have different boundary implementations in the future. */
void Field::bd(double out[])
{
  // loop over the interior cells only
  for(int i=1; i<obj->L-1; i++) { 
    for(int j=1; j<obj->W-1; j++) { 
      for(int k=1; k<obj->H-1; k++) {
        
        int ind = I(X,i,j,k);
        out[ind]   = this->arr[I(ZX,i,j,k-1)] - this->arr[I(ZX,i,j,k+1)]
                   - this->arr[I(XY,i,j-1,k)] + this->arr[I(XY,i,j+1,k)];

        out[ind+1] = this->arr[I(XY,i-1,j,k)] - this->arr[I(XY,i+1,j,k)]
                   - this->arr[I(YZ,i,j,k-1)] + this->arr[I(YZ,i,j,k+1)];

        out[ind+2] = this->arr[I(YZ,i,j-1,k)] - this->arr[I(YZ,i,j+1,k)]
                   - this->arr[I(ZX,i-1,j,k)] + this->arr[I(ZX,i+1,j,k)];
      }
    }
  }
}

