#include "main.h"

Field::Field(int nsize, Grid *test)
{
  obj = test;
  L = obj->L;
  W = obj->W;
  H = obj->H;
  GL = obj->GL;
  GW = obj->GW;
  GH = obj->GH;

  NSIZE = nsize;
  int BUFF_SIZE = (L+2*NGUARD)*(W+2*NGUARD)*(H+2*NGUARD);
  ARR_SIZE = nsize*BUFF_SIZE; 
  arr = new double[ARR_SIZE]; // for pressure or velocity
  cout << "Field object is being created" << endl;  
}

Field::~Field(void)
{
  cout << "Field object is being deleted" << endl;
}

// w refers to x or y or z components
// w is 0 for scalars
int Field::I(int w, int i, int j, int k)
{
  return (i*W*H*NSIZE + j*H*NSIZE + k*NSIZE + w);
}

void Field::initialize(int ic, int* coords, int* numprocs_dir)
{
  switch(ic)
  {
    deafult:
    case TGV : 
      cout << "Initializing Taylor-Green vortex problem" << endl;
      taylor_green_vortex(coords, numprocs_dir);
      break;
  }
}

void Field::taylor_green_vortex(int* coords, int* numprocs_dir)
{
  // intialize in the interior cells only
  for(int i=NGUARD; i < L+NGUARD; i++) { 
    for(int j=NGUARD; j < W+NGUARD; j++) { 
      for(int k=NGUARD; k < H+NGUARD; k++) {

        int ind = I(X,i,j,k); 
        this->arr[ind] = 0.0; // for pressure

        if (NSIZE>1) // if velocity
        {
          // We need to subtract NGUARD here, since i=NGUARD is corresponding to x=0
          // This will make the xCoord, the left endpoint of the cell
          // We have to add the offset associated with the coords of the processor.
          double U0 = 1.0;
          double xCoord = (i - NGUARD)*obj->dx + coords[0]*(2.0*M_PI/numprocs_dir[0]);
          double yCoord = (j - NGUARD)*obj->dx + coords[1]*(2.0*M_PI/numprocs_dir[1]);
          double zCoord = (k - NGUARD)*obj->dx + coords[2]*(2.0*M_PI/numprocs_dir[2]);
          this->arr[ind] = U0 * cos(xCoord)*sin(yCoord)*sin(zCoord);
          this->arr[ind+1] = -U0/2.0 * sin(xCoord)*cos(yCoord)*sin(zCoord);
          this->arr[ind+2] = -U0/2.0 * sin(xCoord)*sin(yCoord)*cos(zCoord);
        }
      }
    }
  }

}

void Field::update_bdry(int bdry)
{
  switch(bdry) 
  {
    default:
    case PERIODIC :
       cout << "Periodic boundary condition" << endl;
       periodic_bdry();
       break;
  }
}

void Field::periodic_bdry()
{
  int w = 0; //FIXME: hardcoded for scalar array update
  // Update the x-bdry
  for (int i = 0; i < NGUARD; ++i) {
    for (int j = 0; j < GW; ++j) { 
      for (int k = 0; k < GH; ++k) {

        // FIXME: add weight to the constructor
        // lower bdry update
        int ind_guard = I(w, i, j, k);
        int ind_int = I(w,i+L,j,k);

        // this updates x, y and z components for vector arrays
        for (int ii=0; ii<NSIZE; ++ii)
          this->arr[ind_guard+ii] = this->arr[ind_int+ii];

        // upper bdry update
        int iu = GL - NGUARD + i;
        ind_guard = I(w, iu, j, k);
        ind_int = I(w,iu-L,j,k);
     
        // this updates x, y and z components for vector arrays
        for (int ii=0; ii<NSIZE; ++ii)
          this->arr[ind_guard+ii] = this->arr[ind_int+ii];        
  
      }
    }
  }

  // Update the y-bdry
  for (int i = 0; i < GL; ++i) {
    for (int j = 0; j < NGUARD; ++j) { 
      for (int k = 0; k < GH; ++k) {
       
        // lower bdry update
        int ind_guard = I(w, i, j, k);
        int ind_int = I(w,i,j+W,k);

        // this updates x, y and z components for vector arrays
        for (int ii=0; ii<NSIZE; ++ii)
          this->arr[ind_guard+ii] = this->arr[ind_int+ii];

        // upper bdry update
        int ju = GW - NGUARD + j;
        ind_guard = I(w, i, ju, k);
        ind_int = I(w,i,ju-W,k);

        // this updates x, y and z components for vector arrays
        for (int ii=0; ii<NSIZE; ++ii)
          this->arr[ind_guard+ii] = this->arr[ind_int+ii];

      }
    }
  }


  // Update the z-bdry
  for (int i = 0; i < 0; ++i) {
    for (int j = 0; j < GW; ++j) { 
      for (int k = 0; k < NGUARD; ++k) {

        // lower bdry update
        int ind_guard = I(w, i, j, k);
        int ind_int = I(w,i,j,k+H);

        // this updates x, y and z components for vector arrays
        for (int ii=0; ii<NSIZE; ++ii)
          this->arr[ind_guard+ii] = this->arr[ind_int+ii];

        // upper bdry update
        int ku = GH - NGUARD + k;
        ind_guard = I(w, i, j, ku);
        ind_int = I(w,i,j,ku-H);

        // this updates x, y and z components for vector arrays
        for (int ii=0; ii<NSIZE; ++ii)
          this->arr[ind_guard+ii] = this->arr[ind_int+ii];
      }
    }
  }

}

/*Boundary operator of a two chain. This is the curl operator. Since the * operator between 2-chain and 1-chain doesn't do anything, we are really using this as *del.i 
NOTE: This method is for updating the interior cells only. Having a separate boundary update makes it easier to have different boundary implementations in the future. */
void Field::bd(double out[])
{
  // loop over the interior cells only
  for(int i=NGUARD; i < L+NGUARD; i++) { 
    for(int j=NGUARD; j < W+NGUARD; j++) { 
      for(int k=NGUARD; k < H+NGUARD; k++) {
        
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

// Function for calculating divergence ( which is the boundary map of the 1-chain.) Outputs to a 0 chain of dimension L,W,H
void Field::bd10(double out[])
{
  for(int i=NGUARD; i<L+NGUARD; i++) {
    for(int j=NGUARD; j<W+NGUARD; j++) {
      for(int k=NGUARD; k<H+NGUARD; k++) {
 
        out[I(X,i,j,k)] = this->arr[I(X,i-1,j,k)] - this->arr[I(X,i+1,j,k)] 
                        + this->arr[I(Y,i,j-1,k)] - this->arr[I(Y,i,j+1,k)] 
                        + this->arr[I(Z,i,j,k-1)] - this->arr[I(Z,i,j,k+1)];
      }
    }
  }
}

//FIXME: this is a temporary replication of the above function for non-linear term
// Function for calculating divergence ( which is the boundary map of the 1-chain.) Outputs to a 0 chain of dimension L,W,H
void Field::bd10_NL(double* deg1ch, double out[])
{
  for(int i=NGUARD; i<L+NGUARD; i++) {
    for(int j=NGUARD; j<W+NGUARD; j++) {
      for(int k=NGUARD; k<H+NGUARD; k++) {

        out[I(X,i,j,k)] = deg1ch[I(X,i-1,j,k)] - deg1ch[I(X,i+1,j,k)]
                        + deg1ch[I(Y,i,j-1,k)] - deg1ch[I(Y,i,j+1,k)]
                        + deg1ch[I(Z,i,j,k-1)] - deg1ch[I(Z,i,j,k+1)];
      }
    }
  }
}


// Function for calculating gradient ( which is the coboundary map of the 0-chain.) Outputs to a 1 chain of dimension L,W,H
void Field::d01(double out[])
{
  for(int i=NGUARD; i < L+NGUARD; i++) {
    for(int j=NGUARD; j < W+NGUARD; j++) {
      for(int k=NGUARD; k < H+NGUARD; k++) {

        int ind = I(X,i,j,k);
        // PRAO: pay attention to the indices here to make sure that they stay consistent
        out[ind] = -this->arr[I(0,i-1,j,k)] + this->arr[I(0,i+1,j,k)];
        out[ind+1] = -this->arr[I(0,i,j-1,k)] + this->arr[I(0,i,j+1,k)];
        out[ind+2] = -this->arr[I(0,i,j,k-1)] + this->arr[I(0,i,j,k+1)];
      }
    }
  }
}

// Add Laplacian for 1-chain. nu is the viscosity. negative (-nu) is required in order to make the Laplacian term negative semidefinite.
//
void Field::laplacian(double out[])
{
  double factor = -nu/(4.0*obj->dx*obj->dx);
  int center,up,down,left,right,front,back;

  for(int i=NGUARD; i<L+NGUARD; i++) {
    for(int j=NGUARD; j<W+NGUARD; j++) {
      for(int k=NGUARD; k<H+NGUARD; k++) {

        center= I(0,i,j,k);
        front = I(0,i-2,j,k);
        back = I(0,i+2,j,k);
        left = I(0,i,j-2,k);
        right = I(0,i,j+2,k);
        down = I(0,i,j,k-2);
        up = I(0,i,j,k+2);

        //PRAO: should l's upper limit be obj->NSIZE?
        for(int l = 0; l < NSIZE; l++) {
          out[center+l] += factor*(   6.0*this->arr[center+l] 
                                    - this->arr[front+l] - this->arr[back+l] 
                                    - this->arr[right+l] - this->arr[left+l] 
                                    - this->arr[down+l] - this->arr[up+l]   );
        }
      }
    }
  }
}

