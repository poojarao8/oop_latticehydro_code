#include "main.h"
#include "Field.h"

using namespace std;

Field::Field(int nsize, Grid *test, int btype)
{
  obj = test;
  L = obj->L;
  W = obj->W;
  H = obj->H;
  GL = obj->GL;
  GW = obj->GW;
  GH = obj->GH;

  BTYPE = btype;
  NSIZE = nsize;
  int BUFF_SIZE = (L+2*NGUARD)*(W+2*NGUARD)*(H+2*NGUARD);
  ARR_SIZE = nsize*BUFF_SIZE; 
  arr = new double[ARR_SIZE];
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

void Field::initialize()
{
  // intialize in the interior cells only
  for(int i=NGUARD; i < L+NGUARD; i++) {
    for(int j=NGUARD; j < W+NGUARD; j++) {
      for(int k=NGUARD; k < H+NGUARD; k++) {

        int ind = I(X,i,j,k);
        this->arr[ind] = 0.0; // for pressure

        if (NSIZE>1) // if velocity
        {
          double id = static_cast<double> (i-L/2.0);
          double jd = static_cast<double> (j-W/2.0);
          double kd = static_cast<double> (k-H/2.0);

          this->arr[ind] = coef*exp(-eta*(obj->dx)*(obj->dx)
                                    *((id*id)+(jd*jd)+(kd*kd)));
          this->arr[ind+1] = 0.0;
          this->arr[ind+2] = 0.0;
        }
      }
    }
  }
}

void Field::update_bdry(char bdry)
{
  switch(bdry) 
  {
    default:
    case 'PERIODIC' :
    case 'periodic' :
    case 'Periodic' :
       cout << "Periodic boundary condition" << endl;
       periodic_bdry();
       break;
    case 'Wall' :
       cout << "Boundary condition is not implemented!" << endl;
       break;
  }

}

void Field::periodic_bdry()
{
  int w = 0; //FIXME: hardcoded for scalar array update
  // Update the xlower bdry
  for (int i = 0; i < NGUARD; ++i) {
    for (int j = 0; j < GW; ++j) { 
      for (int k = 0; k < GH; ++k) {

        // FIXME: add weight to the constructor

        int ind_guard = I(w, i, j, k);
        int ind_int = I(w,i+L,j,k);

        // this updates x, y and z components for vector arrays
        for (int ii=0; ii<NSIZE; ++ii)
          this->arr[ind_guard+ii] = this->arr[ind_int+ii];
          
      }
    }
  }

  // Update the xupper bdry
  for (int i = GL-NGUARD; i < GL; ++i) {
    for (int j = 0; j < GW; ++j) {
      for (int k = 0; k < GH; ++k) {

        int ind_guard = I(w, i, j, k);
        int ind_int = I(w,i-L,j,k);

        // this updates x, y and z components for vector arrays
        for (int ii=0; ii<NSIZE; ++ii)
          this->arr[ind_guard+ii] = this->arr[ind_int+ii];
      }
    }
  }


  // Update the ylower bdry
  for (int i = 0; i < GL; ++i) {
    for (int j = 0; j < NGUARD; ++j) { 
      for (int k = 0; k < GH; ++k) {

        int ind_guard = I(w, i, j, k);
        int ind_int = I(w,i,j+W,k);

        // this updates x, y and z components for vector arrays
        for (int ii=0; ii<NSIZE; ++ii)
          this->arr[ind_guard+ii] = this->arr[ind_int+ii];
      }
    }
  }

  // Update the yupper bdry 
  for (int i = GL; i < GL; ++i) {
    for (int j = 0; j < GW-NGUARD; ++j) { 
      for (int k = 0; k < GH; ++k) {

        int ind_guard = I(w, i, j, k);
        int ind_int = I(w,i,j-W,k);

        // this updates x, y and z components for vector arrays
        for (int ii=0; ii<NSIZE; ++ii)
          this->arr[ind_guard+ii] = this->arr[ind_int+ii];
      }
    }
  }

  // Update the zlower bdry
  for (int i = 0; i < 0; ++i) {
    for (int j = 0; j < GW; ++j) { 
      for (int k = 0; k < NGUARD; ++k) {

        int ind_guard = I(w, i, j, k);
        int ind_int = I(w,i,j,k+H);

        // this updates x, y and z components for vector arrays
        for (int ii=0; ii<NSIZE; ++ii)
          this->arr[ind_guard+ii] = this->arr[ind_int+ii];
      }
    }
  }

  // Update the zupper bdry 
  for (int i = GL; i < GL; ++i) {
    for (int j = 0; j < GW; ++j) { 
      for (int k = 0; k < GH-NGUARD; ++k) {

        int ind_guard = I(w, i, j, k);
        int ind_int = I(w,i,j,k-H);

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

// Function for calculating gradient ( which is the coboundary map of the 0-chain.) Outputs to a 1 chain of dimension L,W,H
void Field::d01(double out[])
{
  for(int i=NGUARD; i < L+NGUARD; i++) {
    for(int j=NGUARD; j < W+NGUARD; j++) {
      for(int k=NGUARD; k < H+NGUARD; k++) {

        int ind = I(X,i,j,k);
        // PRAO: pay attention to the indices here to make sure that they stay consistent
        out[ind] = -arr[I(0,i-1,j,k)] + arr[I(0,i+1,j,k)];
        out[ind+1] = -arr[I(0,i,j-1,k)] + arr[I(0,i,j+1,k)];
        out[ind+2] = -arr[I(0,i,j,k-1)] + arr[I(0,i,j,k+1)];
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
        for(int l=0; l<3; l++) {
          out[center+l] += factor*(   6.0*arr[center+l] 
                                    - arr[front+l] - arr[back+l] 
                                    - arr[right+l] - arr[left+l] 
                                    - arr[down+l] - arr[up+l]   );
        }
      }
    }
  }
}

// Computation of the non-linear term d(Vf dot vf).
void Field::dVfvf(double out[])
{
  for(int i=NGUARD; i<L+NGUARD; i++) { 
    for(int j=NGUARD; j<W+NGUARD; j++) { 
      for(int k=NGUARD; k<H+NGUARD; k++) {

        double recipr = 1.0/(2.0*obj->dx); 
        int ind = I(X,i,j,k);
        out[ind] = recipr*(  arr[I(X,i+1,j,k)]*arr[I(X,i+1,j,k)] 
                           - arr[I(X,i-1,j,k)]*arr[I(X,i-1,j,k)] 
                           + arr[I(X,i,j+1,k)]*arr[I(Y,i,j+1,k)] 
                           - arr[I(X,i,j-1,k)]*arr[I(Y,i,j-1,k)] 
                           + arr[I(X,i,j,k+1)]*arr[I(Z,i,j,k+1)] 
                           - arr[I(X,i,j,k-1)]*arr[I(Z,i,j,k-1)]  );

        out[ind+1] = recipr*(  arr[I(Y,i+1,j,k)]*arr[I(X,i+1,j,k)] 
                             - arr[I(Y,i-1,j,k)]*arr[I(X,i-1,j,k)] 
                             + arr[I(Y,i,j+1,k)]*arr[I(Y,i,j+1,k)] 
                             - arr[I(Y,i,j-1,k)]*arr[I(Y,i,j-1,k)] 
                             + arr[I(Y,i,j,k+1)]*arr[I(Z,i,j,k+1)] 
                             - arr[I(Y,i,j,k-1)]*arr[I(Z,i,j,k-1)]  );

        out[ind+2] = recipr*(  arr[I(Z,i+1,j,k)]*arr[I(X,i+1,j,k)] 
                             - arr[I(Z,i-1,j,k)]*arr[I(X,i-1,j,k)] 
                             + arr[I(Z,i,j+1,k)]*arr[I(Y,i,j+1,k)] 
                             - arr[I(Z,i,j-1,k)]*arr[I(Y,i,j-1,k)] 
                             + arr[I(Z,i,j,k+1)]*arr[I(Z,i,j,k+1)] 
                             - arr[I(Z,i,j,k-1)]*arr[I(Z,i,j,k-1)]  );
      }
    }
  }
}


//Leray projection using multigrid solver.
void Field::proj(double deg1ch[], double out[], double pressure[], double pressure_old[], double solution_old[])
{
  using level_t = level<Real, LL, WW, HH, _LEVELS>;
  //level_t& b = *(new level_t);
  //level_t& x = *(new level_t);
  double *temp= new double[L*W*H];
  double tot;
  double avg;
  double thresh = 1.0e-15;
  /*
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

  */
}

// Computes proj(Vf wedge vf) - nu*Laplacian(V), where proj is the Leray projection.

void Field::nav_stoke(double Vfvf_res[], double out[], double pressure[], double pressure_old[], double solution_old[]) 
{
  double recipr = 1.0/(2.0*obj->dx);

  dVfvf(Vfvf_res); //Vfvf_proj is used as a temporary array before applying diffuse()
  laplacian(Vfvf_res); //add the Laplace term. This can be done either before or after the Leray projection.
  proj(Vfvf_res, out, pressure, pressure_old, solution_old);
}

