// to compile, use gcc -ocfhydro lattice_hydro_coarse_fine_mult.cpp -std=c++17 -lm -lstdc++ -fopenmp
// Needs cxxopts.hpp library in the include subdirectory
// Also gcc/5.3.1 or above.
// On some machines, you may have to do "module load gcc/5.3.1" or similar. Both in compilation and also running the compiled code.
// Once compiled, the code should be run as "./diffused -h 0.01 -n 0.001 -e 0.125"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <omp.h>
#include "include/cxxopts.hpp"
#include "solver.h"

#define YZ 0
#define ZX 1
#define XY 2 //XY, YZ and ZX are the oriented planes (i.e. 2-symplex) at each lattice point.
#define X 0
#define Y 1
#define Z 2  //X, Y and Z are the oriented segments (i.e. 1-symplex) at each lattice point.
#define OUTNAME "Dfd_" //Output file for computed value. Python can read this by numpy package with  np.fromfile command.
#define OUTNAMEC "Coarse_"
#define L 64
#define W 64
#define H 64//Length, width and height of the fine grid. Coarse grid will be half.
#define MU 10.0

using namespace std;


const int LL=L/2;
const int WW=W/2;
const int HH=H/2;
const int ARR_SIZE=3*L*W*H;
const int COARSE_ARR_SIZE=3*LL*WW*HH;
const int OCTALSIZE=LL*WW*HH;
const double factor=4.0*OCTALSIZE;

using level_t = level<Real, LL, WW, HH, _LEVELS>;


cxxopts::ParseResult
parse(int argc, char* argv[])
{
  try
  {
    cxxopts::Options options(argv[0], " -testing cxxopts");
    options
    .positional_help("[optional args]")
    .show_positional_help();

    options
    .allow_unrecognised_options()
    .add_options()
    ("l, length", "length of one side of the domain.", cxxopts::value<double>(), "l")
    ("c, coeff", "coefficient of the initialization function.", cxxopts::value<double>(), "c")
    ("n, nu", "nu value (viscosity) to use", cxxopts::value<double>(), "nu")
    ("e, eta","Eta value (a parameter in the initial condition)", cxxopts::value<double>(), "eta")
    ("p, print","Print every p-th step.", cxxopts::value<int>(),"print")
    ("f, frequency","How often apply coarse-fine.", cxxopts::value<int>(),"freq")
    ("s, silent", "do not print calculations on screen")
    ;



    // options.parse_positional({"input", "output", "positional"});

    auto result = options.parse(argc, argv);

    // std::cout << "Arguments remain = " << argc << std::endl;

    return result;

  } catch (const cxxopts::OptionException& e)
  {
    std::cout << "error parsing options: " << e.what() << std::endl;
    exit(1);
  }
}

// Function that converts coordinates on 3-torus to the array index.
inline int I(int w, int x, int y, int z)
{
  return (x*W*H*3 + y*H*3 +z*3+ w);
}

inline int I(int w, int x, int y, int z, int length, int width, int height)
{
  return (x*width*height*3 + y*height*3 +z*3+ w);
}


inline int J(int x, int y, int z)
{
  return (x*W*H + y*H +z);
}


// Boundary operator of a two chain. This is the curl operator. Since the * operator between 2-chain and 1-chain don't do anything, we are really using this as *del.
// Rather than making the function return an array pointer, it stores the result in the out[] array. This was done to avoid having to make new arrays and having to deal with garbage collection.
void bd(double in_arr[], double out[])
{
  for(int i=0; i<L; i++) {
    for(int j=0; j<W; j++) {
      for(int k=0; k<H; k++) {
        out[I(X,i,j,k)]= in_arr[I(ZX,i,j,(k!=0) ? (k-1) : (H-1))]-in_arr[I(ZX,i,j,(k!=(H-1)) ? (k+1) : 0)]-in_arr[I(XY,i,(j!=0) ? (j-1) : (W-1),k)]+in_arr[I(XY,i,(j!=(W-1)) ? (j+1) : 0,k)];
        out[I(Y,i,j,k)]= in_arr[I(XY,(i!=0) ? (i-1) : (L-1),j,k)]-in_arr[I(XY,(i!=(L-1)) ? (i+1) : 0,j,k)]-in_arr[I(YZ,i,j,(k!=0) ? (k-1) : (H-1))]+in_arr[I(YZ,i,j,(k!=(H-1)) ? (k+1) : 0)];
        out[I(Z,i,j,k)]= in_arr[I(YZ,i,(j!=0) ? (j-1) : (W-1),k)]-in_arr[I(YZ,i,(j!=(W-1)) ? (j+1) : 0,k)]-in_arr[I(ZX,(i!=0) ? (i-1) : (L-1),j,k)]+in_arr[I(ZX,(i!=(L-1)) ? (i+1) : 0,j,k)];
      }
    }
  }
}

void bdC(double in_arr[], double out[]) //boundary operator for coarse
{
  for(int i=0; i<LL; i++) {
    for(int j=0; j<WW; j++) {
      for(int k=0; k<HH; k++) {
        out[I(X,i,j,k,LL,WW,HH)]= in_arr[I(ZX,i,j,(k!=0) ? (k-1) : (HH-1),LL,WW,HH)]-in_arr[I(ZX,i,j,(k!=(HH-1)) ? (k+1) : 0,LL,WW,HH)]-in_arr[I(XY,i,(j!=0) ? (j-1) : (WW-1),k,LL,WW,HH)]+in_arr[I(XY,i,(j!=(WW-1)) ? (j+1) : 0,k,LL,WW,HH)];
        out[I(Y,i,j,k,LL,WW,HH)]= in_arr[I(XY,(i!=0) ? (i-1) : (LL-1),j,k,LL,WW,HH)]-in_arr[I(XY,(i!=(LL-1)) ? (i+1) : 0,j,k,LL,WW,HH)]-in_arr[I(YZ,i,j,(k!=0) ? (k-1) : (HH-1),LL,WW,HH)]+in_arr[I(YZ,i,j,(k!=(HH-1)) ? (k+1) : 0,LL,WW,HH)];
        out[I(Z,i,j,k,LL,WW,HH)]= in_arr[I(YZ,i,(j!=0) ? (j-1) : (WW-1),k,LL,WW,HH)]-in_arr[I(YZ,i,(j!=(WW-1)) ? (j+1) : 0,k,LL,WW,HH)]-in_arr[I(ZX,(i!=0) ? (i-1) : (LL-1),j,k,LL,WW,HH)]+in_arr[I(ZX,(i!=(LL-1)) ? (i+1) : 0,j,k,LL,WW,HH)];
      }
    }
  }
}


//Make a 1 chain to have mean value 0
void meanzero(double arr[]){
  double meanval,tot;
  for (int l=0; l<3; l++) {

    for(int b1=0; l<2; l++) {for(int b2=0; l<2; l++) {for(int b3=0; l<2; l++) {
      tot=0.0;
      for(int i=0; i<LL; i++) for(int j=0; j<WW; j++) for(int k=0; k<HH; k++) tot+=arr[I(l,i*2+b1,j*2+b2,k*2+b3)];
      meanval = tot/OCTALSIZE;
      for(int i=0; i<LL; i++) for(int j=0; j<WW; j++) for(int k=0; k<HH; k++) arr[I(l,i*2+b1,j*2+b2,k*2+b3)] -= meanval;
    }}}
  }
}




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


// Function for calculating divergence ( which is the boundary map of the 1-chain.) Outputs to a 0 chain of dimension L,W,H
void bd10(double arr[], double out[]){
  for(int i=0; i<L; i++) {
    for(int j=0; j<W; j++) {
      for(int k=0; k<H; k++) {
        out[J(i,j,k)] = arr[I(X,(i!=0) ? (i-1) : (L-1),j,k)] -arr[I(X,(i!=(L-1)) ? (i+1) : 0,j,k)]+arr[I(Y,i,(j!=0) ? (j-1) : (W-1),k)]-arr[I(Y,i,(j!=(W-1)) ? (j+1) : 0,k)]+ arr[I(Z,i,j,(k!=0) ? (k-1) : (H-1))] -arr[I(Z,i,j,(k!=(H-1)) ? (k+1) : 0)];
      }
    }
  }
}

void bd10C(double arr[], double out[]){
  for(int i=0; i<LL; i++) {
    for(int j=0; j<WW; j++) {
      for(int k=0; k<HH; k++) {
        out[J(i,j,k)] = arr[I(X,(i!=0) ? (i-1) : (LL-1),j,k)] -arr[I(X,(i!=(LL-1)) ? (i+1) : 0,j,k)]+arr[I(Y,i,(j!=0) ? (j-1) : (WW-1),k)]-arr[I(Y,i,(j!=(WW-1)) ? (j+1) : 0,k)]+ arr[I(Z,i,j,(k!=0) ? (k-1) : (HH-1))] -arr[I(Z,i,j,(k!=(HH-1)) ? (k+1) : 0)];
      }
    }
  }
}


// Function for calculating gradient ( which is the coboundary map of the 0-chain.) Outputs to a 1 chain of dimension L,W,H
void d01(double arr[], double out[])
{
  for(int i=0; i<L; i++) {
    for(int j=0; j<W; j++) {
      for(int k=0; k<H; k++) {
        out[I(X,i,j,k)] = -arr[J((i!=0) ? (i-1) : (L-1),j,k)] + arr[J((i!=(L-1)) ? (i+1) : 0,j,k)];
        out[I(Y,i,j,k)] = -arr[J(i,(j!=0) ? (j-1) : (W-1),k)] + arr[J(i,(j!=(W-1)) ? (j+1) : 0,k)];
        out[I(Z,i,j,k)] = -arr[J(i,j,(k!=0) ? (k-1) : (H-1))] + arr[J(i,j,(k!=(H-1)) ? (k+1) : 0)];
      }
    }
  }
}

// Add Laplacian for 1-chain. nu is the viscosity. negative (-nu) is required in order to make the Laplacian term negative semidefinite.
void addLaplacian(double arr[], double out[], double nu, double h)
{
  double factor = -nu/(4.0*h*h);
  int center,up,down,left,right,front,back;
  for(int i=0; i<L; i++) {
    for(int j=0; j<W; j++) {
      for(int k=0; k<H; k++) {
        center= I(0,i,j,k);
        front= I(0, (i-2+L)%L,j,k);
        back= I(0, (i+2)%L,j,k);
        left= I(0,i, (j-2+W)%W,k);
        right= I(0,i, (j+2)%W,k);
        up= I(0,i,j, (k+2)%H);
        down= I(0,i,j, (k-2+H)%H);
        for(int l=0; l<3; l++) {
          out[center+l]+=factor*(6.0*arr[center+l] -arr[front+l] -arr[back+l] -arr[right+l]  -arr[left+l] -arr[down+l]  -arr[up+l]);
        }
      }
    }
  }
}





// Splits a 0 chain into 8 parts.
void split8(double arr[], double **out) {

  int octal;
  for(int b1=0; b1<2; b1++) {for(int b2=0; b2<2; b2++) {for(int b3=0; b3<2; b3++) {
    for(int i=0; i<LL; i++) {
      for(int j=0; j<WW; j++) {
        for(int k=0; k<HH; k++) {
          octal=b1*4+b2*2+b3;
          out[octal][i*WW*HH+j*HH+k]=arr[(i*2+b1)*W*H+(j*2+b2)*H+k*2+b3];
        }
      }
    }
  }}}
}

//Function to undo split8
void unsplit8(double **arr, double out[]) {

  int octal;
  for(int b1=0; b1<2; b1++) {for(int b2=0; b2<2; b2++) {for(int b3=0; b3<2; b3++) {
    for(int i=0; i<LL; i++) {
      for(int j=0; j<WW; j++) {
        for(int k=0; k<HH; k++) {
          octal=b1*4+b2*2+b3;
          out[(i*2+b1)*W*H+(j*2+b2)*H+k*2+b3]=arr[octal][i*WW*HH+j*HH+k];
        }
      }
    }
  }}}
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

/*
if (abs(tot)>1.0e-17) {
  cout<<"total "<<tot<<" is too large. Subtracting avg."<<endl;
  avg= tot/(LL*WW*HH);
  tot=0.0;
  for(int i=0; i<LL; i++) {
    for(int j=0; j<WW; j++) {
      for(int k=0; k<HH; k++) {
        
        b.dat.at(i,j,k)=pressure[(i*2+b1)*W*H+(j*2+b2)*H+k*2+b3]-avg;
        tot += b.dat.at(i,j,k);
        
      }
    }
  }
}
*/

if (abs(tot)>1.0e-17) {thresh= abs(tot)*10;}
else  {thresh= 1.0e-15;}

cout << "octal: " << octal <<endl;
cout << "total: " <<tot <<endl;
  K=Slash(b, x, 15, 500, thresh);
  cout<<"Multigrid solver applied "<<K<<" times"<<endl;
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

int read_file(double arr[], int timestamp)
  {
    char fname[50];
    sprintf(fname,"%st%d.raw",OUTNAME,timestamp);
    ifstream in(fname, ios::out | ios::binary);
    if(!in) {
      cout << "Cannot open file." << endl;
      return 0;
    }

    for(int i=0; i<ARR_SIZE; i++) in.read((char *) &arr[i], sizeof(double));

    in.close();
    cout << "File "<< fname << " read" << endl;
    return 1;
  }


int write2file(double outdata[], int outdatasize, double h, double nu, double eta)
  {
    char fname[50];
    sprintf(fname,"%sN%uh%.3fnu%.3feta%.3f.raw",OUTNAMEC,L,h,nu,eta);
    ofstream out(fname, ios::out | ios::binary);
    if(!out) {
      cout << "Cannot open file." << endl;
      return 0;
    }

    for(int i=0; i<outdatasize; i++) {
      for(int j=0; j<16; j++) {
        out.write((char *) &outdata[i*16+j], sizeof(double));
      }
    }

    out.close();
    cout << "File "<< fname << " created" << endl;
    return 1;
  }

  int write2fileC(double outdata[], int outdatasize, double h, double nu, double eta)
  {
    char fname[50];
    sprintf(fname,"%sN%uh%.3fnu%.3feta%.3f.raw",OUTNAME,L,h,nu,eta);
    ofstream out(fname, ios::out | ios::binary);
    if(!out) {
      cout << "Cannot open file." << endl;
      return 0;
    }

    for(int i=0; i<outdatasize; i++) {
      for(int j=0; j<16; j++) {
        out.write((char *) &outdata[i*16+j], sizeof(double));
      }
    }

    out.close();
    cout << "File "<< fname << " created" << endl;
    return 1;
  }

// Computes proj(Vf wedge vf) - nu*Laplacian(V), where proj is the Leray projection.
  inline void nav_stoke (double h,double nu,double V[], double Vfvf_res[], double out[], double pressure[], double pressure_old[], double solution_old[]) {
      double recipr = 1.0/(2.0*h);

      dVfvf(recipr,V, Vfvf_res); //Vfvf_proj is used as a temporary array before applying diffuse()
      addLaplacian(V,Vfvf_res,nu,h); //add the Laplace term. This can be done either before or after the Leray projection.
      proj(Vfvf_res, out, pressure, pressure_old, solution_old);
    }

    // L2 norm of a vector
    double norm(double V[]) {
      double tot=0.0;
      for (int i=0; i<ARR_SIZE; i++) tot+= V[i]*V[i];
      return sqrt(tot);
    }

    // maxnorm of a vector
    double maxnorm(double V[]) {
      double maxval=0.0;
      double normsqval;
      for (int i=0; i<(ARR_SIZE-2); i+=3) {
        normsqval = V[i]*V[i]+V[i+1]*V[i+1]+V[i+2]*V[i+2];
        if(normsqval>maxval) maxval=normsqval;
      }
      return sqrt(maxval);
    }

    // maxnorm of a vector
    double maxnormC(double V[]) {
      double maxval=0.0;
      double normsqval;
      for (int i=0; i<(COARSE_ARR_SIZE-2); i+=3) {
        normsqval = V[i]*V[i]+V[i+1]*V[i+1]+V[i+2]*V[i+2];
        if(normsqval>maxval) maxval=normsqval;
      }
      return sqrt(maxval);
    }


    // L2 norm of the gradient
    void gradient_norm(double V[],double *gradnorm, double *gradsymm, double *gradskew ) {
      double Dxx,Dxy,Dxz,Dyx,Dyy,Dyz,Dzx,Dzy,Dzz;
      *gradsymm = 0.0;
      *gradskew = 0.0;
      for(int i=0; i<L; i++) {
        for(int j=0; j<W; j++) {
          for(int k=0; k<H; k++) {
            Dxx=V[I(X,(i!=(L-1)) ? (i+1) : 0,j,k)] - V[I(X,(i!=0) ? (i-1) : (L-1),j,k)];
            Dxy=V[I(X,i,(j!=(W-1)) ? (j+1) : 0,k)] - V[I(X,i,(j!=0) ? (j-1) : (W-1),k)];
            Dxz=V[I(X,i,j,(k!=(H-1)) ? (k+1) : 0)] - V[I(X,i,j,(k!=0) ? (k-1) : (H-1))];
            Dyx=V[I(Y,(i!=(L-1)) ? (i+1) : 0,j,k)] - V[I(Y,(i!=0) ? (i-1) : (L-1),j,k)];
            Dyy=V[I(Y,i,(j!=(W-1)) ? (j+1) : 0,k)] - V[I(Y,i,(j!=0) ? (j-1) : (W-1),k)];
            Dyz=V[I(Y,i,j,(k!=(H-1)) ? (k+1) : 0)] - V[I(Y,i,j,(k!=0) ? (k-1) : (H-1))];
            Dzx=V[I(Z,(i!=(L-1)) ? (i+1) : 0,j,k)] - V[I(Z,(i!=0) ? (i-1) : (L-1),j,k)];
            Dzy=V[I(Z,i,(j!=(W-1)) ? (j+1) : 0,k)] - V[I(Z,i,(j!=0) ? (j-1) : (W-1),k)];
            Dzz=V[I(Z,i,j,(k!=(H-1)) ? (k+1) : 0)] - V[I(Z,i,j,(k!=0) ? (k-1) : (H-1))];
            *gradsymm += Dxx*Dxx + Dyy*Dyy + Dzz*Dzz;
            *gradskew += Dxy*Dxy+Dxz*Dxz+Dyx*Dyx+Dyz*Dyz+Dzx*Dzx+Dzy*Dzy;
          }
        }
      }
      *gradnorm = sqrt(*gradsymm+*gradskew);
      *gradsymm = sqrt(*gradsymm);
      *gradskew = sqrt(*gradskew);

    }

    // L2 norm of the gradient for the coarse grid
    void gradient_normC(double V[],double *gradnorm, double *gradsymm, double *gradskew ) {
      double Dxx,Dxy,Dxz,Dyx,Dyy,Dyz,Dzx,Dzy,Dzz;
      *gradsymm = 0.0;
      *gradskew = 0.0;
      for(int i=0; i<LL; i++) {
        for(int j=0; j<WW; j++) {
          for(int k=0; k<HH; k++) {
            Dxx=V[I(X,(i!=(L-1)) ? (i+1) : 0,j,k)] - V[I(X,(i!=0) ? (i-1) : (L-1),j,k)];
            Dxy=V[I(X,i,(j!=(W-1)) ? (j+1) : 0,k)] - V[I(X,i,(j!=0) ? (j-1) : (W-1),k)];
            Dxz=V[I(X,i,j,(k!=(H-1)) ? (k+1) : 0)] - V[I(X,i,j,(k!=0) ? (k-1) : (H-1))];
            Dyx=V[I(Y,(i!=(L-1)) ? (i+1) : 0,j,k)] - V[I(Y,(i!=0) ? (i-1) : (L-1),j,k)];
            Dyy=V[I(Y,i,(j!=(W-1)) ? (j+1) : 0,k)] - V[I(Y,i,(j!=0) ? (j-1) : (W-1),k)];
            Dyz=V[I(Y,i,j,(k!=(H-1)) ? (k+1) : 0)] - V[I(Y,i,j,(k!=0) ? (k-1) : (H-1))];
            Dzx=V[I(Z,(i!=(L-1)) ? (i+1) : 0,j,k)] - V[I(Z,(i!=0) ? (i-1) : (L-1),j,k)];
            Dzy=V[I(Z,i,(j!=(W-1)) ? (j+1) : 0,k)] - V[I(Z,i,(j!=0) ? (j-1) : (W-1),k)];
            Dzz=V[I(Z,i,j,(k!=(H-1)) ? (k+1) : 0)] - V[I(Z,i,j,(k!=0) ? (k-1) : (H-1))];
            *gradsymm += Dxx*Dxx + Dyy*Dyy + Dzz*Dzz;
            *gradskew += Dxy*Dxy+Dxz*Dxz+Dyx*Dyx+Dyz*Dyz+Dzx*Dzx+Dzy*Dzy;
          }
        }
      }
      *gradnorm = sqrt(*gradsymm+*gradskew);
      *gradsymm = sqrt(*gradsymm);
      *gradskew = sqrt(*gradskew);

    }

    // Percentage of vectors whose norm is above given thresholds.
    void pct_vectors(double V[], double threshold1, double threshold2, double *pct1, double *pct2) {
      unsigned int count1=0;
      unsigned int count2=0;
      double normsqval;
      double th1=threshold1*threshold1;
      double th2=threshold2*threshold2;
      for (int i=0; i<(ARR_SIZE-2); i+=3) {
        normsqval = V[i]*V[i]+V[i+1]*V[i+1]+V[i+2]*V[i+2];
        if(normsqval>th1) {
          count1++;
          if(normsqval>th2) {
            count2++;
          }
        }
      }
      *pct1 = (static_cast<double> (count1))/(static_cast<double> (L*W*H));
      *pct2 = (static_cast<double> (count2))/(static_cast<double> (L*W*H));
    }


    // maxnorm of a scalar field
    double maxnormscal(double scal[]) {
      double maxval=0.0;
      double absval;
      for (int i=0; i<(L*W*H); i++) {
        absval = abs(scal[i]);
        if(absval>maxval) maxval=absval;
      }
      return maxval;
    }

void fine2coarse(double arr[], double out[], int fL, int fW, int fH)
{
int cL= fL/2; //fL,fW,fH are length, width and height of the coarse grid.
int cW= fW/2; //cL,cW,cH are length, width and height of the fine grid.
int cH= fH/2;

      for(int l=0; l< cL/2; l++) {
         for(int m=0; m< cW/2; m++) {
            for(int n=0; n< cH/2; n++){
out[I(0,2*l,2*m,2*n,cL,cW,cH)]=(arr[I(0,4*l,((4*m-2)>=0)?(4*m-2):(4*m-2+fW),((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+2*arr[I(0,4*l,((4*m-2)>=0)?(4*m-2):(4*m-2+fW),4*n,fL,fW,fH)]+arr[I(0,4*l,((4*m-2)>=0)?(4*m-2):(4*m-2+fW),4*n+2,fL,fW,fH)]+2*arr[I(0,4*l,4*m,((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+4*arr[I(0,4*l,4*m,4*n,fL,fW,fH)]+2*arr[I(0,4*l,4*m,4*n+2,fL,fW,fH)]+arr[I(0,4*l,4*m+2,((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+2*arr[I(0,4*l,4*m+2,4*n,fL,fW,fH)]+arr[I(0,4*l,4*m+2,4*n+2,fL,fW,fH)])/4;
out[I(0,2*l,2*m,2*n+1,cL,cW,cH)]=(arr[I(0,4*l,((4*m-2)>=0)?(4*m-2):(4*m-2+fW),4*n+1,fL,fW,fH)]+arr[I(0,4*l,((4*m-2)>=0)?(4*m-2):(4*m-2+fW),4*n+3,fL,fW,fH)]+2*arr[I(0,4*l,4*m,4*n+1,fL,fW,fH)]+2*arr[I(0,4*l,4*m,4*n+3,fL,fW,fH)]+arr[I(0,4*l,4*m+2,4*n+1,fL,fW,fH)]+arr[I(0,4*l,4*m+2,4*n+3,fL,fW,fH)])/2;
out[I(0,2*l,2*m+1,2*n,cL,cW,cH)]=(arr[I(0,4*l,4*m+1,((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+2*arr[I(0,4*l,4*m+1,4*n,fL,fW,fH)]+arr[I(0,4*l,4*m+1,4*n+2,fL,fW,fH)]+arr[I(0,4*l,4*m+3,((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+2*arr[I(0,4*l,4*m+3,4*n,fL,fW,fH)]+arr[I(0,4*l,4*m+3,4*n+2,fL,fW,fH)])/2;
out[I(0,2*l,2*m+1,2*n+1,cL,cW,cH)]=arr[I(0,4*l,4*m+1,4*n+1,fL,fW,fH)]+arr[I(0,4*l,4*m+1,4*n+3,fL,fW,fH)]+arr[I(0,4*l,4*m+3,4*n+1,fL,fW,fH)]+arr[I(0,4*l,4*m+3,4*n+3,fL,fW,fH)];
out[I(0,2*l+1,2*m,2*n,cL,cW,cH)]=(arr[I(0,4*l+1,((4*m-2)>=0)?(4*m-2):(4*m-2+fW),((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+2*arr[I(0,4*l+1,((4*m-2)>=0)?(4*m-2):(4*m-2+fW),4*n,fL,fW,fH)]+arr[I(0,4*l+1,((4*m-2)>=0)?(4*m-2):(4*m-2+fW),4*n+2,fL,fW,fH)]+2*arr[I(0,4*l+1,4*m,((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+4*arr[I(0,4*l+1,4*m,4*n,fL,fW,fH)]+2*arr[I(0,4*l+1,4*m,4*n+2,fL,fW,fH)]+arr[I(0,4*l+1,4*m+2,((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+2*arr[I(0,4*l+1,4*m+2,4*n,fL,fW,fH)]+arr[I(0,4*l+1,4*m+2,4*n+2,fL,fW,fH)]+arr[I(0,4*l+3,((4*m-2)>=0)?(4*m-2):(4*m-2+fW),((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+2*arr[I(0,4*l+3,((4*m-2)>=0)?(4*m-2):(4*m-2+fW),4*n,fL,fW,fH)]+arr[I(0,4*l+3,((4*m-2)>=0)?(4*m-2):(4*m-2+fW),4*n+2,fL,fW,fH)]+2*arr[I(0,4*l+3,4*m,((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+4*arr[I(0,4*l+3,4*m,4*n,fL,fW,fH)]+2*arr[I(0,4*l+3,4*m,4*n+2,fL,fW,fH)]+arr[I(0,4*l+3,4*m+2,((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+2*arr[I(0,4*l+3,4*m+2,4*n,fL,fW,fH)]+arr[I(0,4*l+3,4*m+2,4*n+2,fL,fW,fH)])/8;
out[I(0,2*l+1,2*m,2*n+1,cL,cW,cH)]=(arr[I(0,4*l+1,((4*m-2)>=0)?(4*m-2):(4*m-2+fW),4*n+1,fL,fW,fH)]+arr[I(0,4*l+1,((4*m-2)>=0)?(4*m-2):(4*m-2+fW),4*n+3,fL,fW,fH)]+2*arr[I(0,4*l+1,4*m,4*n+1,fL,fW,fH)]+2*arr[I(0,4*l+1,4*m,4*n+3,fL,fW,fH)]+arr[I(0,4*l+1,4*m+2,4*n+1,fL,fW,fH)]+arr[I(0,4*l+1,4*m+2,4*n+3,fL,fW,fH)]+arr[I(0,4*l+3,((4*m-2)>=0)?(4*m-2):(4*m-2+fW),4*n+1,fL,fW,fH)]+arr[I(0,4*l+3,((4*m-2)>=0)?(4*m-2):(4*m-2+fW),4*n+3,fL,fW,fH)]+2*arr[I(0,4*l+3,4*m,4*n+1,fL,fW,fH)]+2*arr[I(0,4*l+3,4*m,4*n+3,fL,fW,fH)]+arr[I(0,4*l+3,4*m+2,4*n+1,fL,fW,fH)]+arr[I(0,4*l+3,4*m+2,4*n+3,fL,fW,fH)])/4;
out[I(0,2*l+1,2*m+1,2*n,cL,cW,cH)]=(arr[I(0,4*l+1,4*m+1,((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+2*arr[I(0,4*l+1,4*m+1,4*n,fL,fW,fH)]+arr[I(0,4*l+1,4*m+1,4*n+2,fL,fW,fH)]+arr[I(0,4*l+1,4*m+3,((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+2*arr[I(0,4*l+1,4*m+3,4*n,fL,fW,fH)]+arr[I(0,4*l+1,4*m+3,4*n+2,fL,fW,fH)]+arr[I(0,4*l+3,4*m+1,((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+2*arr[I(0,4*l+3,4*m+1,4*n,fL,fW,fH)]+arr[I(0,4*l+3,4*m+1,4*n+2,fL,fW,fH)]+arr[I(0,4*l+3,4*m+3,((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+2*arr[I(0,4*l+3,4*m+3,4*n,fL,fW,fH)]+arr[I(0,4*l+3,4*m+3,4*n+2,fL,fW,fH)])/4;
out[I(0,2*l+1,2*m+1,2*n+1,cL,cW,cH)]=(arr[I(0,4*l+1,4*m+1,4*n+1,fL,fW,fH)]+arr[I(0,4*l+1,4*m+1,4*n+3,fL,fW,fH)]+arr[I(0,4*l+1,4*m+3,4*n+1,fL,fW,fH)]+arr[I(0,4*l+1,4*m+3,4*n+3,fL,fW,fH)]+arr[I(0,4*l+3,4*m+1,4*n+1,fL,fW,fH)]+arr[I(0,4*l+3,4*m+1,4*n+3,fL,fW,fH)]+arr[I(0,4*l+3,4*m+3,4*n+1,fL,fW,fH)]+arr[I(0,4*l+3,4*m+3,4*n+3,fL,fW,fH)])/2;
out[I(1,2*l,2*m,2*n,cL,cW,cH)]=(arr[I(1,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),4*m,((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+2*arr[I(1,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),4*m,4*n,fL,fW,fH)]+arr[I(1,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),4*m,4*n+2,fL,fW,fH)]+2*arr[I(1,4*l,4*m,((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+4*arr[I(1,4*l,4*m,4*n,fL,fW,fH)]+2*arr[I(1,4*l,4*m,4*n+2,fL,fW,fH)]+arr[I(1,4*l+2,4*m,((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+2*arr[I(1,4*l+2,4*m,4*n,fL,fW,fH)]+arr[I(1,4*l+2,4*m,4*n+2,fL,fW,fH)])/4;
out[I(1,2*l,2*m,2*n+1,cL,cW,cH)]=(arr[I(1,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),4*m,4*n+1,fL,fW,fH)]+arr[I(1,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),4*m,4*n+3,fL,fW,fH)]+2*arr[I(1,4*l,4*m,4*n+1,fL,fW,fH)]+2*arr[I(1,4*l,4*m,4*n+3,fL,fW,fH)]+arr[I(1,4*l+2,4*m,4*n+1,fL,fW,fH)]+arr[I(1,4*l+2,4*m,4*n+3,fL,fW,fH)])/2;
out[I(1,2*l,2*m+1,2*n,cL,cW,cH)]=(arr[I(1,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),4*m+1,((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+2*arr[I(1,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),4*m+1,4*n,fL,fW,fH)]+arr[I(1,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),4*m+1,4*n+2,fL,fW,fH)]+arr[I(1,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),4*m+3,((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+2*arr[I(1,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),4*m+3,4*n,fL,fW,fH)]+arr[I(1,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),4*m+3,4*n+2,fL,fW,fH)]+2*arr[I(1,4*l,4*m+1,((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+4*arr[I(1,4*l,4*m+1,4*n,fL,fW,fH)]+2*arr[I(1,4*l,4*m+1,4*n+2,fL,fW,fH)]+2*arr[I(1,4*l,4*m+3,((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+4*arr[I(1,4*l,4*m+3,4*n,fL,fW,fH)]+2*arr[I(1,4*l,4*m+3,4*n+2,fL,fW,fH)]+arr[I(1,4*l+2,4*m+1,((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+2*arr[I(1,4*l+2,4*m+1,4*n,fL,fW,fH)]+arr[I(1,4*l+2,4*m+1,4*n+2,fL,fW,fH)]+arr[I(1,4*l+2,4*m+3,((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+2*arr[I(1,4*l+2,4*m+3,4*n,fL,fW,fH)]+arr[I(1,4*l+2,4*m+3,4*n+2,fL,fW,fH)])/8;
out[I(1,2*l,2*m+1,2*n+1,cL,cW,cH)]=(arr[I(1,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),4*m+1,4*n+1,fL,fW,fH)]+arr[I(1,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),4*m+1,4*n+3,fL,fW,fH)]+arr[I(1,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),4*m+3,4*n+1,fL,fW,fH)]+arr[I(1,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),4*m+3,4*n+3,fL,fW,fH)]+2*arr[I(1,4*l,4*m+1,4*n+1,fL,fW,fH)]+2*arr[I(1,4*l,4*m+1,4*n+3,fL,fW,fH)]+2*arr[I(1,4*l,4*m+3,4*n+1,fL,fW,fH)]+2*arr[I(1,4*l,4*m+3,4*n+3,fL,fW,fH)]+arr[I(1,4*l+2,4*m+1,4*n+1,fL,fW,fH)]+arr[I(1,4*l+2,4*m+1,4*n+3,fL,fW,fH)]+arr[I(1,4*l+2,4*m+3,4*n+1,fL,fW,fH)]+arr[I(1,4*l+2,4*m+3,4*n+3,fL,fW,fH)])/4;
out[I(1,2*l+1,2*m,2*n,cL,cW,cH)]=(arr[I(1,4*l+1,4*m,((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+2*arr[I(1,4*l+1,4*m,4*n,fL,fW,fH)]+arr[I(1,4*l+1,4*m,4*n+2,fL,fW,fH)]+arr[I(1,4*l+3,4*m,((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+2*arr[I(1,4*l+3,4*m,4*n,fL,fW,fH)]+arr[I(1,4*l+3,4*m,4*n+2,fL,fW,fH)])/2;
out[I(1,2*l+1,2*m,2*n+1,cL,cW,cH)]=arr[I(1,4*l+1,4*m,4*n+1,fL,fW,fH)]+arr[I(1,4*l+1,4*m,4*n+3,fL,fW,fH)]+arr[I(1,4*l+3,4*m,4*n+1,fL,fW,fH)]+arr[I(1,4*l+3,4*m,4*n+3,fL,fW,fH)];
out[I(1,2*l+1,2*m+1,2*n,cL,cW,cH)]=(arr[I(1,4*l+1,4*m+1,((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+2*arr[I(1,4*l+1,4*m+1,4*n,fL,fW,fH)]+arr[I(1,4*l+1,4*m+1,4*n+2,fL,fW,fH)]+arr[I(1,4*l+1,4*m+3,((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+2*arr[I(1,4*l+1,4*m+3,4*n,fL,fW,fH)]+arr[I(1,4*l+1,4*m+3,4*n+2,fL,fW,fH)]+arr[I(1,4*l+3,4*m+1,((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+2*arr[I(1,4*l+3,4*m+1,4*n,fL,fW,fH)]+arr[I(1,4*l+3,4*m+1,4*n+2,fL,fW,fH)]+arr[I(1,4*l+3,4*m+3,((4*n-2)>=0)?(4*n-2):(4*n-2+fH),fL,fW,fH)]+2*arr[I(1,4*l+3,4*m+3,4*n,fL,fW,fH)]+arr[I(1,4*l+3,4*m+3,4*n+2,fL,fW,fH)])/4;
out[I(1,2*l+1,2*m+1,2*n+1,cL,cW,cH)]=(arr[I(1,4*l+1,4*m+1,4*n+1,fL,fW,fH)]+arr[I(1,4*l+1,4*m+1,4*n+3,fL,fW,fH)]+arr[I(1,4*l+1,4*m+3,4*n+1,fL,fW,fH)]+arr[I(1,4*l+1,4*m+3,4*n+3,fL,fW,fH)]+arr[I(1,4*l+3,4*m+1,4*n+1,fL,fW,fH)]+arr[I(1,4*l+3,4*m+1,4*n+3,fL,fW,fH)]+arr[I(1,4*l+3,4*m+3,4*n+1,fL,fW,fH)]+arr[I(1,4*l+3,4*m+3,4*n+3,fL,fW,fH)])/2;
out[I(2,2*l,2*m,2*n,cL,cW,cH)]=(arr[I(2,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),((4*m-2)>=0)?(4*m-2):(4*m-2+fW),4*n,fL,fW,fH)]+2*arr[I(2,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),4*m,4*n,fL,fW,fH)]+arr[I(2,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),4*m+2,4*n,fL,fW,fH)]+2*arr[I(2,4*l,((4*m-2)>=0)?(4*m-2):(4*m-2+fW),4*n,fL,fW,fH)]+4*arr[I(2,4*l,4*m,4*n,fL,fW,fH)]+2*arr[I(2,4*l,4*m+2,4*n,fL,fW,fH)]+arr[I(2,4*l+2,((4*m-2)>=0)?(4*m-2):(4*m-2+fW),4*n,fL,fW,fH)]+2*arr[I(2,4*l+2,4*m,4*n,fL,fW,fH)]+arr[I(2,4*l+2,4*m+2,4*n,fL,fW,fH)])/4;
out[I(2,2*l,2*m,2*n+1,cL,cW,cH)]=(arr[I(2,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),((4*m-2)>=0)?(4*m-2):(4*m-2+fW),4*n+1,fL,fW,fH)]+arr[I(2,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),((4*m-2)>=0)?(4*m-2):(4*m-2+fW),4*n+3,fL,fW,fH)]+2*arr[I(2,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),4*m,4*n+1,fL,fW,fH)]+2*arr[I(2,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),4*m,4*n+3,fL,fW,fH)]+arr[I(2,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),4*m+2,4*n+1,fL,fW,fH)]+arr[I(2,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),4*m+2,4*n+3,fL,fW,fH)]+2*arr[I(2,4*l,((4*m-2)>=0)?(4*m-2):(4*m-2+fW),4*n+1,fL,fW,fH)]+2*arr[I(2,4*l,((4*m-2)>=0)?(4*m-2):(4*m-2+fW),4*n+3,fL,fW,fH)]+4*arr[I(2,4*l,4*m,4*n+1,fL,fW,fH)]+4*arr[I(2,4*l,4*m,4*n+3,fL,fW,fH)]+2*arr[I(2,4*l,4*m+2,4*n+1,fL,fW,fH)]+2*arr[I(2,4*l,4*m+2,4*n+3,fL,fW,fH)]+arr[I(2,4*l+2,((4*m-2)>=0)?(4*m-2):(4*m-2+fW),4*n+1,fL,fW,fH)]+arr[I(2,4*l+2,((4*m-2)>=0)?(4*m-2):(4*m-2+fW),4*n+3,fL,fW,fH)]+2*arr[I(2,4*l+2,4*m,4*n+1,fL,fW,fH)]+2*arr[I(2,4*l+2,4*m,4*n+3,fL,fW,fH)]+arr[I(2,4*l+2,4*m+2,4*n+1,fL,fW,fH)]+arr[I(2,4*l+2,4*m+2,4*n+3,fL,fW,fH)])/8;
out[I(2,2*l,2*m+1,2*n,cL,cW,cH)]=(arr[I(2,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),4*m+1,4*n,fL,fW,fH)]+arr[I(2,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),4*m+3,4*n,fL,fW,fH)]+2*arr[I(2,4*l,4*m+1,4*n,fL,fW,fH)]+2*arr[I(2,4*l,4*m+3,4*n,fL,fW,fH)]+arr[I(2,4*l+2,4*m+1,4*n,fL,fW,fH)]+arr[I(2,4*l+2,4*m+3,4*n,fL,fW,fH)])/2;
out[I(2,2*l,2*m+1,2*n+1,cL,cW,cH)]=(arr[I(2,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),4*m+1,4*n+1,fL,fW,fH)]+arr[I(2,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),4*m+1,4*n+3,fL,fW,fH)]+arr[I(2,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),4*m+3,4*n+1,fL,fW,fH)]+arr[I(2,((4*l-2)>=0)?(4*l-2):(4*l-2+fL),4*m+3,4*n+3,fL,fW,fH)]+2*arr[I(2,4*l,4*m+1,4*n+1,fL,fW,fH)]+2*arr[I(2,4*l,4*m+1,4*n+3,fL,fW,fH)]+2*arr[I(2,4*l,4*m+3,4*n+1,fL,fW,fH)]+2*arr[I(2,4*l,4*m+3,4*n+3,fL,fW,fH)]+arr[I(2,4*l+2,4*m+1,4*n+1,fL,fW,fH)]+arr[I(2,4*l+2,4*m+1,4*n+3,fL,fW,fH)]+arr[I(2,4*l+2,4*m+3,4*n+1,fL,fW,fH)]+arr[I(2,4*l+2,4*m+3,4*n+3,fL,fW,fH)])/4;
out[I(2,2*l+1,2*m,2*n,cL,cW,cH)]=(arr[I(2,4*l+1,((4*m-2)>=0)?(4*m-2):(4*m-2+fW),4*n,fL,fW,fH)]+2*arr[I(2,4*l+1,4*m,4*n,fL,fW,fH)]+arr[I(2,4*l+1,4*m+2,4*n,fL,fW,fH)]+arr[I(2,4*l+3,((4*m-2)>=0)?(4*m-2):(4*m-2+fW),4*n,fL,fW,fH)]+2*arr[I(2,4*l+3,4*m,4*n,fL,fW,fH)]+arr[I(2,4*l+3,4*m+2,4*n,fL,fW,fH)])/2;
out[I(2,2*l+1,2*m,2*n+1,cL,cW,cH)]=(arr[I(2,4*l+1,((4*m-2)>=0)?(4*m-2):(4*m-2+fW),4*n+1,fL,fW,fH)]+arr[I(2,4*l+1,((4*m-2)>=0)?(4*m-2):(4*m-2+fW),4*n+3,fL,fW,fH)]+2*arr[I(2,4*l+1,4*m,4*n+1,fL,fW,fH)]+2*arr[I(2,4*l+1,4*m,4*n+3,fL,fW,fH)]+arr[I(2,4*l+1,4*m+2,4*n+1,fL,fW,fH)]+arr[I(2,4*l+1,4*m+2,4*n+3,fL,fW,fH)]+arr[I(2,4*l+3,((4*m-2)>=0)?(4*m-2):(4*m-2+fW),4*n+1,fL,fW,fH)]+arr[I(2,4*l+3,((4*m-2)>=0)?(4*m-2):(4*m-2+fW),4*n+3,fL,fW,fH)]+2*arr[I(2,4*l+3,4*m,4*n+1,fL,fW,fH)]+2*arr[I(2,4*l+3,4*m,4*n+3,fL,fW,fH)]+arr[I(2,4*l+3,4*m+2,4*n+1,fL,fW,fH)]+arr[I(2,4*l+3,4*m+2,4*n+3,fL,fW,fH)])/4;
out[I(2,2*l+1,2*m+1,2*n,cL,cW,cH)]=arr[I(2,4*l+1,4*m+1,4*n,fL,fW,fH)]+arr[I(2,4*l+1,4*m+3,4*n,fL,fW,fH)]+arr[I(2,4*l+3,4*m+1,4*n,fL,fW,fH)]+arr[I(2,4*l+3,4*m+3,4*n,fL,fW,fH)];
out[I(2,2*l+1,2*m+1,2*n+1,cL,cW,cH)]=(arr[I(2,4*l+1,4*m+1,4*n+1,fL,fW,fH)]+arr[I(2,4*l+1,4*m+1,4*n+3,fL,fW,fH)]+arr[I(2,4*l+1,4*m+3,4*n+1,fL,fW,fH)]+arr[I(2,4*l+1,4*m+3,4*n+3,fL,fW,fH)]+arr[I(2,4*l+3,4*m+1,4*n+1,fL,fW,fH)]+arr[I(2,4*l+3,4*m+1,4*n+3,fL,fW,fH)]+arr[I(2,4*l+3,4*m+3,4*n+1,fL,fW,fH)]+arr[I(2,4*l+3,4*m+3,4*n+3,fL,fW,fH)])/2;
}
}
}
}


void coarse2fine(double arr[], double out[], int cL, int cW, int cH)
{
int fL= cL*2; //fL,fW,fH are length, width and height of the coarse grid.
int fW= cW*2; //cL,cW,cH are length, width and height of the fine grid.
int fH= cH*2;

      for(int l=0; l< cL/2; l++) {
         for(int m=0; m< cW/2; m++) {
            for(int n=0; n< cH/2; n++){
out[I(0,4*l,4*m,4*n,fL,fW,fH)]=(arr[I(0,2*l,2*m,2*n,cL,cW,cH)])/4;
out[I(0,4*l,4*m,4*n+1,fL,fW,fH)]=(arr[I(0,2*l,2*m,2*n+1,cL,cW,cH)])/4;
out[I(0,4*l,4*m,4*n+2,fL,fW,fH)]=(arr[I(0,2*l,2*m,2*n,cL,cW,cH)]+arr[I(0,2*l,2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/8;
out[I(0,4*l,4*m,4*n+3,fL,fW,fH)]=(arr[I(0,2*l,2*m,2*n+1,cL,cW,cH)])/4;
out[I(0,4*l,4*m+1,4*n,fL,fW,fH)]=(arr[I(0,2*l,2*m+1,2*n,cL,cW,cH)])/4;
out[I(0,4*l,4*m+1,4*n+1,fL,fW,fH)]=(arr[I(0,2*l,2*m+1,2*n+1,cL,cW,cH)])/4;
out[I(0,4*l,4*m+1,4*n+2,fL,fW,fH)]=(arr[I(0,2*l,2*m+1,2*n,cL,cW,cH)]+arr[I(0,2*l,2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/8;
out[I(0,4*l,4*m+1,4*n+3,fL,fW,fH)]=(arr[I(0,2*l,2*m+1,2*n+1,cL,cW,cH)])/4;
out[I(0,4*l,4*m+2,4*n,fL,fW,fH)]=(arr[I(0,2*l,2*m,2*n,cL,cW,cH)]+arr[I(0,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)])/8;
out[I(0,4*l,4*m+2,4*n+1,fL,fW,fH)]=(arr[I(0,2*l,2*m,2*n+1,cL,cW,cH)]+arr[I(0,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)])/8;
out[I(0,4*l,4*m+2,4*n+2,fL,fW,fH)]=(arr[I(0,2*l,2*m,2*n,cL,cW,cH)]+arr[I(0,2*l,2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(0,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)]+arr[I(0,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/16;
out[I(0,4*l,4*m+2,4*n+3,fL,fW,fH)]=(arr[I(0,2*l,2*m,2*n+1,cL,cW,cH)]+arr[I(0,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)])/8;
out[I(0,4*l,4*m+3,4*n,fL,fW,fH)]=(arr[I(0,2*l,2*m+1,2*n,cL,cW,cH)])/4;
out[I(0,4*l,4*m+3,4*n+1,fL,fW,fH)]=(arr[I(0,2*l,2*m+1,2*n+1,cL,cW,cH)])/4;
out[I(0,4*l,4*m+3,4*n+2,fL,fW,fH)]=(arr[I(0,2*l,2*m+1,2*n,cL,cW,cH)]+arr[I(0,2*l,2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/8;
out[I(0,4*l,4*m+3,4*n+3,fL,fW,fH)]=(arr[I(0,2*l,2*m+1,2*n+1,cL,cW,cH)])/4;
out[I(0,4*l+1,4*m,4*n,fL,fW,fH)]=(arr[I(0,((2*l-1)>=0)?(2*l-1):(2*l-1+cL),2*m,2*n,cL,cW,cH)]+3*arr[I(0,2*l+1,2*m,2*n,cL,cW,cH)])/16;
out[I(0,4*l+1,4*m,4*n+1,fL,fW,fH)]=(arr[I(0,((2*l-1)>=0)?(2*l-1):(2*l-1+cL),2*m,2*n+1,cL,cW,cH)]+3*arr[I(0,2*l+1,2*m,2*n+1,cL,cW,cH)])/16;
out[I(0,4*l+1,4*m,4*n+2,fL,fW,fH)]=(arr[I(0,((2*l-1)>=0)?(2*l-1):(2*l-1+cL),2*m,2*n,cL,cW,cH)]+arr[I(0,((2*l-1)>=0)?(2*l-1):(2*l-1+cL),2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+3*arr[I(0,2*l+1,2*m,2*n,cL,cW,cH)]+3*arr[I(0,2*l+1,2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/32;
out[I(0,4*l+1,4*m,4*n+3,fL,fW,fH)]=(arr[I(0,((2*l-1)>=0)?(2*l-1):(2*l-1+cL),2*m,2*n+1,cL,cW,cH)]+3*arr[I(0,2*l+1,2*m,2*n+1,cL,cW,cH)])/16;
out[I(0,4*l+1,4*m+1,4*n,fL,fW,fH)]=(arr[I(0,((2*l-1)>=0)?(2*l-1):(2*l-1+cL),2*m+1,2*n,cL,cW,cH)]+3*arr[I(0,2*l+1,2*m+1,2*n,cL,cW,cH)])/16;
out[I(0,4*l+1,4*m+1,4*n+1,fL,fW,fH)]=(arr[I(0,((2*l-1)>=0)?(2*l-1):(2*l-1+cL),2*m+1,2*n+1,cL,cW,cH)]+3*arr[I(0,2*l+1,2*m+1,2*n+1,cL,cW,cH)])/16;
out[I(0,4*l+1,4*m+1,4*n+2,fL,fW,fH)]=(arr[I(0,((2*l-1)>=0)?(2*l-1):(2*l-1+cL),2*m+1,2*n,cL,cW,cH)]+arr[I(0,((2*l-1)>=0)?(2*l-1):(2*l-1+cL),2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+3*arr[I(0,2*l+1,2*m+1,2*n,cL,cW,cH)]+3*arr[I(0,2*l+1,2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/32;
out[I(0,4*l+1,4*m+1,4*n+3,fL,fW,fH)]=(arr[I(0,((2*l-1)>=0)?(2*l-1):(2*l-1+cL),2*m+1,2*n+1,cL,cW,cH)]+3*arr[I(0,2*l+1,2*m+1,2*n+1,cL,cW,cH)])/16;
out[I(0,4*l+1,4*m+2,4*n,fL,fW,fH)]=(arr[I(0,((2*l-1)>=0)?(2*l-1):(2*l-1+cL),2*m,2*n,cL,cW,cH)]+arr[I(0,((2*l-1)>=0)?(2*l-1):(2*l-1+cL),((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)]+3*arr[I(0,2*l+1,2*m,2*n,cL,cW,cH)]+3*arr[I(0,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)])/32;
out[I(0,4*l+1,4*m+2,4*n+1,fL,fW,fH)]=(arr[I(0,((2*l-1)>=0)?(2*l-1):(2*l-1+cL),2*m,2*n+1,cL,cW,cH)]+arr[I(0,((2*l-1)>=0)?(2*l-1):(2*l-1+cL),((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)]+3*arr[I(0,2*l+1,2*m,2*n+1,cL,cW,cH)]+3*arr[I(0,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)])/32;
out[I(0,4*l+1,4*m+2,4*n+2,fL,fW,fH)]=(arr[I(0,((2*l-1)>=0)?(2*l-1):(2*l-1+cL),2*m,2*n,cL,cW,cH)]+arr[I(0,((2*l-1)>=0)?(2*l-1):(2*l-1+cL),2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(0,((2*l-1)>=0)?(2*l-1):(2*l-1+cL),((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)]+arr[I(0,((2*l-1)>=0)?(2*l-1):(2*l-1+cL),((2*m+2)<cW)?(2*m+2):(2*m+2-cW),((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+3*arr[I(0,2*l+1,2*m,2*n,cL,cW,cH)]+3*arr[I(0,2*l+1,2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+3*arr[I(0,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)]+3*arr[I(0,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/64;
out[I(0,4*l+1,4*m+2,4*n+3,fL,fW,fH)]=(arr[I(0,((2*l-1)>=0)?(2*l-1):(2*l-1+cL),2*m,2*n+1,cL,cW,cH)]+arr[I(0,((2*l-1)>=0)?(2*l-1):(2*l-1+cL),((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)]+3*arr[I(0,2*l+1,2*m,2*n+1,cL,cW,cH)]+3*arr[I(0,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)])/32;
out[I(0,4*l+1,4*m+3,4*n,fL,fW,fH)]=(arr[I(0,((2*l-1)>=0)?(2*l-1):(2*l-1+cL),2*m+1,2*n,cL,cW,cH)]+3*arr[I(0,2*l+1,2*m+1,2*n,cL,cW,cH)])/16;
out[I(0,4*l+1,4*m+3,4*n+1,fL,fW,fH)]=(arr[I(0,((2*l-1)>=0)?(2*l-1):(2*l-1+cL),2*m+1,2*n+1,cL,cW,cH)]+3*arr[I(0,2*l+1,2*m+1,2*n+1,cL,cW,cH)])/16;
out[I(0,4*l+1,4*m+3,4*n+2,fL,fW,fH)]=(arr[I(0,((2*l-1)>=0)?(2*l-1):(2*l-1+cL),2*m+1,2*n,cL,cW,cH)]+arr[I(0,((2*l-1)>=0)?(2*l-1):(2*l-1+cL),2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+3*arr[I(0,2*l+1,2*m+1,2*n,cL,cW,cH)]+3*arr[I(0,2*l+1,2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/32;
out[I(0,4*l+1,4*m+3,4*n+3,fL,fW,fH)]=(arr[I(0,((2*l-1)>=0)?(2*l-1):(2*l-1+cL),2*m+1,2*n+1,cL,cW,cH)]+3*arr[I(0,2*l+1,2*m+1,2*n+1,cL,cW,cH)])/16;
out[I(0,4*l+2,4*m,4*n,fL,fW,fH)]=(arr[I(0,2*l,2*m,2*n,cL,cW,cH)]+arr[I(0,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,2*n,cL,cW,cH)])/8;
out[I(0,4*l+2,4*m,4*n+1,fL,fW,fH)]=(arr[I(0,2*l,2*m,2*n+1,cL,cW,cH)]+arr[I(0,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,2*n+1,cL,cW,cH)])/8;
out[I(0,4*l+2,4*m,4*n+2,fL,fW,fH)]=(arr[I(0,2*l,2*m,2*n,cL,cW,cH)]+arr[I(0,2*l,2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(0,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,2*n,cL,cW,cH)]+arr[I(0,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/16;
out[I(0,4*l+2,4*m,4*n+3,fL,fW,fH)]=(arr[I(0,2*l,2*m,2*n+1,cL,cW,cH)]+arr[I(0,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,2*n+1,cL,cW,cH)])/8;
out[I(0,4*l+2,4*m+1,4*n,fL,fW,fH)]=(arr[I(0,2*l,2*m+1,2*n,cL,cW,cH)]+arr[I(0,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,2*n,cL,cW,cH)])/8;
out[I(0,4*l+2,4*m+1,4*n+1,fL,fW,fH)]=(arr[I(0,2*l,2*m+1,2*n+1,cL,cW,cH)]+arr[I(0,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,2*n+1,cL,cW,cH)])/8;
out[I(0,4*l+2,4*m+1,4*n+2,fL,fW,fH)]=(arr[I(0,2*l,2*m+1,2*n,cL,cW,cH)]+arr[I(0,2*l,2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(0,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,2*n,cL,cW,cH)]+arr[I(0,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/16;
out[I(0,4*l+2,4*m+1,4*n+3,fL,fW,fH)]=(arr[I(0,2*l,2*m+1,2*n+1,cL,cW,cH)]+arr[I(0,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,2*n+1,cL,cW,cH)])/8;
out[I(0,4*l+2,4*m+2,4*n,fL,fW,fH)]=(arr[I(0,2*l,2*m,2*n,cL,cW,cH)]+arr[I(0,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)]+arr[I(0,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,2*n,cL,cW,cH)]+arr[I(0,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)])/16;
out[I(0,4*l+2,4*m+2,4*n+1,fL,fW,fH)]=(arr[I(0,2*l,2*m,2*n+1,cL,cW,cH)]+arr[I(0,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)]+arr[I(0,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,2*n+1,cL,cW,cH)]+arr[I(0,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)])/16;
out[I(0,4*l+2,4*m+2,4*n+2,fL,fW,fH)]=(arr[I(0,2*l,2*m,2*n,cL,cW,cH)]+arr[I(0,2*l,2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(0,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)]+arr[I(0,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(0,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,2*n,cL,cW,cH)]+arr[I(0,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(0,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)]+arr[I(0,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),((2*m+2)<cW)?(2*m+2):(2*m+2-cW),((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/32;
out[I(0,4*l+2,4*m+2,4*n+3,fL,fW,fH)]=(arr[I(0,2*l,2*m,2*n+1,cL,cW,cH)]+arr[I(0,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)]+arr[I(0,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,2*n+1,cL,cW,cH)]+arr[I(0,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)])/16;
out[I(0,4*l+2,4*m+3,4*n,fL,fW,fH)]=(arr[I(0,2*l,2*m+1,2*n,cL,cW,cH)]+arr[I(0,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,2*n,cL,cW,cH)])/8;
out[I(0,4*l+2,4*m+3,4*n+1,fL,fW,fH)]=(arr[I(0,2*l,2*m+1,2*n+1,cL,cW,cH)]+arr[I(0,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,2*n+1,cL,cW,cH)])/8;
out[I(0,4*l+2,4*m+3,4*n+2,fL,fW,fH)]=(arr[I(0,2*l,2*m+1,2*n,cL,cW,cH)]+arr[I(0,2*l,2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(0,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,2*n,cL,cW,cH)]+arr[I(0,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/16;
out[I(0,4*l+2,4*m+3,4*n+3,fL,fW,fH)]=(arr[I(0,2*l,2*m+1,2*n+1,cL,cW,cH)]+arr[I(0,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,2*n+1,cL,cW,cH)])/8;
out[I(0,4*l+3,4*m,4*n,fL,fW,fH)]=(3*arr[I(0,2*l+1,2*m,2*n,cL,cW,cH)]+arr[I(0,((2*l+3)<cL)?(2*l+3):(2*l+3-cL),2*m,2*n,cL,cW,cH)])/16;
out[I(0,4*l+3,4*m,4*n+1,fL,fW,fH)]=(3*arr[I(0,2*l+1,2*m,2*n+1,cL,cW,cH)]+arr[I(0,((2*l+3)<cL)?(2*l+3):(2*l+3-cL),2*m,2*n+1,cL,cW,cH)])/16;
out[I(0,4*l+3,4*m,4*n+2,fL,fW,fH)]=(3*arr[I(0,2*l+1,2*m,2*n,cL,cW,cH)]+3*arr[I(0,2*l+1,2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(0,((2*l+3)<cL)?(2*l+3):(2*l+3-cL),2*m,2*n,cL,cW,cH)]+arr[I(0,((2*l+3)<cL)?(2*l+3):(2*l+3-cL),2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/32;
out[I(0,4*l+3,4*m,4*n+3,fL,fW,fH)]=(3*arr[I(0,2*l+1,2*m,2*n+1,cL,cW,cH)]+arr[I(0,((2*l+3)<cL)?(2*l+3):(2*l+3-cL),2*m,2*n+1,cL,cW,cH)])/16;
out[I(0,4*l+3,4*m+1,4*n,fL,fW,fH)]=(3*arr[I(0,2*l+1,2*m+1,2*n,cL,cW,cH)]+arr[I(0,((2*l+3)<cL)?(2*l+3):(2*l+3-cL),2*m+1,2*n,cL,cW,cH)])/16;
out[I(0,4*l+3,4*m+1,4*n+1,fL,fW,fH)]=(3*arr[I(0,2*l+1,2*m+1,2*n+1,cL,cW,cH)]+arr[I(0,((2*l+3)<cL)?(2*l+3):(2*l+3-cL),2*m+1,2*n+1,cL,cW,cH)])/16;
out[I(0,4*l+3,4*m+1,4*n+2,fL,fW,fH)]=(3*arr[I(0,2*l+1,2*m+1,2*n,cL,cW,cH)]+3*arr[I(0,2*l+1,2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(0,((2*l+3)<cL)?(2*l+3):(2*l+3-cL),2*m+1,2*n,cL,cW,cH)]+arr[I(0,((2*l+3)<cL)?(2*l+3):(2*l+3-cL),2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/32;
out[I(0,4*l+3,4*m+1,4*n+3,fL,fW,fH)]=(3*arr[I(0,2*l+1,2*m+1,2*n+1,cL,cW,cH)]+arr[I(0,((2*l+3)<cL)?(2*l+3):(2*l+3-cL),2*m+1,2*n+1,cL,cW,cH)])/16;
out[I(0,4*l+3,4*m+2,4*n,fL,fW,fH)]=(3*arr[I(0,2*l+1,2*m,2*n,cL,cW,cH)]+3*arr[I(0,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)]+arr[I(0,((2*l+3)<cL)?(2*l+3):(2*l+3-cL),2*m,2*n,cL,cW,cH)]+arr[I(0,((2*l+3)<cL)?(2*l+3):(2*l+3-cL),((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)])/32;
out[I(0,4*l+3,4*m+2,4*n+1,fL,fW,fH)]=(3*arr[I(0,2*l+1,2*m,2*n+1,cL,cW,cH)]+3*arr[I(0,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)]+arr[I(0,((2*l+3)<cL)?(2*l+3):(2*l+3-cL),2*m,2*n+1,cL,cW,cH)]+arr[I(0,((2*l+3)<cL)?(2*l+3):(2*l+3-cL),((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)])/32;
out[I(0,4*l+3,4*m+2,4*n+2,fL,fW,fH)]=(3*arr[I(0,2*l+1,2*m,2*n,cL,cW,cH)]+3*arr[I(0,2*l+1,2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+3*arr[I(0,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)]+3*arr[I(0,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(0,((2*l+3)<cL)?(2*l+3):(2*l+3-cL),2*m,2*n,cL,cW,cH)]+arr[I(0,((2*l+3)<cL)?(2*l+3):(2*l+3-cL),2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(0,((2*l+3)<cL)?(2*l+3):(2*l+3-cL),((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)]+arr[I(0,((2*l+3)<cL)?(2*l+3):(2*l+3-cL),((2*m+2)<cW)?(2*m+2):(2*m+2-cW),((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/64;
out[I(0,4*l+3,4*m+2,4*n+3,fL,fW,fH)]=(3*arr[I(0,2*l+1,2*m,2*n+1,cL,cW,cH)]+3*arr[I(0,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)]+arr[I(0,((2*l+3)<cL)?(2*l+3):(2*l+3-cL),2*m,2*n+1,cL,cW,cH)]+arr[I(0,((2*l+3)<cL)?(2*l+3):(2*l+3-cL),((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)])/32;
out[I(0,4*l+3,4*m+3,4*n,fL,fW,fH)]=(3*arr[I(0,2*l+1,2*m+1,2*n,cL,cW,cH)]+arr[I(0,((2*l+3)<cL)?(2*l+3):(2*l+3-cL),2*m+1,2*n,cL,cW,cH)])/16;
out[I(0,4*l+3,4*m+3,4*n+1,fL,fW,fH)]=(3*arr[I(0,2*l+1,2*m+1,2*n+1,cL,cW,cH)]+arr[I(0,((2*l+3)<cL)?(2*l+3):(2*l+3-cL),2*m+1,2*n+1,cL,cW,cH)])/16;
out[I(0,4*l+3,4*m+3,4*n+2,fL,fW,fH)]=(3*arr[I(0,2*l+1,2*m+1,2*n,cL,cW,cH)]+3*arr[I(0,2*l+1,2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(0,((2*l+3)<cL)?(2*l+3):(2*l+3-cL),2*m+1,2*n,cL,cW,cH)]+arr[I(0,((2*l+3)<cL)?(2*l+3):(2*l+3-cL),2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/32;
out[I(0,4*l+3,4*m+3,4*n+3,fL,fW,fH)]=(3*arr[I(0,2*l+1,2*m+1,2*n+1,cL,cW,cH)]+arr[I(0,((2*l+3)<cL)?(2*l+3):(2*l+3-cL),2*m+1,2*n+1,cL,cW,cH)])/16;
out[I(1,4*l,4*m,4*n,fL,fW,fH)]=(arr[I(1,2*l,2*m,2*n,cL,cW,cH)])/4;
out[I(1,4*l,4*m,4*n+1,fL,fW,fH)]=(arr[I(1,2*l,2*m,2*n+1,cL,cW,cH)])/4;
out[I(1,4*l,4*m,4*n+2,fL,fW,fH)]=(arr[I(1,2*l,2*m,2*n,cL,cW,cH)]+arr[I(1,2*l,2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/8;
out[I(1,4*l,4*m,4*n+3,fL,fW,fH)]=(arr[I(1,2*l,2*m,2*n+1,cL,cW,cH)])/4;
out[I(1,4*l,4*m+1,4*n,fL,fW,fH)]=(arr[I(1,2*l,((2*m-1)>=0)?(2*m-1):(2*m-1+cW),2*n,cL,cW,cH)]+3*arr[I(1,2*l,2*m+1,2*n,cL,cW,cH)])/16;
out[I(1,4*l,4*m+1,4*n+1,fL,fW,fH)]=(arr[I(1,2*l,((2*m-1)>=0)?(2*m-1):(2*m-1+cW),2*n+1,cL,cW,cH)]+3*arr[I(1,2*l,2*m+1,2*n+1,cL,cW,cH)])/16;
out[I(1,4*l,4*m+1,4*n+2,fL,fW,fH)]=(arr[I(1,2*l,((2*m-1)>=0)?(2*m-1):(2*m-1+cW),2*n,cL,cW,cH)]+arr[I(1,2*l,((2*m-1)>=0)?(2*m-1):(2*m-1+cW),((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+3*arr[I(1,2*l,2*m+1,2*n,cL,cW,cH)]+3*arr[I(1,2*l,2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/32;
out[I(1,4*l,4*m+1,4*n+3,fL,fW,fH)]=(arr[I(1,2*l,((2*m-1)>=0)?(2*m-1):(2*m-1+cW),2*n+1,cL,cW,cH)]+3*arr[I(1,2*l,2*m+1,2*n+1,cL,cW,cH)])/16;
out[I(1,4*l,4*m+2,4*n,fL,fW,fH)]=(arr[I(1,2*l,2*m,2*n,cL,cW,cH)]+arr[I(1,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)])/8;
out[I(1,4*l,4*m+2,4*n+1,fL,fW,fH)]=(arr[I(1,2*l,2*m,2*n+1,cL,cW,cH)]+arr[I(1,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)])/8;
out[I(1,4*l,4*m+2,4*n+2,fL,fW,fH)]=(arr[I(1,2*l,2*m,2*n,cL,cW,cH)]+arr[I(1,2*l,2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(1,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)]+arr[I(1,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/16;
out[I(1,4*l,4*m+2,4*n+3,fL,fW,fH)]=(arr[I(1,2*l,2*m,2*n+1,cL,cW,cH)]+arr[I(1,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)])/8;
out[I(1,4*l,4*m+3,4*n,fL,fW,fH)]=(3*arr[I(1,2*l,2*m+1,2*n,cL,cW,cH)]+arr[I(1,2*l,((2*m+3)<cW)?(2*m+3):(2*m+3-cW),2*n,cL,cW,cH)])/16;
out[I(1,4*l,4*m+3,4*n+1,fL,fW,fH)]=(3*arr[I(1,2*l,2*m+1,2*n+1,cL,cW,cH)]+arr[I(1,2*l,((2*m+3)<cW)?(2*m+3):(2*m+3-cW),2*n+1,cL,cW,cH)])/16;
out[I(1,4*l,4*m+3,4*n+2,fL,fW,fH)]=(3*arr[I(1,2*l,2*m+1,2*n,cL,cW,cH)]+3*arr[I(1,2*l,2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(1,2*l,((2*m+3)<cW)?(2*m+3):(2*m+3-cW),2*n,cL,cW,cH)]+arr[I(1,2*l,((2*m+3)<cW)?(2*m+3):(2*m+3-cW),((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/32;
out[I(1,4*l,4*m+3,4*n+3,fL,fW,fH)]=(3*arr[I(1,2*l,2*m+1,2*n+1,cL,cW,cH)]+arr[I(1,2*l,((2*m+3)<cW)?(2*m+3):(2*m+3-cW),2*n+1,cL,cW,cH)])/16;
out[I(1,4*l+1,4*m,4*n,fL,fW,fH)]=(arr[I(1,2*l+1,2*m,2*n,cL,cW,cH)])/4;
out[I(1,4*l+1,4*m,4*n+1,fL,fW,fH)]=(arr[I(1,2*l+1,2*m,2*n+1,cL,cW,cH)])/4;
out[I(1,4*l+1,4*m,4*n+2,fL,fW,fH)]=(arr[I(1,2*l+1,2*m,2*n,cL,cW,cH)]+arr[I(1,2*l+1,2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/8;
out[I(1,4*l+1,4*m,4*n+3,fL,fW,fH)]=(arr[I(1,2*l+1,2*m,2*n+1,cL,cW,cH)])/4;
out[I(1,4*l+1,4*m+1,4*n,fL,fW,fH)]=(arr[I(1,2*l+1,((2*m-1)>=0)?(2*m-1):(2*m-1+cW),2*n,cL,cW,cH)]+3*arr[I(1,2*l+1,2*m+1,2*n,cL,cW,cH)])/16;
out[I(1,4*l+1,4*m+1,4*n+1,fL,fW,fH)]=(arr[I(1,2*l+1,((2*m-1)>=0)?(2*m-1):(2*m-1+cW),2*n+1,cL,cW,cH)]+3*arr[I(1,2*l+1,2*m+1,2*n+1,cL,cW,cH)])/16;
out[I(1,4*l+1,4*m+1,4*n+2,fL,fW,fH)]=(arr[I(1,2*l+1,((2*m-1)>=0)?(2*m-1):(2*m-1+cW),2*n,cL,cW,cH)]+arr[I(1,2*l+1,((2*m-1)>=0)?(2*m-1):(2*m-1+cW),((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+3*arr[I(1,2*l+1,2*m+1,2*n,cL,cW,cH)]+3*arr[I(1,2*l+1,2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/32;
out[I(1,4*l+1,4*m+1,4*n+3,fL,fW,fH)]=(arr[I(1,2*l+1,((2*m-1)>=0)?(2*m-1):(2*m-1+cW),2*n+1,cL,cW,cH)]+3*arr[I(1,2*l+1,2*m+1,2*n+1,cL,cW,cH)])/16;
out[I(1,4*l+1,4*m+2,4*n,fL,fW,fH)]=(arr[I(1,2*l+1,2*m,2*n,cL,cW,cH)]+arr[I(1,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)])/8;
out[I(1,4*l+1,4*m+2,4*n+1,fL,fW,fH)]=(arr[I(1,2*l+1,2*m,2*n+1,cL,cW,cH)]+arr[I(1,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)])/8;
out[I(1,4*l+1,4*m+2,4*n+2,fL,fW,fH)]=(arr[I(1,2*l+1,2*m,2*n,cL,cW,cH)]+arr[I(1,2*l+1,2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(1,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)]+arr[I(1,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/16;
out[I(1,4*l+1,4*m+2,4*n+3,fL,fW,fH)]=(arr[I(1,2*l+1,2*m,2*n+1,cL,cW,cH)]+arr[I(1,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)])/8;
out[I(1,4*l+1,4*m+3,4*n,fL,fW,fH)]=(3*arr[I(1,2*l+1,2*m+1,2*n,cL,cW,cH)]+arr[I(1,2*l+1,((2*m+3)<cW)?(2*m+3):(2*m+3-cW),2*n,cL,cW,cH)])/16;
out[I(1,4*l+1,4*m+3,4*n+1,fL,fW,fH)]=(3*arr[I(1,2*l+1,2*m+1,2*n+1,cL,cW,cH)]+arr[I(1,2*l+1,((2*m+3)<cW)?(2*m+3):(2*m+3-cW),2*n+1,cL,cW,cH)])/16;
out[I(1,4*l+1,4*m+3,4*n+2,fL,fW,fH)]=(3*arr[I(1,2*l+1,2*m+1,2*n,cL,cW,cH)]+3*arr[I(1,2*l+1,2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(1,2*l+1,((2*m+3)<cW)?(2*m+3):(2*m+3-cW),2*n,cL,cW,cH)]+arr[I(1,2*l+1,((2*m+3)<cW)?(2*m+3):(2*m+3-cW),((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/32;
out[I(1,4*l+1,4*m+3,4*n+3,fL,fW,fH)]=(3*arr[I(1,2*l+1,2*m+1,2*n+1,cL,cW,cH)]+arr[I(1,2*l+1,((2*m+3)<cW)?(2*m+3):(2*m+3-cW),2*n+1,cL,cW,cH)])/16;
out[I(1,4*l+2,4*m,4*n,fL,fW,fH)]=(arr[I(1,2*l,2*m,2*n,cL,cW,cH)]+arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,2*n,cL,cW,cH)])/8;
out[I(1,4*l+2,4*m,4*n+1,fL,fW,fH)]=(arr[I(1,2*l,2*m,2*n+1,cL,cW,cH)]+arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,2*n+1,cL,cW,cH)])/8;
out[I(1,4*l+2,4*m,4*n+2,fL,fW,fH)]=(arr[I(1,2*l,2*m,2*n,cL,cW,cH)]+arr[I(1,2*l,2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,2*n,cL,cW,cH)]+arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/16;
out[I(1,4*l+2,4*m,4*n+3,fL,fW,fH)]=(arr[I(1,2*l,2*m,2*n+1,cL,cW,cH)]+arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,2*n+1,cL,cW,cH)])/8;
out[I(1,4*l+2,4*m+1,4*n,fL,fW,fH)]=(arr[I(1,2*l,((2*m-1)>=0)?(2*m-1):(2*m-1+cW),2*n,cL,cW,cH)]+3*arr[I(1,2*l,2*m+1,2*n,cL,cW,cH)]+arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),((2*m-1)>=0)?(2*m-1):(2*m-1+cW),2*n,cL,cW,cH)]+3*arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,2*n,cL,cW,cH)])/32;
out[I(1,4*l+2,4*m+1,4*n+1,fL,fW,fH)]=(arr[I(1,2*l,((2*m-1)>=0)?(2*m-1):(2*m-1+cW),2*n+1,cL,cW,cH)]+3*arr[I(1,2*l,2*m+1,2*n+1,cL,cW,cH)]+arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),((2*m-1)>=0)?(2*m-1):(2*m-1+cW),2*n+1,cL,cW,cH)]+3*arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,2*n+1,cL,cW,cH)])/32;
out[I(1,4*l+2,4*m+1,4*n+2,fL,fW,fH)]=(arr[I(1,2*l,((2*m-1)>=0)?(2*m-1):(2*m-1+cW),2*n,cL,cW,cH)]+arr[I(1,2*l,((2*m-1)>=0)?(2*m-1):(2*m-1+cW),((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+3*arr[I(1,2*l,2*m+1,2*n,cL,cW,cH)]+3*arr[I(1,2*l,2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),((2*m-1)>=0)?(2*m-1):(2*m-1+cW),2*n,cL,cW,cH)]+arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),((2*m-1)>=0)?(2*m-1):(2*m-1+cW),((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+3*arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,2*n,cL,cW,cH)]+3*arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/64;
out[I(1,4*l+2,4*m+1,4*n+3,fL,fW,fH)]=(arr[I(1,2*l,((2*m-1)>=0)?(2*m-1):(2*m-1+cW),2*n+1,cL,cW,cH)]+3*arr[I(1,2*l,2*m+1,2*n+1,cL,cW,cH)]+arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),((2*m-1)>=0)?(2*m-1):(2*m-1+cW),2*n+1,cL,cW,cH)]+3*arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,2*n+1,cL,cW,cH)])/32;
out[I(1,4*l+2,4*m+2,4*n,fL,fW,fH)]=(arr[I(1,2*l,2*m,2*n,cL,cW,cH)]+arr[I(1,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)]+arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,2*n,cL,cW,cH)]+arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)])/16;
out[I(1,4*l+2,4*m+2,4*n+1,fL,fW,fH)]=(arr[I(1,2*l,2*m,2*n+1,cL,cW,cH)]+arr[I(1,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)]+arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,2*n+1,cL,cW,cH)]+arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)])/16;
out[I(1,4*l+2,4*m+2,4*n+2,fL,fW,fH)]=(arr[I(1,2*l,2*m,2*n,cL,cW,cH)]+arr[I(1,2*l,2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(1,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)]+arr[I(1,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,2*n,cL,cW,cH)]+arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)]+arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),((2*m+2)<cW)?(2*m+2):(2*m+2-cW),((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/32;
out[I(1,4*l+2,4*m+2,4*n+3,fL,fW,fH)]=(arr[I(1,2*l,2*m,2*n+1,cL,cW,cH)]+arr[I(1,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)]+arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,2*n+1,cL,cW,cH)]+arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)])/16;
out[I(1,4*l+2,4*m+3,4*n,fL,fW,fH)]=(3*arr[I(1,2*l,2*m+1,2*n,cL,cW,cH)]+arr[I(1,2*l,((2*m+3)<cW)?(2*m+3):(2*m+3-cW),2*n,cL,cW,cH)]+3*arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,2*n,cL,cW,cH)]+arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),((2*m+3)<cW)?(2*m+3):(2*m+3-cW),2*n,cL,cW,cH)])/32;
out[I(1,4*l+2,4*m+3,4*n+1,fL,fW,fH)]=(3*arr[I(1,2*l,2*m+1,2*n+1,cL,cW,cH)]+arr[I(1,2*l,((2*m+3)<cW)?(2*m+3):(2*m+3-cW),2*n+1,cL,cW,cH)]+3*arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,2*n+1,cL,cW,cH)]+arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),((2*m+3)<cW)?(2*m+3):(2*m+3-cW),2*n+1,cL,cW,cH)])/32;
out[I(1,4*l+2,4*m+3,4*n+2,fL,fW,fH)]=(3*arr[I(1,2*l,2*m+1,2*n,cL,cW,cH)]+3*arr[I(1,2*l,2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(1,2*l,((2*m+3)<cW)?(2*m+3):(2*m+3-cW),2*n,cL,cW,cH)]+arr[I(1,2*l,((2*m+3)<cW)?(2*m+3):(2*m+3-cW),((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+3*arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,2*n,cL,cW,cH)]+3*arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),((2*m+3)<cW)?(2*m+3):(2*m+3-cW),2*n,cL,cW,cH)]+arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),((2*m+3)<cW)?(2*m+3):(2*m+3-cW),((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/64;
out[I(1,4*l+2,4*m+3,4*n+3,fL,fW,fH)]=(3*arr[I(1,2*l,2*m+1,2*n+1,cL,cW,cH)]+arr[I(1,2*l,((2*m+3)<cW)?(2*m+3):(2*m+3-cW),2*n+1,cL,cW,cH)]+3*arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,2*n+1,cL,cW,cH)]+arr[I(1,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),((2*m+3)<cW)?(2*m+3):(2*m+3-cW),2*n+1,cL,cW,cH)])/32;
out[I(1,4*l+3,4*m,4*n,fL,fW,fH)]=(arr[I(1,2*l+1,2*m,2*n,cL,cW,cH)])/4;
out[I(1,4*l+3,4*m,4*n+1,fL,fW,fH)]=(arr[I(1,2*l+1,2*m,2*n+1,cL,cW,cH)])/4;
out[I(1,4*l+3,4*m,4*n+2,fL,fW,fH)]=(arr[I(1,2*l+1,2*m,2*n,cL,cW,cH)]+arr[I(1,2*l+1,2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/8;
out[I(1,4*l+3,4*m,4*n+3,fL,fW,fH)]=(arr[I(1,2*l+1,2*m,2*n+1,cL,cW,cH)])/4;
out[I(1,4*l+3,4*m+1,4*n,fL,fW,fH)]=(arr[I(1,2*l+1,((2*m-1)>=0)?(2*m-1):(2*m-1+cW),2*n,cL,cW,cH)]+3*arr[I(1,2*l+1,2*m+1,2*n,cL,cW,cH)])/16;
out[I(1,4*l+3,4*m+1,4*n+1,fL,fW,fH)]=(arr[I(1,2*l+1,((2*m-1)>=0)?(2*m-1):(2*m-1+cW),2*n+1,cL,cW,cH)]+3*arr[I(1,2*l+1,2*m+1,2*n+1,cL,cW,cH)])/16;
out[I(1,4*l+3,4*m+1,4*n+2,fL,fW,fH)]=(arr[I(1,2*l+1,((2*m-1)>=0)?(2*m-1):(2*m-1+cW),2*n,cL,cW,cH)]+arr[I(1,2*l+1,((2*m-1)>=0)?(2*m-1):(2*m-1+cW),((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+3*arr[I(1,2*l+1,2*m+1,2*n,cL,cW,cH)]+3*arr[I(1,2*l+1,2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/32;
out[I(1,4*l+3,4*m+1,4*n+3,fL,fW,fH)]=(arr[I(1,2*l+1,((2*m-1)>=0)?(2*m-1):(2*m-1+cW),2*n+1,cL,cW,cH)]+3*arr[I(1,2*l+1,2*m+1,2*n+1,cL,cW,cH)])/16;
out[I(1,4*l+3,4*m+2,4*n,fL,fW,fH)]=(arr[I(1,2*l+1,2*m,2*n,cL,cW,cH)]+arr[I(1,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)])/8;
out[I(1,4*l+3,4*m+2,4*n+1,fL,fW,fH)]=(arr[I(1,2*l+1,2*m,2*n+1,cL,cW,cH)]+arr[I(1,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)])/8;
out[I(1,4*l+3,4*m+2,4*n+2,fL,fW,fH)]=(arr[I(1,2*l+1,2*m,2*n,cL,cW,cH)]+arr[I(1,2*l+1,2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(1,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)]+arr[I(1,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/16;
out[I(1,4*l+3,4*m+2,4*n+3,fL,fW,fH)]=(arr[I(1,2*l+1,2*m,2*n+1,cL,cW,cH)]+arr[I(1,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)])/8;
out[I(1,4*l+3,4*m+3,4*n,fL,fW,fH)]=(3*arr[I(1,2*l+1,2*m+1,2*n,cL,cW,cH)]+arr[I(1,2*l+1,((2*m+3)<cW)?(2*m+3):(2*m+3-cW),2*n,cL,cW,cH)])/16;
out[I(1,4*l+3,4*m+3,4*n+1,fL,fW,fH)]=(3*arr[I(1,2*l+1,2*m+1,2*n+1,cL,cW,cH)]+arr[I(1,2*l+1,((2*m+3)<cW)?(2*m+3):(2*m+3-cW),2*n+1,cL,cW,cH)])/16;
out[I(1,4*l+3,4*m+3,4*n+2,fL,fW,fH)]=(3*arr[I(1,2*l+1,2*m+1,2*n,cL,cW,cH)]+3*arr[I(1,2*l+1,2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(1,2*l+1,((2*m+3)<cW)?(2*m+3):(2*m+3-cW),2*n,cL,cW,cH)]+arr[I(1,2*l+1,((2*m+3)<cW)?(2*m+3):(2*m+3-cW),((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/32;
out[I(1,4*l+3,4*m+3,4*n+3,fL,fW,fH)]=(3*arr[I(1,2*l+1,2*m+1,2*n+1,cL,cW,cH)]+arr[I(1,2*l+1,((2*m+3)<cW)?(2*m+3):(2*m+3-cW),2*n+1,cL,cW,cH)])/16;
out[I(2,4*l,4*m,4*n,fL,fW,fH)]=(arr[I(2,2*l,2*m,2*n,cL,cW,cH)])/4;
out[I(2,4*l,4*m,4*n+1,fL,fW,fH)]=(arr[I(2,2*l,2*m,((2*n-1)>=0)?(2*n-1):(2*n-1+cH),cL,cW,cH)]+3*arr[I(2,2*l,2*m,2*n+1,cL,cW,cH)])/16;
out[I(2,4*l,4*m,4*n+2,fL,fW,fH)]=(arr[I(2,2*l,2*m,2*n,cL,cW,cH)]+arr[I(2,2*l,2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/8;
out[I(2,4*l,4*m,4*n+3,fL,fW,fH)]=(3*arr[I(2,2*l,2*m,2*n+1,cL,cW,cH)]+arr[I(2,2*l,2*m,((2*n+3)<cH)?(2*n+3):(2*n+3-cH),cL,cW,cH)])/16;
out[I(2,4*l,4*m+1,4*n,fL,fW,fH)]=(arr[I(2,2*l,2*m+1,2*n,cL,cW,cH)])/4;
out[I(2,4*l,4*m+1,4*n+1,fL,fW,fH)]=(arr[I(2,2*l,2*m+1,((2*n-1)>=0)?(2*n-1):(2*n-1+cH),cL,cW,cH)]+3*arr[I(2,2*l,2*m+1,2*n+1,cL,cW,cH)])/16;
out[I(2,4*l,4*m+1,4*n+2,fL,fW,fH)]=(arr[I(2,2*l,2*m+1,2*n,cL,cW,cH)]+arr[I(2,2*l,2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/8;
out[I(2,4*l,4*m+1,4*n+3,fL,fW,fH)]=(3*arr[I(2,2*l,2*m+1,2*n+1,cL,cW,cH)]+arr[I(2,2*l,2*m+1,((2*n+3)<cH)?(2*n+3):(2*n+3-cH),cL,cW,cH)])/16;
out[I(2,4*l,4*m+2,4*n,fL,fW,fH)]=(arr[I(2,2*l,2*m,2*n,cL,cW,cH)]+arr[I(2,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)])/8;
out[I(2,4*l,4*m+2,4*n+1,fL,fW,fH)]=(arr[I(2,2*l,2*m,((2*n-1)>=0)?(2*n-1):(2*n-1+cH),cL,cW,cH)]+3*arr[I(2,2*l,2*m,2*n+1,cL,cW,cH)]+arr[I(2,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),((2*n-1)>=0)?(2*n-1):(2*n-1+cH),cL,cW,cH)]+3*arr[I(2,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)])/32;
out[I(2,4*l,4*m+2,4*n+2,fL,fW,fH)]=(arr[I(2,2*l,2*m,2*n,cL,cW,cH)]+arr[I(2,2*l,2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(2,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)]+arr[I(2,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/16;
out[I(2,4*l,4*m+2,4*n+3,fL,fW,fH)]=(3*arr[I(2,2*l,2*m,2*n+1,cL,cW,cH)]+arr[I(2,2*l,2*m,((2*n+3)<cH)?(2*n+3):(2*n+3-cH),cL,cW,cH)]+3*arr[I(2,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)]+arr[I(2,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),((2*n+3)<cH)?(2*n+3):(2*n+3-cH),cL,cW,cH)])/32;
out[I(2,4*l,4*m+3,4*n,fL,fW,fH)]=(arr[I(2,2*l,2*m+1,2*n,cL,cW,cH)])/4;
out[I(2,4*l,4*m+3,4*n+1,fL,fW,fH)]=(arr[I(2,2*l,2*m+1,((2*n-1)>=0)?(2*n-1):(2*n-1+cH),cL,cW,cH)]+3*arr[I(2,2*l,2*m+1,2*n+1,cL,cW,cH)])/16;
out[I(2,4*l,4*m+3,4*n+2,fL,fW,fH)]=(arr[I(2,2*l,2*m+1,2*n,cL,cW,cH)]+arr[I(2,2*l,2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/8;
out[I(2,4*l,4*m+3,4*n+3,fL,fW,fH)]=(3*arr[I(2,2*l,2*m+1,2*n+1,cL,cW,cH)]+arr[I(2,2*l,2*m+1,((2*n+3)<cH)?(2*n+3):(2*n+3-cH),cL,cW,cH)])/16;
out[I(2,4*l+1,4*m,4*n,fL,fW,fH)]=(arr[I(2,2*l+1,2*m,2*n,cL,cW,cH)])/4;
out[I(2,4*l+1,4*m,4*n+1,fL,fW,fH)]=(arr[I(2,2*l+1,2*m,((2*n-1)>=0)?(2*n-1):(2*n-1+cH),cL,cW,cH)]+3*arr[I(2,2*l+1,2*m,2*n+1,cL,cW,cH)])/16;
out[I(2,4*l+1,4*m,4*n+2,fL,fW,fH)]=(arr[I(2,2*l+1,2*m,2*n,cL,cW,cH)]+arr[I(2,2*l+1,2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/8;
out[I(2,4*l+1,4*m,4*n+3,fL,fW,fH)]=(3*arr[I(2,2*l+1,2*m,2*n+1,cL,cW,cH)]+arr[I(2,2*l+1,2*m,((2*n+3)<cH)?(2*n+3):(2*n+3-cH),cL,cW,cH)])/16;
out[I(2,4*l+1,4*m+1,4*n,fL,fW,fH)]=(arr[I(2,2*l+1,2*m+1,2*n,cL,cW,cH)])/4;
out[I(2,4*l+1,4*m+1,4*n+1,fL,fW,fH)]=(arr[I(2,2*l+1,2*m+1,((2*n-1)>=0)?(2*n-1):(2*n-1+cH),cL,cW,cH)]+3*arr[I(2,2*l+1,2*m+1,2*n+1,cL,cW,cH)])/16;
out[I(2,4*l+1,4*m+1,4*n+2,fL,fW,fH)]=(arr[I(2,2*l+1,2*m+1,2*n,cL,cW,cH)]+arr[I(2,2*l+1,2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/8;
out[I(2,4*l+1,4*m+1,4*n+3,fL,fW,fH)]=(3*arr[I(2,2*l+1,2*m+1,2*n+1,cL,cW,cH)]+arr[I(2,2*l+1,2*m+1,((2*n+3)<cH)?(2*n+3):(2*n+3-cH),cL,cW,cH)])/16;
out[I(2,4*l+1,4*m+2,4*n,fL,fW,fH)]=(arr[I(2,2*l+1,2*m,2*n,cL,cW,cH)]+arr[I(2,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)])/8;
out[I(2,4*l+1,4*m+2,4*n+1,fL,fW,fH)]=(arr[I(2,2*l+1,2*m,((2*n-1)>=0)?(2*n-1):(2*n-1+cH),cL,cW,cH)]+3*arr[I(2,2*l+1,2*m,2*n+1,cL,cW,cH)]+arr[I(2,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),((2*n-1)>=0)?(2*n-1):(2*n-1+cH),cL,cW,cH)]+3*arr[I(2,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)])/32;
out[I(2,4*l+1,4*m+2,4*n+2,fL,fW,fH)]=(arr[I(2,2*l+1,2*m,2*n,cL,cW,cH)]+arr[I(2,2*l+1,2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(2,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)]+arr[I(2,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/16;
out[I(2,4*l+1,4*m+2,4*n+3,fL,fW,fH)]=(3*arr[I(2,2*l+1,2*m,2*n+1,cL,cW,cH)]+arr[I(2,2*l+1,2*m,((2*n+3)<cH)?(2*n+3):(2*n+3-cH),cL,cW,cH)]+3*arr[I(2,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)]+arr[I(2,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),((2*n+3)<cH)?(2*n+3):(2*n+3-cH),cL,cW,cH)])/32;
out[I(2,4*l+1,4*m+3,4*n,fL,fW,fH)]=(arr[I(2,2*l+1,2*m+1,2*n,cL,cW,cH)])/4;
out[I(2,4*l+1,4*m+3,4*n+1,fL,fW,fH)]=(arr[I(2,2*l+1,2*m+1,((2*n-1)>=0)?(2*n-1):(2*n-1+cH),cL,cW,cH)]+3*arr[I(2,2*l+1,2*m+1,2*n+1,cL,cW,cH)])/16;
out[I(2,4*l+1,4*m+3,4*n+2,fL,fW,fH)]=(arr[I(2,2*l+1,2*m+1,2*n,cL,cW,cH)]+arr[I(2,2*l+1,2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/8;
out[I(2,4*l+1,4*m+3,4*n+3,fL,fW,fH)]=(3*arr[I(2,2*l+1,2*m+1,2*n+1,cL,cW,cH)]+arr[I(2,2*l+1,2*m+1,((2*n+3)<cH)?(2*n+3):(2*n+3-cH),cL,cW,cH)])/16;
out[I(2,4*l+2,4*m,4*n,fL,fW,fH)]=(arr[I(2,2*l,2*m,2*n,cL,cW,cH)]+arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,2*n,cL,cW,cH)])/8;
out[I(2,4*l+2,4*m,4*n+1,fL,fW,fH)]=(arr[I(2,2*l,2*m,((2*n-1)>=0)?(2*n-1):(2*n-1+cH),cL,cW,cH)]+3*arr[I(2,2*l,2*m,2*n+1,cL,cW,cH)]+arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,((2*n-1)>=0)?(2*n-1):(2*n-1+cH),cL,cW,cH)]+3*arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,2*n+1,cL,cW,cH)])/32;
out[I(2,4*l+2,4*m,4*n+2,fL,fW,fH)]=(arr[I(2,2*l,2*m,2*n,cL,cW,cH)]+arr[I(2,2*l,2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,2*n,cL,cW,cH)]+arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/16;
out[I(2,4*l+2,4*m,4*n+3,fL,fW,fH)]=(3*arr[I(2,2*l,2*m,2*n+1,cL,cW,cH)]+arr[I(2,2*l,2*m,((2*n+3)<cH)?(2*n+3):(2*n+3-cH),cL,cW,cH)]+3*arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,2*n+1,cL,cW,cH)]+arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,((2*n+3)<cH)?(2*n+3):(2*n+3-cH),cL,cW,cH)])/32;
out[I(2,4*l+2,4*m+1,4*n,fL,fW,fH)]=(arr[I(2,2*l,2*m+1,2*n,cL,cW,cH)]+arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,2*n,cL,cW,cH)])/8;
out[I(2,4*l+2,4*m+1,4*n+1,fL,fW,fH)]=(arr[I(2,2*l,2*m+1,((2*n-1)>=0)?(2*n-1):(2*n-1+cH),cL,cW,cH)]+3*arr[I(2,2*l,2*m+1,2*n+1,cL,cW,cH)]+arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,((2*n-1)>=0)?(2*n-1):(2*n-1+cH),cL,cW,cH)]+3*arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,2*n+1,cL,cW,cH)])/32;
out[I(2,4*l+2,4*m+1,4*n+2,fL,fW,fH)]=(arr[I(2,2*l,2*m+1,2*n,cL,cW,cH)]+arr[I(2,2*l,2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,2*n,cL,cW,cH)]+arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/16;
out[I(2,4*l+2,4*m+1,4*n+3,fL,fW,fH)]=(3*arr[I(2,2*l,2*m+1,2*n+1,cL,cW,cH)]+arr[I(2,2*l,2*m+1,((2*n+3)<cH)?(2*n+3):(2*n+3-cH),cL,cW,cH)]+3*arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,2*n+1,cL,cW,cH)]+arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,((2*n+3)<cH)?(2*n+3):(2*n+3-cH),cL,cW,cH)])/32;
out[I(2,4*l+2,4*m+2,4*n,fL,fW,fH)]=(arr[I(2,2*l,2*m,2*n,cL,cW,cH)]+arr[I(2,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)]+arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,2*n,cL,cW,cH)]+arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)])/16;
out[I(2,4*l+2,4*m+2,4*n+1,fL,fW,fH)]=(arr[I(2,2*l,2*m,((2*n-1)>=0)?(2*n-1):(2*n-1+cH),cL,cW,cH)]+3*arr[I(2,2*l,2*m,2*n+1,cL,cW,cH)]+arr[I(2,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),((2*n-1)>=0)?(2*n-1):(2*n-1+cH),cL,cW,cH)]+3*arr[I(2,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)]+arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,((2*n-1)>=0)?(2*n-1):(2*n-1+cH),cL,cW,cH)]+3*arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,2*n+1,cL,cW,cH)]+arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),((2*m+2)<cW)?(2*m+2):(2*m+2-cW),((2*n-1)>=0)?(2*n-1):(2*n-1+cH),cL,cW,cH)]+3*arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)])/64;
out[I(2,4*l+2,4*m+2,4*n+2,fL,fW,fH)]=(arr[I(2,2*l,2*m,2*n,cL,cW,cH)]+arr[I(2,2*l,2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(2,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)]+arr[I(2,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,2*n,cL,cW,cH)]+arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)]+arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),((2*m+2)<cW)?(2*m+2):(2*m+2-cW),((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/32;
out[I(2,4*l+2,4*m+2,4*n+3,fL,fW,fH)]=(3*arr[I(2,2*l,2*m,2*n+1,cL,cW,cH)]+arr[I(2,2*l,2*m,((2*n+3)<cH)?(2*n+3):(2*n+3-cH),cL,cW,cH)]+3*arr[I(2,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)]+arr[I(2,2*l,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),((2*n+3)<cH)?(2*n+3):(2*n+3-cH),cL,cW,cH)]+3*arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,2*n+1,cL,cW,cH)]+arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m,((2*n+3)<cH)?(2*n+3):(2*n+3-cH),cL,cW,cH)]+3*arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)]+arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),((2*m+2)<cW)?(2*m+2):(2*m+2-cW),((2*n+3)<cH)?(2*n+3):(2*n+3-cH),cL,cW,cH)])/64;
out[I(2,4*l+2,4*m+3,4*n,fL,fW,fH)]=(arr[I(2,2*l,2*m+1,2*n,cL,cW,cH)]+arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,2*n,cL,cW,cH)])/8;
out[I(2,4*l+2,4*m+3,4*n+1,fL,fW,fH)]=(arr[I(2,2*l,2*m+1,((2*n-1)>=0)?(2*n-1):(2*n-1+cH),cL,cW,cH)]+3*arr[I(2,2*l,2*m+1,2*n+1,cL,cW,cH)]+arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,((2*n-1)>=0)?(2*n-1):(2*n-1+cH),cL,cW,cH)]+3*arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,2*n+1,cL,cW,cH)])/32;
out[I(2,4*l+2,4*m+3,4*n+2,fL,fW,fH)]=(arr[I(2,2*l,2*m+1,2*n,cL,cW,cH)]+arr[I(2,2*l,2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,2*n,cL,cW,cH)]+arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/16;
out[I(2,4*l+2,4*m+3,4*n+3,fL,fW,fH)]=(3*arr[I(2,2*l,2*m+1,2*n+1,cL,cW,cH)]+arr[I(2,2*l,2*m+1,((2*n+3)<cH)?(2*n+3):(2*n+3-cH),cL,cW,cH)]+3*arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,2*n+1,cL,cW,cH)]+arr[I(2,((2*l+2)<cL)?(2*l+2):(2*l+2-cL),2*m+1,((2*n+3)<cH)?(2*n+3):(2*n+3-cH),cL,cW,cH)])/32;
out[I(2,4*l+3,4*m,4*n,fL,fW,fH)]=(arr[I(2,2*l+1,2*m,2*n,cL,cW,cH)])/4;
out[I(2,4*l+3,4*m,4*n+1,fL,fW,fH)]=(arr[I(2,2*l+1,2*m,((2*n-1)>=0)?(2*n-1):(2*n-1+cH),cL,cW,cH)]+3*arr[I(2,2*l+1,2*m,2*n+1,cL,cW,cH)])/16;
out[I(2,4*l+3,4*m,4*n+2,fL,fW,fH)]=(arr[I(2,2*l+1,2*m,2*n,cL,cW,cH)]+arr[I(2,2*l+1,2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/8;
out[I(2,4*l+3,4*m,4*n+3,fL,fW,fH)]=(3*arr[I(2,2*l+1,2*m,2*n+1,cL,cW,cH)]+arr[I(2,2*l+1,2*m,((2*n+3)<cH)?(2*n+3):(2*n+3-cH),cL,cW,cH)])/16;
out[I(2,4*l+3,4*m+1,4*n,fL,fW,fH)]=(arr[I(2,2*l+1,2*m+1,2*n,cL,cW,cH)])/4;
out[I(2,4*l+3,4*m+1,4*n+1,fL,fW,fH)]=(arr[I(2,2*l+1,2*m+1,((2*n-1)>=0)?(2*n-1):(2*n-1+cH),cL,cW,cH)]+3*arr[I(2,2*l+1,2*m+1,2*n+1,cL,cW,cH)])/16;
out[I(2,4*l+3,4*m+1,4*n+2,fL,fW,fH)]=(arr[I(2,2*l+1,2*m+1,2*n,cL,cW,cH)]+arr[I(2,2*l+1,2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/8;
out[I(2,4*l+3,4*m+1,4*n+3,fL,fW,fH)]=(3*arr[I(2,2*l+1,2*m+1,2*n+1,cL,cW,cH)]+arr[I(2,2*l+1,2*m+1,((2*n+3)<cH)?(2*n+3):(2*n+3-cH),cL,cW,cH)])/16;
out[I(2,4*l+3,4*m+2,4*n,fL,fW,fH)]=(arr[I(2,2*l+1,2*m,2*n,cL,cW,cH)]+arr[I(2,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)])/8;
out[I(2,4*l+3,4*m+2,4*n+1,fL,fW,fH)]=(arr[I(2,2*l+1,2*m,((2*n-1)>=0)?(2*n-1):(2*n-1+cH),cL,cW,cH)]+3*arr[I(2,2*l+1,2*m,2*n+1,cL,cW,cH)]+arr[I(2,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),((2*n-1)>=0)?(2*n-1):(2*n-1+cH),cL,cW,cH)]+3*arr[I(2,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)])/32;
out[I(2,4*l+3,4*m+2,4*n+2,fL,fW,fH)]=(arr[I(2,2*l+1,2*m,2*n,cL,cW,cH)]+arr[I(2,2*l+1,2*m,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)]+arr[I(2,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n,cL,cW,cH)]+arr[I(2,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/16;
out[I(2,4*l+3,4*m+2,4*n+3,fL,fW,fH)]=(3*arr[I(2,2*l+1,2*m,2*n+1,cL,cW,cH)]+arr[I(2,2*l+1,2*m,((2*n+3)<cH)?(2*n+3):(2*n+3-cH),cL,cW,cH)]+3*arr[I(2,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),2*n+1,cL,cW,cH)]+arr[I(2,2*l+1,((2*m+2)<cW)?(2*m+2):(2*m+2-cW),((2*n+3)<cH)?(2*n+3):(2*n+3-cH),cL,cW,cH)])/32;
out[I(2,4*l+3,4*m+3,4*n,fL,fW,fH)]=(arr[I(2,2*l+1,2*m+1,2*n,cL,cW,cH)])/4;
out[I(2,4*l+3,4*m+3,4*n+1,fL,fW,fH)]=(arr[I(2,2*l+1,2*m+1,((2*n-1)>=0)?(2*n-1):(2*n-1+cH),cL,cW,cH)]+3*arr[I(2,2*l+1,2*m+1,2*n+1,cL,cW,cH)])/16;
out[I(2,4*l+3,4*m+3,4*n+2,fL,fW,fH)]=(arr[I(2,2*l+1,2*m+1,2*n,cL,cW,cH)]+arr[I(2,2*l+1,2*m+1,((2*n+2)<cH)?(2*n+2):(2*n+2-cH),cL,cW,cH)])/8;
out[I(2,4*l+3,4*m+3,4*n+3,fL,fW,fH)]=(3*arr[I(2,2*l+1,2*m+1,2*n+1,cL,cW,cH)]+arr[I(2,2*l+1,2*m+1,((2*n+3)<cH)?(2*n+3):(2*n+3-cH),cL,cW,cH)])/16;
}
}
}
}


int main(int argc, char *argv[])
    {
      int outdatasize;
      //values to report: L,W,H,mu,eta,nu,h,t,VnormL2,VnormLinf,gradVnormL2,pct1,pct2,divVLinf
      double *outdata = new double[10001*16];
      double VnormL2,maximumV,pct1,pct2,divVLinf,gradnorm, gradsymm, gradskew;
      double maximumVC, gradnormC, gradsymmC, gradskewC;
      double *V = new double[ARR_SIZE]; //Velocity on fine grid
	    double *Vc = new double[COARSE_ARR_SIZE]; // coarse grid velocity
	    double *tempC = new double[COARSE_ARR_SIZE]; //coarse grid temp arr.
      double *Laplace_V = new double[ARR_SIZE];
      double *Vfvf_res = new double[ARR_SIZE];
      double *temp = new double[ARR_SIZE];
      double *R= new double[ARR_SIZE];
      double *Total= new double[ARR_SIZE];
      double *scal= new double[L*W*H];
      double *pressure_old= new double[L*W*H];
      double *solution_old= new double[L*W*H];
      double id, jd, kd;
      double threshold1,threshold2;
      double dt,t;
      double nu=0.01;
      double h=0.01;
      double eta=5.0;
      double phys_length = 6.283185307179586;
      double coeff=MU;
      int pr=1;
      int freq=1;
      int basetime=0;
      auto result = parse(argc, argv);
      auto arguments = result.arguments();
      bool silent=false;


      if (result.count("c")) coeff = result["c"].as<double>();
      if (result.count("l")) phys_length = result["l"].as<double>();
      if (result.count("nu")) nu = result["nu"].as<double>();
      if (result.count("eta")) eta = result["eta"].as<double>();
      if (result.count("silent")) silent = result["silent"].as<bool>();
      if (result.count("print")) pr=result["print"].as<int>();
      if (result.count("freq")) freq=result["freq"].as<int>();
      if (!silent) {
        std::cout<<"Using values : phys_length "<<phys_length<<", nu "<<nu<<", eta "<<eta<<", smoothing frequency:"<<freq<<std::endl;
      }

      h = phys_length/L;

//Testing
char plane[] ="xyz";
for(int i=0; i<COARSE_ARR_SIZE; i++) Vc[i]=0.0;

for(int i=0; i<(L*W*H); i++) scal[i]=0.0;
for(int i=0; i<(L*W*H); i++) pressure_old[i]=0.0;
for(int i=0; i<(L*W*H); i++) solution_old[i]=0.0;

/*
Vc[I(0,0,0,1,LL,WW,HH)]=1.0; //Put simplex X at 0,0,0 in the coarse grid
Vc[I(0,1,1,0,LL,WW,HH)]=1.0; //Put simplex X at 0,0,0 in the coarse grid

coarse2fine(Vc,V,WW,LL,HH);
cout<<"Testing. coarse2fine of x_(0,0,0) is:"<<endl;
for(int l=0; l<3; l++){
	for(int i=0; i<L; i++) {
		for(int j=0; j<W; j++) {
			for(int k=0; k<H; k++){
                if (V[I(l,i,j,k)]!=0.0) {
                    cout<<plane[l]<<"_("<<i<<" ,"<<j<<" ,"<<k<<") :"<<V[I(l,i,j,k)]<<endl;
                }
            }
        }
    }
}
*/




      outdatasize=0;
      t=0.0;
      dt=0.01;
      if (2.0*nu*dt>1.0) dt=1.0/(2.0*nu); //Needed for stability at high viscosity (low Reynolds number).
      cout<<"dt: "<<dt<<endl;
      for (int i=0; i<COARSE_ARR_SIZE; i++) tempC[i]=0.0;

      // initialize a divergence free vector field.

      for(int i=0; i<LL; i++) {
        for(int j=0; j<WW; j++) {
          for(int k=0; k<HH; k++) {
            id= static_cast<double> (i-LL/2);
            jd= static_cast<double> (j-WW/2);
            kd= static_cast<double> (k-HH/2);
            tempC[I(X,i,j,k,LL,WW,HH)] = coeff*exp(-eta*h*h*((id*id)+(jd*jd)+(kd*kd)));
          }
        }
        // cout<< temp[I(X,i,16,16)] << endl;
      }

      bdC(tempC,Vc);   // This makes Vc a divergence free vector field since div(curl (vector)) = 0.



	  coarse2fine(Vc,V,WW,LL,HH); //make V the fine version of Vc

cout<<"L_inf norm of div V is: "<<maxnormscal(scal)<<endl;



      cout<<"Lattice Size: ("<< L << "," << W << "," << "," << H << ")" << endl;
      cout<<"h: " << h << endl;
      cout<<"Viscosity: " << nu << endl;

      threshold1=0.01;
      for (int rep=0; rep<10001; rep++) {


        maximumV=maxnorm(V);
        if (maximumV!=maximumV) { //if maximumV is nan.
          cout<<"Blowup encountered."<<endl;
          break;
        }
        if (dt> 0.05/maximumV) {
          dt=0.05/maximumV;
          if (dt<0.0000001) {
            cout<<"Blowup encountered."<<endl;
            break;
          }
        }

        if (maximumV<1e-6) {
          cout<<"Decays to zero."<<endl;
          break;
        }


        if(rep%pr==0) {

          VnormL2=norm(V);
          gradient_norm(V,&gradnorm,&gradsymm,&gradskew);
          threshold2= 0.25*maximumV;
          pct_vectors(V,threshold1,threshold2,&pct1,&pct2);
          bd10(V,scal);
          divVLinf=maxnormscal(scal);

          //print
          if (!silent) {
            cout<<"h,nu: "<<h <<" , "<<nu<<endl;
            cout << "time: " <<t <<endl;
            cout << "step: " <<rep <<endl;
            cout << "velocity L-inf norm: "<< maximumV <<endl;

            cout << "velocity L2 norm: "<< VnormL2 <<endl;

            cout << "grad(V) L2 norm: "<<gradnorm <<endl;
            cout << "L2 norm of symmetric part of grad(V): "<<gradsymm<<endl;
            cout << "L2 norm of skew-symmetric part of grad(V): "<<gradskew<<endl;


            cout << "Percent of velocities with norm greater than "<<threshold1<<" :"<<pct1<<endl;
            cout << "Percent of velocities with norm greater than "<<threshold2<<" :"<<pct2<<endl;


            cout << "divergence L-inf norm: "<< divVLinf <<endl;
            cout << endl;
          }

          //values to report: L,W,H,mu,eta,nu,h,t,VNormL2,VnormLinf,gradnorm,gradsymm,gradskew,pct1,pct2,divVLinf
          outdata[16*outdatasize+0]= L;
          outdata[16*outdatasize+1]= W;
          outdata[16*outdatasize+2]= H;
          outdata[16*outdatasize+3]= MU;
          outdata[16*outdatasize+4]= eta;
          outdata[16*outdatasize+5]= nu;
          outdata[16*outdatasize+6]= h;
          outdata[16*outdatasize+7]= t;
          outdata[16*outdatasize+8]= VnormL2;
          outdata[16*outdatasize+9]= maximumV;
          outdata[16*outdatasize+10]= gradnorm;
          outdata[16*outdatasize+11]= gradsymm;
          outdata[16*outdatasize+12]= gradskew;
          outdata[16*outdatasize+13]= pct1;
          outdata[16*outdatasize+14]= pct2;
          outdata[16*outdatasize+15]= divVLinf;
          outdatasize++;
        }

		// Runge Kutta
        nav_stoke(h,nu,V, Vfvf_res, R, scal, pressure_old, solution_old); //scal is pressure.
        for(int i = 0; i<ARR_SIZE; i++) temp[i] = V[i] + R[i]*dt*0.5;
        for(int i = 0; i<ARR_SIZE; i++) Total[i]=R[i]; //separate loop to make compiler to use SIMD
        nav_stoke(h,nu,temp, Vfvf_res, R, scal, pressure_old, solution_old);
        for(int i = 0; i<ARR_SIZE; i++) temp[i] = V[i] + R[i]*dt*0.5;
        for(int i = 0; i<ARR_SIZE; i++) Total[i]+=2*R[i];
        nav_stoke(h,nu,temp, Vfvf_res,  R, scal, pressure_old,  solution_old);
        for(int i = 0; i<ARR_SIZE; i++) temp[i] = V[i] + R[i]*dt;
        for(int i = 0; i<ARR_SIZE; i++) Total[i]+=2*R[i];
        nav_stoke(h,nu,temp, Vfvf_res, R, scal, pressure_old,  solution_old);
        for(int i = 0; i<ARR_SIZE; i++) Total[i]+=R[i];
        for(int i = 0; i<ARR_SIZE; i++) V[i] += Total[i]*(dt/6);
        t+= dt;

	if (rep%freq==0) {
        fine2coarse(V,Vc,L,W,H);
		coarse2fine(Vc,V,LL,WW,HH);}

   // print out maxnormC and gradient every 50th step
   if (rep%10==0) {
        maximumVC = maxnorm(Vc);
        gradient_norm(Vc,&gradnormC,&gradsymmC,&gradskewC);

 	cout << "Coarse scale velocity L-inf norm: "<< maximumVC <<endl;

   	cout << "Coarse scale grad(V) L2 norm: "<<gradnormC <<endl;
    	cout << "Coarse scale L2 norm of symmetric part of grad(V): "<<gradsymmC <<endl;
    	cout << "Coarse scale L2 norm of skew-symmetric part of grad(V): "<<gradskewC <<endl<<endl;
      }

      /* adding code for writing visit files
      char outfile[100];
      output_file(rep, outfile);
      visit_plot(outfile, L, W, H, pressure_old, V);
      */

      }


      //write2file(outdata,outdatasize,h,nu,eta);

      return 1;
    }
