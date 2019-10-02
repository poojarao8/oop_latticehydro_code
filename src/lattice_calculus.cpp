#include "main.h"

// Boundary operator of a two chain. This is the curl operator. Since the * operator between 2-chain and 1-chain don't do anything, we are really using this as *del.
// // Rather than making the function return an array pointer, it stores the result in the out[] array. This was done to avoid having to make new arrays and having to deal with garbage collection.
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
//
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

