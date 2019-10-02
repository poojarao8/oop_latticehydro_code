#include "main.h"

// L2 norm of a vector
double norm(double V[], int arr_size)
{
  double tot=0.0;
  for (int i=0; i<arr_size; i++) tot+= V[i]*V[i];
  return sqrt(tot);
}

// maxnorm of a vector
double maxnorm(double V[], int arr_size) 
{
  double maxval=0.0;
  double normsqval;
  for (int i=0; i<(arr_size-2); i+=3) 
  {
    normsqval = V[i]*V[i]+V[i+1]*V[i+1]+V[i+2]*V[i+2];
    if(normsqval>maxval) maxval=normsqval;
  }
  return sqrt(maxval);
}


// L2 norm of the gradient
void gradient_norm(double V[], int l, int w, int h, double *gradnorm, double *gradsymm, double *gradskew) 
{
  double Dxx,Dxy,Dxz,Dyx,Dyy,Dyz,Dzx,Dzy,Dzz;
  *gradsymm = 0.0;
  *gradskew = 0.0;
  for(int i=0; i<l; i++) {
    for(int j=0; j<w; j++) {
      for(int k=0; k<h; k++) {
        Dxx = V[I(X,(i!=(l-1)) ? (i+1) : 0,j,k)] - V[I(X,(i!=0) ? (i-1) : (l-1),j,k)];
        Dxy = V[I(X,i,(j!=(w-1)) ? (j+1) : 0,k)] - V[I(X,i,(j!=0) ? (j-1) : (w-1),k)];
        Dxz = V[I(X,i,j,(k!=(h-1)) ? (k+1) : 0)] - V[I(X,i,j,(k!=0) ? (k-1) : (h-1))];
        Dyx = V[I(Y,(i!=(l-1)) ? (i+1) : 0,j,k)] - V[I(Y,(i!=0) ? (i-1) : (l-1),j,k)];
        Dyy = V[I(Y,i,(j!=(w-1)) ? (j+1) : 0,k)] - V[I(Y,i,(j!=0) ? (j-1) : (w-1),k)];
        Dyz = V[I(Y,i,j,(k!=(h-1)) ? (k+1) : 0)] - V[I(Y,i,j,(k!=0) ? (k-1) : (h-1))];
        Dzx = V[I(Z,(i!=(l-1)) ? (i+1) : 0,j,k)] - V[I(Z,(i!=0) ? (i-1) : (l-1),j,k)];
        Dzy = V[I(Z,i,(j!=(w-1)) ? (j+1) : 0,k)] - V[I(Z,i,(j!=0) ? (j-1) : (w-1),k)];
        Dzz = V[I(Z,i,j,(k!=(h-1)) ? (k+1) : 0)] - V[I(Z,i,j,(k!=0) ? (k-1) : (h-1))];
        *gradsymm += Dxx*Dxx + Dyy*Dyy + Dzz*Dzz;
        *gradskew += Dxy*Dxy+Dxz*Dxz+Dyx*Dyx+Dyz*Dyz+Dzx*Dzx+Dzy*Dzy;
      }
    }
  }
  *gradnorm = sqrt(*gradsymm+*gradskew);
  *gradsymm = sqrt(*gradsymm);
  *gradskew = sqrt(*gradskew);

}

//PRAO: added to plot grad v
// L2 norm of the gradient
void gradient_norm_domain(double V[], int l, int w, int h, double *gradnormdomain)
{
  double Dxx,Dxy,Dxz,Dyx,Dyy,Dyz,Dzx,Dzy,Dzz;
  for(int i=0; i<l; i++) {
    for(int j=0; j<w; j++) {
      for(int k=0; k<h; k++) {
        Dxx = V[I(X,(i!=(l-1)) ? (i+1) : 0,j,k)] - V[I(X,(i!=0) ? (i-1) : (l-1),j,k)];
        Dxy = V[I(X,i,(j!=(w-1)) ? (j+1) : 0,k)] - V[I(X,i,(j!=0) ? (j-1) : (w-1),k)];
        Dxz = V[I(X,i,j,(k!=(h-1)) ? (k+1) : 0)] - V[I(X,i,j,(k!=0) ? (k-1) : (h-1))];
        Dyx = V[I(Y,(i!=(l-1)) ? (i+1) : 0,j,k)] - V[I(Y,(i!=0) ? (i-1) : (l-1),j,k)];
        Dyy = V[I(Y,i,(j!=(w-1)) ? (j+1) : 0,k)] - V[I(Y,i,(j!=0) ? (j-1) : (w-1),k)];
        Dyz = V[I(Y,i,j,(k!=(h-1)) ? (k+1) : 0)] - V[I(Y,i,j,(k!=0) ? (k-1) : (h-1))];
        Dzx = V[I(Z,(i!=(l-1)) ? (i+1) : 0,j,k)] - V[I(Z,(i!=0) ? (i-1) : (l-1),j,k)];
        Dzy = V[I(Z,i,(j!=(w-1)) ? (j+1) : 0,k)] - V[I(Z,i,(j!=0) ? (j-1) : (w-1),k)];
        Dzz = V[I(Z,i,j,(k!=(h-1)) ? (k+1) : 0)] - V[I(Z,i,j,(k!=0) ? (k-1) : (h-1))];
        gradnormdomain[J(i,j,k,l,w,h)] = Dxx*Dxx + Dyy*Dyy + Dzz*Dzz +
          Dxy*Dxy+Dxz*Dxz+Dyx*Dyx+Dyz*Dyz+Dzx*Dzx+Dzy*Dzy;
      }
    }
  }
}


double grad_ratio_norm(double V[], int l, int w, int h) 
{
  double Rxn,Ryn,Rzn,Rxd,Ryd,Rzd;
  double Vx0,Vy0,Vz0;
  double Vxx1,Vxy1,Vxz1,Vyx1,Vyy1,Vyz1,Vzx1,Vzy1,Vzz1;
  double gradratio, gradratio_max=0.0;

  for(int i=0; i<l; i++) {
    for(int j=0; j<w; j++) {
      for(int k=0; k<h; k++) {
        Vx0 = V[I(X,i,j,k)];
        Vy0 = V[I(Y,i,j,k)];
        Vz0 = V[I(Z,i,j,k)];

        Vxx1 = V[I(X,(i!=(l-1)) ? (i+1) : 0,j,k)];
        Vxy1 = V[I(X,i,(j!=(w-1)) ? (j+1) : 0,k)];
        Vxz1 = V[I(X,i,j,(k!=(h-1)) ? (k+1) : 0)];
        Vyx1 = V[I(Y,(i!=(l-1)) ? (i+1) : 0,j,k)];
        Vyy1 = V[I(Y,i,(j!=(w-1)) ? (j+1) : 0,k)];
        Vyz1 = V[I(Y,i,j,(k!=(h-1)) ? (k+1) : 0)];
        Vzx1 = V[I(Z,(i!=(l-1)) ? (i+1) : 0,j,k)];
        Vzy1 = V[I(Z,i,(j!=(w-1)) ? (j+1) : 0,k)];
        Vzz1 = V[I(Z,i,j,(k!=(h-1)) ? (k+1) : 0)];

        Rxn = (Vxx1-Vx0)*(Vxx1-Vx0)+(Vyx1-Vy0)*(Vyx1-Vy0)+(Vzx1-Vz0)*(Vzx1-Vz0);
        Rxd = Vxx1*Vxx1+Vyx1*Vyx1+Vzx1*Vzx1+Vx0*Vx0+Vy0*Vy0+Vz0*Vz0;
        Ryn = (Vxy1-Vx0)*(Vxy1-Vx0)+(Vyy1-Vy0)*(Vyy1-Vy0)+(Vzy1-Vz0)*(Vzy1-Vz0);
        Ryd = Vxy1*Vxy1+Vyy1*Vyy1+Vzy1*Vzy1+Vx0*Vx0+Vy0*Vy0+Vz0*Vz0;
        Rzn = (Vxz1-Vx0)*(Vxz1-Vx0)+(Vyz1-Vy0)*(Vyz1-Vy0)+(Vzz1-Vz0)*(Vzz1-Vz0);
        Rzd = Vxz1*Vxz1+Vyz1*Vyz1+Vzz1*Vzz1+Vx0*Vx0+Vy0*Vy0+Vz0*Vz0;

        gradratio = ((Rxn!=0) ? Rxn/Rxd : 0)+((Ryn!=0) ? Ryn/Ryd : 0)+((Rzn!=0) ? Rzn/Rzd : 0);
        if (gradratio_max<gradratio) gradratio_max=gradratio;
      }
    }
  }
  return sqrt(gradratio_max);
}


// Percentage of vectors whose norm is above given thresholds.
void pct_vectors(double V[], int arr_size, int l, int w, int h, double threshold1, double threshold2, double *pct1, double *pct2) 
{
  unsigned int count1=0;
  unsigned int count2=0;
  double normsqval;
  double th1=threshold1*threshold1;
  double th2=threshold2*threshold2;

  for (int i=0; i<(arr_size-2); i+=3) 
  {
    normsqval = V[i]*V[i]+V[i+1]*V[i+1]+V[i+2]*V[i+2];
    if(normsqval>th1) 
    {
      count1++;
      if(normsqval>th2) 
      {
        count2++;
      }
    }
  }
  *pct1 = (static_cast<double> (count1))/(static_cast<double> (l*w*h));
  *pct2 = (static_cast<double> (count2))/(static_cast<double> (l*w*h));
}

// maxnorm of a scalar field
double maxnormscal(double scal[], int scal_size) 
{
  double maxval=0.0;
  double absval;

  for (int i=0; i<scal_size; i++) 
  {
    absval = abs(scal[i]);
    if(absval>maxval) maxval=absval;
  }
  return maxval;
}

// max element of an array
double maxElement(double* arr, int arr_size)
{
  double largest = arr[0];
  for(int i = 1;i < arr_size; i++) 
    if(largest < arr[i]) largest = arr[i];
}
