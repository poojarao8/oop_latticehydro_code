#ifndef Field_h
#define Field_h

#include "Grid.h"

class Grid;

class Field
{
  private:
    // grid info 
    int L, W, H;
    int GL, GW, GH;
    int NSIZE; // 1 or 3 depending on scalar or vector
    int ARR_SIZE;

    // bdry update and indexing function
    void periodic_bdry();

    // lattice calculus functions
    void bd(double*); //curl
    void bd10(double*); // divergence
    void d01(double*); // gradient
    void laplacian(double*); 

  public:
    Field(int, Grid *); // constructor
    ~Field(); // destructor
    Grid *obj;
    double* arr; 
    int I(int, int, int, int); // indexing function
    void initialize();  
    void update_bdry(char); // apply boundary on guard cells
};

#endif
