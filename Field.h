#ifndef Field_h
#define Field_h

#include "Grid.h"

class Grid;

class Field
{
  public:
    Field(int, Grid *); // constructor
    ~Field(); // destructor
    int I(int, int, int, int);
    double* arr; 
    void initialize();  
    void bd(double out[]);
  private:
    int NSIZE; // 1 or 3 depending on scalar or vector
    int ARR_SIZE;
    Grid *obj;
};

#endif
