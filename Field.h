#ifndef Field_h
#define Field_h

#include "Grid.h"

class Grid;

class Field
{
  public:
    Field(int, Grid *, int); // constructor
    ~Field(); // destructor
    int I(int, int, int, int);
    double* arr; 
    void initialize();  
    void bd(double out[]); //curl
    void update_guard(); // apply boundary on guard cells
  private:
    int BTYPE;
    int NSIZE; // 1 or 3 depending on scalar or vector
    int ARR_SIZE;
    Grid *obj;
};

#endif
