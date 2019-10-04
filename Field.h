#ifndef Field_h
#define Field_h

#include "Grid.h"

class Grid;

class Field
{

  public:
    void initialize();
    Field(int, Grid *); // constructor
    ~Field(); // destructor
    int I(int, int, int, int, Grid &);
  private:
    int NSIZE; // 1 or 3 depending on scalar or vector
    Grid *obj;
    int ARR_SIZE;
};

#endif
