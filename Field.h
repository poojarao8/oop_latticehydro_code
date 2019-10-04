#include "Grid.h"

class Field
{

  public:
    void initialize(Grid &);
    Field(int, Grid &); // constructor
    ~Field(); // destructor

  private:
    int ARR_SIZE;
};

