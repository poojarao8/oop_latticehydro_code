#ifndef Field_h
#define Field_h


class Grid;

class Field
{
  private:
    // grid info 
    int L, W, H;
    int GL, GW, GH;
    int NSIZE; // 1 or 3 depending on scalar or vector

    // bdry update and indexing function
    void periodic_bdry();
    void taylor_green_vortex(int*, int*);
  public:
    Field(int, Grid *); // constructor
    ~Field(); // destructor
    Grid *obj;
    int ARR_SIZE;
    double* arr; // for velocity or pressure 
    int I(int, int, int, int); // indexing function
    void initialize(int, int*, int*);  
    void update_bdry(int); // apply boundary on guard cells
    // lattice calculus functions
    void bd(double*); //curl
    void bd10(double*); // divergence
    void bd10_NL(double*, double*); // divergence
    void d01(double*); // gradient
    void laplacian(double*);


};

#endif
