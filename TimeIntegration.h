#ifndef TimeIntegration_h
#define TimeIntegration_h

#include "TimeIntegration.h"

class Field;

class TimeIntegration
{
  private:
    double T;
    double dt;
    int L, W, H;
    //constexpr int LL, WW, HH;
    Field *vel;
    Field *pres;    
    // arrays for storing field values
    double* Vfvf_res; // for non-linear term
    double* pressure_old; // pressure at previous iteration 
    double* pressure_lap; // pressure laplacian at previous iteration
    double* pressure_lap_old; // div of pressure from prev itr

    void dVfvf(double*);
    void proj(double*);
    void nav_stoke(double*); 

  public:
    TimeIntegration(Field*, Field*, double); // constructor
    ~TimeIntegration(); // destructor
    void time_stepping(char);
    void runge_kutta4();
    void forward_euler();
    void pressure_solve(MPI_Comm, double*, double* );
};

#endif
