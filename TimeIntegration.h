#ifndef TimeIntegration_h
#define TimeIntegration_h

#include "TimeIntegration.h"

class Field;

class TimeIntegration
{
  private:
    double dt;
    int L, W, H;
    Field *vobj;
    Field *sobj;    

  public:
    TimeIntegration(Field*, Field*, double, char); // constructor
    ~TimeIntegration(); // destructor
    void dVfvf(double out[]); // non-linear term
    void time_stepping(char);
    void runge_kutta4();
};

#endif
