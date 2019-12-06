#ifndef TimeIntegration_h
#define TimeIntegration_h

class Field;

class TimeIntegration
{
  private:
    double T;
    double dt;
    int L, W, H;
    Field *vel;
    Field *pres;    
    // arrays for storing field values

    void dVfvf(double*);
    void proj(double*);
    void nav_stoke(double*); 

  public:
    TimeIntegration(Field*, Field*, double); // constructor
    ~TimeIntegration(); // destructor
    int I(int w, int i, int j, int k, int);
    void bd10(double*, double*);
    void time_stepping(char);
    void runge_kutta4();
    void forward_euler();
    void pressure_solve(MPI_Comm, double*, double*);
};

#endif
