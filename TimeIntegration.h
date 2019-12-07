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
    int *coords;
    int* igrid; // interior grid size on each proc    
    int* ggrid; // global interior grid size on the whole domain

    void dVfvf(double*);
    void proj(double*);
    void nav_stoke(double*); 

  public:
    TimeIntegration(Field*, Field*, double, MPI_Comm, int*, int*, int*); // constructor
    ~TimeIntegration(); // destructor
    MPI_Comm mpi_comm; //MPI communictaor needed in HYPRE
    int I(int w, int i, int j, int k, int);
    void bd10(double*, double*);
    void time_stepping(char);
    void runge_kutta4();
    void forward_euler();
    void pressure_solve(double*, double*);
};

#endif
