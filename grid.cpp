#include <iostream>


/*
Grid Class
*/


using namespace std;

class Grid
{
  public:
   int L;
   int W;
   int H;
   int VEC_ARR_SIZE;
   int SCL_ARR_SIZE;
   Grid(int, int, int); // constructor declaration
    ~Grid(); // destructor declaration
   
  private:
    double dx;
};

// member function declarations including the constructor

Grid::Grid(int l, int w, int h)
{
  L = l;
  W = w;
  H = h;
  dx = 2*3.14/L;
  SCL_ARR_SIZE = L*W*H;
  VEC_ARR_SIZE = 3*L*W*H;

  cout << "Object is being created" << endl;
  cout << "dx = " << dx << endl;  
  cout << "SCL_ARR_SIZE = " << SCL_ARR_SIZE << endl;
  cout << "VEC_ARR_SIZE = " << VEC_ARR_SIZE << endl;
}

Grid::~Grid(void)
{
  cout << "Grid object is being deleted" << endl;
}

class Field
{
 
  // how to specify that its a vel or pressure data or vec or scalar data?
  public:
    void initialize(Grid& obj);       
    Field(); //constructor
    ~Field(); //destructor

};

Field::Field(void)
{
  cout << "Grid object is being created" << endl;
}

Field::~Field(void)
{
  cout << "Field object is being deleted" << endl;
}

void Field::initialize(Grid& obj)
{
  for(int i=0; i<obj.L; i++) {
    for(int j=0; j<obj.W; j++) {
      for(int k=0; k<obj.H; k++) {

        double id = static_cast<double> (i-obj.L/2);
        double jd = static_cast<double> (j-obj.W/2);
        double kd = static_cast<double> (k-obj.H/2);

       //vel[I(X,i,j,k,LL,WW,HH)] = coef*exp(-eta*(h/2.0)*(h/2.0)*((id*id)+(jd*jd)+(kd*kd)));
      }
    }
  }
 
}

int main()
{
  Grid coarse_grid(32, 32, 32);
  Grid fine_grid(64, 64, 64);

  Field pressure;

  return 0;
}
