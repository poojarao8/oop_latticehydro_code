#include <iostream>


/*
Grid Class
*/


using namespace std;

class Grid
{
  public:
   int VEC_ARR_SIZE;
   int SCL_ARR_SIZE;
   Grid(int, int, int); // constructor declaration
    ~Grid(); // destructor declaration

  private:
    double dx;
};

// memmber function declarations including the constructor

Grid::Grid(int L, int W, int H)
{
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
  cout << "Object is being deleted" << endl;
}


int main()
{
  Grid coarse_grid(32, 32, 32);
  Grid fine_grid(64, 64, 64);

  return 0;
}
