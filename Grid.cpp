#include "main.h"

Grid::Grid(int l, int w, int h)
{
  L = l;
  W = w;
  H = h;

  grid_pts = L*W*H;
  dx = 2*3.14/L;

  GL = L + 2*NGUARD; // grid containing the guard cells
  GW = W + 2*NGUARD;
  GH = H + 2*NGUARD; 

  cout << "Grid object is being created" << endl;
}

Grid::~Grid(void)
{
  cout << "Grid object is being deleted" << endl;
}
