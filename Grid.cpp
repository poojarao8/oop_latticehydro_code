#include <iostream>
#include "main.h"
#include "Grid.h" 

using namespace std;

Grid::Grid(int l, int w, int h)
{
  L = l;
  W = w;
  H = h;

  grid_pts = L*W*H;
  dx = 2*3.14/L;

  GL = L + 2*NGUARD;
  GW = W + 2*NGUARD;
  GH = H + 2*NGUARD; 

  cout << "Object is being created" << endl;
  cout << "dx = " << dx << endl;  
}

Grid::~Grid(void)
{
  cout << "Grid object is being deleted" << endl;
}
