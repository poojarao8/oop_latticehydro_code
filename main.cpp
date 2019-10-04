#include "main.h"

using namespace std;

int main()
{
  Grid coarse_grid(32, 32, 32);
  Grid fine_grid(64, 64, 64);

  Field pressure(1, coarse_grid);

  return 0;
}
