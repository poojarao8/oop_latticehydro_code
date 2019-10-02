#include "main.h"

// Function that converts coordinates on 3-torus to the array index.
int I(int w, int x, int y, int z)
{
  return (x*W*H*3 + y*H*3 +z*3+ w);
}

int I(int w, int x, int y, int z, int length, int width, int height)
{
  return (x*width*height*3 + y*height*3 +z*3+ w);
}


//PRAO: added this for scalar field indexing
int J(int x, int y, int z, int l, int w, int h)
{
  return (x*w*h + y*h + z);
}

int J(int x, int y, int z)
{
  return (x*W*H + y*H +z);
}

