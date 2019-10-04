#ifndef Grid_h
#define Grid_h

class Grid
{
  public:
   int L;
   int W;
   int H;
   int grid_pts;
   double dx;
   Grid(int, int, int); // constructor declaration
    ~Grid(); // destructor declaration
   
};

#endif
