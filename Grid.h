#ifndef Grid_h
#define Grid_h

class Grid
{
  public:
   int L;
   int W;
   int H;
   int grid_pts;
   Grid(int, int, int); // constructor declaration
    ~Grid(); // destructor declaration
   
  private:
   double dx;
};

#endif
