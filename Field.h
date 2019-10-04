class Field
{

  public:
    void initialize();
    Field(int, Grid&); // constructor
    ~Field(); // destructor
    Grid& obj;
  private:
    int ARR_SIZE;
};

