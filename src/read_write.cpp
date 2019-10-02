#include <iostream>
#include "main.h"

using namespace std;

int read_file(double arr[], int timestamp)
  {
    char fname[50];
    sprintf(fname,"%st%d.raw",OUTNAME,timestamp);
    ifstream in(fname, ios::out | ios::binary);
    if(!in) {
      cout << "Cannot open file." << endl;
      return 0;
    }

    for(int i=0; i<ARR_SIZE; i++) in.read((char *) &arr[i], sizeof(double));

    in.close();
    cout << "File "<< fname << " read" << endl;
    return 1;
  }

int write2file(double outdata[], int outdatasize, double nu, double eta, int cfonoff, double L2ratio, double gradVratio)
  {
    char fname[50];
    sprintf(fname,"%sN%unu%.3feta%.3fc%ib%.2fg%.2f.raw",OUTNAME,L,nu,eta,cfonoff,L2ratio,gradVratio);
    ofstream out(fname, ios::out | ios::binary);
    if(!out) {
      cout << "Cannot open file." << endl;
      return 0;
    }

    for(int i=0; i<outdatasize; i++) {
      for(int j=0; j<16; j++) {
        out.write((char *) &outdata[i*16+j], sizeof(double));
      }
    }

    out.close();
    cout << "File "<< fname << " created" << endl;
    return 1;
  }


int write2fileC(double outdata[], int outdatasize, double h, double nu, double eta)
  {
    char fname[50];
    sprintf(fname,"%sN%uh%.3fnu%.3feta%.3f.raw",OUTNAMEC,L,h,nu,eta);
    ofstream out(fname, ios::out | ios::binary);
    if(!out) {
      cout << "Cannot open file." << endl;
      return 0;
    }

    for(int i=0; i<outdatasize; i++) {
      for(int j=0; j<16; j++) {
        out.write((char *) &outdata[i*16+j], sizeof(double));
      }
    }

    out.close();
    cout << "File "<< fname << " created" << endl;
    return 1;
  }
