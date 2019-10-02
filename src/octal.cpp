void split8(double arr[], double **out) {

  int octal;
  for(int b1=0; b1<2; b1++) {for(int b2=0; b2<2; b2++) {for(int b3=0; b3<2; b3++) {
    for(int i=0; i<LL; i++) {
      for(int j=0; j<WW; j++) {
        for(int k=0; k<HH; k++) {
          octal=b1*4+b2*2+b3;
          out[octal][i*WW*HH+j*HH+k]=arr[(i*2+b1)*W*H+(j*2+b2)*H+k*2+b3];
        }
      }
    }
  }}}
}

//Function to undo split8
void unsplit8(double **arr, double out[]) {

  int octal;
  for(int b1=0; b1<2; b1++) {for(int b2=0; b2<2; b2++) {for(int b3=0; b3<2; b3++) {
    for(int i=0; i<LL; i++) {
      for(int j=0; j<WW; j++) {
        for(int k=0; k<HH; k++) {
          octal=b1*4+b2*2+b3;
          out[(i*2+b1)*W*H+(j*2+b2)*H+k*2+b3]=arr[octal][i*WW*HH+j*HH+k];
        }
      }
    }
  }}}
}

