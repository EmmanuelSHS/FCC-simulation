#ifndef PAIR_H
#define PAIR_H

#include "memory.h"

class PairCF {
public:
  PairCF();
  ~PairCF();

  char *fname;
  double **x, L, *gr;
  int *nr, nmax;
  double delr, r;
  int natom;

  Memory *memory;
private:
  int readxyz();// return 0 if programming is right, else >1
  void compute();
  void output();
};

#endif
