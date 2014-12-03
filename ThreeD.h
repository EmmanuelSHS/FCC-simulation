#ifndef ThreeD_H
#define ThreeD_H

#include"stdio.h"
#include"stdlib.h"
#include"memory.h"
#include"random.h"

class ThreeD{
public:
  ThreeD();
  ~ThreeD();

  int nstep, istep, ifreq;
  int nx, ny, nz, nucell, natom;
  double time, dt, hdt;
  double sigma, epsilon;
  double m, inv_m, hm;
  double d, rcut2, length, dr2, inv_dr2;
  double ke, pe;
  double tau, t0, t;
  double corcut, cotau;

  double *v_ave, *dx;
  double **x, **v, **f;
  Memory *memo;
  RanPark *ran;

  void init();
  void MD_loop();

private:
  void force();
  void output();
  void input();

  FILE *fp1, *fp2;
};

#endif
