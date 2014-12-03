#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "ThreeD.h"
#include "random.h"
#include "memory.h"
#include "math.h"

#define KB 1.
#define MAXLINE 5120

ThreeD::ThreeD()
{
  nstep = istep = 0;
  nx = ny = nz = nucell = natom = 0;
  ifreq = 0;
  time = 0.;
  sigma = epsilon = 0.;
  m = inv_m = hm = 0.;
  dt = hdt = 0.;
  d = rcut2 = length = dr2 = inv_dr2 = 0.;
  corcut = cotau = 0.;
  ke = pe = 0.;
  x = v = f = NULL;
  memo = new Memory();
  ran = new RanPark(1234);
  tau = 0.;
  t0 = t = 0.;
  v_ave = NULL;
  dx = NULL;
}

ThreeD::~ThreeD()
{
  if (x) memo -> destroy(x);
  if (v) memo -> destroy(v);
  if (f) memo -> destroy(f);
  delete ran;
  delete memo;
  if (v_ave) delete v_ave;
  if (dx) delete dx;
}

void ThreeD::init()
{
/*
  nstep = 100;
  nx = ny = nz = 3;
  ifreq = 10;
  sigma = epsilon = 1.;
  d = 1.5;
  rcut2 = 1.25 * sigma * sigma;
  m = 1.;
  t0 = 100.;

  dt = 0.0005;
  hdt = 0.5 * dt;
  tau = 1000. * dt;
*/

  input();
  tau = cotau * dt;
  rcut2 = corcut * corcut * sigma * sigma;
  hdt = 0.5 * dt;

/*
  nstep = 100;
  nx = ny = nz = 3;
  ifreq = 10;
  sigma = epsilon = 1.;
  d = 1.5;
  rcut2 = 1.25 * sigma * sigma;
  m = 1.;
  t0 = 100.;

  dt = 0.0005;
  hdt = 0.5 * dt;
  tau = 1000. * dt;
*/

  nucell = 4;
  natom = nx * ny * nz * nucell;
  istep = 0;
  length = double(nx) * d;
  inv_m = 1./m;
  hm = 0.5 * m;
  
  memo -> create(x, natom, 3, "x");
  memo -> create(v, natom, 3, "v");
  memo -> create(f, natom, 3, "f");

  t = 0.;
  time = 0.;

  ke = pe = 0.;

  v_ave = new double [3];
  dx = new double [3];
  
  // initialize atomic position
  int i = 0;
  for (int ix = 0; ix < nx; ++ix)
  for (int iy = 0; iy < ny; ++iy)
  for (int iz = 0; iz < nz; ++iz){
    x[i][0] = d * double(ix);
    x[i][1] = d * double(iy);
    x[i][2] = d * double(iz);
    ++i;

    x[i][0] = d * (double(ix) + 0.5);
    x[i][1] = d * (double(iy) + 0.5);
    x[i][2] = d * double(iz);
    ++i;

    x[i][0] = d * double(ix);
    x[i][1] = d * (double(iy) + 0.5);
    x[i][2] = d * (double(iz) + 0.5);
    ++i;

    x[i][0] = d * (double(ix) + 0.5);
    x[i][1] = d * double(iy);
    x[i][2] = d * (double(iz) + 0.5);
    ++i;
  }

  // initialize the velocity
  for (int k = 0; k < 3; ++k) v_ave[k] = 0.;

  for (int i = 0; i < natom; ++i)
  for (int k = 0; k < 3; ++k){
    v[i][k] = ran -> gaussian() - 0.5;
    v_ave[k] += v[i][k];
  }

  for (int k = 0; k < 3; ++k) v_ave[k] /= double(natom);

  for (int i = 0; i < natom; ++i)
  for (int k = 0; k < 3; ++k){
    v[i][k] -= v_ave[k];
    ke += v[i][k] * v[i][k];
  }
  ke *= hm;
  t = ke / (1.5 * double(natom) * KB);

  ke = 0.;
  double t0tot = sqrt(t0 / t);
  for (int i = 0; i < natom; ++i)
  for (int k = 0; k < 3; ++k){
    v[i][k] *= t0tot;
    ke += v[i][k] * v[i][k];
  }
  ke *= hm;
  t = ke / (1.5 * double(natom) * KB);

  force();

return;
}

void ThreeD::MD_loop()
{
  fp1 = fopen( "log.dat", "w");
  fp2 = fopen( "atomcfg.xyz", "w");
  fprintf(fp1, "# istep time pe ke temperature\n");

  for (int i = 0; i < nstep; ++i){
    if ( istep % ifreq == 0) output();
    ++istep;

    // velocity verlet
    for ( int i = 0; i < natom; ++i)
    for ( int k = 0; k < 3; ++k)  v[i][k] += f[i][k] * hdt * inv_m; // the force at time 0 is calculated at init() -> force()

    for ( int i = 0; i < natom; ++i)
    for ( int k = 0; k < 3; ++k)  x[i][k] += v[i][k] * dt;

    force();

    ke = 0.;

    for ( int i = 0; i < natom; ++i)
    for ( int k = 0; k < 3; ++k){
      v[i][k] += f[i][k] * inv_m * hdt;
      ke += v[i][k] * v[i][k];
    }

    ke *= hm;
    t = ke / (1.5 * double(natom) * KB);

    time += dt;
  }
  output();

  fclose(fp1);
  fclose(fp2);

return;
}

void ThreeD::force()
{
  pe = 0.;

  double m_tau = m / tau;
  double w = sqrt(24. * m * KB * t0 / (tau * dt));

  for ( int i = 0; i < natom; ++i)
  for ( int k = 0; k < 3; ++k) f[i][k] = 0.;

  for ( int i = 0; i < natom; ++i){
  for ( int j = i + 1; j < natom; ++j){
    for (int k = 0; k < 3; ++k){
      dx[k] = x[j][k] - x[i][k];
        
      // periodical boundary condition
      while ( dx[k] > 0.5 * length) dx[k] -= length;
      while ( dx[k] < -0.5 * length) dx[k] += length;
    }

    dr2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

    if ( dr2 < rcut2){
      inv_dr2 = 1./(dr2);
//      double s_d = sigma * inv_dr; // rule: if *, no underline, divided by, with        
      double s_d2 = sigma * sigma / dr2;
      double s_d6 = s_d2 * s_d2 * s_d2;
      double s_d12 = s_d6 * s_d6;
      pe += s_d12 - s_d6;

      for (int k = 0; k < 3; ++k){
        f[i][k] += (-12. * s_d12 + 6. * s_d6) * dx[k] * 4. * epsilon * inv_dr2; 
        f[j][k] -= (-12. * s_d12 + 6. * s_d6) * dx[k] * 4. * epsilon * inv_dr2; 
      }
    }
  }
    for ( int k = 0; k < 3; ++k) f[i][k] += -m_tau * v[i][k] + (ran -> uniform() - 0.5) * w;
  }
  pe *= 4. * epsilon;
return;
}

void ThreeD::input()
{
  int len = 0;
  char str[MAXLINE];
  FILE *fpara = NULL;

  fpara = fopen("parameters.txt", "r");
  
  while (feof(fpara) == 0){
    fgets(str, MAXLINE, fpara);
    len = strlen(str);
    char *p = strtok(str, " \n\t\r\0");// at the end of txt file there exist a '\0'

    if (strcmp(p, "nx") == 0){
      p = strtok(NULL, " \n\t\r\0");
      nx = atoi(p);
    }

    if (strcmp(p, "ny") == 0){
      p = strtok(NULL, " \n\t\r\0");
      ny = atoi(p);
    }

    if (strcmp(p, "nz") == 0){
      p = strtok(NULL, " \n\t\r\0");
      nz = atoi(p);
    }

    if (strcmp(p, "ifreq") == 0){
      p = strtok(NULL, " \n\t\r\0");
      ifreq = atoi(p);
    }

    if (strcmp(p, "nstep") == 0){
      p = strtok(NULL, " \n\t\r\0");
      nstep = atoi(p);
    }

    if (strcmp(p, "sigma") == 0){
      p = strtok(NULL, " \n\t\r\0");
      sigma = double(atof(p));
    }

    if (strcmp(p, "epsilon") == 0){
      p = strtok(NULL, " \n\t\r\0");
      epsilon = double(atof(p));
    }

    if (strcmp(p, "rcut") == 0){
      p = strtok(NULL, " \n\t\r\0");
      corcut = double(atof(p));
    }

    if (strcmp(p, "tau") == 0){
      p = strtok(NULL, " \n\t\r\0");
      cotau = double(atof(p));
    }

    if (strcmp(p, "m") == 0){
      p = strtok(NULL, " \n\t\r\0");
      m = double(atof(p));
    }

    if (strcmp(p, "d") == 0){
      p = strtok(NULL, " \n\t\r\0");
      d = double(atof(p));
    }

    if (strcmp(p, "t0") == 0){
      p = strtok(NULL, " \n\t\r\0");
      t0 = double(atof(p));
    }

    if (strcmp(p, "dt") == 0){
      p = strtok(NULL, " \n\t\r\0");
      dt = double(atof(p));
    }

  }

fclose(fpara);
}

void ThreeD::output()
{
  // thermal info
  fprintf(fp1, "%d %lg %lg %lg %lg\n", istep, time, pe, ke, t);

  // configuration
  fprintf(fp2, "%d\n", natom);
  fprintf(fp2, "Step: %d, PE = %lg, KE = %lg\n", istep, pe, ke);
  for (int i = 0; i < natom; ++i) fprintf(fp2, "H %lg %lg %lg\n", x[i][0], x[i][1], x[i][2]);
}
