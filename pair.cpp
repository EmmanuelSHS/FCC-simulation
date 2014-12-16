#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "memory.h"
#include "math.h"
#include "pair.h"

#define MAXLINE 256
#define DIM 3

PairCF::PairCF()
{
  x = NULL;
  gr = NULL;
  nr = NULL;
  fname = NULL;
  natom = 0; // later input needed
  memory = new Memory();

  char str[MAXLINE], *ptr;
  while (1){
    printf("\n Please input your xyz file name: ");
    fgets(str, MAXLINE, stdin);
    ptr = strtok(str, " \n\t\r\f");
    if (ptr) break;
  }
  fname = new char [strlen(ptr) + 1];
  strcpy(fname, ptr);

  L = 0.;
  while (1){
    printf("\n Please input your box length: ");
    fgets(str, MAXLINE, stdin);
    char *ptr = strtok(str, " \n\t\r\f");
    if (ptr) break;
  }
  L = atof(ptr);
  if (L < 0.) exit(0);

  delr = -1.;
  while (delr <= 0.){
    printf("\n Please input your step size for g(r): ");
    fgets(str, MAXLINE, stdin);
    char *ptr = strtok(str, " \n\t\r\f");
    if (ptr) delr = atof(ptr);
  }

  nmax = 0;
  while (nmax <= 0){
    printf("\n Please input your # of bins for g(r): ");
    fgets(str, MAXLINE, stdin);
    char *ptr = strtok(str, " \n\t\r\f");
    if (ptr) nmax = atoi(ptr);
  }

  double LMax = double(nmax) *delr;
  if (LMax > 0.5 * L){
    nmax = (0.5 * L / delr);
  }
  // read xyz
  readxyz();
  compute();
  output();

return;
}

PairCF::~PairCF()
{
  if (x) memory -> destroy(x);
  if (nr);
  if (gr);
  if (fname);

  delete memory;

return;
}

int PairCF::readxyz()
{
  FILE *fp = fopen(fname, "r");
  if (fp == NULL) return 1;

  char str[MAXLINE], *ptr;
  // read # of atoms
  fgets(str, MAXLINE, fp);
  ptr = strtok(str, " \n\r\t\f");
  if (ptr == NULL) return 2;
  natom = atoi(ptr);
  if (natom < 1){
    fclose(fp);
    return 3;
  } // to close fp before return

  memory -> create(x, natom, 3, "x");

  // comment line
  fgets(str, MAXLINE, fp);

  // read all the position data
  for (int i = 0; i < natom; ++i){
    fgets(str, MAXLINE, fp);
    ptr = strtok(str, " \n\t\r\f");
    if (ptr == NULL) return 4;
    for (int j = 0; j < 3; ++j){
      ptr = strtok(NULL, " \n\t\r\f");
      if (ptr == NULL) return 5;
      x[i][j] = atof(ptr);
    }
  }

  fclose(fp);

return 0;
}

void PairCF::compute()
{
  memory -> create(nr, nmax, "nr");
  memory -> create(gr, nmax, "gr");

  for (int i = 0; i < nmax; ++i) nr[i] = 0;

  double hL = 0.5 * L, Rij[DIM];
  // compute n(r)
  for (int i = 0; i < natom; ++i){
  for (int j = i+1; j < natom; ++j){
    double r2 = 0.;
    for (int idim = 0; idim < DIM; ++idim){
      Rij[idim] = x[j][idim] - x[i][idim];
      while (Rij[idim] > hL)  Rij[idim] -= L;
      while (Rij[idim] <-hL)  Rij[idim] += L;
      r2 += Rij[idim] * Rij[idim];
    }

    r = sqrt(r2);
    int ibin = int(r / delr + 0.5);// get a inv_delr
    if (ibin < nmax) nr[ibin] += 2;
  }
  }

  double factor = 1./((16. * atan(1.)) * delr * double(natom) * double(natom) / (L*L*L));

  for (int i = 1; i < nmax; ++i){ // i = 0 will get a not a number
    r = double(i) * delr;
    double inv_r2 = 1./(r * r);
    gr[i] = nr[i] * factor * inv_r2;
  }

}

void PairCF::output()
{
  char str[MAXLINE], *ptr;
  char *fout;

  while (1){
    printf("\n Please name your output file: ");
    fgets(str, MAXLINE, stdin);
    ptr = strtok(str, " \n\t\r\f");
    if (ptr) break;
  }
 
  fout = new char [strlen(ptr) + 1];
  strcpy(fout, ptr);

  FILE *fp = fopen(fout, "w");
  fprintf(fp, "# r  g(r)\n");
  for (int i = 1; i < nmax; ++i){
    r = double(i) * r;
    fprintf(fp, "%lg %lg \n", r, gr[i]);
  }

  fclose(fp);

}
