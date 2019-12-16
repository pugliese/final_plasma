#ifndef C_H
#define C_H

#include "math.h"

struct Particles {
  int n;
  double* x;
  double* v;
  double* f;
  double mass;
  double charge;

  double energy_pot;
  double kinetic;
};

double energia_cinetica(struct Particles *parts);
int save_lammpstrj(char *filename, struct Particles *parts, double L, int append);
int load_lammpstrj(char *filename, struct Particles *parts, double* L, double rcut);
int liberar(struct Particles *parts);

#endif
