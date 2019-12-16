#ifndef INTERACCION_H
#define INTERACCION_H

#include "math.h"
#include "celda.h"

struct Interaccion {
  double* density;
  double* potential;
  double* field;
  double* background;

  int M;
  double dx;
  double L;
};

double create_neutralizing_background(struct Interaccion *electric, double Q, int M, double L);
int calculate_density(struct Particles *parts, struct Interaccion *electric);
int calculate_potential(struct Interaccion *electric);
int calculate_field(struct Interaccion *electric);
double forces(struct Particles *parts, struct Interaccion *electric);
int save_vector(char* filename, double* vector, int M, int append);
int liberar_electric(struct Interaccion* electric);

#endif
