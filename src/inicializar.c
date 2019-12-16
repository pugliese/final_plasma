#include "general.h"
#include "celda.h"
#include "inicializar.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

double set_box(struct Particles *parts, double L){
  double dL = L/parts->n;
  for(int i = 0; i < parts->n; i++){
    parts->x[i] = dL*(i+0.5);
  }
  return dL;
}

double set_perturbation(struct Particles *parts, double L, int k, double f){
  double dL = L/parts->n, K = 2*M_PI*k/L, x_o, d = f*dL;
  for(int i = 0; i < parts->n; i++){
    x_o = dL*(i+0.5);
    parts->x[i] = x_o + d*cos(K*x_o);
  }
  return dL;
}

double set_p(struct Particles *parts, double T){
  double sigma = sqrt(T/parts->mass);
  for(int k = 0; k < parts->n; k++){
    parts->v[k] = boltzmann(sigma);
  }
  double res = 0;
  for(int k = 0; k < parts->n; k++){
    res = res + parts->v[k]*parts->v[k]*parts->mass;
  }
  parts->kinetic = res/(2*parts->n);
  return res/parts->n;
}

double set_p_square(struct Particles *parts, double T){
  double sigma = sqrt(3*T/parts->mass);
  for(int k = 0; k < parts->n; k++){
    parts->v[k] = (2*uniform() - 1)*sigma;
  }
  double res = 0;
  for(int k = 0; k < parts->n; k++){
    res = res + parts->v[k]*parts->v[k]*parts->mass;
  }
  parts->kinetic = res/(2*parts->n);
  return res/parts->n;
}

double set_p_pulse(struct Particles *parts, double v){
  parts->v[0] = v;
  for(int k = 1; k < parts->n; k++){
    parts->v[k] = 0;
  }
  return 0;
}
