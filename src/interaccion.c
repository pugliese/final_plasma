#include "general.h"
#include "fft.h"
#include "celda.h"
#include "interaccion.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

double create_neutralizing_background(struct Interaccion *electric, double Q, int M, double L){
  electric->M = M;
  electric->dx = L/M;
  electric->L = L;
  electric->density = malloc(M*sizeof(double));
  electric->potential = malloc(M*sizeof(double));
  electric->field = malloc(M*sizeof(double));
  electric->background = malloc(M*sizeof(double));
  double charge = -Q/L;
  for(int m = 0; m < M; m++)  electric->background[m] = charge;
  return charge;
}

int calculate_density(struct Particles *parts, struct Interaccion *electric){
  int m, i, M = electric->M, N = parts->n;
  double dx = electric->dx, f;
  for(m = 0; m < M; m++)  electric->density[m] = electric->background[m];
  for(i = 0; i < N; i++){
    f = parts->x[i]/dx;
    m = (int) floor(f);
    f = f - m;
    electric->density[(m+1)%M] += parts->charge*f/dx;
    electric->density[m] += parts->charge*(1 - f)/dx;
  }
  return 0;
}

int calculate_potential(struct Interaccion *electric){
  int m, M = electric->M;
  double *imaginary = malloc(M*sizeof(double));
  double K, K2, dx_2 = electric->dx*0.5;
  for(m = 0; m < M; m++){
    electric->potential[m] = electric->density[m]/M;
    imaginary[m] = 0;
  }
  Fft_transform(electric->potential, imaginary, M);
  for(m = 1; m < M; m++){
    K = sin(M_PI*m/M)/dx_2;
    K2 = K*K;
    electric->potential[m] = electric->potential[m]/K2;
    imaginary[m] = imaginary[m]/K2;
  }
  Fft_inverseTransform(electric->potential, imaginary, M);
  free(imaginary);
  return 0;
}

int calculate_field(struct Interaccion *electric){
  int m, M = electric->M;
  double dx_2 = 2*electric->dx;
  for(m = 0; m < M; m++){
    electric->field[m] = (electric->potential[(m+M-1)%M] - electric->potential[(m+1)%M])/(dx_2);
  }
  return 0;
}

double forces(struct Particles *parts, struct Interaccion *electric){
  double pot = 0, dx = electric->dx, f, E, V, P;
  int i, M = electric->M, m;
  for (i = 0; i < parts->n; i++){
    f = parts->x[i]/dx;
    m = (int) floor(f);
    f = f - m;
    E = (1-f)*electric->field[m] + f*electric->field[(m+1)%M];
    V = (1-f)*electric->potential[m] + f*electric->potential[(m+1)%M];
    parts->f[i] = parts->charge*E;
    pot += parts->charge*V;
  }
  return pot;
}

int save_vector(char* filename, double* vector, int M, int append){
  FILE *f;
  if (append) f = fopen(filename, "a");
  else f = fopen(filename, "w");
  fprintf(f, "%d", append);
  for(int m = 0; m < M; m++)  fprintf(f, " %lf", vector[m]);
  fprintf(f, "\n");
  fclose(f);
  return 0;
}

int liberar_electric(struct Interaccion* electric){
  free(electric->density);
  free(electric->field);
  free(electric->potential);
  free(electric->background);

  return 0;
}
