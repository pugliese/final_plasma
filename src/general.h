#ifndef GENERAL_H
#define GENERAL_H

#include "math.h"

int elegir_part(int N);
double uniform();
double normaldist();
double boltzmann(double sigma);
double max(double a, double b);
double min(double a, double b);
int shuffle_array(int *array, int n);
int delta_celdas(double* q1, double* q2, double l, int *delta_idxs, double *delta);
double norma2(double* v);
double max_vec(double* v, int n);
int idx_max_norm(double* v, int n);

#endif
