#include "general.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

int elegir_part(int N){
  int i = rand();
  while (i > (RAND_MAX/N)*N){
    i = rand();
  }
  i = i % N;
  return i;
}

double uniform(){
  double p = ((double) rand())/((double) RAND_MAX);
  return p;
}

double normaldist()
{
  double x = 0,y,r2,c;
  r2 = 0.0;
  while(r2 > 1.0 || r2 == 0.0){
    x = 2.0*((double) rand() / (double)RAND_MAX) - 1.0;
    y = 2.0*((double) rand() / (double)RAND_MAX) - 1.0;
    r2 = x*x + y*y;
  }
  c = sqrt(-2.0*log(r2)/r2);
  return x*c;
}

double boltzmann(double sigma){
  return normaldist()*sigma;
}

int delta_celdas(double* q1, double* q2, double l, int *delta_idxs, double *delta){
  for(int k = 0; k < 3; k++)  delta[k] = q1[k] - q2[k] - l*delta_idxs[k];
  return 0;
}

double norma2(double* v){
  double res = 0;
  for(int k = 0; k < 3; k++)  res += v[k]*v[k];
  return res;
}

double max(double a, double b){
  if (a>b){
    return a;
  }else{
    return b;
  }
}

double max_vec(double* v, int n){
  double max = v[0];
  for(int i = 1; i < n; i++){
    if(max < v[i]) max = v[i];
  }
  return max;
}

int idx_max_norm(double* v, int n){
  int idx_max = 0;
  double rsq_max = norma2(v);
  for(int i = 1; i < n; i++){
    if(rsq_max < norma2(v+3*i)){
      idx_max = i;
      rsq_max = norma2(v+3*i);
    }
  }
  return idx_max;
}

double min(double a, double b){
  if (a<b){
    return a;
  }else{
    return b;
  }
}

int shuffle_array(int *array, int n){
  if (n > 1){
    int i, t;
    for (i = 0; i < n - 1; i++){
      size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
      t = array[j];
      array[j] = array[i];
      array[i] = t;
    }
  }
  return 0;
}
