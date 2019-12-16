#include "celda.h"
#include "general.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

double energia_cinetica(struct Particles *parts){
  double Ekin = 0;
  for(int i = 0; i < parts->n; i++){
    Ekin += parts->v[i]*parts->v[i];
  }
  Ekin = 0.5*parts->mass*Ekin/parts->n;
  parts->kinetic = Ekin;
  return Ekin;
}
/*
int save_lammpstrj(char *filename, struct Particles *parts, double L, int append){
  FILE *f;
  if (append) f = fopen(filename, "a");
  else f = fopen(filename, "w");
	int m, M = parts->M, i, t;
	double l = parts->l;
	fprintf(f, "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp\n", append, parts->n);
	for(int l = 0; l < 3; l++){
		fprintf(f, "0 %f\n", L);
	}
  fprintf(f, "ITEM: ATOMS id type x vx \n");
	for(m = 0; m < M; m++){
		i = parts->primero[m];
		while(i!=-1){
      fprintf(f, "%d %d %lf %lf\n", i, t, parts->x[i], parts->v[i]);
			i = parts->siguiente[i];
		}
	}
  fclose(f);
  return 0;
}

int load_lammpstrj(char *filename, struct Particles *parts, double* L, double rcut){
  FILE *f = fopen(filename, "r");
  char buffer[255];
  int id;
  for(int l = 0; l < 3; l++){
    fgets(buffer, 255, f);
  }
  id = fscanf(f, "%d\n", &parts->n);
  for(int l = 0; l < 2; l++){
    fgets(buffer, 255, f);
  }
  id = fscanf(f, "0 %lf\n", L);
  for(int l = 0; l < 2; l++){
    fgets(buffer, 255, f);
  }
  parts->x = (double *) malloc(parts->n*sizeof(double));
  parts->v = (double *) malloc(parts->n*sizeof(double));
	for(int l = 0; l < parts->n; l++){
    id = fscanf(f, "%d %d %lf %lf \n", &id, &id, parts->x+l, parts->v+l);
	}
  fclose(f);
  return 0;
}
*/
int liberar(struct Particles *parts){
  free(parts->x);
  free(parts->v);
  free(parts->f);

  return 0;
}
