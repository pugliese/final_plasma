#include "general.h"
#include "fft.h"
#include "interaccion.h"
#include "inicializar.h"
#include "avanzar.h"
#include "celda.h"
#include "termostato.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

int main(int argc, char *argv[]){
  srand(time(NULL));
  char filename[255];
  char opcion = 'g';
  if (argc >= 2){
    sscanf(argv[1], "%c\n", &opcion);
  }

  struct Particles parts;
  parts.mass = 1;
  parts.charge = -1;
  struct Interaccion electric;

  struct Thermostat_Langevin thermostat;
  thermostat.gamma = 0.5;

  if(opcion == '1'){
    char filename_density[255];
    char filename_potential[255];
    sprintf(filename, "data/point_charge/cinetica.txt");
    sprintf(filename_potential, "data/point_charge/potencial.txt");
    sprintf(filename_density, "data/point_charge/densidad.txt");
    FILE *fp = fopen(filename, "w");

    double l = 1, To = 0;
    thermostat.T = To;

    parts.n = 4096/2;
    parts.x = (double *) malloc(parts.n*sizeof(double));
    parts.v = (double *) malloc(parts.n*sizeof(double));
    parts.f = (double *) malloc(parts.n*sizeof(double));

    int M = 1024;
    double L = l*parts.n, Q = 10;

    create_neutralizing_background(&electric, parts.charge*parts.n + Q, M, L);
    electric.background[M/2] += Q/electric.dx;

    set_box(&parts, L);
    set_p(&parts, To);
    calculate_density(&parts, &electric);
    calculate_potential(&electric);
    calculate_field(&electric);

    parts.energy_pot = forces(&parts, &electric)/parts.n;

    double h = 0.001;
    int N_steps = 100001, t;
    int segs = time(NULL);
    for(t = 0; t < N_steps; t++){
      printf("Progreso: %d%%\r", (100*t)/N_steps);
      fflush(stdout);
      if(t % 100 == 0){
        save_vector(filename_density, electric.density, electric.M, t);
        save_vector(filename_potential, electric.potential, electric.M, t);
      }
      fprintf(fp, "%f %f\n", parts.kinetic, parts.energy_pot);
      velocity_verlet_NVT(&parts, &electric, h, &thermostat);
    }
    printf("%f %f %f\n", electric.density[0], electric.potential[0], electric.field[0]);
    segs = time(NULL) - segs;
    printf("Progreso: 100%%\nTardo %dmin, %dsegs\n", segs/60, segs%60);
    fclose(fp);
  }

  if(opcion == 'p'){
    char filename_density[255];
    char filename_potential[255];
    int k = 1;
    if(argc > 2)  sscanf(argv[2], "%d", &k);
    sprintf(filename, "data/perturbacion/cinetica_k=%d.txt", k);
    sprintf(filename_potential, "data/perturbacion/potencial_k=%d.txt", k);
    sprintf(filename_density, "data/perturbacion/densidad_k=%d.txt", k);
    FILE *fp = fopen(filename, "w");

    double l = 1, To = 0;

    parts.n = 4096/2;
    parts.x = (double *) malloc(parts.n*sizeof(double));
    parts.v = (double *) malloc(parts.n*sizeof(double));
    parts.f = (double *) malloc(parts.n*sizeof(double));

    int M = 1024;
    double L = l*parts.n;

    create_neutralizing_background(&electric, parts.charge*parts.n, M, L);

    set_perturbation(&parts, L, k, 0.1);
    set_p(&parts, To);
    calculate_density(&parts, &electric);
    calculate_potential(&electric);
    calculate_field(&electric);

    parts.energy_pot = forces(&parts, &electric)/parts.n;

    double h = 0.001;
    int N_steps = 10001, t;
    int segs = time(NULL);
    for(t = 0; t < N_steps; t++){
      printf("Progreso: %d%%\r", (100*t)/N_steps);
      fflush(stdout);
      if(t % 10 == 0){
        save_vector(filename_density, electric.density, electric.M, t);
        save_vector(filename_potential, electric.potential, electric.M, t);
      }
      fprintf(fp, "%f %f\n", parts.kinetic, parts.energy_pot);
      velocity_verlet(&parts, &electric, h);
    }
    printf("%f %f %f\n", electric.density[0], electric.potential[0], electric.field[0]);
    segs = time(NULL) - segs;
    printf("Progreso: 100%%\nTardo %dmin, %dsegs\n", segs/60, segs%60);
    fclose(fp);
  }

  liberar(&parts);
  liberar_electric(&electric);

  return 0;
}
