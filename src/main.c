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
    double l = 1, To = 0;
    if(argc > 2)  sscanf(argv[2], "%lf", &To);
    sprintf(filename_potential, "data/point_charge/potencial_medio_T=%.1f.txt", To);
    sprintf(filename_density, "data/point_charge/densidad_medio_T=%.1f.txt", To);
    FILE *fp = fopen(filename, "w");

    thermostat.T = To;

    parts.n = 256;
    parts.x = (double *) malloc(parts.n*sizeof(double));
    parts.v = (double *) malloc(parts.n*sizeof(double));
    parts.f = (double *) malloc(parts.n*sizeof(double));

    int M = parts.n;
    double L = l*parts.n, Q = 10;

    double* mean_potential = (double *) malloc(M*sizeof(double));
    double* mean_density = (double *) malloc(M*sizeof(double));
    for(int m = 0; m < M; m++){
      mean_potential[m] = 0;
      mean_density[m] = 0;
    }

    create_neutralizing_background(&electric, parts.charge*parts.n + Q, M, L);
    electric.background[M/2] += Q/electric.dx;

    set_box(&parts, L);
    set_p(&parts, To);
    calculate_density(&parts, &electric);
    calculate_potential(&electric);
    calculate_field(&electric);

    parts.energy_pot = forces(&parts, &electric)/parts.n;

    double h = 0.001;
    int N_steps = 10001000, t, N_skip = 1000;
    int segs = time(NULL);
    for(t = 0; t < N_steps; t++){
      //printf("Progreso: %d%%\r", (100*t)/N_steps);
      fflush(stdout);
      /*
      if(t % 50 == 0){
        save_vector(filename_density, electric.density, electric.M, t);
        save_vector(filename_potential, electric.potential, electric.M, t);
      }
      fprintf(fp, "%f %f\n", parts.kinetic, parts.energy_pot);
      */
      velocity_verlet_NVT(&parts, &electric, h, &thermostat);
      if (t >= N_skip){
        for(int m = 0; m < M; m++){
          mean_potential[m] += electric.potential[m]/(N_steps-N_skip);
          mean_density[m] += electric.density[m]/(N_steps-N_skip);
        }
      }
    }
    segs = time(NULL) - segs;
    printf("Progreso: 100%%\nTardo %dmin, %dsegs\n", segs/60, segs%60);
    save_vector(filename_density, mean_density, electric.M, 0);
    save_vector(filename_potential, mean_potential, electric.M, 0);
    fclose(fp);
  }

  if(opcion == '2'){
    char filename_density[255];
    char filename_potential[255];
    char filename_phase_x[255];
    char filename_phase_v[255];
    sprintf(filename, "data/two_streams/v5/energia.txt");
    sprintf(filename_potential, "data/two_streams/v5/potencial.txt");
    sprintf(filename_density, "data/two_streams/v5/densidad.txt");
    sprintf(filename_phase_x, "data/two_streams/v5/fases_x.txt");
    sprintf(filename_phase_v, "data/two_streams/v5/fases_v.txt");
    FILE *fp = fopen(filename, "w");


    parts.n = 4096;
    parts.x = (double *) malloc(parts.n*sizeof(double));
    parts.v = (double *) malloc(parts.n*sizeof(double));
    parts.f = (double *) malloc(parts.n*sizeof(double));

    double l = 1, v = parts.n/20;
    int M = 1024;
    double L = l*parts.n;

    create_neutralizing_background(&electric, parts.charge*parts.n, M, L);

    set_perturbation(&parts, L, 2, 0.03);
    set_two_streams(&parts, v);
    calculate_density(&parts, &electric);
    calculate_potential(&electric);
    calculate_field(&electric);

    parts.energy_pot = forces(&parts, &electric)/parts.n;

    double h = 0.01;
    int N_steps = 10001, t;
    int segs = time(NULL);
    for(t = 0; t < N_steps; t++){
      printf("Progreso: %d%%\r", (100*t)/N_steps);
      fflush(stdout);
      if(t % 10 == 0){
        save_vector(filename_density, electric.density, electric.M, t);
        save_vector(filename_potential, electric.potential, electric.M, t);
        save_vector(filename_phase_x, parts.x, parts.n, t);
        save_vector(filename_phase_v, parts.v, parts.n, t);
        fprintf(fp, "%f %f\n", parts.kinetic, parts.energy_pot);
      }
      velocity_verlet(&parts, &electric, h);
    }
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
