#include "general.h"
#include "interaccion.h"
#include "termostato.h"
#include "avanzar.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

int first_step(struct Particles *parts, double h){
  for(int i = 0; i < parts->n; i++){
    parts->x[i] += h*(parts->v[i] + 0.5*h*parts->f[i]/parts->mass);
    parts->v[i] += 0.5*h*parts->f[i]/parts->mass;
  }
  return 0;
}

int apply_PBC(struct Particles *parts, double L){
  for(int i = 0; i < parts->n; i++){
    if (parts->x[i] > L){
      parts->x[i] -= L;
    }else{
      if (parts->x[i] < 0){
        parts->x[i] += L;
      }
    }
  }
  return 0;
}

int second_step(struct Particles *parts, double h){
  for(int i = 0; i < parts->n; i++){
    parts->v[i] += 0.5*h*parts->f[i]/parts->mass;
  }
  return 0;
}

int velocity_verlet(struct Particles *parts, struct Interaccion *electric, double h){
  first_step(parts, h);
  apply_PBC(parts, electric->L);
  calculate_density(parts, electric);
  calculate_potential(electric);
  calculate_field(electric);
  parts->energy_pot = forces(parts, electric)/parts->n;
  second_step(parts, h);
  energia_cinetica(parts);
  return 0;
}

int velocity_verlet_NVT(struct Particles *parts, struct Interaccion *electric, double h, struct Thermostat_Langevin *thermostat){
  first_step(parts, h);
  apply_PBC(parts, electric->L);
  calculate_density(parts, electric);
  calculate_potential(electric);
  calculate_field(electric);
  parts->energy_pot = forces(parts, electric)/parts->n;
  add_forces_langevin(parts, thermostat, h);
  second_step(parts, h);
  energia_cinetica(parts);
  return 0;
}

int rescale_velocities(struct Particles *parts, double factor){
  for(int i = 0; i < parts->n; i++) parts->v[i] *= factor;
  return 0;
}
