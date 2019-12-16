#include "general.h"
#include "celda.h"
#include "interaccion.h"
#include "termostato.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

double add_forces_langevin(struct Particles *parts, struct Thermostat_Langevin *thermostat, double h){
  double sigma = sqrt(2*parts->mass*thermostat->gamma*thermostat->T/h), T_act = 0;
  int i;
  for (i = 0; i < parts->n; i++){
    parts->f[i] -= thermostat->gamma*parts->v[i] + boltzmann(sigma);
    T_act += parts->v[i]*parts->v[i];
  }
  T_act = T_act/parts->n;
  return T_act;
}
