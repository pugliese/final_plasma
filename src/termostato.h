#ifndef TERMOSTATO_H
#define TERMOSTATO_H

#include "math.h"
#include "celda.h"
#include "general.h"

struct Thermostat_Langevin{
  double gamma;
  double T;
};

double add_forces_langevin(struct Particles *parts, struct Thermostat_Langevin *thermostat, double h);

#endif
