#ifndef AVANZAR_H
#define AVANZAR_H

#include "math.h"
#include "termostato.h"

int first_step(struct Particles *parts, double h);
int second_step(struct Particles *parts, double h);
int velocity_verlet(struct Particles *parts, struct Interaccion *lj, double h);
int velocity_verlet_NVT(struct Particles *parts, struct Interaccion *electric, double h, struct Thermostat_Langevin *thermostat);
int rescale_velocities(struct Particles *parts, double factor);

#endif
