#ifndef INICIALIZAR_H
#define INICIALIZAR_H

#include "math.h"

double set_box(struct Particles *parts, double L);
double set_perturbation(struct Particles *parts, double L, int k, double f);
double set_p(struct Particles *parts, double T);
double set_p_square(struct Particles *parts, double T);
double set_p_pulse(struct Particles *parts, double v);

#endif
