#ifndef VELOCITY_H
#define VELOCITY_H

#include "types.h"

/* Build MB velocity lookup table (23 temps from 20K to 240K) */
void velocity_table_init(VelocityTable *vt, double mmass, double mole,
                         double k_B);

/* Sample a speed from the table for given temperature */
double velocity_sample(const VelocityTable *vt, double T_local,
                       unsigned int *seed);

/* Compute outgoing speed with sticking coefficient */
double velocity_sticking(double vr, double v_old_sq, double alpha);

#endif /* VELOCITY_H */
