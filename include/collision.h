#ifndef COLLISION_H
#define COLLISION_H

#include "types.h"

/* Check if collision should occur (accumulative probability) */
int collision_check(double *real_col_p, double col_threshold,
                    double ds, double local_density, double cross_section);

/* Execute elastic hard-sphere collision in CM frame.
 * v = test particle velocity (modified in place)
 * v2 = collision partner velocity */
void collision_execute(Vec3 *v, Vec3 v2);

/* Sample collision partner speed from polynomial CDF using bisection */
double collision_sample_speed(const PolyFit *pf, double rand_val, int nn);

#endif /* COLLISION_H */
