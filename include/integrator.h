#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "types.h"

/* Result codes for a Verlet step */
typedef enum {
    VERLET_OK = 0,
    VERLET_LANDED,      /* Ri <= Radii */
    VERLET_ESCAPED       /* Ri >= max_height */
} VerletStatus;

/*
 * Störmer-Verlet time step.
 * Updates p->pos and p->vel in place.
 * Returns old acceleration in oa (for next step).
 * ds_out = displacement magnitude this step (for collision probability).
 */
VerletStatus verlet_step(Particle *p, const PhysicsParams *params,
                         double h, Vec3 *oa, double *ds_out);

#endif /* INTEGRATOR_H */
