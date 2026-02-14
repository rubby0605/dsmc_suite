#ifndef PARTICLE_H
#define PARTICLE_H

#include "types.h"

/* Initialize particle at surface with random position */
void particle_init_surface(Particle *p, double Radii);

/* Update grid indices from current position */
void particle_update_indices(Particle *p, double Radii, int ftd, int ffd,
                             double unih, double alpha2, int use_exp_grid);

#endif /* PARTICLE_H */
