#ifndef GRID_H
#define GRID_H

#include "types.h"

/* Create/free spherical grid */
SphericalGrid *grid_create(int ftd, int ffd, int hd);
void grid_free(SphericalGrid *g);

/* Flat index into 3D array [ht][td][fd] (0-based) */
static inline int grid_index(int ht, int td, int fd, int ftd, int ffd) {
    return ht * ftd * ffd + td * ffd + fd;
}

/* Convert spherical coords to grid indices (0-based) */
void grid_spherical_to_index(double thi, double fi, int ftd, int ffd,
                             int *td, int *fd);

/* Compute cell volumes for uniform height grid */
void grid_compute_volume_uniform(SphericalGrid *g, double Radii, double unih);

/* Compute cell volumes for exponential height grid */
void grid_compute_volume_exp(SphericalGrid *g, double Radii, double unih,
                             double alpha2);

/* Get height for exponential grid (returns radius at top of cell) */
double grid_exp_height(int ht, double alpha2, double unih, double Radii);

/* Convert radius to exponential grid height index (0-based) */
int grid_exp_ht_index(double Ri, double Radii, double unih, double alpha2);

#endif /* GRID_H */
