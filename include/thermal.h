#ifndef THERMAL_H
#define THERMAL_H

#include "types.h"

/* Thermal model parameters */
typedef struct {
    double density;       /* regolith density [kg/m^3] */
    double conductivity;  /* thermal conductivity [W/m/K] */
    double Cp;            /* specific heat [J/kg/K] */
    double albedo;        /* surface albedo */
    double distance_AU;   /* heliocentric distance [AU] */
    double depth;         /* max depth [m] */
    int n_substep;        /* substeps per rotation step */
    double z_angle;       /* obliquity [rad] */
    /* Antoine equation parameters for H2O vapor pressure */
    double p_mat[5];
    double LH_mat[5];
    double tran_unit_p;   /* Torr to Pa */
    double tran_unit_LH;  /* cal to J */
} ThermalParams;

/* Initialize default Ceres thermal parameters */
void thermal_params_default(ThermalParams *tp);

/* Run the thermal model, write surface temperature to grid->temperature.
 * Flux[ftd][ffd] = solar flux at each surface cell.
 * If flux is NULL, compute simple cos(theta)*sin(phi) flux. */
void thermal_solve(SphericalGrid *grid, const ThermalParams *tp,
                   const double *flux, int nthreads);

#endif /* THERMAL_H */
