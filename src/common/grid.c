#include "grid.h"
#include "constants.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

SphericalGrid *grid_create(int ftd, int ffd, int hd) {
    SphericalGrid *g = calloc(1, sizeof(SphericalGrid));
    g->ftd = ftd;
    g->ffd = ffd;
    g->hd = hd;
    g->density = calloc((size_t)hd * ftd * ffd, sizeof(double));
    g->temperature = calloc((size_t)ftd * ffd, sizeof(double));
    g->volume = calloc((size_t)hd * ftd * ffd, sizeof(double));
    g->horizon = calloc((size_t)ftd * ffd, sizeof(double));
    return g;
}

void grid_free(SphericalGrid *g) {
    if (!g) return;
    free(g->density);
    free(g->temperature);
    free(g->volume);
    free(g->horizon);
    free(g);
}

void grid_spherical_to_index(double thi, double fi, int ftd, int ffd,
                             int *td, int *fd) {
    /* thi in [-pi/2, pi/2], fi in [0, 2*pi) */
    /* MATLAB: td = ceil((thi+pi/2)/pi * ftd), fd = ceil(fi/2/pi * ffd) */
    *td = (int)ceil((thi + M_PI / 2.0) / M_PI * ftd) - 1; /* 0-based */
    *fd = (int)ceil(fi / (2.0 * M_PI) * ffd) - 1;          /* 0-based */

    if (*td < 0) *td = 0;
    if (*td >= ftd) *td = ftd - 1;
    if (*fd < 0) *fd = 0;
    if (*fd >= ffd) *fd = ffd - 1;
}

void grid_compute_volume_uniform(SphericalGrid *g, double Radii, double unih) {
    double dfi = 2.0 * M_PI / g->ffd;
    double oRi = Radii;
    double Ri = Radii;

    for (int ht = 0; ht < g->hd; ht++) {
        Ri = oRi + unih;
        double oth = -M_PI / 2.0;
        for (int td = 0; td < g->ftd; td++) {
            double thi = (td + 1.0) / g->ftd * M_PI - M_PI / 2.0;
            double vol = fabs(oRi * oRi * oRi + Ri * Ri * Ri) / 6.0 *
                         dfi * fabs(sin(oth) - sin(thi));
            for (int fd = 0; fd < g->ffd; fd++) {
                g->volume[grid_index(ht, td, fd, g->ftd, g->ffd)] = vol;
            }
            oth = thi;
        }
        oRi = Ri;
    }
}

void grid_compute_volume_exp(SphericalGrid *g, double Radii, double unih,
                             double alpha2) {
    double dfi = 2.0 * M_PI / g->ffd;
    double oRi = Radii;
    double Ri = Radii;

    for (int ht = 0; ht < g->hd; ht++) {
        Ri = oRi + unih * pow(1.0 + alpha2, ht);
        double oth = -M_PI / 2.0;
        for (int td = 0; td < g->ftd; td++) {
            if (td > 0) {
                oth = (double)td / g->ftd * M_PI - M_PI / 2.0;
            }
            double thi = (td + 1.0) / g->ftd * M_PI - M_PI / 2.0;
            double vol = fabs(oRi * oRi * oRi + Ri * Ri * Ri) / 6.0 *
                         dfi * fabs(sin(oth) - sin(thi));
            for (int fd = 0; fd < g->ffd; fd++) {
                g->volume[grid_index(ht, td, fd, g->ftd, g->ffd)] = vol;
            }
        }
        oRi = Ri;
    }
}

double grid_exp_height(int ht, double alpha2, double unih, double Radii) {
    /* h_line[ht] = ((1+alpha2)^(ht+1) - 1) / alpha2 * unih + Radii */
    return (pow(1.0 + alpha2, ht + 1) - 1.0) / alpha2 * unih + Radii;
}

int grid_exp_ht_index(double Ri, double Radii, double unih, double alpha2) {
    /* Inverse of exponential height: solve for ht */
    /* ht = ceil(log(alpha2*(Ri-Radii)/unih + 1) / log(1+alpha2)) - 1 */
    double arg = alpha2 * (Ri - Radii) / unih + 1.0;
    if (arg <= 1.0) return 0;
    int ht = (int)ceil(log(arg) / log(1.0 + alpha2)) - 1;
    if (ht < 0) ht = 0;
    return ht;
}
