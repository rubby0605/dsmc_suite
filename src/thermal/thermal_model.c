#include "thermal.h"
#include "constants.h"
#include "grid.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

/*
 * 1D spherical thermal diffusion model with sublimation.
 * Converted from "m0-1 thermal model.m"
 *
 * For each surface cell (td, fd):
 *   - Solve 1D heat equation with Crank-Nicolson scheme
 *   - Surface BC: solar flux - sigma*T^4 - sublimation heat loss = 0
 *   - Bottom BC: zero heat flux (insulating)
 *   - Use bisection for nonlinear surface temperature
 */

void thermal_params_default(ThermalParams *tp) {
    tp->density = 730.0;         /* kg/m^3 (dust) */
    tp->conductivity = 6e-4;    /* W/m/K */
    tp->Cp = 1300.0;            /* J/kg/K */
    tp->albedo = 0.02;
    tp->distance_AU = 2.6;
    tp->depth = 0.05;           /* m */
    tp->n_substep = 5;
    tp->z_angle = 12.3 / 180.0 * M_PI;

    /* Antoine equation for H2O vapor pressure */
    tp->p_mat[0] = -2445.5646;
    tp->p_mat[1] = -6.757169;
    tp->p_mat[2] = -0.01677006;
    tp->p_mat[3] = 1.20514e-5;
    tp->p_mat[4] = 8.2312;
    tp->tran_unit_p = 133.32237;  /* Torr to Pa */

    /* Latent heat parameters */
    tp->LH_mat[0] = 0.0;
    tp->LH_mat[1] = 12420.0;
    tp->LH_mat[2] = -4.8;
    tp->LH_mat[3] = 0.0;
    tp->LH_mat[4] = 0.0;
    tp->tran_unit_LH = 4.2;      /* cal to J */
}

/* Evaluate Antoine vapor pressure at temperature T */
static double antoine_pressure(const ThermalParams *tp, double T) {
    double log10_p = tp->p_mat[0] / T + tp->p_mat[1] +
                     tp->p_mat[2] * T + tp->p_mat[3] * T * T +
                     tp->p_mat[4] * log10(T);
    return tp->tran_unit_p * pow(10.0, log10_p);
}

/* Evaluate latent heat at temperature T */
static double latent_heat(const ThermalParams *tp, double T) {
    double LH = tp->LH_mat[0] / T + tp->LH_mat[1] +
                tp->LH_mat[2] * T + tp->LH_mat[3] * T * T +
                tp->LH_mat[4] * log10(T);
    return tp->tran_unit_LH * LH;
}

/* Surface energy balance function for bisection:
 * f(T) = F + k/dr*T_below - k/dr*T - sigma*T^4 - sublimation */
static double surface_balance(double T, double F, double k_over_dr,
                              double oT, const ThermalParams *tp) {
    double mmol = CONST_MMOL_H2O;
    double p = antoine_pressure(tp, T);
    double LH = latent_heat(tp, T);
    double sublim = LH * sqrt(mmol / (2.0 * M_PI * CONST_K_B * T)) * p;
    return F + k_over_dr * oT - k_over_dr * T - CONST_SIGMA * T * T * T * T - sublim;
}

/*
 * Thomas algorithm for tridiagonal system:
 *   a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1] = d[i]
 */
static void thomas_solve(double *a, double *b, double *c, double *d,
                         double *x, int n) {
    /* Forward sweep */
    double *cp = malloc(n * sizeof(double));
    double *dp = malloc(n * sizeof(double));

    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];

    for (int i = 1; i < n; i++) {
        double denom = b[i] - a[i] * cp[i - 1];
        cp[i] = c[i] / denom;
        dp[i] = (d[i] - a[i] * dp[i - 1]) / denom;
    }

    /* Back substitution */
    x[n - 1] = dp[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = dp[i] - cp[i] * x[i + 1];
    }

    free(cp);
    free(dp);
}

void thermal_solve(SphericalGrid *grid, const ThermalParams *tp,
                   const double *flux, int nthreads) {
    int ftd = grid->ftd;
    int ffd = grid->ffd;

    double P = CONST_P_CERES;
    int n_sub = tp->n_substep;
    double dt = P / ffd / n_sub;
    double kk = tp->conductivity;
    double lo = tp->density;
    double Cp = tp->Cp;

    double courant = 0.1;
    double dr = sqrt(courant * 2.0 * dt * kk / lo / Cp);

    int h_depth = (int)ceil(tp->depth / dr);
    if (h_depth < 3) h_depth = 3;

    printf("Thermal model: dr=%.6f m, depth points=%d\n", dr, h_depth);

    double R = kk * dt / lo / 2.0 / Cp / (dr * dr);
    double k_over_dr = kk / dr;

    /* Compute solar flux if not provided */
    double *local_flux = NULL;
    if (!flux) {
        double De = 1.0;
        double Si = CONST_SOLAR_1AU * (De / tp->distance_AU) *
                    (De / tp->distance_AU);
        double dF = Si * (1.0 - tp->albedo);

        local_flux = calloc((size_t)ftd * ffd, sizeof(double));
        for (int td = 0; td < ftd; td++) {
            for (int fd = 0; fd < ffd; fd++) {
                double thi = (td + 0.5) / ftd * M_PI - M_PI / 2.0;
                double fi = (fd + 0.5) / ffd * 2.0 * M_PI;
                if (fi >= M_PI) {
                    local_flux[td * ffd + fd] = 0.0;
                } else {
                    local_flux[td * ffd + fd] = dF * cos(thi) * sin(fi);
                }
            }
        }
        flux = local_flux;
    }

    /* Build Crank-Nicolson matrix components.
     * For each (td, fd) we solve the same tridiagonal system,
     * but the surface BC changes. We use Thomas algorithm directly
     * rather than inverting a matrix like MATLAB does with inv(a)*b.
     *
     * Interior: implicit CN scheme
     *   -R*T(i-1,n+1) + (1+2R)*T(i,n+1) - R*T(i+1,n+1)
     *     = R*T(i-1,n) + (1-2R)*T(i,n) + R*T(i+1,n)
     *
     * Bottom BC: zero flux  => T(h) = T(h-1)
     * Top: set by bisection then fixed
     */

    /* Initialize temperature field */
    /* T[td * ffd * h_depth + fd * h_depth + depth_idx] */
    int T_size = ftd * ffd * h_depth;
    double *T = malloc(T_size * sizeof(double));
    for (int i = 0; i < T_size; i++) T[i] = 220.0;

    /* Better initial guess based on flux */
    for (int td = 0; td < ftd; td++) {
        for (int fd = 0; fd < ffd; fd++) {
            double F = fabs(flux[td * ffd + fd]);
            double T_init = (F > 0) ? pow(F / CONST_SIGMA, 0.25) : 100.0;
            for (int d = 0; d < h_depth; d++) {
                T[td * ffd * h_depth + fd * h_depth + d] = T_init;
            }
        }
    }

    /* Run thermal diffusion */
    printf("Running thermal diffusion...\n");

    #pragma omp parallel for collapse(2) num_threads(nthreads) schedule(dynamic)
    for (int td = 0; td < ftd; td++) {
        for (int turn = 0; turn < 1; turn++) {  /* single revolution */
            /* Working arrays for this (td, *) column */
            double *qT = malloc(h_depth * sizeof(double));
            double *a_lo = malloc(h_depth * sizeof(double));
            double *a_di = malloc(h_depth * sizeof(double));
            double *a_up = malloc(h_depth * sizeof(double));
            double *rhs = malloc(h_depth * sizeof(double));

            for (int fd = 0; fd < ffd; fd++) {
                /* Load current depth profile */
                int base = td * ffd * h_depth + fd * h_depth;
                for (int d = 0; d < h_depth; d++) {
                    qT[d] = T[base + d];
                }

                double F = fabs(flux[td * ffd + fd]);

                for (int hh = 0; hh < n_sub; hh++) {
                    double oT = qT[1];  /* temperature just below surface */

                    /* Bisection for surface temperature */
                    double T0 = 50.0, T1 = 300.0;
                    double fT0 = surface_balance(T0, F, k_over_dr, oT, tp);
                    double fT1 = surface_balance(T1, F, k_over_dr, oT, tp);

                    /* Expand range if needed */
                    for (int m = 0; m < 30; m++) {
                        if (T1 >= 10000.0) break;
                        if (fT0 * fT1 > 0) {
                            T1 *= 5.0;
                            fT1 = surface_balance(T1, F, k_over_dr, oT, tp);
                            continue;
                        }
                        double T2 = (T0 + T1) / 2.0;
                        double fT2 = surface_balance(T2, F, k_over_dr, oT, tp);
                        if (fT0 * fT2 <= 0) {
                            T1 = T2; fT1 = fT2;
                        } else {
                            T0 = T2; fT0 = fT2;
                        }
                        if (fabs(T0 - T1) <= 0.01) break;
                    }

                    qT[0] = (T0 + T1) / 2.0;

                    /* Apply Crank-Nicolson for interior + bottom BC */
                    /* Build RHS: b * qT_old */
                    rhs[0] = qT[0];  /* surface is fixed by bisection */
                    for (int i = 1; i < h_depth - 1; i++) {
                        rhs[i] = R * qT[i - 1] + (1.0 - 2.0 * R) * qT[i] +
                                 R * qT[i + 1];
                    }
                    /* Bottom: zero flux => T(h) = T(h-1) */
                    rhs[h_depth - 1] = 0.0;

                    /* Build tridiagonal for implicit side */
                    a_lo[0] = 0.0;
                    a_di[0] = 1.0;
                    a_up[0] = 0.0;
                    for (int i = 1; i < h_depth - 1; i++) {
                        a_lo[i] = -R;
                        a_di[i] = 1.0 + 2.0 * R;
                        a_up[i] = -R;
                    }
                    a_lo[h_depth - 1] = -1.0 / dr;
                    a_di[h_depth - 1] = 1.0 / dr;
                    a_up[h_depth - 1] = 0.0;

                    thomas_solve(a_lo, a_di, a_up, rhs, qT, h_depth);
                }

                /* Store result */
                for (int d = 0; d < h_depth; d++) {
                    T[base + d] = qT[d];
                }
            }

            free(qT);
            free(a_lo);
            free(a_di);
            free(a_up);
            free(rhs);
        }
    }

    /* Extract surface temperature */
    for (int td = 0; td < ftd; td++) {
        for (int fd = 0; fd < ffd; fd++) {
            grid->temperature[td * ffd + fd] =
                T[td * ffd * h_depth + fd * h_depth];
        }
    }

    printf("Thermal model done.\n");

    free(T);
    free(local_flux);
}
