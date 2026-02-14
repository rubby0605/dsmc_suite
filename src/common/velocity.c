#include "velocity.h"
#include "constants.h"
#include <math.h>
#include <stdlib.h>

void velocity_table_init(VelocityTable *vt, double mmass, double mole,
                         double k_B) {
    /*
     * Build lookup table: 23 temperature points from 20K to 240K.
     * For each T, sample MB distribution into a speed array.
     *
     * MATLAB logic:
     *   for MIV_T = linspace(20,240,23)
     *     for vv = linspace(0,1000,101)
     *       aa = 4*pi*vv^2 * sqrt((1/(2*pi)*mmass/mole/k/T)^3)
     *            * exp(-mmass/mole*3*vv^2/(2*k*T))
     *       aa = aa * 5e5
     *       fill aa copies of vv into table
     */
    double dT = (240.0 - 20.0) / (VEL_TABLE_NTEMP - 1);
    double m = mmass / mole;  /* molecular mass */

    for (int i = 0; i < VEL_TABLE_NTEMP; i++) {
        double T = 20.0 + i * dT;
        vt->temps[i] = T;
        int n_v = 0;

        for (int j = 0; j <= 100; j++) {
            double vv = j * 10.0;  /* 0 to 1000 in steps of 10 */
            double factor = 1.0 / (2.0 * M_PI) * m / (k_B * T);
            double aa = 4.0 * M_PI * vv * vv
                        * pow(factor, 1.5)
                        * exp(-m * 3.0 * vv * vv / (2.0 * k_B * T));
            aa *= 5e5;
            int count = (int)floor(aa);
            if (count == 0) continue;

            /* Fill 'count' entries with this speed */
            for (int k = 0; k < count && n_v < VEL_TABLE_NVEL; k++) {
                vt->speeds[i][n_v] = vv;
                n_v++;
            }
        }
        if (n_v == 0) n_v = 1;  /* safety */
        vt->count[i] = n_v;
    }
}

double velocity_sample(const VelocityTable *vt, double T_local,
                       unsigned int *seed) {
    /*
     * MATLAB logic:
     *   if T >= 220: use index 21 (0-based: 20)
     *   else: index = ceil(T/10) - 2 (MATLAB 1-based)
     * In C (0-based): index = (int)(T/10) - 2, clamped
     */
    int idx;
    if (T_local >= 220.0) {
        idx = 20;  /* corresponds to MATLAB index 21 */
    } else {
        idx = (int)ceil(T_local / 10.0) - 2;
    }
    if (idx < 0) idx = 0;
    if (idx >= VEL_TABLE_NTEMP) idx = VEL_TABLE_NTEMP - 1;

    int r = rand_r(seed) % vt->count[idx];
    return vt->speeds[idx][r];
}

double velocity_sticking(double vr, double v_old_sq, double alpha) {
    /* vout = sqrt(alpha * vr^2 + (1-alpha) * v_old^2) */
    return sqrt(alpha * vr * vr + (1.0 - alpha) * v_old_sq);
}
