#include "collision.h"
#include <math.h>

int collision_check(double *real_col_p, double col_threshold,
                    double ds, double local_density, double cross_section) {
    /*
     * Accumulate collision probability:
     *   real_col_p += ds * density * cross_section
     *   if (col_threshold < 1 - exp(-real_col_p)) -> collision
     */
    *real_col_p += ds * local_density * cross_section;
    return (col_threshold < 1.0 - exp(-(*real_col_p)));
}

void collision_execute(Vec3 *v, Vec3 v2) {
    /*
     * Elastic hard-sphere collision in center-of-mass frame.
     * Equal mass assumption (same molecule).
     *
     * MATLAB:
     *   del_v = v - v2
     *   vcm = v + v2
     *   v = vcm + 0.5 * del_v    (= (v + v2)/2 + (v - v2)/2 = v)
     *
     * Wait, the MATLAB code actually does:
     *   del_v(i) = v(i) - v2(i)
     *   vcm(i) = v(i) + v2(i)
     *   v(i) = vcm(i) + 0.5 * del_v(i)
     *
     * This equals: v(i) = v(i) + v2(i) + 0.5*(v(i) - v2(i))
     *            = 1.5*v(i) + 0.5*v2(i)
     *
     * This is the MATLAB formula as written. We replicate it exactly.
     */
    Vec3 del_v = vec3_sub(*v, v2);
    Vec3 vcm = vec3_add(*v, v2);

    v->x = vcm.x + 0.5 * del_v.x;
    v->y = vcm.y + 0.5 * del_v.y;
    v->z = vcm.z + 0.5 * del_v.z;
}

double collision_sample_speed(const PolyFit *pf, double rand_val, int nn) {
    /*
     * Sample from polynomial CDF using bisection.
     * MATLAB:
     *   df = polyval(p2, max_v);
     *   ranking = rand * (nn - df) + df;
     *   then bisect to find v where polyval(p2, v) == ranking
     */
    if (!pf->valid) return 0.0;

    double min_v = pf->min_v;
    double max_v = pf->max_v;

    /* Evaluate polynomial at endpoints */
    double val_at_max = 0.0;
    for (int i = 0; i <= POLY_DEGREE; i++) {
        val_at_max = val_at_max * max_v + pf->coeffs[i];
    }

    double ranking = rand_val * ((double)nn - val_at_max) + val_at_max;

    /* Bisection */
    double v0 = min_v, v1 = max_v;

    double f0 = 0.0;
    {
        double x = v0;
        for (int i = 0; i <= POLY_DEGREE; i++) f0 = f0 * x + pf->coeffs[i];
        f0 -= ranking;
    }

    for (int iter = 0; iter < 100; iter++) {
        double vm = (v0 + v1) * 0.5;
        double fm = 0.0;
        {
            double x = vm;
            for (int i = 0; i <= POLY_DEGREE; i++) fm = fm * x + pf->coeffs[i];
            fm -= ranking;
        }

        double f1 = 0.0;
        {
            double x = v1;
            for (int i = 0; i <= POLY_DEGREE; i++) f1 = f1 * x + pf->coeffs[i];
            f1 -= ranking;
        }

        if (fm * f1 < 0) {
            f0 = fm;
            v0 = vm;
        } else {
            v1 = vm;
        }

        if (fabs(v1 - v0) < 0.01) break;
    }

    return (v0 + v1) * 0.5;
}
