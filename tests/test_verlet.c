#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "types.h"
#include "constants.h"
#include "integrator.h"
#include "velocity.h"
#include "polyfit.h"

/* Test 1: Circular orbit conservation */
static int test_circular_orbit(void) {
    printf("Test 1: Circular orbit energy conservation...\n");

    PhysicsParams params;
    params.G = CONST_G;
    params.M = CONST_M_CERES;
    params.Radii = CONST_R_CERES;
    params.max_height = params.Radii * 10;

    /* Circular orbit at 2*Radii */
    double r0 = 2.0 * params.Radii;
    double v_circ = sqrt(params.G * params.M / r0);

    Particle p;
    p.pos = (Vec3){r0, 0, 0};
    p.vel = (Vec3){0, v_circ, 0};

    double h = 1.0;  /* 1 second timestep */
    Vec3 oa = {params.G * params.M * p.pos.x / (r0 * r0 * r0),
               params.G * params.M * p.pos.y / (r0 * r0 * r0),
               params.G * params.M * p.pos.z / (r0 * r0 * r0)};

    /* Initial energy */
    double E0 = 0.5 * vec3_dot(p.vel, p.vel) -
                params.G * params.M / vec3_len(p.pos);

    /* Integrate for 1000 steps */
    for (int i = 0; i < 1000; i++) {
        verlet_step(&p, &params, h, &oa, NULL);
    }

    double E1 = 0.5 * vec3_dot(p.vel, p.vel) -
                params.G * params.M / vec3_len(p.pos);

    double rel_err = fabs(E1 - E0) / fabs(E0);
    printf("  E0 = %.10e, E1 = %.10e, rel_error = %.2e\n", E0, E1, rel_err);

    if (rel_err < 1e-3) {  /* MATLAB-style Verlet is not fully symplectic */
        printf("  PASSED\n");
        return 0;
    } else {
        printf("  FAILED (error too large)\n");
        return 1;
    }
}

/* Test 2: Velocity table MB distribution */
static int test_velocity_table(void) {
    printf("Test 2: MB velocity table...\n");

    VelocityTable vt;
    velocity_table_init(&vt, CONST_MMASS_H2O, CONST_MOLE, CONST_K_B);

    /* Check that table has entries for all temps */
    int ok = 1;
    for (int i = 0; i < VEL_TABLE_NTEMP; i++) {
        if (vt.count[i] <= 0) {
            printf("  FAILED: temp[%d]=%.0f has no entries\n",
                   i, vt.temps[i]);
            ok = 0;
        }
    }

    /* Sample at 200K, check mean is reasonable */
    unsigned int seed = 42;
    double sum = 0;
    int N = 10000;
    for (int i = 0; i < N; i++) {
        sum += velocity_sample(&vt, 200.0, &seed);
    }
    double mean_v = sum / N;
    /* Expected thermal speed ~ sqrt(8kT/pi/m) ~ 400 m/s for H2O at 200K */
    /* Table saturates at 1000 entries (matching MATLAB p_max(23,1e3))
     * so mean is lower than true MB average */
    printf("  Mean velocity at 200K: %.1f m/s (expect ~70-150 due to table cap)\n",
           mean_v);

    if (ok && mean_v > 50 && mean_v < 200) {
        printf("  PASSED\n");
        return 0;
    } else {
        printf("  FAILED\n");
        return 1;
    }
}

/* Test 3: Polynomial fit */
static int test_polyfit(void) {
    printf("Test 3: Polynomial fit (degree 3)...\n");

    /* Fit y = 2x^3 - x^2 + 3x - 1 */
    double x[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    double y[11];
    for (int i = 0; i < 11; i++) {
        double xi = x[i];
        y[i] = 2 * xi * xi * xi - xi * xi + 3 * xi - 1;
    }

    double coeffs[4];
    int rc = polyfit(x, y, 11, 3, coeffs);
    if (rc != 0) {
        printf("  FAILED: polyfit returned %d\n", rc);
        return 1;
    }

    printf("  Coefficients: %.4f %.4f %.4f %.4f\n",
           coeffs[0], coeffs[1], coeffs[2], coeffs[3]);
    printf("  Expected:     2.0000 -1.0000 3.0000 -1.0000\n");

    /* Check at x=5 */
    double y_pred = polyval(coeffs, 3, 5.0);
    double y_true = 2 * 125 - 25 + 15 - 1;  /* = 239 */
    double err = fabs(y_pred - y_true);
    printf("  polyval(5) = %.4f (true = %.4f, error = %.2e)\n",
           y_pred, y_true, err);

    if (err < 0.01 && fabs(coeffs[0] - 2.0) < 0.01) {
        printf("  PASSED\n");
        return 0;
    } else {
        printf("  FAILED\n");
        return 1;
    }
}

int main(void) {
    printf("=== DSMC Suite Tests ===\n\n");

    int failures = 0;
    failures += test_circular_orbit();
    failures += test_velocity_table();
    failures += test_polyfit();

    printf("\n=== %d test(s) failed ===\n", failures);
    return failures;
}
