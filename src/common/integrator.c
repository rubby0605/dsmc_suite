#include "integrator.h"
#include <math.h>

VerletStatus verlet_step(Particle *p, const PhysicsParams *params,
                         double h, Vec3 *oa, double *ds_out) {
    /*
     * Störmer-Verlet integration from MATLAB lines 250-256:
     *   ns(i) = s(i) + v(i)*h - G*M*s(i)/(Ri^3)*(1.5*h^2) - oa(i)*h^2/6
     *   v(i)  = v(i) - G*M*s(i)/(Ri^3)*h
     *   oa(i) = G*M*s(i)/(Ri^3)
     */
    double GM = params->G * params->M;
    double Ri = vec3_len(p->pos);
    double Ri3 = Ri * Ri * Ri;
    double h2 = h * h;

    Vec3 ns;
    double gm_over_r3 = GM / Ri3;

    /* Position update */
    ns.x = p->pos.x + p->vel.x * h - gm_over_r3 * p->pos.x * (1.5 * h2)
           - oa->x * (h2 / 6.0);
    ns.y = p->pos.y + p->vel.y * h - gm_over_r3 * p->pos.y * (1.5 * h2)
           - oa->y * (h2 / 6.0);
    ns.z = p->pos.z + p->vel.z * h - gm_over_r3 * p->pos.z * (1.5 * h2)
           - oa->z * (h2 / 6.0);

    /* Displacement for collision probability */
    if (ds_out) {
        double dx = ns.x - p->pos.x;
        double dy = ns.y - p->pos.y;
        double dz = ns.z - p->pos.z;
        *ds_out = sqrt(dx * dx + dy * dy + dz * dz);
    }

    /* Velocity update */
    p->vel.x -= gm_over_r3 * p->pos.x * h;
    p->vel.y -= gm_over_r3 * p->pos.y * h;
    p->vel.z -= gm_over_r3 * p->pos.z * h;

    /* Update old acceleration */
    oa->x = gm_over_r3 * p->pos.x;
    oa->y = gm_over_r3 * p->pos.y;
    oa->z = gm_over_r3 * p->pos.z;

    /* Update position */
    p->pos = ns;

    /* Check boundaries */
    double Ri_new = vec3_len(p->pos);
    if (Ri_new <= params->Radii) return VERLET_LANDED;
    if (Ri_new >= params->max_height) return VERLET_ESCAPED;

    return VERLET_OK;
}
