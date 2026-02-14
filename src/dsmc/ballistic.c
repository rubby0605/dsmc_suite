#include "types.h"
#include "constants.h"
#include "grid.h"
#include "velocity.h"
#include "integrator.h"
#include "fileio.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

/*
 * Basic DSMC ballistic motion simulation.
 * Converted from m2_ballistic_motion_correct0227.m
 *
 * Each particle:
 *   1. Starts at random surface location
 *   2. If T < 100K, wait (frozen) until rotation brings to warm side
 *   3. Sample MB velocity, add rotational velocity
 *   4. Integrate orbit with Störmer-Verlet
 *   5. Photodissociation decay on dayside (fd <= ffd/2)
 *   6. On landing: check polar deposit, apply rotation time offset, repeat
 *   7. Record velocity for collision model, accumulate density
 */

/* Convert position to spherical angles using atan2 (fixed quadrant) */
static void pos_to_angles(Vec3 s, double Ri, double *thi, double *fi) {
    *thi = asin(s.z / Ri);
    *fi = atan2(s.x, s.y);
    if (*fi < 0) *fi += 2.0 * M_PI;
}

void run_ballistic(int totaln, const char *temp_file, const char *output_prefix,
                   int nthreads) {
    /* Read temperature */
    int ftd, ffd;
    double *T_surface = read_temperature_dat(temp_file, &ftd, &ffd);
    if (!T_surface) {
        fprintf(stderr, "Failed to read temperature file\n");
        return;
    }

    /* Also build time-shifted temperature (for m2 style with nt index) */
    /* For the basic version we use surface T directly (test_ballistic style) */

    /* Physics parameters */
    PhysicsParams params;
    params.G = CONST_G;
    params.M = CONST_M_CERES;
    params.k_B = CONST_K_B;
    params.mmass = CONST_MMASS_H2O;
    params.mole = CONST_MOLE;
    params.Radii = CONST_R_CERES;
    params.ve = sqrt(2.0 * params.G * params.M / params.Radii);
    params.lifetime = CONST_LIFETIME;
    params.P = CONST_P_CERES;
    params.alpha = 0.5;
    params.dt_step = 5.0;
    params.dt_rot = params.P / ffd;
    params.use_exp_grid = 0;
    params.unih = params.Radii * 0.1;
    params.alpha2 = 0.0;

    int hd = 100;
    params.max_height = params.unih * hd + params.Radii;

    /* Create grid */
    SphericalGrid *grid = grid_create(ftd, ffd, hd);

    /* Build velocity table */
    VelocityTable vt;
    velocity_table_init(&vt, params.mmass, params.mole, params.k_B);

    /* Velocity record for collision model (phase 2) */
    int max_v_records = 100 * ftd * hd;
    double *v_record = calloc(max_v_records * 3, sizeof(double));
    int *v_td = calloc(max_v_records, sizeof(int));
    int *v_fd = calloc(max_v_records, sizeof(int));
    int *v_ht = calloc(max_v_records, sizeof(int));
    int *count_nn = calloc((size_t)ftd * hd, sizeof(int));

    /* ttime output: [5][totaln] */
    double *ttime = calloc(5 * totaln, sizeof(double));

    int total_records = 0;
    int total_escaped = 0;

    double h = params.dt_step;

    printf("Starting ballistic simulation: %d particles\n", totaln);
    printf("  ftd=%d, ffd=%d, hd=%d\n", ftd, ffd, hd);
    printf("  Escape velocity = %.2f m/s\n", params.ve);
    printf("  Lifetime = %.1f s\n", params.lifetime);

    #pragma omp parallel num_threads(nthreads)
    {
        unsigned int seed = (unsigned int)(omp_get_thread_num() * 12345 + 67890);

        #pragma omp for schedule(dynamic, 100)
        for (int n = 0; n < totaln; n++) {
            int ou = 0;
            int rot = 0;
            double time = 0.0;
            double W = 1.0;
            double solar_time = 0.0;
            double flight_time = 0.0;
            double del_time = 0.0;

            /* Random initial position on surface */
            double thi = (double)rand_r(&seed) / RAND_MAX * M_PI - M_PI / 2.0;
            double fi = (double)rand_r(&seed) / RAND_MAX * 2.0 * M_PI;
            Vec3 s = {
                params.Radii * cos(thi) * sin(fi),
                params.Radii * cos(thi) * cos(fi),
                params.Radii * sin(thi)
            };

            ttime[0 * totaln + n] = thi;

            Vec3 v = {0, 0, 0};

            while (W >= 0.01) {
                double KF = 2.0 * M_PI * params.Radii * fabs(cos(thi)) / params.P;
                int td, fd;
                grid_spherical_to_index(thi, fi, ftd, ffd, &td, &fd);

                double local_T = T_surface[td * ffd + fd];

                /* Frozen on cold surface: wait for rotation */
                if (local_T < 100.0) {
                    int found_warm = 0;
                    for (int add = 0; add < ffd; add++) {
                        rot++;
                        time += params.dt_rot;

                        int check_fd = (fd + add + 1) % ffd;
                        if (T_surface[td * ffd + check_fd] >= 100.0) {
                            fd = check_fd;
                            fi = (fd + 0.5) / ffd * 2.0 * M_PI;
                            found_warm = 1;
                            break;
                        }
                    }
                    if (!found_warm) {
                        /* Permanently frozen at poles */
                        #pragma omp atomic
                        grid->horizon[td * ffd + fd] += W;
                        break;
                    }
                    local_T = T_surface[td * ffd + fd];
                }

                /* Sample velocity from MB distribution */
                double vr = velocity_sample(&vt, local_T, &seed);
                double v_old_sq = v.x * v.x + v.y * v.y + v.z * v.z;
                double vout = velocity_sticking(vr, v_old_sq, params.alpha);

                /* Random launch direction */
                double flyd = (double)rand_r(&seed) / RAND_MAX * M_PI + thi - M_PI / 2.0;
                double xi = (double)rand_r(&seed) / RAND_MAX * M_PI + fi - M_PI / 2.0;

                v.x = cos(flyd) * sin(xi) * vout;
                v.y = cos(flyd) * cos(xi) * vout;
                v.z = sin(flyd) * vout;

                /* Add rotational velocity */
                v.x += KF * sin(fi + M_PI / 2.0);
                v.y += KF * cos(fi + M_PI / 2.0);

                /* Check escape */
                if (vec3_len(v) >= params.ve) { ou = 1; break; }

                /* Initialize Verlet */
                double Ri = vec3_len(s);
                Vec3 oa = {
                    params.G * params.M * s.x / (Ri * Ri * Ri),
                    params.G * params.M * s.y / (Ri * Ri * Ri),
                    params.G * params.M * s.z / (Ri * Ri * Ri)
                };

                Particle p;
                p.pos = s;
                p.vel = v;

                /* Orbit integration loop */
                for (int ii = 0; ii < 10000000; ii++) {
                    flight_time += h;
                    time += h;

                    double ds;
                    VerletStatus status = verlet_step(&p, &params, h, &oa, &ds);

                    if (status == VERLET_LANDED) break;
                    if (status == VERLET_ESCAPED) { ou = 1; break; }

                    Ri = vec3_len(p.pos);
                    pos_to_angles(p.pos, Ri, &thi, &fi);
                    grid_spherical_to_index(thi, fi, ftd, ffd, &td, &fd);

                    /* Photodissociation on dayside */
                    if (fd < ffd / 2) {
                        W *= exp(-h / params.lifetime);
                        solar_time += h;
                    }

                    if (ii == 0) continue;

                    /* Height index */
                    int ht = (int)ceil((Ri - params.Radii) / params.unih) - 1;
                    if (ht < 0) ht = 0;
                    if (ht >= hd) continue;

                    /* Record velocity (for later collision model) */
                    if (count_nn[td * hd + ht] < 200) {
                        int ridx;
                        #pragma omp critical
                        {
                            if (total_records < max_v_records &&
                                count_nn[td * hd + ht] < 200) {
                                ridx = total_records++;
                                v_record[ridx * 3] = p.vel.x;
                                v_record[ridx * 3 + 1] = p.vel.y;
                                v_record[ridx * 3 + 2] = p.vel.z;
                                v_td[ridx] = td;
                                v_fd[ridx] = fd;
                                v_ht[ridx] = ht;
                                count_nn[td * hd + ht]++;
                            }
                        }
                    }

                    /* Accumulate density */
                    #pragma omp atomic
                    grid->density[grid_index(ht, td, fd, ftd, ffd)] += W;
                }

                s = p.pos;
                v = p.vel;

                if (ou) break;

                /* Rotation time offset for landed particle */
                Ri = vec3_len(s);
                pos_to_angles(s, Ri, &thi, &fi);
                grid_spherical_to_index(thi, fi, ftd, ffd, &td, &fd);

                del_time += h * 10000000; /* approximate */

                /* Polar deposit */
                if (td <= 1 || td >= ftd - 2) {
                    #pragma omp atomic
                    grid->horizon[td * ffd + fd] += 0.05 * W;
                    W *= 0.95;
                }

                /* Rotation accounting */
                if (del_time >= params.dt_rot) {
                    int del_n = (int)(del_time / params.dt_rot);
                    rot += del_n;
                    fd = (fd - del_n % ffd + ffd) % ffd;
                    del_time -= params.dt_rot * del_n;
                }
            }

            if (ou) {
                #pragma omp atomic
                total_escaped++;
                continue;
            }

            ttime[1 * totaln + n] = thi;
            ttime[2 * totaln + n] = time;
            ttime[3 * totaln + n] = W;
            ttime[4 * totaln + n] = flight_time;

            if (n % 1000 == 0) {
                printf("  Particle %d/%d done\n", n, totaln);
            }
        }
    }

    printf("Simulation complete. Escaped: %d/%d\n", total_escaped, totaln);

    /* Write outputs */
    char fname[256];
    snprintf(fname, sizeof(fname), "%s_den.dat", output_prefix);
    write_density_dat(fname, grid->density, ftd, ffd, hd);

    snprintf(fname, sizeof(fname), "%s_horizon.dat", output_prefix);
    write_horizon_dat(fname, grid->horizon, ftd, ffd);

    snprintf(fname, sizeof(fname), "%s_v_record.dat", output_prefix);
    write_v_record_dat(fname, v_record, v_td, v_fd, v_ht, total_records);

    snprintf(fname, sizeof(fname), "%s_ttime.dat", output_prefix);
    write_ttime_dat(fname, ttime, 5, totaln);

    /* Cleanup */
    grid_free(grid);
    free(v_record);
    free(v_td);
    free(v_fd);
    free(v_ht);
    free(count_nn);
    free(ttime);
    free(T_surface);
}
