#include "types.h"
#include "constants.h"
#include "grid.h"
#include "velocity.h"
#include "integrator.h"
#include "collision.h"
#include "polyfit.h"
#include "fileio.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

/*
 * Advanced DSMC with exponential grid + collisions + polynomial velocity fit.
 * Converted from test_ballistic.m
 *
 * Two phases:
 *   Phase 1: Collisionless ballistic run, record v_record
 *   Phase 2: Build polynomial CDF fits from v_record, run with collisions
 */

static void pos_to_angles(Vec3 s, double Ri, double *thi, double *fi) {
    *thi = asin(s.z / Ri);
    *fi = atan2(s.x, s.y);
    if (*fi < 0) *fi += 2.0 * M_PI;
}

/* Build polynomial fits from velocity records */
static void build_poly_fits(PolyFit *poly_mat, int *count_num,
                            double *new_v_record, int ftd, int ffd, int hd,
                            int max_per_cell, int nthreads) {
    int num_poly = POLY_DEGREE;

    #pragma omp parallel for collapse(3) num_threads(nthreads) schedule(dynamic)
    for (int td = 0; td < ftd; td++) {
        for (int fd = 0; fd < ffd; fd++) {
            for (int ht = 0; ht < hd; ht++) {
                int idx3 = td * ffd * hd + fd * hd + ht;
                int nn = count_num[idx3];

                if (nn <= 8) {
                    poly_mat[idx3].valid = 0;
                    continue;
                }

                /* Extract and sort velocities */
                double *v_dist = malloc(nn * sizeof(double));
                for (int i = 0; i < nn; i++) {
                    v_dist[i] = new_v_record[idx3 * max_per_cell + i];
                }

                /* Bubble sort (small arrays) */
                for (int i = 0; i < nn; i++) {
                    for (int j = i + 1; j < nn; j++) {
                        if (v_dist[j] < v_dist[i]) {
                            double tmp = v_dist[j];
                            v_dist[j] = v_dist[i];
                            v_dist[i] = tmp;
                        }
                    }
                }

                double min_v = v_dist[0];
                double max_v = v_dist[nn - 1];

                /* Build histogram and cumulative */
                int nbins = (nn + 1) / 2;
                if (nbins < 4) nbins = 4;
                double *xx = malloc(nbins * sizeof(double));
                double *yy = calloc(nbins, sizeof(double));
                double *yy2 = calloc(nbins, sizeof(double));

                double dv = (max_v - min_v) / nbins;
                for (int i = 0; i < nbins; i++) {
                    xx[i] = min_v + (i + 0.5) * dv;
                }

                /* Histogram */
                for (int i = 0; i < nn; i++) {
                    int bin = (int)((v_dist[i] - min_v) / dv);
                    if (bin >= nbins) bin = nbins - 1;
                    if (bin < 0) bin = 0;
                    yy[bin]++;
                }

                /* Reverse cumulative */
                for (int i = 0; i < nbins; i++) {
                    yy2[i] = 0;
                    for (int j = i; j < nbins; j++) {
                        yy2[i] += yy[j];
                    }
                }

                /* Polynomial fit */
                double coeffs[POLY_DEGREE + 1];
                if (polyfit(xx, yy2, nbins, num_poly, coeffs) == 0) {
                    for (int i = 0; i <= num_poly; i++) {
                        poly_mat[idx3].coeffs[i] = coeffs[i];
                    }
                    poly_mat[idx3].min_v = min_v;
                    poly_mat[idx3].max_v = max_v;
                    poly_mat[idx3].valid = 1;
                } else {
                    poly_mat[idx3].valid = 0;
                }

                free(v_dist);
                free(xx);
                free(yy);
                free(yy2);
            }
        }
    }
}

void run_advanced(int totaln, const char *temp_file,
                  const char *prev_den_file, const char *prev_v_record_file,
                  const char *output_prefix, int nthreads) {
    /* Read temperature (surface only for test_ballistic style) */
    int ftd, ffd;
    double *T_surface = read_temperature_dat(temp_file, &ftd, &ffd);
    if (!T_surface) {
        fprintf(stderr, "Failed to read temperature file\n");
        return;
    }

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
    params.dt_step = 0.5;  /* smaller for collision accuracy */
    params.dt_rot = params.P / ffd;
    params.cross_section = CONST_CROSS_SECTION;
    params.use_exp_grid = 1;

    /* Exponential grid setup */
    int hd = 50;
    params.max_height = 3e6 * sqrt(2.0);
    params.alpha2 = 0.04;
    params.unih = (params.max_height - params.Radii) /
                  (pow(1.0 + params.alpha2, hd) - 1.0) * params.alpha2;

    double h = params.dt_step;

    /* Build velocity table */
    VelocityTable vt;
    velocity_table_init(&vt, params.mmass, params.mole, params.k_B);

    /* Read previous density field (from Phase 1) */
    int prev_ftd, prev_ffd, prev_hd;
    double *copy_atm = NULL;
    if (prev_den_file) {
        copy_atm = read_density_dat(prev_den_file, &prev_ftd, &prev_ffd,
                                    &prev_hd, 1);
    }

    /* Read or build velocity records */
    int max_per_cell = 50;
    (void)0;  /* v_rec dimensions match ftd/ffd */

    /* count_num[td][fd][ht] */
    int *count_num = calloc((size_t)ftd * ffd * hd, sizeof(int));
    /* new_v_record[td][fd][ht][sample] -- flat */
    double *new_v_record = calloc((size_t)ftd * ffd * hd * max_per_cell,
                                  sizeof(double));

    if (prev_v_record_file) {
        /* Read velocity records from previous phase */
        double *v_rec_data;
        int *td_list, *fd_list, *ht_list;
        int nrecords;
        if (read_v_record_dat(prev_v_record_file, &v_rec_data,
                              &td_list, &fd_list, &ht_list, &nrecords) == 0) {
            printf("Read %d velocity records\n", nrecords);

            for (int i = 0; i < nrecords; i++) {
                int td = td_list[i];
                int fd = fd_list[i];
                int ht = ht_list[i];
                if (td < 0 || td >= ftd || fd < 0 || fd >= ffd ||
                    ht < 0 || ht >= hd) continue;

                int idx3 = td * ffd * hd + fd * hd + ht;
                if (count_num[idx3] < max_per_cell) {
                    double speed = sqrt(v_rec_data[i * 3] * v_rec_data[i * 3] +
                                        v_rec_data[i * 3 + 1] * v_rec_data[i * 3 + 1] +
                                        v_rec_data[i * 3 + 2] * v_rec_data[i * 3 + 2]);
                    new_v_record[idx3 * max_per_cell + count_num[idx3]] = speed;
                    count_num[idx3]++;
                }
            }

            free(v_rec_data);
            free(td_list);
            free(fd_list);
            free(ht_list);
        }
    }

    /* Build polynomial fits */
    printf("Building polynomial velocity fits...\n");
    PolyFit *poly_mat = calloc((size_t)ftd * ffd * hd, sizeof(PolyFit));
    build_poly_fits(poly_mat, count_num, new_v_record, ftd, ffd, hd,
                    max_per_cell, nthreads);

    /* Count valid fits */
    int valid_fits = 0;
    for (int i = 0; i < ftd * ffd * hd; i++) {
        if (poly_mat[i].valid) valid_fits++;
    }
    printf("Valid polynomial fits: %d / %d\n", valid_fits, ftd * ffd * hd);

    /* Phase 2: Run with collisions */
    SphericalGrid *grid = grid_create(ftd, ffd, hd);
    grid_compute_volume_exp(grid, params.Radii, params.unih, params.alpha2);

    int total_escaped = 0;
    int total_collisions = 0;

    double *ttime = calloc(9 * totaln, sizeof(double));
    int *bomb_line = calloc(totaln, sizeof(int));

    printf("Starting advanced DSMC: %d particles with collisions\n", totaln);

    #pragma omp parallel num_threads(nthreads)
    {
        unsigned int seed = (unsigned int)(omp_get_thread_num() * 54321 + 98765);

        #pragma omp for schedule(dynamic, 100) reduction(+:total_collisions)
        for (int n = 0; n < totaln; n++) {
            int ou = 0;
            int rot = 0;
            int num_col = 0;
            double real_col_p = 0.0;
            double col_p = (double)rand_r(&seed) / RAND_MAX;
            double time = 0.0;
            double W = 1.0;
            double solar_time = 0.0;
            double flight_time = 0.0;
            double rest_time = 0.0;
            double del_time = 0.0;

            double thi = (double)rand_r(&seed) / RAND_MAX * M_PI - M_PI / 2.0;
            double fi = (double)rand_r(&seed) / RAND_MAX * 2.0 * M_PI;
            Vec3 s = {
                params.Radii * cos(thi) * sin(fi),
                params.Radii * cos(thi) * cos(fi),
                params.Radii * sin(thi)
            };

            ttime[0 * totaln + n] = thi;
            ttime[5 * totaln + n] = fi;

            int td, fd;
            grid_spherical_to_index(thi, fi, ftd, ffd, &td, &fd);
            ttime[6 * totaln + n] = T_surface[td * ffd + fd];

            Vec3 v = {0, 0, 0};

            while (W >= 0.01) {
                double KF = 2.0 * M_PI * params.Radii * fabs(cos(thi)) / params.P;
                grid_spherical_to_index(thi, fi, ftd, ffd, &td, &fd);

                double local_T = T_surface[td * ffd + fd];

                /* Frozen surface handling */
                if (local_T < 100.0) {
                    #pragma omp atomic
                    grid->horizon[td * ffd + fd] += W;

                    int found_warm = 0;
                    for (int add = 1; add <= ffd; add++) {
                        rot++;
                        time += params.dt_rot;
                        rest_time += params.dt_rot;

                        if (T_surface[td * ffd + (add % ffd)] >= 100.0) {
                            double rot_fi = (double)(add - fd) / ffd * 2.0 * M_PI;
                            double cs = cos(rot_fi), sn = sin(rot_fi);
                            Vec3 ns = {
                                s.x * cs - s.y * sn,
                                s.x * sn + s.y * cs,
                                s.z
                            };
                            s = ns;
                            fd = add % ffd;
                            fi = (fd + 0.5) / ffd * 2.0 * M_PI;
                            found_warm = 1;
                            break;
                        }
                    }
                    if (!found_warm) {
                        #pragma omp atomic
                        grid->horizon[td * ffd + fd] += W;
                        break;
                    }
                    local_T = T_surface[td * ffd + fd];
                }

                /* Sample velocity */
                double vr = velocity_sample(&vt, local_T, &seed);
                double v_old_sq = v.x * v.x + v.y * v.y + v.z * v.z;
                double vout = velocity_sticking(vr, v_old_sq, params.alpha);

                double flyd = (double)rand_r(&seed) / RAND_MAX * M_PI + thi - M_PI / 2.0;
                double xi = (double)rand_r(&seed) / RAND_MAX * M_PI + fi - M_PI / 2.0;

                v.x = cos(flyd) * sin(xi) * vout;
                v.y = cos(flyd) * cos(xi) * vout;
                v.z = sin(flyd) * vout;

                v.x += KF * sin(fi + M_PI / 2.0);
                v.y += KF * cos(fi + M_PI / 2.0);

                if (vec3_len(v) >= params.ve) { ou = 1; break; }

                double Ri = vec3_len(s);
                Vec3 oa = {
                    params.G * params.M * s.x / (Ri * Ri * Ri),
                    params.G * params.M * s.y / (Ri * Ri * Ri),
                    params.G * params.M * s.z / (Ri * Ri * Ri)
                };

                Particle p;
                p.pos = s;
                p.vel = v;

                /* Orbit integration with collisions */
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

                    /* Photodissociation */
                    if (fd < ffd / 2) {
                        W *= exp(-h / params.lifetime);
                        solar_time += h;
                    }

                    if (ii == 0) continue;

                    /* Exponential height index */
                    int ht = grid_exp_ht_index(Ri, params.Radii,
                                               params.unih, params.alpha2);
                    if (ht < 0) ht = 0;
                    if (ht >= hd) continue;

                    /* --- Collision check --- */
                    if (copy_atm && ht < hd / 3) {
                        int atm_idx = grid_index(ht, td, fd, ftd, ffd);
                        real_col_p += ds * copy_atm[atm_idx] *
                                      params.cross_section;

                        int poly_idx = td * ffd * hd + fd * hd + ht;
                        if (col_p < 1.0 - exp(-real_col_p) &&
                            poly_mat[poly_idx].valid) {
                            real_col_p = 0.0;
                            num_col++;

                            /* Sample collision partner velocity */
                            double rand_val = (double)rand_r(&seed) / RAND_MAX;
                            double vv2 = collision_sample_speed(
                                &poly_mat[poly_idx], rand_val,
                                count_num[poly_idx]);

                            /* Random direction */
                            double thi2 = M_PI * (double)rand_r(&seed) / RAND_MAX - M_PI / 2.0;
                            double fi2 = 2.0 * M_PI * (double)rand_r(&seed) / RAND_MAX;
                            Vec3 v2 = {
                                vv2 * cos(thi2) * cos(fi2),
                                vv2 * cos(thi2) * sin(fi2),
                                vv2 * sin(thi2)
                            };

                            collision_execute(&p.vel, v2);
                            continue;
                        }
                    }

                    /* Accumulate density */
                    #pragma omp atomic
                    grid->density[grid_index(ht, td, fd, ftd, ffd)] += W;
                }

                s = p.pos;
                v = p.vel;

                if (ou) break;

                Ri = vec3_len(s);
                pos_to_angles(s, Ri, &thi, &fi);
                grid_spherical_to_index(thi, fi, ftd, ffd, &td, &fd);

                del_time += h * 10000000;

                if (td <= 1 || td >= ftd - 2) {
                    #pragma omp atomic
                    grid->horizon[td * ffd + fd] += 0.05 * W;
                    W *= 0.95;
                }

                if (del_time >= params.dt_rot) {
                    int del_n = (int)(del_time / params.dt_rot);
                    rot += del_n;
                    fd = (fd - del_n % ffd + ffd) % ffd;
                    del_time -= params.dt_rot * del_n;
                }
            }

            bomb_line[n] = num_col;
            total_collisions += num_col;

            if (ou) {
                #pragma omp atomic
                total_escaped++;
                continue;
            }

            ttime[1 * totaln + n] = thi;
            ttime[2 * totaln + n] = time;
            ttime[3 * totaln + n] = W;
            ttime[4 * totaln + n] = flight_time;
            ttime[7 * totaln + n] = rest_time;
            ttime[8 * totaln + n] = (double)num_col;

            if (n % 1000 == 0) {
                printf("  Particle %d/%d (collisions so far: %d)\n",
                       n, totaln, num_col);
            }
        }
    }

    printf("Advanced DSMC complete. Escaped: %d/%d, Total collisions: %d\n",
           total_escaped, totaln, total_collisions);

    /* Write outputs */
    char fname[256];
    snprintf(fname, sizeof(fname), "%s_den.dat", output_prefix);
    write_density_dat(fname, grid->density, ftd, ffd, hd);

    snprintf(fname, sizeof(fname), "%s_horizon.dat", output_prefix);
    write_horizon_dat(fname, grid->horizon, ftd, ffd);

    snprintf(fname, sizeof(fname), "%s_ttime.dat", output_prefix);
    write_ttime_dat(fname, ttime, 9, totaln);

    /* Cleanup */
    grid_free(grid);
    free(poly_mat);
    free(count_num);
    free(new_v_record);
    free(ttime);
    free(bomb_line);
    free(T_surface);
    free(copy_atm);
}
