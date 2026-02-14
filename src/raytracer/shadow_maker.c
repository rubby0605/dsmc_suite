#include "raytracer.h"
#include "constants.h"
#include "fileio.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

/*
 * MC shadow/flux computation.
 * Converted from 1_Shadow_maker(NObug0720).m
 *
 * For each face at each time step:
 *   1. Compute solar angle → base flux = dot(solar, normal)
 *   2. If face is illuminated, MC sample points on face
 *   3. For each sample, ray-trace to check shadow by neighbor faces
 *   4. Flux = (1 - shadow_fraction) * base_flux
 */

Mesh *mesh_load(const char *vertex_file, const char *face_file,
                int npoint, int nface) {
    Mesh *m = calloc(1, sizeof(Mesh));
    m->npoint = npoint;
    m->nface = nface;
    m->spot = calloc(npoint, sizeof(double[3]));
    m->facet = calloc(nface, sizeof(int[3]));
    m->face_normal = calloc(nface, sizeof(double[3]));
    m->face_area = calloc(nface, sizeof(double));
    m->loc_xyz = calloc(nface, sizeof(double[3]));
    m->dir = calloc(nface, sizeof(double[2]));

    /* Read vertices */
    FILE *fp = fopen(vertex_file, "r");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", vertex_file); return m; }
    for (int i = 0; i < npoint; i++) {
        fscanf(fp, "%lf %lf %lf", &m->spot[i][0], &m->spot[i][1],
               &m->spot[i][2]);
    }
    fclose(fp);

    /* Read faces (1-based indices in file) */
    fp = fopen(face_file, "r");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", face_file); return m; }
    for (int i = 0; i < nface; i++) {
        fscanf(fp, "%d %d %d", &m->facet[i][0], &m->facet[i][1],
               &m->facet[i][2]);
        /* Convert to 0-based */
        m->facet[i][0]--;
        m->facet[i][1]--;
        m->facet[i][2]--;
    }
    fclose(fp);

    return m;
}

void mesh_free(Mesh *m) {
    if (!m) return;
    free(m->spot);
    free(m->facet);
    free(m->face_normal);
    free(m->face_area);
    free(m->loc_xyz);
    free(m->dir);
    free(m);
}

void mesh_compute_normals(Mesh *m) {
    /* Normalize vertices to unit sphere */
    for (int i = 0; i < m->npoint; i++) {
        double len = sqrt(m->spot[i][0] * m->spot[i][0] +
                          m->spot[i][1] * m->spot[i][1] +
                          m->spot[i][2] * m->spot[i][2]);
        if (len > 0) {
            m->spot[i][0] /= len;
            m->spot[i][1] /= len;
            m->spot[i][2] /= len;
        }
    }

    for (int nf = 0; nf < m->nface; nf++) {
        Vec3 aa = {m->spot[m->facet[nf][0]][0], m->spot[m->facet[nf][0]][1],
                   m->spot[m->facet[nf][0]][2]};
        Vec3 bb = {m->spot[m->facet[nf][1]][0], m->spot[m->facet[nf][1]][1],
                   m->spot[m->facet[nf][1]][2]};
        Vec3 cc = {m->spot[m->facet[nf][2]][0], m->spot[m->facet[nf][2]][1],
                   m->spot[m->facet[nf][2]][2]};

        /* Centroid */
        m->loc_xyz[nf][0] = (aa.x + bb.x + cc.x) / 3.0;
        m->loc_xyz[nf][1] = (aa.y + bb.y + cc.y) / 3.0;
        m->loc_xyz[nf][2] = (aa.z + bb.z + cc.z) / 3.0;

        /* Edge vectors */
        Vec3 c1 = vec3_sub(bb, aa);
        Vec3 c2 = vec3_sub(cc, aa);

        /* Area (half cross product magnitude) */
        m->face_area[nf] = vec3_dot(c1, c2) / 2.0;

        /* Normal = cross(c2_hat, c1_hat) (outward convention from MATLAB) */
        Vec3 c2n = vec3_normalize(c2);
        Vec3 c1n = vec3_normalize(c1);
        Vec3 c3 = vec3_cross(c2n, c1n);
        Vec3 normal = vec3_normalize(c3);

        m->face_normal[nf][0] = normal.x;
        m->face_normal[nf][1] = normal.y;
        m->face_normal[nf][2] = normal.z;

        /* Direction angles */
        double thi = asin(normal.z);
        double fi = atan2(normal.x, normal.y);
        if (fi < 0) fi += 2.0 * M_PI;
        m->dir[nf][0] = thi;
        m->dir[nf][1] = fi;
    }
}

void shadow_find_neighbors(const Mesh *m, double za, int tt,
                           int **count_out, int ***neighbors_out,
                           int max_neighbors) {
    double initial_sc = 0.89;

    /* Allocate: count[tt * nface], neighbors[tt * nface][max_neighbors] */
    int total = tt * m->nface;
    *count_out = calloc(total, sizeof(int));
    *neighbors_out = calloc(total, sizeof(int *));
    for (int i = 0; i < total; i++) {
        (*neighbors_out)[i] = calloc(max_neighbors, sizeof(int));
    }

    for (int t = 0; t < tt; t++) {
        double fi = ((t + 0.5) / tt) * 2.0 * M_PI;
        Vec3 solar = {cos(-za) * sin(fi), cos(-za) * cos(fi), sin(-za)};

        for (int nf = 0; nf < m->nface; nf++) {
            int jj = 0;
            Vec3 s2 = {m->loc_xyz[nf][0], m->loc_xyz[nf][1], m->loc_xyz[nf][2]};

            for (int nnf = 0; nnf < m->nface && jj < max_neighbors; nnf++) {
                if (nnf == nf) continue;
                Vec3 s1 = {m->loc_xyz[nnf][0], m->loc_xyz[nnf][1],
                           m->loc_xyz[nnf][2]};
                Vec3 line_s = vec3_normalize(vec3_sub(s1, s2));
                if (vec3_dot(line_s, solar) >= initial_sc) {
                    (*neighbors_out)[t * m->nface + nf][jj] = nnf;
                    jj++;
                }
            }
            (*count_out)[t * m->nface + nf] = jj;
        }
        printf("  Neighbor search: time step %d/%d\n", t + 1, tt);
    }
}

Vec3 random_point_on_triangle(Vec3 a, Vec3 b, Vec3 c, unsigned int *seed) {
    /*
     * Uniform random point on triangle using the MATLAB method:
     * Handles the barycentric coordinate sampling with the
     * parallelogram folding approach from the original code.
     */
    Vec3 c1 = vec3_sub(a, b);
    Vec3 c2 = vec3_sub(c, b);
    double cosb = vec3_dot(vec3_normalize(c1), vec3_normalize(c2));

    /* Ensure proper orientation */
    if (cosb < 0) {
        Vec3 tmp = a; a = b; b = tmp;
        c1 = vec3_sub(a, b);
        c2 = vec3_sub(c, b);
        cosb = vec3_dot(vec3_normalize(c1), vec3_normalize(c2));
    }

    Vec3 c3 = vec3_scale(c2, cosb);
    Vec3 c4 = vec3_sub(c1, c3);
    Vec3 c5 = vec3_sub(c2, c3);

    double r1 = (double)rand_r(seed) / RAND_MAX;
    double r2 = (double)rand_r(seed) / RAND_MAX;
    double r3 = r1 - cosb;

    Vec3 pot;
    if (r3 > 0 && (r2 / (1.0 - r1) > vec3_len(c2) / vec3_len(c5))) {
        pot = vec3_add(vec3_add(vec3_scale(c2, 1.0 - r3),
                                vec3_scale(c4, 1.0 - r2)), b);
    } else if (r3 < 0 && (r2 / r1 > vec3_len(c2) / vec3_len(c3))) {
        pot = vec3_add(vec3_add(vec3_scale(c2, -r3),
                                vec3_scale(c4, 1.0 - r2)), b);
    } else {
        pot = vec3_add(vec3_add(vec3_scale(c2, r1),
                                vec3_scale(c4, r2)), b);
    }

    return pot;
}

int ray_triangle_intersect(Vec3 origin, Vec3 dir, Vec3 a, Vec3 b, Vec3 c) {
    /*
     * Check if ray from origin in direction dir hits triangle (a,b,c).
     * Uses the area method from the MATLAB code.
     */
    /* Face normal from vertices */
    Vec3 e1 = vec3_sub(b, a);
    Vec3 e2 = vec3_sub(c, a);
    Vec3 n = vec3_cross(e1, e2);
    double area = vec3_len(n);

    if (area < 1e-30) return 0;

    /* Plane equation: dot(normal, point) = dot(normal, a) */
    double ndir = vec3_dot(n, dir);
    if (fabs(ndir) < 1e-30) return 0;

    double zeropoint = n.x * a.x + n.y * a.y + n.z * a.z;
    double multi = (zeropoint - (n.x * origin.x + n.y * origin.y +
                                 n.z * origin.z)) / ndir;

    if (multi <= 0) return 0;  /* Behind ray */

    Vec3 endc = {
        origin.x + dir.x * multi,
        origin.y + dir.y * multi,
        origin.z + dir.z * multi
    };

    /* Check if point is inside triangle using sub-area test */
    Vec3 ea = vec3_sub(a, endc);
    Vec3 eb = vec3_sub(b, endc);
    Vec3 ec = vec3_sub(c, endc);

    double area2 = vec3_len(vec3_cross(ea, eb)) +
                   vec3_len(vec3_cross(ea, ec)) +
                   vec3_len(vec3_cross(ec, eb));

    return (area2 - area <= 0.001);
}

void flux_compute(const Mesh *m, double za, int tt, int Nmonte,
                  const int *neighbor_count, const int **neighbors,
                  double *F, int nthreads) {
    #pragma omp parallel for collapse(2) num_threads(nthreads) schedule(dynamic)
    for (int t = 0; t < tt; t++) {
        for (int nf = 0; nf < m->nface; nf++) {
            unsigned int seed = (unsigned int)(t * m->nface + nf + 42);

            double fi = ((t + 0.5) / tt) * 2.0 * M_PI + M_PI;
            Vec3 solar = {cos(-za) * sin(fi), cos(-za) * cos(fi), sin(-za)};

            /* Base flux = dot(solar, face_normal) */
            Vec3 fn = {m->face_normal[nf][0], m->face_normal[nf][1],
                       m->face_normal[nf][2]};
            double c6 = vec3_dot(solar, fn);

            if (c6 < 0) {
                F[nf * tt + t] = 0.0;
                continue;
            }

            int nc = neighbor_count[t * m->nface + nf];
            if (nc == 0) {
                F[nf * tt + t] = c6;
                continue;
            }

            /* MC shadow testing */
            Vec3 aa = {m->spot[m->facet[nf][0]][0], m->spot[m->facet[nf][0]][1],
                       m->spot[m->facet[nf][0]][2]};
            Vec3 bb = {m->spot[m->facet[nf][1]][0], m->spot[m->facet[nf][1]][1],
                       m->spot[m->facet[nf][1]][2]};
            Vec3 cc = {m->spot[m->facet[nf][2]][0], m->spot[m->facet[nf][2]][1],
                       m->spot[m->facet[nf][2]][2]};

            int shadow_count = 0;
            for (int mc = 0; mc < Nmonte; mc++) {
                Vec3 pot = random_point_on_triangle(aa, bb, cc, &seed);
                int shadowed = 0;

                for (int j = 0; j < nc; j++) {
                    int nnf = neighbors[t * m->nface + nf][j];

                    Vec3 da = {m->spot[m->facet[nnf][0]][0],
                               m->spot[m->facet[nnf][0]][1],
                               m->spot[m->facet[nnf][0]][2]};
                    Vec3 db = {m->spot[m->facet[nnf][1]][0],
                               m->spot[m->facet[nnf][1]][1],
                               m->spot[m->facet[nnf][1]][2]};
                    Vec3 dc = {m->spot[m->facet[nnf][2]][0],
                               m->spot[m->facet[nnf][2]][1],
                               m->spot[m->facet[nnf][2]][2]};

                    if (ray_triangle_intersect(pot, solar, da, db, dc)) {
                        shadowed = 1;
                        break;
                    }
                }

                if (shadowed) shadow_count++;
            }

            F[nf * tt + t] = (1.0 - (double)shadow_count / Nmonte) * c6;
        }
    }
}

void run_raytracer(const char *mesh_vertex_file, const char *mesh_face_file,
                   int npoint, int nface, double za_deg, int tt, int Nmonte,
                   const char *output_file, int nthreads) {
    printf("Loading mesh: %d vertices, %d faces\n", npoint, nface);
    Mesh *m = mesh_load(mesh_vertex_file, mesh_face_file, npoint, nface);
    mesh_compute_normals(m);

    double za = za_deg / 180.0 * M_PI;

    printf("Finding shadow neighbors...\n");
    int max_neighbors = 10;
    int *neighbor_count;
    int **neighbor_list;
    shadow_find_neighbors(m, za, tt, &neighbor_count, &neighbor_list,
                          max_neighbors);

    printf("Computing flux with %d MC samples...\n", Nmonte);
    double *F = calloc((size_t)nface * tt, sizeof(double));
    flux_compute(m, za, tt, Nmonte, neighbor_count,
                 (const int **)neighbor_list, F, nthreads);

    /* Write output */
    write_flux_dat(output_file, F, nface, tt);
    printf("Flux written to %s\n", output_file);

    /* Cleanup */
    for (int i = 0; i < tt * nface; i++) {
        free(neighbor_list[i]);
    }
    free(neighbor_list);
    free(neighbor_count);
    free(F);
    mesh_free(m);
}
