#ifndef RAYTRACER_H
#define RAYTRACER_H

#include "types.h"

/* Load asteroid mesh from file (vertices + faces).
 * Format: first line = npoint nface, then vertices, then faces. */
Mesh *mesh_load(const char *vertex_file, const char *face_file,
                int npoint, int nface);

void mesh_free(Mesh *m);

/* Compute face normals, areas, centroids, directions */
void mesh_compute_normals(Mesh *m);

/* Find shadow neighbor candidates for each face at each time step.
 * neighbors[tt][nface][max_neighbors], count[tt][nface] */
void shadow_find_neighbors(const Mesh *m, double za, int tt,
                           int **count_out, int ***neighbors_out,
                           int max_neighbors);

/* Compute flux with MC shadow testing.
 * F[nface * tt], row-major.
 * Nmonte = number of MC samples per face per time. */
void flux_compute(const Mesh *m, double za, int tt, int Nmonte,
                  const int *neighbor_count, const int **neighbors,
                  double *F, int nthreads);

/* Random point on triangle (a, b, c) */
Vec3 random_point_on_triangle(Vec3 a, Vec3 b, Vec3 c, unsigned int *seed);

/* Ray-triangle intersection test.
 * Returns 1 if ray from origin in direction dir hits triangle (a,b,c). */
int ray_triangle_intersect(Vec3 origin, Vec3 dir, Vec3 a, Vec3 b, Vec3 c);

/* Run full shadow maker pipeline */
void run_raytracer(const char *mesh_vertex_file, const char *mesh_face_file,
                   int npoint, int nface, double za_deg, int tt, int Nmonte,
                   const char *output_file, int nthreads);

#endif /* RAYTRACER_H */
