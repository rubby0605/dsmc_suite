#ifndef RAYTRACER_H
#define RAYTRACER_H

#include "types.h"

/* Load asteroid mesh from file (vertices + faces). */
Mesh *mesh_load(const char *vertex_file, const char *face_file,
                int npoint, int nface);

void mesh_free(Mesh *m);

/* Compute face normals, areas, centroids, directions */
void mesh_compute_normals(Mesh *m);

/* Random point on triangle (a, b, c) — standard barycentric method */
Vec3 random_point_on_triangle(Vec3 a, Vec3 b, Vec3 c, unsigned int *seed);

/* Run full shadow maker pipeline (uses Embree BVH internally) */
void run_raytracer(const char *mesh_vertex_file, const char *mesh_face_file,
                   int npoint, int nface, double za_deg, int tt, int Nmonte,
                   const char *output_file, int nthreads);

#endif /* RAYTRACER_H */
