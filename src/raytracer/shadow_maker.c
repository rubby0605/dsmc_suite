#include "raytracer.h"
#include "constants.h"
#include "fileio.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#include <embree4/rtcore.h>

/*
 * MC shadow/flux computation using Intel Embree for BVH-accelerated
 * ray-triangle intersection.
 *
 * Replaces the old O(n^2) brute-force neighbor search with Embree's
 * BVH which gives O(n log n) build + O(log n) per ray query.
 */

/* ── Mesh I/O (unchanged) ── */

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

    FILE *fp = fopen(vertex_file, "r");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", vertex_file); return m; }
    for (int i = 0; i < npoint; i++) {
        fscanf(fp, "%lf %lf %lf", &m->spot[i][0], &m->spot[i][1],
               &m->spot[i][2]);
    }
    fclose(fp);

    fp = fopen(face_file, "r");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", face_file); return m; }
    for (int i = 0; i < nface; i++) {
        fscanf(fp, "%d %d %d", &m->facet[i][0], &m->facet[i][1],
               &m->facet[i][2]);
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

/* ── Normal computation (unchanged) ── */

void mesh_compute_normals(Mesh *m) {
    for (int nf = 0; nf < m->nface; nf++) {
        Vec3 aa = {m->spot[m->facet[nf][0]][0], m->spot[m->facet[nf][0]][1],
                   m->spot[m->facet[nf][0]][2]};
        Vec3 bb = {m->spot[m->facet[nf][1]][0], m->spot[m->facet[nf][1]][1],
                   m->spot[m->facet[nf][1]][2]};
        Vec3 cc = {m->spot[m->facet[nf][2]][0], m->spot[m->facet[nf][2]][1],
                   m->spot[m->facet[nf][2]][2]};

        m->loc_xyz[nf][0] = (aa.x + bb.x + cc.x) / 3.0;
        m->loc_xyz[nf][1] = (aa.y + bb.y + cc.y) / 3.0;
        m->loc_xyz[nf][2] = (aa.z + bb.z + cc.z) / 3.0;

        Vec3 c1 = vec3_sub(bb, aa);
        Vec3 c2 = vec3_sub(cc, aa);

        m->face_area[nf] = vec3_len(vec3_cross(c1, c2)) / 2.0;

        Vec3 c2n = vec3_normalize(c2);
        Vec3 c1n = vec3_normalize(c1);
        Vec3 c3 = vec3_cross(c2n, c1n);
        Vec3 normal = vec3_normalize(c3);

        /* Ensure normal points outward (away from body center).
         * For a closed mesh roughly centered at origin, the face centroid
         * points outward — flip normal if it disagrees. */
        Vec3 centroid = {m->loc_xyz[nf][0], m->loc_xyz[nf][1], m->loc_xyz[nf][2]};
        if (vec3_dot(normal, centroid) < 0) {
            normal = vec3_scale(normal, -1.0);
        }

        m->face_normal[nf][0] = normal.x;
        m->face_normal[nf][1] = normal.y;
        m->face_normal[nf][2] = normal.z;

        double thi = asin(normal.z);
        double fi = atan2(normal.x, normal.y);
        if (fi < 0) fi += 2.0 * M_PI;
        m->dir[nf][0] = thi;
        m->dir[nf][1] = fi;
    }
}

/* ── Random point on triangle (simplified to standard barycentric) ── */

Vec3 random_point_on_triangle(Vec3 a, Vec3 b, Vec3 c, unsigned int *seed) {
    double r1 = sqrt((double)rand_r(seed) / RAND_MAX);
    double r2 = (double)rand_r(seed) / RAND_MAX;
    /* P = (1-r1)*A + r1*(1-r2)*B + r1*r2*C */
    return vec3_add(vec3_add(
        vec3_scale(a, 1.0 - r1),
        vec3_scale(b, r1 * (1.0 - r2))),
        vec3_scale(c, r1 * r2));
}

/* ── Embree BVH setup ── */

static RTCDevice embree_device = NULL;

static RTCScene embree_build_scene(const Mesh *m) {
    if (!embree_device) {
        embree_device = rtcNewDevice(NULL);
    }

    RTCScene scene = rtcNewScene(embree_device);
    /* Allow concurrent ray queries from multiple threads */
    rtcSetSceneFlags(scene, RTC_SCENE_FLAG_ROBUST);
    rtcSetSceneBuildQuality(scene, RTC_BUILD_QUALITY_HIGH);

    RTCGeometry geom = rtcNewGeometry(embree_device, RTC_GEOMETRY_TYPE_TRIANGLE);

    /* Set vertex buffer (Embree wants float, not double) */
    float *verts = (float *)rtcSetNewGeometryBuffer(
        geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3,
        3 * sizeof(float), m->npoint);
    for (int i = 0; i < m->npoint; i++) {
        verts[3 * i + 0] = (float)m->spot[i][0];
        verts[3 * i + 1] = (float)m->spot[i][1];
        verts[3 * i + 2] = (float)m->spot[i][2];
    }

    /* Set index buffer */
    unsigned *indices = (unsigned *)rtcSetNewGeometryBuffer(
        geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,
        3 * sizeof(unsigned), m->nface);
    for (int i = 0; i < m->nface; i++) {
        indices[3 * i + 0] = (unsigned)m->facet[i][0];
        indices[3 * i + 1] = (unsigned)m->facet[i][1];
        indices[3 * i + 2] = (unsigned)m->facet[i][2];
    }

    rtcCommitGeometry(geom);
    rtcAttachGeometry(scene, geom);
    rtcReleaseGeometry(geom);
    rtcCommitScene(scene);

    return scene;
}

/* Check if a ray from origin toward solar direction hits any face.
 * Returns 1 if occluded (shadowed), 0 if clear.
 * Offset origin along face normal to avoid self-intersection on closed meshes. */
static int embree_occluded(RTCScene scene, Vec3 origin, Vec3 face_normal, Vec3 dir) {
    struct RTCRay ray;
    /* Offset along face normal (outward) to start clearly outside the surface */
    float eps = 1e-3f;
    ray.org_x = (float)(origin.x + eps * face_normal.x);
    ray.org_y = (float)(origin.y + eps * face_normal.y);
    ray.org_z = (float)(origin.z + eps * face_normal.z);
    ray.dir_x = (float)dir.x;
    ray.dir_y = (float)dir.y;
    ray.dir_z = (float)dir.z;
    ray.tnear = 0.0f;
    ray.tfar = 1e30f;
    ray.mask = (unsigned)-1;
    ray.flags = 0;

    struct RTCOccludedArguments args;
    rtcInitOccludedArguments(&args);

    rtcOccluded1(scene, &ray, &args);

    /* If tfar becomes -inf, the ray was occluded */
    return (ray.tfar < 0.0f) ? 1 : 0;
}

/* ── Flux computation with Embree ── */

static void flux_compute_embree(const Mesh *m, RTCScene scene,
                                double za, int tt, int Nmonte,
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

            /* MC shadow testing via Embree */
            Vec3 aa = {m->spot[m->facet[nf][0]][0], m->spot[m->facet[nf][0]][1],
                       m->spot[m->facet[nf][0]][2]};
            Vec3 bb = {m->spot[m->facet[nf][1]][0], m->spot[m->facet[nf][1]][1],
                       m->spot[m->facet[nf][1]][2]};
            Vec3 cc = {m->spot[m->facet[nf][2]][0], m->spot[m->facet[nf][2]][1],
                       m->spot[m->facet[nf][2]][2]};

            int shadow_count = 0;
            for (int mc = 0; mc < Nmonte; mc++) {
                Vec3 pot = random_point_on_triangle(aa, bb, cc, &seed);
                if (embree_occluded(scene, pot, fn, solar)) {
                    shadow_count++;
                }
            }

            F[nf * tt + t] = (1.0 - (double)shadow_count / Nmonte) * c6;
        }
    }
}

/* ── Main entry point ── */

void run_raytracer(const char *mesh_vertex_file, const char *mesh_face_file,
                   int npoint, int nface, double za_deg, int tt, int Nmonte,
                   const char *output_file, int nthreads) {
    printf("Loading mesh: %d vertices, %d faces\n", npoint, nface);
    Mesh *m = mesh_load(mesh_vertex_file, mesh_face_file, npoint, nface);
    mesh_compute_normals(m);

    double za = za_deg / 180.0 * M_PI;

    printf("Building Embree BVH...\n");
    RTCScene scene = embree_build_scene(m);
    printf("BVH built.\n");

    printf("Computing flux with %d MC samples (Embree accelerated)...\n", Nmonte);
    double *F = calloc((size_t)nface * tt, sizeof(double));
    flux_compute_embree(m, scene, za, tt, Nmonte, F, nthreads);

    /* Write output */
    write_flux_dat(output_file, F, nface, tt);
    printf("Flux written to %s\n", output_file);

    /* Cleanup */
    rtcReleaseScene(scene);
    if (embree_device) {
        rtcReleaseDevice(embree_device);
        embree_device = NULL;
    }
    free(F);
    mesh_free(m);
}
