#ifndef TYPES_H
#define TYPES_H

#include <math.h>

/* 3D vector */
typedef struct {
    double x, y, z;
} Vec3;

static inline double vec3_dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline Vec3 vec3_cross(Vec3 a, Vec3 b) {
    return (Vec3){
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

static inline double vec3_len(Vec3 v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

static inline Vec3 vec3_scale(Vec3 v, double s) {
    return (Vec3){v.x * s, v.y * s, v.z * s};
}

static inline Vec3 vec3_add(Vec3 a, Vec3 b) {
    return (Vec3){a.x + b.x, a.y + b.y, a.z + b.z};
}

static inline Vec3 vec3_sub(Vec3 a, Vec3 b) {
    return (Vec3){a.x - b.x, a.y - b.y, a.z - b.z};
}

static inline Vec3 vec3_normalize(Vec3 v) {
    double l = vec3_len(v);
    if (l < 1e-30) return (Vec3){0, 0, 0};
    return vec3_scale(v, 1.0 / l);
}

/* Particle state */
typedef struct {
    Vec3 pos;
    Vec3 vel;
    double weight;
    int td, fd, ht;          /* grid indices (1-based like MATLAB) */
    double flight_time;
    int collision_count;
} Particle;

/* Spherical grid */
typedef struct {
    int ftd;                  /* number of theta bins */
    int ffd;                  /* number of phi bins */
    int hd;                   /* number of height bins */
    double *density;          /* [hd * ftd * ffd] */
    double *temperature;      /* [ftd * ffd] surface temperature */
    double *volume;           /* [hd * ftd * ffd] */
    double *horizon;          /* [ftd * ffd] frozen particles */
} SphericalGrid;

/* Physics parameters */
typedef struct {
    double G, M, k_B, mmass, mole;
    double Radii, ve, lifetime, P;
    double cross_section, alpha;
    double dt_step;           /* Verlet time step */
    double dt_rot;            /* rotation time step = P/ffd */
    double max_height;
    double unih;              /* uniform/base height step */
    double alpha2;            /* exponential grid growth factor */
    int use_exp_grid;         /* 0=uniform, 1=exponential */
} PhysicsParams;

/* Velocity lookup table */
#define VEL_TABLE_NTEMP 23
#define VEL_TABLE_NVEL  1000
typedef struct {
    double temps[VEL_TABLE_NTEMP];   /* temperature points */
    double speeds[VEL_TABLE_NTEMP][VEL_TABLE_NVEL]; /* sampled speeds */
    int count[VEL_TABLE_NTEMP];      /* valid entries per temp */
} VelocityTable;

/* Polynomial fit data for collision velocity sampling */
#define POLY_DEGREE 3
typedef struct {
    double coeffs[POLY_DEGREE + 1];  /* polynomial coefficients */
    double min_v, max_v;
    int valid;                        /* 1 if enough data */
} PolyFit;

/* Mesh for ray tracer */
typedef struct {
    int npoint;
    int nface;
    double (*spot)[3];        /* vertex positions [npoint][3] */
    int (*facet)[3];          /* face vertex indices [nface][3] */
    double (*face_normal)[3]; /* face normals [nface][3] */
    double *face_area;        /* face areas [nface] */
    double (*loc_xyz)[3];     /* face centroids [nface][3] */
    double (*dir)[2];         /* face direction (theta, phi) [nface][2] */
} Mesh;

#endif /* TYPES_H */
