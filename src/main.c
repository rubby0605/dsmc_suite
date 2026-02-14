#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "types.h"
#include "constants.h"
#include "grid.h"
#include "velocity.h"
#include "integrator.h"
#include "collision.h"
#include "polyfit.h"
#include "thermal.h"
#include "raytracer.h"
#include "fileio.h"

/* External function declarations */
extern void run_ballistic(int totaln, const char *temp_file,
                          const char *output_prefix, int nthreads);
extern void run_advanced(int totaln, const char *temp_file,
                         const char *prev_den_file,
                         const char *prev_v_record_file,
                         const char *output_prefix, int nthreads);

static void print_usage(const char *prog) {
    printf("DSMC Suite - Planetary Atmosphere Simulation\n");
    printf("=============================================\n\n");
    printf("Usage: %s <command> [options]\n\n", prog);
    printf("Commands:\n");
    printf("  raytrace  - MC ray tracing shadow/flux computation\n");
    printf("  thermal   - 1D spherical thermal diffusion model\n");
    printf("  dsmc      - Basic ballistic DSMC (no collisions)\n");
    printf("  dsmc-adv  - Advanced DSMC with collisions\n");
    printf("  pipeline  - Run full pipeline sequentially\n");
    printf("\n");
    printf("Options:\n");
    printf("  -n <particles>     Number of test particles (default: 10000)\n");
    printf("  -t <threads>       Number of OpenMP threads (default: auto)\n");
    printf("  -T <temp_file>     Temperature input file\n");
    printf("  -d <den_file>      Density input file (for dsmc-adv)\n");
    printf("  -v <v_record_file> Velocity record file (for dsmc-adv)\n");
    printf("  -o <output_prefix> Output file prefix (default: output)\n");
    printf("  -m <vertex_file>   Mesh vertex file (for raytrace)\n");
    printf("  -f <face_file>     Mesh face file (for raytrace)\n");
    printf("  --npoint <n>       Number of mesh vertices\n");
    printf("  --nface <n>        Number of mesh faces\n");
    printf("  --za <degrees>     Obliquity angle in degrees\n");
    printf("  --tt <timesteps>   Number of time steps\n");
    printf("  --nmonte <n>       MC samples per face (raytrace)\n");
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    const char *command = argv[1];

    /* Default options */
    int totaln = 10000;
    int nthreads = omp_get_max_threads();
    const char *temp_file = "Ceres_T_VP.dat";
    const char *den_file = NULL;
    const char *v_record_file = NULL;
    const char *output_prefix = "output";
    const char *vertex_file = NULL;
    const char *face_file = NULL;
    int npoint = 0, nface_mesh = 0;
    double za_deg = 12.3;
    int tt = 24;
    int nmonte = 50;

    /* Parse options */
    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-n") == 0 && i + 1 < argc)
            totaln = atoi(argv[++i]);
        else if (strcmp(argv[i], "-t") == 0 && i + 1 < argc)
            nthreads = atoi(argv[++i]);
        else if (strcmp(argv[i], "-T") == 0 && i + 1 < argc)
            temp_file = argv[++i];
        else if (strcmp(argv[i], "-d") == 0 && i + 1 < argc)
            den_file = argv[++i];
        else if (strcmp(argv[i], "-v") == 0 && i + 1 < argc)
            v_record_file = argv[++i];
        else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc)
            output_prefix = argv[++i];
        else if (strcmp(argv[i], "-m") == 0 && i + 1 < argc)
            vertex_file = argv[++i];
        else if (strcmp(argv[i], "-f") == 0 && i + 1 < argc)
            face_file = argv[++i];
        else if (strcmp(argv[i], "--npoint") == 0 && i + 1 < argc)
            npoint = atoi(argv[++i]);
        else if (strcmp(argv[i], "--nface") == 0 && i + 1 < argc)
            nface_mesh = atoi(argv[++i]);
        else if (strcmp(argv[i], "--za") == 0 && i + 1 < argc)
            za_deg = atof(argv[++i]);
        else if (strcmp(argv[i], "--tt") == 0 && i + 1 < argc)
            tt = atoi(argv[++i]);
        else if (strcmp(argv[i], "--nmonte") == 0 && i + 1 < argc)
            nmonte = atoi(argv[++i]);
    }

    printf("Using %d OpenMP threads\n", nthreads);

    /* Dispatch command */
    if (strcmp(command, "raytrace") == 0) {
        if (!vertex_file || !face_file || npoint == 0 || nface_mesh == 0) {
            fprintf(stderr, "raytrace requires: -m <vertex> -f <face> "
                    "--npoint <n> --nface <n>\n");
            return 1;
        }
        char flux_file[256];
        snprintf(flux_file, sizeof(flux_file), "%s_flux.dat", output_prefix);
        run_raytracer(vertex_file, face_file, npoint, nface_mesh,
                      za_deg, tt, nmonte, flux_file, nthreads);

    } else if (strcmp(command, "thermal") == 0) {
        /* Run thermal model */
        int ftd = 36, ffd = 72;
        SphericalGrid *grid = grid_create(ftd, ffd, 1);

        ThermalParams tp;
        thermal_params_default(&tp);
        tp.z_angle = za_deg / 180.0 * M_PI;

        thermal_solve(grid, &tp, NULL, nthreads);

        /* Write temperature output */
        char t_file[256];
        snprintf(t_file, sizeof(t_file), "%s_temperature.dat", output_prefix);
        FILE *fp = fopen(t_file, "w");
        if (fp) {
            fprintf(fp, "%d\n", ftd);
            fprintf(fp, "%d\n", ffd);
            for (int td = 0; td < ftd; td++) {
                for (int fd = 0; fd < ffd; fd++) {
                    fprintf(fp, "%8.6f\n", grid->temperature[td * ffd + fd]);
                }
            }
            fclose(fp);
            printf("Temperature written to %s\n", t_file);
        }
        grid_free(grid);

    } else if (strcmp(command, "dsmc") == 0) {
        run_ballistic(totaln, temp_file, output_prefix, nthreads);

    } else if (strcmp(command, "dsmc-adv") == 0) {
        run_advanced(totaln, temp_file, den_file, v_record_file,
                     output_prefix, nthreads);

    } else if (strcmp(command, "pipeline") == 0) {
        printf("=== FULL PIPELINE ===\n\n");

        /* Step 1: Ray trace (if mesh provided) */
        char flux_file[256];
        snprintf(flux_file, sizeof(flux_file), "%s_flux.dat", output_prefix);
        if (vertex_file && face_file && npoint > 0 && nface_mesh > 0) {
            printf("--- Step 1: Ray Tracing ---\n");
            run_raytracer(vertex_file, face_file, npoint, nface_mesh,
                          za_deg, tt, nmonte, flux_file, nthreads);
        } else {
            printf("--- Step 1: Skipped (no mesh provided) ---\n");
        }

        /* Step 2: Thermal model */
        printf("\n--- Step 2: Thermal Model ---\n");
        int ftd = 36, ffd = 72;
        SphericalGrid *grid = grid_create(ftd, ffd, 1);
        ThermalParams tp;
        thermal_params_default(&tp);
        thermal_solve(grid, &tp, NULL, nthreads);

        char t_file[256];
        snprintf(t_file, sizeof(t_file), "%s_temperature.dat", output_prefix);
        FILE *fp = fopen(t_file, "w");
        if (fp) {
            fprintf(fp, "%d\n", ftd);
            fprintf(fp, "%d\n", ffd);
            for (int td = 0; td < ftd; td++) {
                for (int fd = 0; fd < ffd; fd++) {
                    fprintf(fp, "%8.6f\n", grid->temperature[td * ffd + fd]);
                }
            }
            fclose(fp);
        }
        grid_free(grid);

        /* Step 3: Basic DSMC */
        printf("\n--- Step 3: Basic DSMC ---\n");
        char dsmc_prefix[256];
        snprintf(dsmc_prefix, sizeof(dsmc_prefix), "%s_phase1", output_prefix);
        run_ballistic(totaln, t_file, dsmc_prefix, nthreads);

        /* Step 4: Advanced DSMC */
        printf("\n--- Step 4: Advanced DSMC ---\n");
        char den_path[256], vr_path[256], adv_prefix[256];
        snprintf(den_path, sizeof(den_path), "%s_den.dat", dsmc_prefix);
        snprintf(vr_path, sizeof(vr_path), "%s_v_record.dat", dsmc_prefix);
        snprintf(adv_prefix, sizeof(adv_prefix), "%s_phase2", output_prefix);
        run_advanced(totaln, t_file, den_path, vr_path, adv_prefix, nthreads);

        printf("\n=== PIPELINE COMPLETE ===\n");

    } else {
        fprintf(stderr, "Unknown command: %s\n", command);
        print_usage(argv[0]);
        return 1;
    }

    return 0;
}
