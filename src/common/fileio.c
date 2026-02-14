#include "fileio.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double *read_density_dat(const char *filename, int *ftd, int *ffd, int *hd,
                         int has_hd) {
    FILE *fp = fopen(filename, "r");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", filename); return NULL; }

    fscanf(fp, "%d", ftd);
    fscanf(fp, "%d", ffd);
    if (has_hd) {
        fscanf(fp, "%d", hd);
    }

    int total = (*hd) * (*ftd) * (*ffd);
    double *data = calloc(total, sizeof(double));

    for (int i = 0; i < total; i++) {
        fscanf(fp, "%lf", &data[i]);
    }
    fclose(fp);
    return data;
}

int write_density_dat(const char *filename, const double *den,
                      int ftd, int ffd, int hd) {
    FILE *fp = fopen(filename, "w");
    if (!fp) return -1;

    fprintf(fp, "%d\n", ftd);
    fprintf(fp, "%d\n", ffd);
    fprintf(fp, "%d\n", hd);

    for (int ht = 0; ht < hd; ht++) {
        for (int td = 0; td < ftd; td++) {
            for (int fd = 0; fd < ffd; fd++) {
                fprintf(fp, "%9.6f\n", den[ht * ftd * ffd + td * ffd + fd]);
            }
        }
    }
    fclose(fp);
    return 0;
}

double *read_temperature_dat(const char *filename, int *ftd, int *ffd) {
    FILE *fp = fopen(filename, "r");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", filename); return NULL; }

    fscanf(fp, "%d", ftd);
    fscanf(fp, "%d", ffd);

    int total = (*ftd) * (*ffd);
    double *data = calloc(total, sizeof(double));

    for (int i = 0; i < total; i++) {
        fscanf(fp, "%lf", &data[i]);
    }
    fclose(fp);
    return data;
}

double *read_temperature_time_dat(const char *filename, int *ftd, int *ffd) {
    /*
     * Read surface temperature and create time-shifted copies.
     * Input: ftd*ffd values (one snapshot).
     * Output: ffd * ftd * ffd array (ffd time steps, each shifted by 1 column).
     *
     * MATLAB logic from test_ballistic.m:
     *   copyT2 = read ftd*ffd
     *   for ii = 1:ffd
     *       copyT(ii,:,:) = copyT2
     *       copyT2 shift columns by 1
     */
    int f_td, f_fd;
    double *T_surface = read_temperature_dat(filename, &f_td, &f_fd);
    if (!T_surface) return NULL;

    *ftd = f_td;
    *ffd = f_fd;

    /* Allocate ffd * ftd * ffd */
    double *copyT = calloc((size_t)f_fd * f_td * f_fd, sizeof(double));
    /* Temp working copy */
    double *work = calloc((size_t)f_td * f_fd, sizeof(double));
    memcpy(work, T_surface, sizeof(double) * f_td * f_fd);

    for (int t = 0; t < f_fd; t++) {
        /* Copy current state into copyT[t, :, :] */
        for (int td = 0; td < f_td; td++) {
            for (int fd = 0; fd < f_fd; fd++) {
                copyT[t * f_td * f_fd + td * f_fd + fd] = work[td * f_fd + fd];
            }
        }
        /* Shift columns: move col 0 to end, shift rest left */
        for (int td = 0; td < f_td; td++) {
            double tmp = work[td * f_fd];
            for (int fd = 0; fd < f_fd - 1; fd++) {
                work[td * f_fd + fd] = work[td * f_fd + fd + 1];
            }
            work[td * f_fd + f_fd - 1] = tmp;
        }
    }

    free(work);
    free(T_surface);
    return copyT;
}

int write_ttime_dat(const char *filename, const double *ttime,
                    int nrows, int totaln) {
    FILE *fp = fopen(filename, "w");
    if (!fp) return -1;

    for (int r = 0; r < nrows; r++) {
        for (int n = 0; n < totaln; n++) {
            fprintf(fp, "%9.5f\n", ttime[r * totaln + n]);
        }
    }
    fclose(fp);
    return 0;
}

int write_v_record_dat(const char *filename, const double *v_record,
                       const int *td_list, const int *fd_list,
                       const int *ht_list, int nrecords) {
    FILE *fp = fopen(filename, "w");
    if (!fp) return -1;

    fprintf(fp, "%d\n", nrecords);
    for (int i = 0; i < nrecords; i++) {
        fprintf(fp, "%f %f %f\n",
                v_record[i * 3], v_record[i * 3 + 1], v_record[i * 3 + 2]);
        fprintf(fp, "%d\n", td_list[i]);
        fprintf(fp, "%d\n", fd_list[i]);
        fprintf(fp, "%d\n", ht_list[i]);
    }
    fclose(fp);
    return 0;
}

int read_v_record_dat(const char *filename, double **v_record,
                      int **td_list, int **fd_list, int **ht_list,
                      int *nrecords) {
    FILE *fp = fopen(filename, "r");
    if (!fp) return -1;

    fscanf(fp, "%d", nrecords);
    int n = *nrecords;

    *v_record = calloc(n * 3, sizeof(double));
    *td_list = calloc(n, sizeof(int));
    *fd_list = calloc(n, sizeof(int));
    *ht_list = calloc(n, sizeof(int));

    for (int i = 0; i < n; i++) {
        fscanf(fp, "%lf %lf %lf",
               &(*v_record)[i * 3], &(*v_record)[i * 3 + 1],
               &(*v_record)[i * 3 + 2]);
        fscanf(fp, "%d", &(*td_list)[i]);
        fscanf(fp, "%d", &(*fd_list)[i]);
        fscanf(fp, "%d", &(*ht_list)[i]);
    }
    fclose(fp);
    return 0;
}

double *read_flux_dat(const char *filename, int *nface, int *tt) {
    FILE *fp = fopen(filename, "r");
    if (!fp) return NULL;

    fscanf(fp, "%d", nface);
    fscanf(fp, "%d", tt);

    int total = (*nface) * (*tt);
    double *F = calloc(total, sizeof(double));

    for (int nf = 0; nf < *nface; nf++) {
        for (int t = 0; t < *tt; t++) {
            fscanf(fp, "%lf", &F[nf * (*tt) + t]);
        }
    }
    fclose(fp);
    return F;
}

int write_flux_dat(const char *filename, const double *F,
                   int nface, int tt) {
    FILE *fp = fopen(filename, "w");
    if (!fp) return -1;

    fprintf(fp, "%d\n", nface);
    fprintf(fp, "%d\n", tt);

    for (int nf = 0; nf < nface; nf++) {
        for (int t = 0; t < tt; t++) {
            fprintf(fp, "%9.7f\n", F[nf * tt + t]);
        }
    }
    fclose(fp);
    return 0;
}

int write_horizon_dat(const char *filename, const double *horizon,
                      int ftd, int ffd) {
    FILE *fp = fopen(filename, "w");
    if (!fp) return -1;

    fprintf(fp, "%d\n", ftd);
    fprintf(fp, "%d\n", ffd);
    for (int td = 0; td < ftd; td++) {
        for (int fd = 0; fd < ffd; fd++) {
            fprintf(fp, "%9.6f\n", horizon[td * ffd + fd]);
        }
    }
    fclose(fp);
    return 0;
}
