#ifndef FILEIO_H
#define FILEIO_H

#include "types.h"

/* Read density file: header has ftd, ffd [, hd], then hd*ftd*ffd values.
 * If has_hd=1, read 3 header lines; else read 2 (hd from caller). */
double *read_density_dat(const char *filename, int *ftd, int *ffd, int *hd,
                         int has_hd);

/* Write density file with header ftd, ffd, hd */
int write_density_dat(const char *filename, const double *den,
                      int ftd, int ffd, int hd);

/* Read surface temperature file (ftd, ffd header, then ftd*ffd values) */
double *read_temperature_dat(const char *filename, int *ftd, int *ffd);

/* Read temperature with time dimension (ffd copies shifted) */
double *read_temperature_time_dat(const char *filename, int *ftd, int *ffd);

/* Write ttime data */
int write_ttime_dat(const char *filename, const double *ttime,
                    int nrows, int totaln);

/* Write v_record data */
int write_v_record_dat(const char *filename, const double *v_record,
                       const int *td_list, const int *fd_list,
                       const int *ht_list, int nrecords);

/* Read v_record data */
int read_v_record_dat(const char *filename, double **v_record,
                      int **td_list, int **fd_list, int **ht_list,
                      int *nrecords);

/* Read flux file: nface, tt, then nface*tt values */
double *read_flux_dat(const char *filename, int *nface, int *tt);

/* Write flux file */
int write_flux_dat(const char *filename, const double *F,
                   int nface, int tt);

/* Write horizon file */
int write_horizon_dat(const char *filename, const double *horizon,
                      int ftd, int ffd);

#endif /* FILEIO_H */
