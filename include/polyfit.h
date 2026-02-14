#ifndef POLYFIT_H
#define POLYFIT_H

/*
 * Least-squares polynomial fit (no external dependencies).
 * Fits polynomial of given degree to (x, y) data of n points.
 * coeffs[0] = highest degree coefficient (like MATLAB polyfit).
 * Returns 0 on success, -1 on failure.
 */
int polyfit(const double *x, const double *y, int n, int degree,
            double *coeffs);

/* Evaluate polynomial using Horner's method.
 * coeffs[0] = highest degree, coeffs[degree] = constant term. */
double polyval(const double *coeffs, int degree, double x);

#endif /* POLYFIT_H */
