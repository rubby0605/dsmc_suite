#include "polyfit.h"
#include <math.h>
#include <string.h>

/*
 * Least-squares polynomial fit using normal equations.
 * For degree=3, this is a 4x4 system: (A^T A) c = A^T y
 * Solve with Gaussian elimination (small matrix, no pivoting issues).
 */
int polyfit(const double *x, const double *y, int n, int degree,
            double *coeffs) {
    int m = degree + 1;
    if (n < m) return -1;

    /* Build normal equations: (X^T X) c = X^T y */
    /* S[k] = sum(x_i^k) for k=0..2*degree */
    double S[14] = {0};  /* max 2*6+1 = 13, but we cap at degree=6 */
    (void)0;  /* B unused - rhs computed directly below */

    for (int i = 0; i < n; i++) {
        double xp = 1.0;
        for (int k = 0; k <= 2 * degree; k++) {
            S[k] += xp;
            xp *= x[i];
        }
        /* B[j] = sum(x_i^(degree-j) * y_i) -- MATLAB convention */
    }

    /* B: we want coeffs[0] = highest power, so map:
     * coeffs[0]*x^degree + coeffs[1]*x^(degree-1) + ... + coeffs[degree]
     * Internal: a[0] = constant, a[1] = x, ..., a[degree] = x^degree
     * Then remap at end. */
    double a[7] = {0};  /* internal coefficients, low-to-high */
    double rhs[7] = {0};

    for (int i = 0; i < n; i++) {
        double xp = 1.0;
        for (int j = 0; j < m; j++) {
            rhs[j] += xp * y[i];
            xp *= x[i];
        }
    }

    /* Build m x m matrix (normal equations) */
    double mat[7][8];  /* augmented matrix */
    memset(mat, 0, sizeof(mat));

    for (int r = 0; r < m; r++) {
        for (int c = 0; c < m; c++) {
            mat[r][c] = S[r + c];
        }
        mat[r][m] = rhs[r];
    }

    /* Gaussian elimination with partial pivoting */
    for (int col = 0; col < m; col++) {
        /* Find pivot */
        int max_row = col;
        double max_val = fabs(mat[col][col]);
        for (int row = col + 1; row < m; row++) {
            if (fabs(mat[row][col]) > max_val) {
                max_val = fabs(mat[row][col]);
                max_row = row;
            }
        }
        if (max_val < 1e-30) return -1;  /* singular */

        /* Swap rows */
        if (max_row != col) {
            for (int k = 0; k <= m; k++) {
                double tmp = mat[col][k];
                mat[col][k] = mat[max_row][k];
                mat[max_row][k] = tmp;
            }
        }

        /* Eliminate below */
        for (int row = col + 1; row < m; row++) {
            double factor = mat[row][col] / mat[col][col];
            for (int k = col; k <= m; k++) {
                mat[row][k] -= factor * mat[col][k];
            }
        }
    }

    /* Back substitution */
    for (int row = m - 1; row >= 0; row--) {
        a[row] = mat[row][m];
        for (int col = row + 1; col < m; col++) {
            a[row] -= mat[row][col] * a[col];
        }
        a[row] /= mat[row][row];
    }

    /* Reverse to MATLAB convention: coeffs[0] = highest power */
    for (int i = 0; i < m; i++) {
        coeffs[i] = a[degree - i];
    }

    return 0;
}

double polyval(const double *coeffs, int degree, double x) {
    /* Horner's method: coeffs[0]*x^degree + ... + coeffs[degree] */
    double result = coeffs[0];
    for (int i = 1; i <= degree; i++) {
        result = result * x + coeffs[i];
    }
    return result;
}
