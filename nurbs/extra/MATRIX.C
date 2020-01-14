/**
// (C) 2012 CHIV, all rights reserved
//
// \brief resolve the linear system equations
// \author chiv
*/
#include "matrix_error.h"
#include "matrix_api.h"
#include <malloc.h>
#include <memory.h>
#include <math.h>

//#define _LUdecompose_CHIV

#if defined(_DEBUG)
#include <stdio.h>
#include <tchar.h>
#endif
/**
 \brief Gauss-Jordan, calculate the inverse of the Matrix
 \note \par 'row' must = \par 'col', if not, an error [-E_BAD_MATRIX] thrown
 \return -E_OK, success, others failed, defined in file:nurbs_error.h, eg.
         -E_NO_INV, your matrix is singular!
		 -E_NO_RESOURCE, there's no more memory to allocate, or can't open file, etc
 \par out, don't give 0, or an c excepation occurs
*/
int matrix_inverse(const double *dmatrix, int row, int col, double *out)
{
	/** temp vars */
	int    *is, *js, i, j, k, l, u, v, n;
	double  d, p;
	double *matrix;
	/** check the matrix, must be square */
	if (row != col) return -E_BAD_MATRIX;
	n = row;
	is = (int *)malloc(sizeof(int) * row);
	if (!is) return -E_NO_RESOURCE;
	js = (int *)malloc(sizeof(int) * row);
	if (!js) { free(is); return -E_NO_RESOURCE; }
	matrix = (double *)malloc(sizeof(double) * n * n);
	if (!matrix) { free(is);free(js);return -E_NO_RESOURCE; }
	memcpy(matrix, dmatrix, sizeof(double) * n * n);
	/** okay, resources ready */
	for (k = 0; k < n; k++) {
		d = 0.;
		for (i = k; i < n; i++) {
			for (j = k; j < n; j++) {
				l = i * n + j;
				p = fabs(matrix[l]);
				if (p > d) {
					d = p;
					is[k] = i;
					js[k] = j;
				}
			}
		}
		if (d == 0.) {
			free(is);
			free(js);
#if defined(_DEBUG)
			memcpy(out, matrix, sizeof(double) * n * n);
#endif
            free(matrix);
			return -E_NO_INV;
		}
		if (is[k] != k) {
			for (j = 0; j < n; j++) {
				u = k * n + j;
				v = is[k] * n + j;
				p = matrix[u];
				matrix[u] = matrix[v];
				matrix[v] = p;
			}
		}
		if (js[k] != k) {
			for (i = 0; i < n; i++) {
				u = i * n + k;
				v = i * n + js[k];
				p = matrix[u];
				matrix[u] = matrix[v];
				matrix[v] = p;
			}
		}
		l = k * n + k;
		matrix[l] = 1. / matrix[l];
		for (j = 0; j < n; j++) {
			if (j != k) {
				u = k * n + j;
				matrix[u] = matrix[u] * matrix[l];
			}
		}
		for (i = 0; i < n; i++) {
			if (i != k) {
				for (j = 0; j < n; j++) {
					if (j != k) {
						u = i * n + j;
						matrix[u] = matrix[u] - matrix[i * n + k] * matrix[k * n + j];
					}
				}
				u = i * n + k;
				matrix[u] = -matrix[u] * matrix[l];
			}
		}
	}
	for (k = n - 1; k >= 0; k--) {
		if (js[k] != k) {
			for (j = 0; j < n; j++) {
				u = k * n + j;
				v = js[k] * n + j;
				p = matrix[u];
				matrix[u] = matrix[v];
				matrix[v] = p;
			}
		}
		if (is[k] != k) {
			for (i = 0; i < n; i++) {
				u = i * n + k;
				v = i * n + is[k];
				p = matrix[u];
				matrix[u] = matrix[v];
				matrix[v] = p;
			}
		}
	}
	memcpy(out, matrix, sizeof(double) * n * n);
	free(is);
	free(js);
	free(matrix);
	
	return -E_OK;
}
/**
  \brief calculate the product of two matrix!
  \note colA == rowB, or error throwed
  \return -E_OK, success, others failed, defined in file:nurbs_error.h,eg.
          -E_NO_RESOURCE, etc
*/
int matrix_multiply(const double *matrixA, int rowA, int colA, const double *matrixB, int rowB, int colB, double *out)
{
	int     i, j, k;
	double  temp, tempc;
	double *prow;
	/** check the matrix */
	if (colA != rowB) {
		return -E_BAD_MATRIX;
	}
	prow = (double *)malloc(sizeof(double) * colA);
	if (!prow) { return -E_NO_RESOURCE; }
	/** you should know, rowB = colA */
	for (i = 0; i < rowA; i++) {
		memcpy(prow, &matrixA[i * colA], sizeof(double) * colA);
		for (j = 0; j < colB; j++) {
			temp = 0.;
			for (k = 0; k < rowB; k++) {
				tempc = matrixB[k * colB + j];
				temp += tempc * prow[k];
			}
			*out++ = temp;
		}
	}
	free(prow);
	
	return -E_OK;
}
#if !defined(_LUdecompose_CHIV)
int matrix_LUdecompose(const double *matrix_source, int n, double *matrixL, double *matrixU)
{
	int i, j, k, w, v, ll;

	double *matrix = (double *)malloc(sizeof(double) * n * n);
	if (!matrix) {
		return -E_NO_RESOURCE;
	}
	memcpy(matrix, matrix_source, sizeof(double) * n * n);
	
	for (k = 0; k <= n - 2; k++) {
		ll = k * n + k;
		if (fabs(matrix[ll]) + 1.0 == 1.0) {
			return -E_NO_SOLUTION;
		}
		for (i = k + 1; i < n; i++) {
			w = i * n + k;
			matrix[w] /= matrix[ll];
		}
		for (i = k + 1; i < n; i++) {
			w = i * n + k;
			for (j = k + 1; j < n; j++) {
				v = i * n + j;
				matrix[v] -= matrix[w] * matrix[k * n + j];
			}
		}
	}
	for (i = 0; i < n; i++) {
		for (j = 0; j < i; j++) {
			w = i * n + j;
			matrixL[w] = matrix[w];
			matrixU[w] = 0.;
		}
		w = i * n + i;
		matrixL[w] = 1.0;
		matrixU[w] = matrix[w];
		for (j = i + 1; j < n; j++) {
			w = i * n + j;
			matrixL[w] = 0.;
			matrixU[w] = matrix[w];
		}
	}
	free(matrix);

	return -E_OK;
}
#else
int matrix_LUdecompose(const double *matrix_source, int n, double *matrix_l, double *matrix_u)
{
	int i, k, r;
	double d;
	double *matrix = (double *)malloc(sizeof(double) * n * n);
	if (!matrix) {
		return -E_NO_RESOURCE;
	}
	memcpy(matrix, matrix_source, sizeof(double) * n * n);
	memset(matrix_u, 0, sizeof(double) * n * n);
	memset(matrix_l, 0, sizeof(double) * n * n);
	
	for (i = 0; i < n; i++) {
		matrix_u[i] = matrix[i];
	}
	if (matrix_u[0] == 0.) { free(matrix); return -E_NO_SOLUTION; }
	for (i = 1; i < n; i++) {
		matrix_l[i * n] = matrix[i * n] / matrix_u[0];
	}
	for (r = 0; r < n; r++) {
		for (i = r; i < n; i++) {
			d = 0.;
			for (k = 0; k < r; k++) {
				d += matrix_l[r * n + k] * matrix_u[k * n + i];
			}
			matrix_u[r * n + i] = matrix[r * n + i] - d;
			d = 0.;
			for (k = 0; k < r; k++) {
				d += matrix_l[i * n + k] * matrix_u[k * n + r];
			}
			if (matrix_u[r * n + r] == 0.) { free(matrix); return -E_NO_SOLUTION; }
			matrix_l[i * n + r] = (matrix[i * n + r] - d) / matrix_u[r * n + r];
		}
	}
	free(matrix);
	
	return -E_OK;
}
#endif

/**
  \brief solve the linear system equation :
         coefficients * matrix_x = matrix_b
		 the matrix 'coefficients' and 'matrix_b' are known, we calculate
		 the matrix x, the unknown variant
  \return error number, defined in file:nurbs_error.h
*/

int matrix_solve(
	const double *coefficients,
	int crow,
	int ccol,
	double *matrix_x,
	const double *matrix_b,
	int brow,
	int bcol)
{
	double *inverse;
	int     error;

	inverse = (double *)malloc(sizeof(double) * crow * ccol);
	if (!inverse) { return -E_NO_RESOURCE; }
	error = matrix_inverse(coefficients, crow, ccol, inverse);

	if (error != -E_OK) { free(inverse); return error; }
	error = matrix_multiply(inverse, crow, ccol, matrix_b, brow, bcol, matrix_x);
	free(inverse);

	return error;
}
/**
  \brief solve the linear system equation, but we don't calculate the inverse
         of the coefficients-matrix. firstly, we decompose the coefficient-matrix
		 (shorted c-matrix) to two triangular matrix, L and U.we know, solve system
		 equations like this is simple, we can assume:
		 Ly = b, Ux = y, exactly, it's the same as Lx = b
  \return -E_OK, success, others failed, defined in file:nurbs_error.h
*/
int matrix_LUsolve(
	const double *coefficients,
	int crow,
	int ccol,
	double *matrix_x,
	const double *matrix_b,
	int brow,
	int bcol)
{
	double *matrix_l, *matrix_u;
	double *y, d;
	int     error;
	int     k, i, j;
	
	/** calculate the LU factorization */
	if (crow != ccol) { return -E_BAD_MATRIX; }
	matrix_l = (double *)malloc(sizeof(double) * crow * ccol);
	if (!matrix_l) { return -E_NO_RESOURCE; }
	matrix_u = (double *)malloc(sizeof(double) * crow * ccol);
	if (!matrix_u) { free(matrix_l); return -E_NO_RESOURCE; }
	error = matrix_LUdecompose(coefficients, crow, matrix_l, matrix_u);
	if (error != -E_OK) { free(matrix_l); free(matrix_u); return error; }
	/**
	  \brief Ly=b
	         y0 = b0, yi = bi - sum(likyk).(k=1...i-1, i=1,2,...n)
	*/
	y = (double *)malloc(sizeof(double) * brow * bcol);
	if (!y) { free(matrix_l); free(matrix_u); return -E_NO_RESOURCE; }
	/** y0 = b0 */
	for (i = 0; i < bcol; i++) {
		y[i] = matrix_b[i];
	}
	/** next */
	for (j = 0; j < bcol; j++) {
		for (i = 1; i < brow; i++) {
			d = 0.;
			for (k = 0; k < i; k++) {
				d += matrix_l[i * crow + k] * y[k * bcol + j];
			}
			y[i * bcol + j] = matrix_b[i * bcol + j] - d;
		}
	}
	/**
	  \brief Ux=y
	         xn = yn / unn, xi = (yi - sum(uikxk)/uii,i=n-1,n-2,...1
	*/
	for (i = 0; i < bcol; i++) {
		matrix_x[(brow - 1) * bcol + i] = y[(brow - 1) * bcol + i];
	}
	for (j = 0; j < bcol; j++) {
		for (i = brow - 2; i >= 0; i--) {
			d = 0.;
			for (k = i + 1 ; k <= brow; k++) {
				d += matrix_u[i * crow + k] * matrix_x[k * bcol + j];
			}
			matrix_x[i * bcol + j] = (y[i * bcol + j] - d) / matrix_u[i * crow + i];
		}
	}

	free(matrix_l);
	free(matrix_u);
	free(y);
	
	return -E_OK;
}
