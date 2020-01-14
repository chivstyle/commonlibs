/**
// copyright (c) chiv 2011-2012 all rights reseved
// \author chiv
// \date 2012.7.5
*/
#ifndef _C_MATRIX_API_H
#define _C_MATRIX_API_H

#define MATRIX_API

#ifdef __cplusplus
extern "C" {
#endif
MATRIX_API
int matrix_inverse(
	const double *dmatrix, 
	int row, 
	int col, 
	double *out);
MATRIX_API
int matrix_multiply(
	const double *matrixA, 
	int rowA, 
	int colA, 
	const double *matrixB, 
	int rowB, 
	int colB, 
	double *out);
MATRIX_API
int matrix_LUdecompose(
	const double *matrixA, 
	int n, 
	double *Lmatrix, 
	double *Umatrix);
/**
  \brief Ax=b, we solve this equation like:
  A(-1) * A * x = A(-1) * b
  then x = A(-1) * b, the key step of this routine is calculate
  the inverse of the matrix \par coefficients, luckly, we have
  a routine named "matrix_inverse" which could calculate the inverse
  of any real matrix if it has inverse indeed!
  \note crow == ccol, or you'll get an error '-E_BAD_MATRIX' defined in file:matrix_error.h
*/
MATRIX_API
int matrix_solve(
	const double *coefficients,
	int crow,
	int ccol,
	double *matrix_x,
	const double *matrix_b,
	int brow,
	int bcol);
/**
  \brief Ax=b, we solve this equation like:
  A=L*U, LUx=b, LyUx=by
  Ly=x, Ux=y
  firstly, we invoke "matrix_LUdecompose"(doolittle decomposition) to decompose the matrix \par
  coefficients to two triangular matrix L and U. L is lower triangular
  matrix.
  you must find the formula about solving triangular matrix coefficients linear system equation.
  I(chiv) found one! so i did it!
  \note crow == ccol, or you'll get an error '-E_BAD_MATRIX' defined in file:nurbs_error.h
*/
MATRIX_API
int matrix_LUsolve(
	const double *coefficients,
	int crow,
	int ccol,
	double *matrix_x,
	const double *matrix_b,
	int brow,
	int bcol);

#ifdef __cplusplus
}
#endif

#endif
