//
// (C) 2012 chiv
//
#include "interpolation.h"
#include "matrix_error.h"
#include "matrix_api.h"
#include <math.h>
#include <malloc.h>
#include <memory.h>

#if defined(_DEBUG)
#include <stdio.h>
#endif

void global_do_knots_equally(int count, double *knots)
{
	int    i;
	double d = (double)count - 1.;
	
	for (i = 0; i < count; i++) {
		if (knots) {
			*knots++ = (double)i / d;
		}
	}
}
/**
 \par Q         is the value
 \par count     is the number of the Q
 \par degree    is the degree of the curve
 \return        is the count of the vector knots
*/
static
double chord2d(const POINT_2D *Q, int count, double d, int index)
{
	if (index == 0) {
		return 0.;
	}
	else if (index == count - 1) {
		return 1.;
	}
	else {
		return chord2d(Q, count, d, index - 1) + \
			sqrt(pow(Q[index].x - Q[index - 1].x, 2.) + pow(Q[index].y - Q[index - 1].y, 2.)) / d;
	}
}
static
double chord3d(const POINT_3D *Q, int count, double d, int index)
{
	if (index == 0) {
		return 0.;
	}
	else if (index == count - 1) {
		return 1.;
	}
	else {
		return chord3d(Q, count, d, index - 1) + \
			sqrt(pow(Q[index].x - Q[index - 1].x, 2.) + pow(Q[index].y - Q[index - 1].y, 2.) + \
			     pow(Q[index].z - Q[index - 1].z, 2.)) / d;
	}
}
/**
 \brief method - chord, calculate the vector uk
*/
void global_do_knots_chord2d(const POINT_2D *Q, int count, double *knots)
{
	/* calculate the chord */
	double d = 0.;
	int    i;
	for (i = 1; i < count; i++) {
		d += sqrt(pow(Q[i].x - Q[i - 1].x, 2.) + pow(Q[i].y - Q[i - 1].y, 2.));
	}
	for (i = 0; i < count; i++) {
		*knots++ = chord2d(Q, count, d, i);
	}
}

void global_do_knots_chord3d(const POINT_3D *Q, int count, double *knots)
{
	/* calculate the chord */
	double d = 0.;
	int    i;
	for (i = 1; i < count; i++) {
		d += sqrt(pow(Q[i].x - Q[i - 1].x, 2.) + pow(Q[i].y - Q[i - 1].y, 2.) + pow(Q[i].z - Q[i - 1].z, 2.));
	}
	for (i = 0; i < count; i++) {
		*knots++ = chord3d(Q, count, d, i);
	}
}

static
double centripetal2d(const POINT_2D *Q, int count, double d, int index)
{
	if (index == 0) {
		return 0.;
	}
	else if (index == count - 1) {
		return 1.;
	}
	else {
		return centripetal2d(Q, count, d, index - 1) + \
			sqrt(sqrt(pow(Q[index].x - Q[index - 1].x, 2.) + pow(Q[index].y - Q[index - 1].y, 2.))) / d;
	}
}
static
double centripetal3d(const POINT_3D *Q, int count, double d, int index)
{
	if (index == 0) {
		return 0.;
	}
	else if (index == count - 1) {
		return 1.;
	}
	else {
		return centripetal3d(Q, count, d, index - 1) + \
			sqrt(sqrt(pow(Q[index].x - Q[index - 1].x, 2.) + pow(Q[index].y - Q[index - 1].y, 2.) +
			          pow(Q[index].z - Q[index - 1].z, 2.))) / d;
	}
}
/**
 \brief method - centripetal, calculate the vector uk
*/
void global_do_knots_centripetal2d(const POINT_2D *Q, int count, double *knots)
{
	/* calculate the d */
	double d = 0.;
	int    i;
	for (i = 1; i < count; i++) {
		d += sqrt(fabs(Q[i].x - Q[i - 1].x));
	}
	for (i = 0; i < count; i++) {
		*knots++ = centripetal2d(Q, count, d, i);
	}
}
void global_do_knots_centripetal3d(const POINT_3D *Q, int count, double *knots)
{
	/* calculate the d */
	double d = 0.;
	int    i;
	for (i = 1; i < count; i++) {
		d += sqrt(fabs(Q[i].x - Q[i - 1].x));
	}
	for (i = 0; i < count; i++) {
		*knots++ = centripetal3d(Q, count, d, i);
	}
}
/**
 \brief calculate the knots
*/
void global_do_knots(double *knots, int count, const double *knots_vector, int vector_count, int degree)
{
	double sum;
	int    dimon = degree + 1;
	int    i, j;
	/** check error */
	if (vector_count + degree + 1 != count) return;
	
	for (i = 0; i < dimon; i++) {
		knots[i] = 0.;
	}
	for (i = count - dimon; i < count; i++) {
		knots[i] = 1.;
	}
	for (i = dimon; i < count - dimon; i++) {
		sum = 0.;
		for (j = i - degree; j < i; j++) {
			sum += knots_vector[j];
		}
		knots[i] = sum / degree;
	}
}

int do_index(double *values, int count, double cod)
{
	int head = 0;
	int tail = count - 1;
	int mid;
	
	while (head <= tail) {
		mid = head + (tail - head) / 2;
		if (values[mid] - cod > 0.)
			tail = mid - 1;
		else if (values[mid] - cod < 0.)
			head = mid + 1;
		else
			return mid;
	}
	return 0;
}
#if 1
static
void calc_coefficients(
    double *matrix,
    int row_count,
    int col_count,
    const double *knots,
    int count,
    const double *knots_vector,
    int vector_count,
    int degree)
{
    for (int i = 0; i < row_count; ++i) {
        for (int j = 0; j < col_count; ++j) {
            matrix[i * col_count + j] = nurbs(knots,
                                              count,
                                              degree,
                                              j,
                                              vector_count,
                                              knots_vector[i]);
        }
    }
}
#else
/**
 \brief calculat the coefficients
*/
static
void calc_coefficients(
	double *matrix, 
	int row_count, 
	int col_count, 
	const double *knots,
	int count, 
	const double *knots_vector, 
	int vector_count,
	int degree)
{
	int dimon;
	int left;
	int i, j;

	dimon = degree + 1;
	matrix[0] = 1.;
	memset(&matrix[1], 0, sizeof(double) * (col_count - 1));
	matrix[(row_count - 1) * col_count + col_count - 1] = 1.;
	memset(&matrix[(row_count - 1) * col_count], 0, sizeof(double) * (col_count - 1));
	
	for (i = 1; i < degree; i++) {
		for (j = 0; j < dimon; j++) {
			matrix[i * col_count + j] = nurbs(knots,
			                                  count,
											  degree,
											  j,
											  vector_count,
											  knots_vector[i]);
		}
		for (j = dimon; j < col_count; j++) {
			matrix[i * col_count + j] = 0.;
		}
	}
	left = 0;
	for (; i < row_count - 1; i++) {
		left++;
		for (j = 0; j < left; j++) {
			matrix[i * col_count + j] = 0.;
		}
		for (j = left; j < left + dimon; j++) {
			matrix[i * col_count + j] = nurbs(knots,
			                                  count,
											  degree,
											  j,
											  vector_count,
											  knots_vector[i]);
		}
		for (j = left + dimon; j < col_count; j++) {
			matrix[i * col_count + j] = 0.;
		}
	}
}
#endif
int nurbs_calc_control3d(
	const POINT_3D *Q,
	int count,
	int degree,
	const double *uk_vector,
	const double *knots,
	int knot_count,
	POINT_3D *ctrls)
{
	double *coefficients;
	int error = 0;
	/** first calculate the cols, like the nurbs_calc_control2d,
	in fact, you can invoke that function */
	coefficients = (double *)malloc(sizeof(double) * count * count);
	if (!coefficients) { return -E_NO_RESOURCE; }
	calc_coefficients(coefficients, count, count, knots, knot_count, uk_vector, count, degree);
	error = matrix_solve(coefficients, count, count, (double *)ctrls, (const double *)Q, count, 3);
	if (error != -E_OK) {
		free(coefficients);
		return -E_NO_SOLUTION;
	}
	/** free resources */
	free(coefficients);

	return -E_OK;
}

int nurbs_calc_control2d_surface(
	const POINT_2D *Q,
	int count_u,
	int count_v,
	int degree_u,
	int degree_v,
	const double *uk_vector_u,
	const double *uk_vector_v,
	const double *knots_u,
	int knots_count_u,
	const double *knots_v,
	int knots_count_v,
	POINT_2D *ctrls)
{
	double *coefficients_u, *coefficients_v;
	POINT_2D *temp1, *temp2;
	int     error;
	int     i, j;
	/** first calculate the cols, like the nurbs_calc_control2d,
	    in fact, you can invoke that function */
	coefficients_v = (double *)malloc(sizeof(double) * count_v * count_v);
	if (!coefficients_v) { return -E_NO_RESOURCE; }
	temp1 = (POINT_2D *)malloc(sizeof(POINT_2D) * count_v);
	if (!temp1) { free(coefficients_v);return -E_NO_RESOURCE; }
	temp2 = (POINT_2D *)malloc(sizeof(POINT_2D) * count_v);
	if (!temp2) { free(coefficients_v);free(temp1);return -E_NO_RESOURCE; }
	calc_coefficients(coefficients_v, count_v, count_v, knots_v, knots_count_v, uk_vector_v, count_v, degree_v);
	/** calculate the col matrix, it's not real */
	for (i = 0; i < count_u; i++) {
		/** fetch one col */
		for (j = 0; j < count_v; j++) {
			temp1[j] = Q[j * count_u + i];
		}
		error = matrix_LUsolve(coefficients_v, count_v, count_v, (double *)temp2, (const double *)temp1, count_v, 2);
		if (error != -E_OK) {
			free(coefficients_v);free(temp1);free(temp2);
			return -E_NO_SOLUTION;
		}
		for (j = 0; j < count_v; j++) {
			ctrls[j * count_u + i] = temp2[j];
		}
	}
	/** free resources */
	free(coefficients_v);free(temp1);free(temp2);
	/** now, you would get the col matrix, since it's just a temple matrix,
	    we will recalculate it */
	coefficients_u = (double *)malloc(sizeof(double) * count_u * count_u);
	if (!coefficients_u) { return -E_NO_RESOURCE; }
	temp1 = (POINT_2D *)malloc(sizeof(POINT_2D) * count_u);
	if (!temp1) { free(coefficients_u);return -E_NO_RESOURCE; }
	temp2 = (POINT_2D *)malloc(sizeof(POINT_2D) * count_u);
	if (!temp2) { free(coefficients_u);free(temp1);return -E_NO_RESOURCE; }
	calc_coefficients(coefficients_u, count_u, count_u, knots_u, knots_count_u, uk_vector_u, count_u, degree_u);
	/** think so? haha, got it! */
	for (i = 0; i < count_v; i++) {
		/** fetch one col */
		memcpy(temp1, &ctrls[i * count_u], sizeof(POINT_2D) * count_u);
		error = matrix_LUsolve(coefficients_u, count_u, count_u, (double *)temp2, (double *)temp1, count_u, 2);
		if (error != -E_OK) {
			free(coefficients_v);free(temp1);free(temp2);
			return -E_NO_SOLUTION;
		}
		memcpy(&ctrls[i * count_u], temp2, sizeof(POINT_2D) * count_u);
	}
	free(coefficients_u);free(temp1);free(temp2);
	
	return -E_OK;
}
