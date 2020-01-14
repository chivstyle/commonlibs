/**
// (C) 2012 CHIV, all rights reserved
//
//  an implementation of NURBS
//
// version : 1.0
*/
#include "nurbs.h"
#include <errno.h>
#include <float.h>

#if _VERBOSE
#define EBADDEGREE   1
#define EBADCOUNT    2
#endif

int MakeClampedUniformKnotsVector(double *knots, int knots_count, int cvcount, int degree)
{
	int i;
	int kncount = cvcount + degree + 1;
	double d = 1. / (cvcount - degree);
	if (kncount > knots_count) {
		return -1;
	}
	for (i = 0; i <= degree; ++i) {
		knots[i] = 0;
		knots[kncount-1-i] = 1;
	}
	for (i = degree+1; i < kncount-degree-1; ++i) {
		knots[i] = d * (i-degree);
	}

	return 0;
}

int MakeClampedChordKnotsVector(double *knots, int knots_count, int cvcount, int degree)
{
	int i, j;
	int kncount = cvcount + degree + 1;
	double d = 1. / (cvcount - 1);
	if (kncount > knots_count) {
		return -1;
	}
	for (i = 0; i <= degree; ++i) {
		knots[i] = 0;
		knots[kncount-1-i] = 1;
	}
	for (i = degree+1; i < kncount-degree-1; ++i) {
		double sum = 0;
		for (j = i-degree; j < i; ++j) {
			sum += j*d;
		}
		knots[i] = sum / 3;
	}

	return 0;
}

#if 1

double nurbs(const double* U, int count, int p, int i, int count_pt, double u)
{
    if (p == 0) {
        if (i == count_pt-1) {
            if (u >= U[i] && u <= U[i+1]) return 1.; else return 0.;
        } else {
            if (u >= U[i] && u < U[i+1]) return 1.; else return 0.;
        }
    } else {
        double f1 = (u - U[i]) * nurbs(U, count, p-1, i, count_pt, u);
        double d1 = U[i+p] - U[i];
        double f2 = (U[i+p+1] - u) * nurbs(U, count, p-1, i+1, count_pt, u);
        double d2 = U[i+p+1] - U[i+1];
        if (d1 == 0.) d1 = FLT_MIN;
        if (d2 == 0.) d2 = FLT_MIN;
        return f1 / d1 + f2 / d2;
    }
}

#else
/* NURBS basic-function */
double nurbs(
	const double * knots,     /* knots vector */
	int      count,     /* number of knots */
	int      degree,    /* degree */
	int      index,     /* index */
	int      count_pt,
	double   slice
)
{
	double numerator1, numerator2;
	double denominator1, denominator2;

	if (degree == 0) {
		if (index == count_pt - 1) {
			if (knots[index] <= slice && slice <= knots[index + 1]) {
				return 1.0;
			}
			else {
				return 0.0;
			}
		}
		else {
			if (knots[index] <= slice && slice < knots[index + 1]) {
				return 1.0;
			}
			else {
				return 0.0;
			}
		}
	}
	else {
		numerator1 = (slice - knots[index]) * nurbs(knots, count, degree - 1, index, count_pt, slice);
		numerator2 = (knots[index + degree + 1] - slice) * nurbs(knots, count, degree - 1, index + 1, count_pt, slice);
		denominator1 = knots[index + degree] - knots[index];
		denominator2 = knots[index + degree + 1] - knots[index + 1];

		if (denominator1 == 0.0)
			denominator1 = FLT_MIN;
		if (denominator2 == 0.0)
			denominator2 = FLT_MIN;

		return numerator1 / denominator1 + numerator2 / denominator2;
	}
}
#endif

POINT_3D NURBS_3D(
	double * knots,
	int count,
	int degree,
	int index,
	double slice,
	POINT_3D * point,
	int count_pt
)
{
	POINT_3D value_point = {0.0, 0.0, 0.0};
	double   basic;
	int      i;
#if _VERBOSE
	if (count_pt + degree + 1 != count) {
		errno = -EBADCOUNT;
		return value_point;
	}
#endif

	for (i = 0; i < count_pt; i++) {
		basic = nurbs(knots, count, degree, i, count_pt, slice);
		value_point.x += (point + i)->x * basic;
		value_point.y += (point + i)->y * basic;
		value_point.z += (point + i)->z * basic;
	}

	return value_point;
}

POINT_4D NURBS_4D(
	double * knots,
	int count,
	int degree,
	int index,
	double slice,
	POINT_4D * point,
	int count_pt
)
{
	POINT_4D value_point = {0.0, 0.0, 0.0, 1.0};
	double   basic;
	double   sum = 0.0;
	int      i;
#if _VERBOSE
	if (degree + count_pt + 1 != count) {
		errno = -EBADCOUNT;
		return value_point;
	}
#endif

	for (i = 0; i < count_pt; i++) {
		basic = nurbs(knots, count, degree, i, count_pt, slice);
		value_point.x += (point + i)->x * (point + i)->w * basic;
		value_point.y += (point + i)->y * (point + i)->w * basic;
		value_point.z += (point + i)->z * (point + i)->w * basic;
		sum += (point + i)->w * basic;
	}
	value_point.x /= sum;
	value_point.y /= sum;
	value_point.z /= sum;

	return value_point;
}

POINT_2D NURBS_2D(
	double * knots,
	int count,
	int degree,
	int index,
	double slice,
	POINT_2D * point,
	int count_pt
)
{
	POINT_2D value_point = {0.0, 0.0};
	double   basic, sum = 0;
	int      i;
#if _VERBOSE
	if (degree + count_pt + 1 != count) {
		errno = -EBADCOUNT;
		return value_point;
	}
#endif
	for (i = 0; i < count_pt; i++) {
		basic = nurbs(knots, count, degree, i, count_pt, slice);
		value_point.x += (point + i)->x * basic;
		value_point.y += (point + i)->y * basic;
	}

	return value_point;
}
