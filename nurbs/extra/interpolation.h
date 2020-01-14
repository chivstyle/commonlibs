/**
// (C) 2012 CHIV
//
// \brief interpolation of NURBS, algorithms from <<the nurbs book>>
*/
#ifndef _NURBS_INTERPOLATION_H
#define _NURBS_INTERPOLATION_H

#include "nurbs.h"

#ifdef __cplusplus
extern "C" {
#endif
/**
 \brief method
 \note don't/never modify this file if your are not the author!!!
*/
typedef enum {
	F_EQUAL = 0,               /** for 2d application, the z is ignored, so if you want this, use F_*_Z */
	F_CHORD,
	F_CENTRIPETAL,
} F_Method;

void global_do_knots_equally(
	int count, 
	double *uk);

void global_do_knots_chord2d(
	const POINT_2D *Q, 
	int count, 
	double *uk);
void global_do_knots_chord3d(
	const POINT_3D *Q,
	int count,
	double *uk);

void global_do_knots_centripetal2d(
	const POINT_2D *Q, 
	int count, 
	double *uk);
void global_do_knots_centripetal3d(
	const POINT_3D *Q,
	int count,
	double *uk);

void global_do_knots(
	double *knots, 
	int count, 
	const double *uk_vector, 
	int uk_count, 
	int degree);

int nurbs_calc_control3d(
	const POINT_3D *Q,
	int count,
	int degree,
	const double *uk_vector,
	const double *knots,
	int knot_count,
	POINT_3D *ctrls);

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
	POINT_2D *ctrls);

#ifdef __cplusplus
}
#endif

#endif
