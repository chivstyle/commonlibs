/**
// (C) 2012 CHIV, all rights reserved
//
// \brief an implementation of NURBS
// this program is free, you can redistribute it without any permission
// writen or oral.but you should keep it free!
//
// \version : 1.0
*/
#ifndef __NURBS_H
#define __NURBS_H

#define NURBS_API

#ifdef __cplusplus
extern "C" {
#endif

#define _VERBOSE   1
#define _POINT_3D  1
#define _POINT_4D  2

typedef struct POINT_2D_ {
	double x;
	double y;
} POINT_2D;

typedef struct POINT_2Di_ {
	int x;
	int y;
} POINT_2Di;

typedef struct POINT_3D_ {
	double x;
	double y;
	double z;
} POINT_3D;

typedef struct POINT_4D_ {
	double x;
	double y;
	double z;
	double w;
} POINT_4D;

NURBS_API
double nurbs(
	const double * knots,     /* knots vector */
	int      count,     /* number of knots */
	int      degree,    /* degree */
	int      index,     /* index */
	int      count_pt,
	double   slice
);
/*!
  return     -1   knots is not large enough
  return      0   ok
*/
NURBS_API
int MakeClampedUniformKnotsVector(double *knots, int knots_count, int cvcount, int degree);
NURBS_API
int MakePeriodicUniformKnotsVector(double *knots, int knots_count, int cvcount, int degree);

NURBS_API
int MakeClampedChordKnotsVector(double *knots, int knots_count, int cvcount, int degree);
NURBS_API
int MakePeriodicChordKnotsVector(double *knots, int knots_count, int cvcount, int degree);

/* standard */
/**
 \brief *_2D?, *_3D? is Bspline, *_4D? is NURBS
*/
NURBS_API
POINT_2D NURBS_2D(
	double * knots, 
	int count, 
	int degree, 
	int index, 
	double slice, 
	POINT_2D * point, 
	int count_pt);
NURBS_API
POINT_3D NURBS_3D(
	double * knots, 
	int count, 
	int degree, 
	int index, 
	double slice, 
	POINT_3D * point, 
	int count_pt);
/**
 \brief
 the member 'w' of the structure POINT_4D is the omega
*/
NURBS_API
POINT_4D NURBS_4D(
	double * knots, 
	int count, 
	int degree, 
	int index, 
	double slice, 
	POINT_4D * point, 
	int count_pt);

#ifdef __cplusplus
}
#endif

#endif
