/**
// (C) 2012 CHIV, all rights reserved
//
// \brief status of the routines in this lib
// \author chiv
*/
#ifndef _C_NURBS_ERROR_H
#define _C_NURBS_ERROR_H

#define E_OK           0          /** sucess */
#define E_NO_INV       1          /** matrix doesn't have inverse */
#define E_SINGULAR     2          /** matrix is singular */
#define E_NO_RESOURCE  3          /** no enough memory or can't open/create file */
#define E_NO_SOLUTION  4          /** system equation doesn't have solution */
#define E_BAD_MATRIX   5          /** the matrix you input is invalid or isn't proper */
#define E_UNKNOWN      6          /** the unknown result */
#define E_BAD_PARA     7          /** the parameter you give isn't proper */

#endif
