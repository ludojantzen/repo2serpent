/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : vectorprod.c                                   */
/*                                                                           */
/* Created:       2016/10/24 (JLe)                                           */
/* Last modified: 2016/10/24 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Calculates vector product of two 3D vectors                  */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "VectorProd:"

/*****************************************************************************/

void VectorProd(double *v, const double *v1, const double *v2)
{
  /* Calculate components */

  v[0] = v1[1]*v2[2] - v1[2]*v2[1];
  v[1] = v1[2]*v2[0] - v1[0]*v2[2];
  v[2] = v1[0]*v2[1] - v1[1]*v2[0];
  
  /* Check */

  CheckValue(FUNCTION_NAME, "vx", "", v[0], -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "vy", "", v[1], -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "vz", "", v[2], -INFTY, INFTY);
}

/*****************************************************************************/
