/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scalarprod.c                                   */
/*                                                                           */
/* Created:       2016/10/24 (JLe)                                           */
/* Last modified: 2016/10/24 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Calculates scalar product of two 3D vectors                  */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScalarProd:"

/*****************************************************************************/

double ScalarProd(const double *v1, const double *v2)
{
  double d;

  /* Calculate scalar product */
  
  d = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
  CheckValue(FUNCTION_NAME, "d", "", d, -INFTY, INFTY);

  /* Return value */

  return d;
}

/*****************************************************************************/
