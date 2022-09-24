/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : polarangle.c                                   */
/*                                                                           */
/* Created:       2012/06/14 (JLe)                                           */
/* Last modified: 2017/02/01 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Calculates polar angle from coordinates                      */
/*                                                                           */
/* Comments: Return value is between 0 and 2*PI. Return value is zero if     */
/*           vector length is zero.                                          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PolarAngle:"

/*****************************************************************************/

double PolarAngle(double x, double y)
{
  double t;

  /* Check quadrature */

  if ((x == 0.0) && (y > 0.0))
    t = PI/2.0;
  else if ((x == 0.0) && (y < 0.0))
    t = 3.0*PI/2.0;
  else if ((x > 0.0) && (y == 0.0))
    t = 0.0;
  else if ((x < 0.0) && (y == 0.0))
    t = PI;
  else if ((x > 0.0) && (y > 0.0))
    t = atan(y/x);
  else if ((x < 0.0) && (y > 0.0))
    t = atan(y/x) + PI;
  else if ((x < 0.0) && (y < 0.0))
    t = atan(y/x) + PI;
  else if ((x > 0.0) && (y < 0.0))
    t = atan(y/x) + 2.0*PI;
  else
    {
      /* Special case: vector length is zero */

      return 0.0;
    }

  /* Return angle */

  return t;
}

/*****************************************************************************/
