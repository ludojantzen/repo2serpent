/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : cyldis.c                                       */
/*                                                                           */
/* Created:       2015/10/07 (JLe)                                           */
/* Last modified: 2015/10/07 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Calculates distence to right circular cylinder on z-axis     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "CylDis:"

/*****************************************************************************/

double CylDis(double x, double y, double u, double v, double r)
{
  double a, b, c, d0, d;

  /* Calculate constants */
  
  if ((a = u*u + v*v) == 0.0)
    return INFTY;
        
  b = u*x + v*y;      
  c = x*x + y*y - r*r;
  
  if ((d0 = b*b - a*c) < 0.0)
    {
      /* No line-of-sight */
      
      return INFTY;
    }
  else if (c < 0)
    {
      /* Point is inside, only one root possible */
      
      d = -(b - sqrt(d0))/a;
    }
  else if ((d = -(b + sqrt(d0))/a) < 0.0)
    {
      /* Point is outside, line-of-sight, opposite direction */
      
      return INFTY;
    }
  
  /* Return distance */
  
  return d;
}

/*****************************************************************************/
