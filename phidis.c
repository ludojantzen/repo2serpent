/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : phidis.c                                       */
/*                                                                           */
/* Created:       2016/05/26 (JLe)                                           */
/* Last modified: 2016/05/31 (JLe)                                           */
/* Version:       2.1.27                                                     */
/*                                                                           */
/* Description: Calculates distence to angular sector                        */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "PhiDis:"

/*****************************************************************************/

double PhiDis(double x, double y, double u, double v, double phi)
{  
  double k, d;

  /* Check special */
  
  if ((fabs(phi) < 1E-6) || 
      (fabs(phi - PI) < 1E-6) || 
      (fabs(phi - 2.0*PI) < 1E-6))
    {
      if (v == 0.0)
        return INFTY;
      else if ((d = -y/v) < 0.0)
        return INFTY;
      else
        return d;
    }
  else if ((fabs(phi - PI/2.0) < 1E-6) || 
           (fabs(phi - 3.0*PI/2.0) < 1E-6))
    {
      if (u == 0.0)
        return INFTY;
      else if ((d = -x/u) < 0.0)
        return INFTY;
      else
        return d;
    }
  
  /* Check quadrature */
  
  if (phi > 2.0*PI)
    Die(FUNCTION_NAME, "phi = %E", phi - 2.0*PI);
  else if (phi > 3.0*PI/2.0)
    phi = phi - 2.0*PI;
  else if (phi > PI/2.0)
    phi = phi - PI;

  /* Slope */

  k = sin(phi)/cos(phi);

  if ((v - k*u) == 0.0)
    return INFTY;
  else if ((d = (k*x - y)/(v - k*u)) < 0.0)
    return INFTY;
  else
    return d;
}

/*****************************************************************************/
