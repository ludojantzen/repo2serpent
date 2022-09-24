/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : symmetryboundary.c                             */
/*                                                                           */
/* Created:       2014/05/28 (JLe)                                           */
/* Last modified: 2014/06/12 (JLe)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/* Description: Calculates minimum distance to symmetry boundary             */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SymmetryBoundary;"

/*****************************************************************************/

double SymmetryBoundary(long sym, double x, double y, double z,
                        double u, double v, double w)
{
  double t0, rot, params[9], d, min, x1, y1, z1, x2, y2, z2;
  
  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(sym)", DATA_ARRAY, sym);

  /* Get angles */
  
  t0 = RDB[sym + SYMMETRY_THETA0];
  rot = RDB[sym + SYMMETRY_ROT];

  /* Check origin (doesn't work if not zeros) */

  if ((RDB[sym + SYMMETRY_X0] != 0.0) || (RDB[sym + SYMMETRY_Y0]))
    Die(FUNCTION_NAME, "Not centered at origin");

  /* Avoid compiler warning */

  x1 = 0.0;
  y1 = 0.0;
  z1 = 0.0;

  x2 = 0.0;
  y2 = 0.0;
  z2 = 0.0;

  /* Check axis */

  if ((long)RDB[sym + SYMMETRY_AXIS] == 1)
    {
      /* Set points */

      x1 = 10.0;
      y1 = 0.0;
      z1 = 0.0;

      x2 = 0.0;
      y2 = sin(t0);
      z2 = cos(t0);
    }
  else if ((long)RDB[sym + SYMMETRY_AXIS] == 2)
    {
      /* Set points */

      x1 = 0.0;
      y1 = 10.0;
      z1 = 0.0;

      x2 = cos(t0);
      y2 = 0.0;
      z2 = sin(t0);
    }
  else if ((long)RDB[sym + SYMMETRY_AXIS] == 3)
    {
      /* Set points */

      x1 = 0.0;
      y1 = 0.0;
      z1 = 10.0;

      x2 = cos(t0);
      y2 = sin(t0);
      z2 = 0.0;
    }
  else
    Die(FUNCTION_NAME, "Invalid axis");

  /* Reset minimum distance */

  min = INFTY;

  /* Set first two points */

  params[0] = 0.0;
  params[1] = 0.0;
  params[2] = 0.0;

  params[3] = x1;
  params[4] = y1;
  params[5] = z1;
  
  /* Third point for first symmetry line */

  params[6] = x2;
  params[7] = y2;
  params[8] = z2;

  /* Calculate distance */

  d = SurfaceDistance(-1, params, SURF_PLANE, 9, x, y, z, u, v, w, -1);
  
  /* Compare */

  if (d > 0.0)
    min = d;

  /* Check axis */

  if ((long)RDB[sym + SYMMETRY_AXIS] == 1)
    {
      /* Set point */

      x2 = 0.0;
      y2 = sin(t0 + rot);
      z2 = cos(t0 + rot);
    }
  else if ((long)RDB[sym + SYMMETRY_AXIS] == 2)
    {
      /* Set point */

      x2 = cos(t0 + rot);
      y2 = 0.0;
      z2 = sin(t0 + rot);
    }
  else if ((long)RDB[sym + SYMMETRY_AXIS] == 3)
    {
      /* Set point */

      x2 = cos(t0 + rot);
      y2 = sin(t0 + rot);
      z2 = 0.0;
    }
  else
    Die(FUNCTION_NAME, "Invalid axis");

  /* Third point for second symmetry line */

  params[6] = x2;
  params[7] = y2;
  params[8] = z2;

  d = SurfaceDistance(-1, params, SURF_PLANE, 9, x, y, z, u, v, w, -1);

  /* Compare */
  
  if ((d > 0.0) && (d < min))
    min = d;

  /* Return minimum distance */

  return min;
}

/****************************************************************************/
