/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : unisym.c                                       */
/*                                                                           */
/* Created:       2012/05/30 (JLe)                                           */
/* Last modified: 2019/10/12 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Coordinate transformations for universe symmetry             */
/*                                                                           */
/* Comments: - Re-written 30.5.2012 (2.1.6)                                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "UniSym:"

/*****************************************************************************/

void UniSym(long sym, double *x, double *y, double *z, double *u,
            double *v, double *w)
{
  long ax, n;
  double t, x0, y0, z0, u0, v0, w0, t0, rot, rot0, ref;
  double r1, r2, r3, r4, f1, f2, f3, f4;

  /* Check symmetry pointer */

  CheckPointer(FUNCTION_NAME, "(sym)", DATA_ARRAY, sym);

  /* Check coordinates and direction cosines */

  CheckValue(FUNCTION_NAME, "x", "1", *x, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "y", "1", *y, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "z", "1", *z, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "u", "1", *u, -1.00001, 1.00001);
  CheckValue(FUNCTION_NAME, "v", "1", *v, -1.00001, 1.00001);
  CheckValue(FUNCTION_NAME, "w", "1", *w, -1.00001, 1.00001);

  /* Get axis */

  ax = (long)RDB[sym + SYMMETRY_AXIS];

  /* Get angles */

  t0 = RDB[sym + SYMMETRY_THETA0];
  rot0 = RDB[sym + SYMMETRY_ROT];

  /* Avoid compiler warning */

  t = 0.0;

  /* Get polar angle */

  if (ax == 1)
    t = PolarAngle(*y, *z);
  else if (ax == 2)
    t = PolarAngle(*x, *z);
  else if (ax == 3)
    t = PolarAngle(*x, *y);
  else
    Die(FUNCTION_NAME, "Invalid axis");

  /* Adjust polar angle */

  if ((t = t - t0) < 0.0)
    t = t + 2.0*PI;

  /* Calculate sector */

  if ((n = (long)(t/rot0)) == 0)
    return;

  /* Check origin */

  if ((RDB[sym + SYMMETRY_X0] != 0.0) || (RDB[sym + SYMMETRY_Y0]))
    Die(FUNCTION_NAME, "Not centered at origin");

  /***************************************************************************/

  /***** Rotate coordinates and direction cosines ****************************/

  /* Calculate rotation angle */

  rot = -n*rot0;

  /* Check value */

  CheckValue(FUNCTION_NAME, "t1", "", t + rot, -1E-12, rot0);

  /* Make coordinate transformations */

  if (ax == 1)
    {
      /* Rotation around x-axis */

      y0 = *y;
      z0 = *z;
      v0 = *v;
      w0 = *w;

      r1 = cosl(rot);
      r2 = -sinl(rot);
      r3 = sinl(rot);
      r4 = cosl(rot);

      *y = y0*r1 + z0*r2;
      *z = y0*r3 + z0*r4;
      *v = v0*r1 + w0*r2;
      *w = v0*r3 + w0*r4;
    }
  else if (ax == 2)
    {
      /* Rotation around y-axis */

      x0 = *x;
      z0 = *z;
      u0 = *u;
      w0 = *w;

      r1 = cosl(rot);
      r2 = sinl(rot);
      r3 = -sinl(rot);
      r4 = cosl(rot);

      *x = x0*r1 + z0*r2;
      *z = x0*r3 + z0*r4;
      *u = u0*r1 + w0*r2;
      *w = u0*r3 + w0*r4;
    }
  else
    {
      /* Rotation around z-axis */

      x0 = *x;
      y0 = *y;
      u0 = *u;
      v0 = *v;

      r1 = cosl(rot);
      r2 = -sinl(rot);
      r3 = sinl(rot);
      r4 = cosl(rot);

      *x = x0*r1 + y0*r2;
      *y = x0*r3 + y0*r4;
      *u = u0*r1 + v0*r2;
      *v = u0*r3 + v0*r4;
    }

  /* Check */

  CheckValue(FUNCTION_NAME, "dis1", "", r1*r4 - r2*r3 - 1.0, -1E-9, 1E-9);

  /* Check coordinates and direction cosines */

  CheckValue(FUNCTION_NAME, "x", "2", *x, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "y", "2", *y, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "z", "2", *z, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "u", "2", *u, -1.00001, 1.00001);
  CheckValue(FUNCTION_NAME, "v", "2", *v, -1.00001, 1.00001);
  CheckValue(FUNCTION_NAME, "w", "2", *w, -1.00001, 1.00001);

  /* Check boundary condition */

  if ((long)RDB[sym + SYMMETRY_BC] == BC_PERIODIC)
    return;

  /* Check parity */

  if (!(n % 2))
    return;

  /***************************************************************************/

  /***** Reflect coordinates *************************************************/

  /* Avoid compiler warning */

  t = 0.0;

  /* Get polar angle */

  if (ax == 1)
    t = PolarAngle(*y, *z);
  else if (ax == 2)
    t = PolarAngle(*x, *z);
  else if (ax == 3)
    t = PolarAngle(*x, *y);
  else
    Die(FUNCTION_NAME, "Invalid axis");

  /* Adjust polar angle */

  if ((t = t - t0) < 0.0)
    t = t + 2.0*PI;

  /* Calculate reflection angle */

  ref = rot0 - 2.0*t;

  /* Make coordinate transformations */

  if (ax == 1)
    {
      /* Rotation around x-axis */

      y0 = *y;
      z0 = *z;

      f1 = cosl(ref);
      f2 = -sinl(ref);
      f3 = sinl(ref);
      f4 = cosl(ref);

      *y = y0*f1 + z0*f2;
      *z = y0*f3 + z0*f4;
    }
  else if (ax == 2)
    {
      /* Rotation around y-axis */

      x0 = *x;
      z0 = *z;

      f1 = cosl(ref);
      f2 = sinl(ref);
      f3 = -sinl(ref);
      f4 = cosl(ref);

      *x = x0*f1 + z0*f2;
      *z = x0*f3 + z0*f4;
    }
  else
    {
      /* Rotation around z-axis */

      x0 = *x;
      y0 = *y;

      f1 = cosl(ref);
      f2 = -sinl(ref);
      f3 = sinl(ref);
      f4 = cosl(ref);

      *x = x0*f1 + y0*f2;
      *y = x0*f3 + y0*f4;
    }

  /* Check */

  CheckValue(FUNCTION_NAME, "dis2", "", f1*f4 - f2*f3 - 1.0, -1E-9, 1E-9);

  /***************************************************************************/

  /***** Reflect direction cosines *******************************************/

  /* Avoid compiler warning */

  t = 0.0;

  /* Get polar angle */

  if (ax == 1)
    t = PolarAngle(*v, *w);
  else if (ax == 2)
    t = PolarAngle(*u, *w);
  else if (ax == 3)
    t = PolarAngle(*u, *v);
  else
    Die(FUNCTION_NAME, "Invalid axis");

  /* Adjust polar angle */

  if ((t = t - t0) < 0.0)
    t = t + 2.0*PI;

  /* Check value */

  CheckValue(FUNCTION_NAME, "t3", "", t, 0.0, 2.0*PI);

  /* Calculate reflection angle */

  ref = rot0 - 2.0*t;

  /* Make transformations on direction cosines */

  if (ax == 1)
    {
      /* Rotation around x-axis */

      v0 = *v;
      w0 = *w;

      f1 = cosl(ref);
      f2 = -sinl(ref);
      f3 = sinl(ref);
      f4 = cosl(ref);

      *v = v0*f1 + w0*f2;
      *w = v0*f3 + w0*f4;
    }
  else if (ax == 2)
    {
      /* Rotation around y-axis */

      u0 = *u;
      w0 = *w;

      f1 = cosl(ref);
      f2 = sinl(ref);
      f3 = -sinl(ref);
      f4 = cosl(ref);

      *u = u0*f1 + w0*f2;
      *w = u0*f3 + w0*f4;
    }
  else
    {
      /* Rotation around z-axis */

      u0 = *u;
      v0 = *v;

      f1 = cosl(ref);
      f2 = -sinl(ref);
      f3 = sinl(ref);
      f4 = cosl(ref);

      *u = u0*f1 + v0*f2;
      *v = u0*f3 + v0*f4;
    }

  /* Check */

  CheckValue(FUNCTION_NAME, "dis3", "", f1*f4 - f2*f3 - 1.0, -1E-9, 1E-9);

  /* Check coordinates and direction cosines */

  CheckValue(FUNCTION_NAME, "x", "3", *x, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "y", "3", *y, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "z", "3", *z, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "u", "3", *u, -1.00001, 1.00001);
  CheckValue(FUNCTION_NAME, "v", "3", *v, -1.00001, 1.00001);
  CheckValue(FUNCTION_NAME, "w", "3", *w, -1.00001, 1.00001);

  /***************************************************************************/
}

/*****************************************************************************/
