/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : divsubdistence.c                               */
/*                                                                           */
/* Created:       2020/01/15 (JLe)                                           */
/* Last modified: 2020/05/08 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Calculates distance to divisor boundaries.                   */
/*                                                                           */
/* Comments: - Only radial subdivision is currently supported.               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DivSubDistance:"

/*****************************************************************************/

double DivSubDistance(long div, long uni, double x, double y, double z,
                      double u, double v, double w, double t)
{
  long ptr, nx, ny, nz, nrad, nseg, i, j, k;
  double r, phi, r0, r1, dmin, d, px, py, pz, dx, dy, dz;
  double xmin, xmax, ymin, ymax, zmin, zmax;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(div)", DATA_ARRAY, div);

  /* Get bin sizes */

  nx = (long)RDB[div + DIV_NX];
  ny = (long)RDB[div + DIV_NY];
  nz = (long)RDB[div + DIV_NZ];
  nrad = (long)RDB[div + DIV_NRAD];
  nseg = (long)RDB[div + DIV_NSEG];

  /* Check for single-valued */

  if (nx*ny*nz*nrad*nseg == 1)
    return INFTY;

  /* Coordinate transformation to cold conditions */

  if ((ptr = (long)RDB[uni + UNIVERSE_PTR_IFC_FUEP]) > VALID_PTR)
    CoordExpans(ptr, &x, &y, &z, t, 1);

  /* Reset indexes */

  i = 0;
  j = 0;
  k = 0;

  /* Reset minimum distance */

  dmin = INFTY;

  /* Check type */

  if ((nx > 1) || (ny > 1))
    {
      /***********************************************************************/

      /***** Cartesian division **********************************************/

      /* x-coordinate */

      if (nx > 1)
        {
          /* Minimum and maximum coordinates */

          xmin = RDB[div + DIV_XMIN];
          CheckValue(FUNCTION_NAME, "xmin", "", xmin, -INFTY, INFTY);

          xmax = RDB[div + DIV_XMAX];
          CheckValue(FUNCTION_NAME, "xmax", "", xmax, xmin, INFTY);

          /* Calculate pitch */

          px = (xmax - xmin)/((double)nx);
          CheckValue(FUNCTION_NAME, "px", "", px, ZERO, INFTY);

          /* Calculate mesh indexes */

          i = (long)((x - xmin)/px);
          CheckValue(FUNCTION_NAME, "i", "", i, 0, nx - 1);

          /* Check direction */

          if (u != 0.0)
            {
              /* Calculate point distance to mesh wall */

              dx = x - xmin - (double)i*px;
              CheckValue(FUNCTION_NAME, "dx", "", dx, -0.5*px, 0.5*px);

              /* Calculate optical distance to mesh wall */

              if ((d = -dx/u) >= 0.0)
                {
                  if (d < dmin)
                    dmin = d;
                }
              else if ((d = -(dx - px)/u) > 0.0)
                {
                  if (d < dmin)
                    dmin = d;
                }
              else
                Warn(FUNCTION_NAME, "dx = %E, u = %E, d = %E", dx, u, d);
            }
        }

      /* y-coordinate */

      if (ny > 1)
        {
          /* Minimum and maximum coordinates */

          ymin = RDB[div + DIV_YMIN];
          CheckValue(FUNCTION_NAME, "ymin", "", ymin, -INFTY, INFTY);

          ymax = RDB[div + DIV_YMAX];
          CheckValue(FUNCTION_NAME, "ymax", "", ymax, ymin, INFTY);

          /* Calculate pitch */

          py = (ymax - ymin)/((double)ny);
          CheckValue(FUNCTION_NAME, "py", "", py, ZERO, INFTY);

          /* Calculate mesh indexes */

          j = (long)((y - ymin)/py);
          CheckValue(FUNCTION_NAME, "j", "", j, 0, ny - 1);

          /* Check direction */

          if (v != 0.0)
            {
              /* Calculate point distance to mesh wall */

              dy = y - ymin - (double)j*py;
              CheckValue(FUNCTION_NAME, "dy", "", dy, -0.5*py, 0.5*py);

              /* Calculate optical distance to mesh wall */

              if ((d = -dy/v) >= 0.0)
                {
                  if (d < dmin)
                    dmin = d;
                }
              else if ((d = -(dy - py)/v) > 0.0)
                {
                  if (d < dmin)
                    dmin = d;
                }
              else
                Warn(FUNCTION_NAME, "dy = %E, v = %E, d = %E", dy, v, d);
            }
        }

      /* z-coordinate */

      if (nz > 1)
        {
          /* Minimum and maximum coordinates */

          zmin = RDB[div + DIV_ZMIN];
          CheckValue(FUNCTION_NAME, "zmin", "", zmin, -INFTY, INFTY);

          zmax = RDB[div + DIV_ZMAX];
          CheckValue(FUNCTION_NAME, "zmax", "", zmax, zmin, INFTY);

          /* Calculate pitch */

          pz = (zmax - zmin)/((double)nz);
          CheckValue(FUNCTION_NAME, "pz", "", pz, ZERO, INFTY);

          /* Calculate mesh indexes */

          k = (long)((z - zmin)/pz);
          CheckValue(FUNCTION_NAME, "k", "", k, 0, nz - 1);

          /* Check direction */

          if (w != 0.0)
            {
              /* Calculate point distance to mesh wall */

              dz = z - zmin - (double)k*pz;
              CheckValue(FUNCTION_NAME, "dz", "", dz, -0.5*pz, 0.5*pz);

              /* Calculate optical distance to mesh wall */

              if ((d = -dz/w) >= 0.0)
                {
                  if (d < dmin)
                    dmin = d;
                }
              else if ((d = -(dz - pz)/w) > 0.0)
                {
                  if (d < dmin)
                    dmin = d;
                }
              else
                Warn(FUNCTION_NAME, "dz = %E, w = %E, d = %E", dz, w, d);
            }
        }

      /* Return distance */

      return dmin;
    }
  else
    {
      /***********************************************************************/

      /***** Radial division *************************************************/

      /* z-coordinate */

      if (nz > 1)
        {
          /* Minimum and maximum coordinates */

          zmin = RDB[div + DIV_ZMIN];
          CheckValue(FUNCTION_NAME, "zmin", "", zmin, -INFTY, INFTY);

          zmax = RDB[div + DIV_ZMAX];
          CheckValue(FUNCTION_NAME, "zmax", "", zmax, zmin, INFTY);

          /* Calculate pitch */

          pz = (zmax - zmin)/((double)nz);
          CheckValue(FUNCTION_NAME, "pz", "", pz, ZERO, INFTY);

          /* Calculate mesh indexes */

          i = (long)((z - zmin)/pz);
          CheckValue(FUNCTION_NAME, "i", "", i, 0, nz - 1);

          /* Check direction */

          if (w != 0.0)
            {
              /* Calculate point distance to mesh wall */

              dz = z - zmin - (double)i*pz;
              CheckValue(FUNCTION_NAME, "dz", "", dz, -0.5*pz, 0.5*pz);

              /* Calculate optical distance to mesh wall */

              if ((d = -dz/w) >= 0.0)
                {
                  if (d < dmin)
                    dmin = d;
                }
              else if ((d = -(dz - pz)/w) > 0.0)
                {
                  if (d < dmin)
                    dmin = d;
                }
              else
                Warn(FUNCTION_NAME, "dz = %E, w = %E, d = %E", dz, w, d);
            }
        }

      /* Radial coordinate */

      if (nrad > 1)
        {
          /* Calculate square radius */

          r = x*x + y*y;

          /* Get limits */

          r0 = RDB[div + DIV_RMIN];
          r1 = RDB[div + DIV_RMAX];

          CheckValue(FUNCTION_NAME, "r0", "", r0, 0.0, r1);
          CheckValue(FUNCTION_NAME, "r1", "", r1, r0, INFTY);

          /* Calculate portion of volume inside this radius */

          r = (r - r0*r0)/(r1*r1 - r0*r0);

          /* Check mode */

          if (((long)RDB[div + DIV_LIMS_CHECK] == YES) &&
              ((r < 0.0) || (r > 1.0)))
            Error(div, "Radial dimension doesn't cover region");
          else if (r == 1.0)
            j = nrad - 1;
          else if (r < 0.0)
            j = -1;
          else
            j = (long)(r*nrad);

          /* Outer radius */

          r = sqrt(((double)(j + 1))/((double)nrad))*(r1 - r0) + r0;

          /* Calculate distance */

          if ((d = CylDis(x, y, u, v, r)) > 0.0)
            if (d < dmin)
              dmin = d;

          /* Inner radius */

          if ((r = sqrt(((double)j)/((double)nrad))*(r1 - r0) + r0) > 0.0)
            if ((d = CylDis(x, y, u, v, r)) > 0.0)
              if (d < dmin)
                dmin = d;
        }

      /* Angular segment */

      if (nseg > 1)
        {
          /* Calculate angle */

          phi = PolarAngle(x, y);

          /* Add tilt */

          phi = phi + RDB[div + DIV_SEG0];

          /* Normalize */

          phi = phi/(2.0*PI);
          phi = phi - (double)((long)phi);

          /* Check phi */

          CheckValue(FUNCTION_NAME, "phi", "", phi, 0.0, 1.0);

          /* Calculate index */

          k = (long)(phi*nseg);
        }

      /* Return distance */

      return dmin;

      /***********************************************************************/
    }

  /***************************************************************************/

  /* Return infinity */

  return INFTY;
}

/*****************************************************************************/
