/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : surfacenormal.c                                */
/*                                                                           */
/* Created:       2014/02/07 (JLe)                                           */
/* Last modified: 2019/12/22 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Calculates normal vector for surface at (x,y,z)              */
/*                                                                           */
/* Comments: - cylx, cyly and shp have not been tested                       */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SurfaceNormal:"

/*****************************************************************************/

void SurfaceNormal(long surf, double x, double y, double z,
                   double *u, double *v, double *w, long id)
{
  double d, A, B, C, x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3, r, R;
  double u0, v0, w0, L;
  long ptr, type, np, n;

  /* Check surface pointer */

  CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

  /* Do coordinate transformation (use dummies for direction cosines) */

  x1 = 1.0;
  y1 = 0.0;
  z1 = 0.0;

  if ((ptr = (long)RDB[surf + SURFACE_PTR_TRANS]) > VALID_PTR)
    CoordTrans(ptr, &x, &y, &z, &x1, &y1, &z1, id);

  /* Get surface type */

  type = (long)RDB[surf + SURFACE_TYPE];

  /* Pointer to parameter list */

  ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];

  /* Check number of parameters */

#ifdef DEBUG

  if((type != SURF_INF) && (ptr < VALID_PTR))
    Die(FUNCTION_NAME, "Surface %s has no parameters",
        GetText(surf + SURFACE_PTR_NAME));
#endif

  /* Get number of parameters */

  np = (long)RDB[surf + SURFACE_N_PARAMS];
  n = 0;

  /***************************************************************************/

  switch(type)
    {
      /**********************************************************************/

    case SURF_PX:
      {
        /* Put vector */

        *u = 1.0;
        *v = 0.0;
        *w = 0.0;

        /* Break case */

        break;
      }

      /**********************************************************************/

    case SURF_PY:
      {
        /* Put vector */

        *u = 0.0;
        *v = 1.0;
        *w = 0.0;

        /* Break case */

        break;
      }

      /**********************************************************************/

    case SURF_PZ:
      {
        /* Put vector */

        *u = 0.0;
        *v = 0.0;
        *w = 1.0;

        /* Break case */

        break;
      }

      /***********************************************************************/

    case SURF_PLANE:
    case SURF_MPLANE:
      {
        /* Check number of parameters */

        if (np == 9)
          {
            /* Surface defined by three points, read coordinates */

            x1 = RDB[ptr + n++];
            y1 = RDB[ptr + n++];
            z1 = RDB[ptr + n++];
            x2 = RDB[ptr + n++];
            y2 = RDB[ptr + n++];
            z2 = RDB[ptr + n++];
            x3 = RDB[ptr + n++];
            y3 = RDB[ptr + n++];
            z3 = RDB[ptr + n++];

            /* Calculate coefficients */

            A = y2*z3 - y3*z2 - y1*(z3 - z2) + z1*(y3 - y2);
            B = z2*x3 - z3*x2 - z1*(x3 - x2) + x1*(z3 - z2);
            C = x2*y3 - x3*y2 - x1*(y3 - y2) + y1*(x3 - x2);
          }
        else
          {
            /* Surface defined by parameters, reset optional*/

            B = 0.0;
            C = 0.0;

            /* Get parameters */

            A = RDB[ptr + n++];

            if (n < np)
              B = RDB[ptr + n++];

            if (n < np)
              C = RDB[ptr + n++];
          }

        /* Calculate denominator */

        if ((d = sqrt(A*A + B*B + C*C)) == 0.0)
          Die(FUNCTION_NAME, "Error in surface definition");

        /* Calculate components */

        *u = A/d;
        *v = B/d;
        *w = C/d;

        /* Break case */

        break;
      }

      /**********************************************************************/

    case SURF_CYL:
    case SURF_CYLZ:
      {
        /* Shift origin */

        x = x - RDB[ptr];
        y = y - RDB[ptr + 1];
        r = RDB[ptr + 2];

        /* Calculate radius */

        R = sqrt(x*x + y*y);

        /* Check if point is on wall */

        if (fabs(R - r) < 2.0*EXTRAP_L)
          {
            /* Calculate normal vector */

            *u = x/R;
            *v = y/R;
            *w = 0.0;
          }
        else if (np == 5)
          {
            /* Compare to bottom */

            if (fabs(z - RDB[ptr + 3]) < 2.0*EXTRAP_L)
              *w = -1.0;
            else if (fabs(z - RDB[ptr + 4]) < 2.0*EXTRAP_L)
              *w = 1.0;
            else
              Die(FUNCTION_NAME, "Impossible position");

            /* Reset remaining */

            *u = 0.0;
            *v = 0.0;
          }
        else
          Die(FUNCTION_NAME, "Impossible position");

        /* Break case */

        break;
      }

      /**********************************************************************/

    case SURF_CYLX:
      {
        /* Shift origin */

        y = y - RDB[ptr];
        z = z - RDB[ptr + 1];
        r = RDB[ptr + 2];

        /* Calculate radius */

        R = sqrt(y*y + z*z);

        /* Check if point is on wall */

        if (fabs(R - r) < 2.0*EXTRAP_L)
          {
            /* Calculate normal vector */

            *u = 0.0;
            *v = y/R;
            *w = z/R;
          }
        else if (np == 5)
          {
            /* Compare to bottom */

            if (fabs(x - RDB[ptr + 3]) < 2.0*EXTRAP_L)
              *u = -1.0;
            else if (fabs(x - RDB[ptr + 4]) < 2.0*EXTRAP_L)
              *u = 1.0;
            else
              Die(FUNCTION_NAME, "Impossible position");

            /* Reset remaining */

            *v = 0.0;
            *w = 0.0;
          }
        else
          Die(FUNCTION_NAME, "Impossible position");

        /* Break case */

        break;
      }

      /**********************************************************************/

    case SURF_CYLY:
      {
        /* Shift origin */

        x = x - RDB[ptr];
        z = z - RDB[ptr + 1];
        r = RDB[ptr + 2];

        /* Calculate radius */

        R = sqrt(x*x + z*z);

        /* Check if point is on wall */

        if (fabs(R - r) < 2.0*EXTRAP_L)
          {
            /* Calculate normal vector */

            *u = x/R;
            *v = 0.0;
            *w = z/R;
          }
        else if (np == 5)
          {
            /* Compare to bottom */

            if (fabs(y - RDB[ptr + 3]) < 2.0*EXTRAP_L)
              *v = -1.0;
            else if (fabs(y - RDB[ptr + 4]) < 2.0*EXTRAP_L)
              *v = 1.0;
            else
              Die(FUNCTION_NAME, "Impossible position");

            /* Reset remaining */

            *u = 0.0;
            *w = 0.0;
          }
        else
          Die(FUNCTION_NAME, "Impossible position");

        /* Break case */

        break;
      }

      /**********************************************************************/

    case SURF_SPH:
      {
        /* Shift origin */

        x = x - RDB[ptr];
        y = y - RDB[ptr + 1];
        z = z - RDB[ptr + 2];
        r = RDB[ptr + 3];

        /* Calculate radius */

        R = sqrt(x*x + y*y + z*z);

        /* Check if point is on wall */

        if (fabs(R - r) < 2.0*EXTRAP_L)
          {
            /* Calculate normal vector */

            *u = x/R;
            *v = y/R;
            *w = z/R;
          }
        else
          Die(FUNCTION_NAME, "Impossible position");

        /* Break case */

        break;
      }

      /**********************************************************************/

    case SURF_CYLV:
    case SURF_RCC:
      {
        /* Get parameters */

        x0 = RDB[ptr + n++];
        y0 = RDB[ptr + n++];
        z0 = RDB[ptr + n++];
        u0 = RDB[ptr + n++];
        v0 = RDB[ptr + n++];
        w0 = RDB[ptr + n++];

        /* Calculate projected distance */

        L = (x - x0)*u0 + (y - y0)*v0 + (z - z0)*w0;

        /* Calculate vector distance from center line */

        *u = x - x0 - L*u0;
        *v = y - y0 - L*v0;
        *w = z - z0 - L*w0;

        /* Normalize */

        R = sqrt((*u)*(*u) + (*v)*(*v) + (*w)*(*w));

        /* Check if point is on wall */

        if (fabs(R - RDB[ptr + 6]) < 2.0*EXTRAP_L)
          {
            /* Calculate normal vector */

            *u = *u/R;
            *v = *v/R;
            *w = *w/R;

            /* Check */

            CheckValue(FUNCTION_NAME, "cylv normal vector", "",
                       *u*u0 + *v*v0 + *w*w0, -1E-9, 1E-9);
          }
        else if (type == SURF_RCC)
          {
            /* Check distance from top and bottom */

            if (fabs(L) < 2.0*EXTRAP_L)
              {
                /* Vector is along the axis */

                *u = u0;
                *v = v0;
                *w = w0;
              }
            else if (fabs((RDB[ptr + 8] - RDB[ptr + 7]) - L) < 2.0*EXTRAP_L)
              {
                /* Vector is opposite the axis */

                *u = -u0;
                *v = -v0;
                *w = -w0;
              }
            else
              Die(FUNCTION_NAME, "Impossible position");
          }
        else
          Die(FUNCTION_NAME, "Impossible position");

        /* Break case */

        break;
      }

      /**********************************************************************/

    default:
      Die(FUNCTION_NAME, "Unsupported surface type %ld", type);
    }

  /* Check direction cosines */

  CheckValue(FUNCTION_NAME, "Normal vector", "",
             (*u)*(*u) + (*v)*(*v) + (*w)*(*w) - 1.0, -1E-4, 1E-4);

  /***************************************************************************/

}

/*****************************************************************************/
