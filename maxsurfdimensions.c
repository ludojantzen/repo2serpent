/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : maxsurfacedimensions.c                         */
/*                                                                           */
/* Created:       2010/11/02 (JLe)                                           */
/* Last modified: 2019/12/18 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Returns maximum dimensions for each surface type             */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MaxSurfDimensions:"

/*****************************************************************************/

void MaxSurfDimensions(long surf, double *xmin, double *xmax, double *ymin,
                       double *ymax, double *zmin, double *zmax)
{
  long ptr, type, np, n;
  double x, y, z, r, h, l;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

  /* Get surface type */

  type = (long)RDB[surf + SURFACE_TYPE];

  /* Get number of parameters */

  np = (long)RDB[surf + SURFACE_N_PARAMS];

  /* Pointer to parameter list */

  ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];

  /* Check pointer */

  if (type != SURF_INF)
    {
      /* Braces needed */

      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
    }

  /* Get maximum values according to surface type */

  switch(type)
    {
      /***********************************************************************/

    case SURF_CYL:
    case SURF_CYLZ:
    case SURF_SQC:
      {
        x = RDB[ptr];
        y = RDB[ptr + 1];
        r = RDB[ptr + 2];

        if (x + r > *xmax)
          *xmax = x + r;
        if (x - r < *xmin)
          *xmin = x - r;
        if (y + r > *ymax)
          *ymax = y + r;
        if (y - r < *ymin)
          *ymin = y - r;

        /* Check number of parameters */

        if (np == 5)
          {
            /* Test top and bottom */

            if (RDB[ptr + 3] < *zmin)
              *zmin = RDB[ptr + 3];

            if (RDB[ptr + 4] > *zmax)
              *zmax = RDB[ptr + 4];
          }

        break;
      }

      /***********************************************************************/

    case SURF_CYLX:
      {
        y = RDB[ptr];
        z = RDB[ptr + 1];
        r = RDB[ptr + 2];

        if (y + r > *ymax)
          *ymax = y + r;
        if (y - r < *ymin)
          *ymin = y - r;
        if (z + r > *zmax)
          *zmax = z + r;
        if (z - r < *zmin)
          *zmin = z - r;

        /* Check number of parameters */

        if (np == 5)
          {
            /* Test top and bottom */

            if (RDB[ptr + 3] < *xmin)
              *xmin = RDB[ptr + 3];

            if (RDB[ptr + 4] > *xmax)
              *xmax = RDB[ptr + 4];
          }

        break;
      }

      /***********************************************************************/

    case SURF_CYLY:
      {
        x = RDB[ptr];
        z = RDB[ptr + 1];
        r = RDB[ptr + 2];

        if (x + r > *xmax)
          *xmax = x + r;
        if (x - r < *xmin)
          *xmin = x - r;
        if (z + r > *zmax)
          *zmax = z + r;
        if (z - r < *zmin)
          *zmin = z - r;

        /* Check number of parameters */

        if (np == 5)
          {
            /* Test top and bottom */

            if (RDB[ptr + 3] < *ymin)
              *ymin = RDB[ptr + 3];

            if (RDB[ptr + 4] > *ymax)
              *ymax = RDB[ptr + 4];
          }

        break;
      }

      /***********************************************************************/

    case SURF_CROSS:
      {
        x = RDB[ptr];
        y = RDB[ptr + 1];
        r = RDB[ptr + 3];

        if (x + r > *xmax)
          *xmax = x + r;
        if (x - r < *xmin)
          *xmin = x - r;
        if (y + r > *ymax)
          *ymax = y + r;
        if (y - r < *ymin)
          *ymin = y - r;

        break;
      }

      /***********************************************************************/

    case SURF_SPH:
    case SURF_CUBE:
      {
        x = RDB[ptr];
        y = RDB[ptr + 1];
        z = RDB[ptr + 2];
        r = RDB[ptr + 3];

        if (x + r > *xmax)
          *xmax = x + r;
        if (x - r < *xmin)
          *xmin = x - r;

        if (y + r > *ymax)
          *ymax = y + r;
        if (y - r < *ymin)
          *ymin = y - r;

        if (z + r > *zmax)
          *zmax = z + r;
        if (z - r < *zmin)
          *zmin = z - r;

        break;
      }

      /***********************************************************************/

    case SURF_CUBOID:
      {
        if (RDB[ptr] < *xmin)
          *xmin = RDB[ptr];
        if (RDB[ptr + 1] > *xmax)
          *xmax = RDB[ptr + 1];

        if (RDB[ptr + 2] < *ymin)
          *ymin = RDB[ptr + 2];
        if (RDB[ptr + 3] > *ymax)
          *ymax = RDB[ptr + 3];

        if (RDB[ptr + 4] < *zmin)
          *zmin = RDB[ptr + 4];
        if (RDB[ptr + 5] > *zmax)
          *zmax = RDB[ptr + 5];

        break;
      }

      /***********************************************************************/

    case SURF_RECT:
      {
        if (RDB[ptr] < *xmin)
          *xmin = RDB[ptr];
        if (RDB[ptr + 1] > *xmax)
          *xmax = RDB[ptr + 1];

        if (RDB[ptr + 2] < *ymin)
          *ymin = RDB[ptr + 2];
        if (RDB[ptr + 3] > *ymax)
          *ymax = RDB[ptr + 3];

        break;
      }

      /***********************************************************************/

    case SURF_HEXYC:
      {
        x = RDB[ptr];
        y = RDB[ptr + 1];
        r = RDB[ptr + 2];

        if (x + r/SIN60 > *xmax)
          *xmax = x + r/SIN60;
        if (x - r/SIN60 < *xmin)
          *xmin = x - r/SIN60;
        if (y + r > *ymax)
          *ymax = y + r;
        if (y - r < *ymin)
          *ymin = y - r;

        break;
      }

      /***********************************************************************/

    case SURF_HEXXC:
      {
        x = RDB[ptr];
        y = RDB[ptr + 1];
        r = RDB[ptr + 2];

        if (x + r/SIN60 > *ymax)
          *ymax = x + r/SIN60;
        if (x - r/SIN60 < *ymin)
          *ymin = x - r/SIN60;
        if (y + r > *xmax)
          *xmax = y + r;
        if (y - r < *xmin)
          *xmin = y - r;

        break;
      }

      /***********************************************************************/

    case SURF_HEXXAP:
      {
        x = RDB[ptr];
        y = RDB[ptr + 1];

        *xmin = x - RDB[ptr + 2];
        *xmax = x + RDB[ptr + 2];
        *ymin = y - 2.0/SQRT3*RDB[ptr + 2];
        *ymax = y + 2.0/SQRT3*RDB[ptr + 2];

        break;
      }

      /***********************************************************************/

    case SURF_HEXYAP:
      {
        x = RDB[ptr];
        y = RDB[ptr + 1];

        *ymin = y - RDB[ptr + 2];
        *ymax = y + RDB[ptr + 2];
        *xmin = x - 2.0/SQRT3*RDB[ptr + 2];
        *xmax = x + 2.0/SQRT3*RDB[ptr + 2];

        break;
      }

      /***********************************************************************/

    case SURF_HEXYPRISM:
      {
        x = RDB[ptr];
        y = RDB[ptr + 1];
        r = RDB[ptr + 2];

        if (x + r/SIN60 > *xmax)
          *xmax = x + r/SIN60;
        if (x - r/SIN60 < *xmin)
          *xmin = x - r/SIN60;
        if (y + r > *ymax)
          *ymax = y + r;
        if (y - r < *ymin)
          *ymin = y - r;

        if (RDB[ptr + 3] < *zmin)
          *zmin = RDB[ptr + 3];
        if (RDB[ptr + 4] > *zmax)
          *zmax = RDB[ptr + 4];

        break;
      }

      /***********************************************************************/

    case SURF_HEXXPRISM:
      {
        x = RDB[ptr];
        y = RDB[ptr + 1];
        r = RDB[ptr + 2];

        if (x + r/SIN60 > *ymax)
          *ymax = x + r/SIN60;
        if (x - r/SIN60 < *ymin)
          *ymin = x - r/SIN60;
        if (y + r > *xmax)
          *xmax = y + r;
        if (y - r < *xmin)
          *xmin = y - r;

        if (RDB[ptr + 3] < *zmin)
          *zmin = RDB[ptr + 3];
        if (RDB[ptr + 4] > *zmax)
          *zmax = RDB[ptr + 4];

        break;
      }

      /***********************************************************************/

    case SURF_PX:
      {
        x = RDB[ptr];

        if (x > *xmax)
          *xmax = x ;
        if (x < *xmin)
          *xmin = x;
            break;
      }

      /***********************************************************************/

    case SURF_PY:
      {
        y = RDB[ptr];

        if (y > *ymax)
          *ymax = y ;
        if (y < *ymin)
          *ymin = y;
        break;
      }

      /***********************************************************************/

    case SURF_PZ:
      {
        z = RDB[ptr];

        if (z > *zmax)
          *zmax = z ;
        if (z < *zmin)
          *zmin = z;
        break;
      }

      /***********************************************************************/

    case SURF_CKX:
      {
        if (np == 5)
          {
            x = RDB[ptr];
            h = RDB[ptr + 4];

            if (h < 0.0)
              {
                if (x > *xmax)
                  *xmax = x;
              }
            else
              {
                if (x < *xmin)
                  *xmin = x;
              }
          }

        break;
      }

      /***********************************************************************/

    case SURF_CKY:
      {
        if (np == 5)
          {
            y = RDB[ptr + 1];
            h = RDB[ptr + 4];

            if (h < 0.0)
              {
                if (y > *ymax)
                  *ymax = y;
              }
            else
              {
                if (y < *ymin)
                  *ymin = y;
              }
          }

        break;
      }

      /***********************************************************************/

    case SURF_CKZ:
      {
        if (np == 5)
          {
            z = RDB[ptr + 2];
            h = RDB[ptr + 4];

            if (h < 0.0)
              {
                if (z > *zmax)
                  *zmax = z;
              }
            else
              {
                if (z < *zmin)
                  *zmin = z;
              }
          }

        break;
      }

      /***********************************************************************/

    case SURF_CONE:
      {
        z = RDB[ptr + 2];
        h = RDB[ptr + 4];

        if (h > 0.0)
          {
            if (z + h > *zmax)
              *zmax = z + h;
          }
        else
          {
            if (z + h < *zmin)
              *zmin = z + h;
          }

        break;
      }

      /***********************************************************************/

    case SURF_INF:
      {
        *xmin =  INFTY;
        *ymin =  INFTY;
        *zmin =  INFTY;
        *xmax = -INFTY;
        *ymax = -INFTY;
        *zmax = -INFTY;

        break;
      }

      /***********************************************************************/

    case SURF_PAD:
      {
        /* NOTE: Boundaries calculated according to outer radius. */

        x = RDB[ptr];
        y = RDB[ptr + 1];
        r = RDB[ptr + 3];

        if (x + r > *xmax)
          *xmax = x + r;
        if (x - r < *xmin)
          *xmin = x - r;
        if (y + r > *ymax)
          *ymax = y + r;
        if (y - r < *ymin)
          *ymin = y - r;

        break;
      }

      /***********************************************************************/

    case SURF_DODE:
      {
        x = RDB[ptr];
        y = RDB[ptr + 1];
        r = RDB[ptr + 2];

        if (np == 4)
          l = RDB[ptr + 3];
        else
          l = RDB[ptr + 2];

        if (x + r > *xmax)
          *xmax = x + r;
        if (x - r < *xmin)
          *xmin = x - r;
        if (y + l > *ymax)
          *ymax = y + l;
        if (y - l < *ymin)
          *ymin = y - l;

        break;
      }

      /***********************************************************************/

    case SURF_OCTA:
    case SURF_ASTRA:
      {
        x = RDB[ptr];
        y = RDB[ptr + 1];
        r = RDB[ptr + 2];

        if (np == 4)
          l = RDB[ptr + 3];
        else
          l = RDB[ptr + 2];

        if (x + r > *xmax)
          *xmax = x + r;
        if (x - r < *xmin)
          *xmin = x - r;
        if (y + r > *ymax)
          *ymax = y + r;
        if (y - r < *ymin)
          *ymin = y - r;

        break;
      }

      /*******************************************************************/

    case SURF_PLANE:
    case SURF_MPLANE:
    case SURF_QUADRATIC:
    case SURF_PPD:
    case SURF_X:
    case SURF_Y:
    case SURF_Z:
    case SURF_TORX:
    case SURF_TORY:
    case SURF_TORZ:
    case SURF_CYLV:
    case SURF_RCC:
      {
        /* This is tricky */

        break;
      }

      /*******************************************************************/

    case SURF_GCROSS:
      {
        x = RDB[ptr];
        y = RDB[ptr + 1];
        r = RDB[ptr + 2];

        if (x + r > *xmax)
          *xmax = x + r;
        if (x - r < *xmin)
          *xmin = x - r;
        if (y + r > *ymax)
          *ymax = y + r;
        if (y - r < *ymin)
          *ymin = y - r;

        break;
      }

      /*******************************************************************/

    case SURF_BOX:
      {
        /* Loop over 18 points */

        for (n = 0; n < 18; n++)
          {
            x = RDB[ptr++];
            y = RDB[ptr++];
            z = RDB[ptr++];

            if (x > *xmax)
              *xmax = x;
            if (x < *xmin)
              *xmin = x;

            if (y > *ymax)
              *ymax = y;
            if (y < *ymin)
              *ymin = y;

            if (z > *zmax)
              *zmax = z;
            if (z < *zmin)
              *zmin = z;
          }

        break;
      }

      /***********************************************************************/

    case SURF_INVOLUTE:
      {
        x = RDB[ptr++];
        y = RDB[ptr++];
        r = RDB[ptr + 4];

        if (x + r > *xmax)
          *xmax = x + r;
        if (x - r < *xmin)
          *xmin = x - r;
        if (y + r > *ymax)
          *ymax = y + r;
        if (y - r < *ymin)
          *ymin = y - r;

        break;
      }

      /**********************************************************************/

    case SURF_TRIAG:
      {
        x = RDB[ptr++];
        y = RDB[ptr++];
        r = RDB[ptr++];

        /* Check if orientation is given */

        if (np > 3)
          {
            /* Get orientation (W-S-E-N) */

            n = (long)RDB[ptr++];
          }
        else
          {
            /* Default is north */

            n = 4;
          }

        /* Check orientation */

        if (n == 1)
          {
            /* West */

            if (x - 2.0*r < *xmin)
              *xmin = x - 2.0*r;
            if (x + r > *xmax)
              *xmax = x + r;
           if (y - 3.0/SQRT3*r < *ymin)
              *ymin = y - 3.0/SQRT3*r;
            if (y + 3.0/SQRT3*r > *ymax)
              *ymax = y + 3.0/SQRT3*r;
          }
        else if (n == 2)
          {
            /* South */

            if (y - 2.0*r < *ymin)
              *ymin = y - 2.0*r;
            if (y + r > *ymax)
              *ymax = y + r;
            if (x - 3.0/SQRT3*r < *xmin)
              *xmin = x - 3.0/SQRT3*r;
            if (x + 3.0/SQRT3*r > *xmax)
              *xmax = x + 3.0/SQRT3*r;
          }
        else if (n == 3)
          {
            /* East */

            if (x + 2.0*r > *xmax)
              *xmax = x + 2.0*r;
            if (x - r < *xmin)
              *xmin = x - r;
            if (y + 3.0/SQRT3*r > *ymax)
              *ymax = y + 3.0/SQRT3*r;
            if (y - 3.0/SQRT3*r < *ymin)
              *ymin = y - 3.0/SQRT3*r;
          }
        else if (n == 4)
          {
            /* North */

            if (y + 2.0*r > *ymax)
              *ymax = y + 2.0*r;
            if (y - r < *ymin)
              *ymin = y - r;
            if (x + 3.0/SQRT3*r > *xmax)
              *xmax = x + 3.0/SQRT3*r;
            if (x - 3.0/SQRT3*r < *xmin)
              *xmin = x - 3.0/SQRT3*r;
          }
        else
          Die(FUNCTION_NAME, "Horrible error");

        /* Check number of parameters */

        if (np == 6)
          {
            /* Test top and bottom */

            if (RDB[ptr] < *zmin)
              *zmin = RDB[ptr];

            if (RDB[ptr + 1] > *zmax)
              *zmax = RDB[ptr + 1];
          }
        else if (np == 5)
          Error(surf, "Invalid number of surface parameters");

        break;
      }

      /**********************************************************************/

    case SURF_USER:
      {
        /* User-defined surface, call special routine */

        UserSurf(1, np, &RDB[ptr], xmin, xmax, ymin, ymax, zmin, zmax, NULL,
                 NULL, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0);

        break;
      }

      /**********************************************************************/

    default:
      {
        Die(FUNCTION_NAME, "Invalid surface type %ld (surface %s)\n",
            type, GetText(surf + SURFACE_PTR_NAME));
      }

      /**********************************************************************/
    }
}

/*****************************************************************************/
