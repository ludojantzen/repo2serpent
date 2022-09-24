/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : nearestmeshboundary.c                          */
/*                                                                           */
/* Created:       2014/03/01 (JLe)                                           */
/* Last modified: 2019/12/21 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description:  Calculates distance nearest mesh boundary                   */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "NearestMeshBoundary:"

/*****************************************************************************/

double NearestMeshBoundary(long msh, double x, double y, double z,
                           double u, double v, double w, long *fail)
{
  long ptr, loc0, idx, n0, n1, n2, i, j, k, n, type;
  double min0, max0, min1, max1, min2, max2, min, l, p0, p1, p2, dx, dy, dz;
  double r, r0, r1, z0, z1, phi, pitch, x0, y0;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

  /* Reset minimum distance */

  min = INFTY;

  /* Get dimensions */

  min0 = RDB[msh + MESH_MIN0];
  max0 = RDB[msh + MESH_MAX0];
  min1 = RDB[msh + MESH_MIN1];
  max1 = RDB[msh + MESH_MAX1];
  min2 = RDB[msh + MESH_MIN2];
  max2 = RDB[msh + MESH_MAX2];

  /***************************************************************************/

  /***** Distance to outer boundary ******************************************/

  /* Reset index */

  idx = -1;

  /* Get mesh type */

  type = (long)RDB[msh + MESH_TYPE];

  /* Get mesh index */

  if (((type != MESH_TYPE_CYLINDRICAL) && (type != MESH_TYPE_ICYL) &&
       (type != MESH_TYPE_HEXX) && (type != MESH_TYPE_HEXY)) &&
      ((idx = MeshIndex(msh, x, y, z, -1.0)) < 0))
    {
      /* Check boundaries to be sure */

#ifdef DEBUG

      if ((x > min0) && (x < max0) && (y > min1) && (y < max1) &&
          (z > min2) && (z < max2))
        Warn(FUNCTION_NAME, "Point is inside mesh");

#endif

      /* Check if point is outside x-boundaries */

      if ((l = x - min0) < 0.0)
        {
          /* Below lower */

          if (u != 0.0)
            if ((l = -l/u) >= 0.0)
              if (l < min)
                min = l;
        }
      else if ((l = x - max0) > 0.0)
        {
          /* Above upper */

          if (u != 0.0)
            if ((l = -l/u) >= 0.0)
              if (l < min)
                min = l;
        }

      /* Check if point is outside y-boundaries */

      if ((l = y - min1) < 0.0)
        {
          /* Below lower */

          if (v != 0.0)
            if ((l = -l/v) >= 0.0)
              if (l < min)
                min = l;
        }
      else if ((l = y - max1) > 0.0)
        {
          /* Above upper */

          if (v != 0.0)
            if ((l = -l/v) >= 0.0)
              if (l < min)
                min = l;
        }

      /* Check if point is outside z-boundaries */

      if ((l = z - min2) < 0.0)
        {
          /* Below lower */

          if (w != 0.0)
            if ((l = -l/w) >= 0.0)
              if (l < min)
                min = l;
        }
      else if ((l = z - max2) > 0.0)
        {
          /* Above upper */

          if (w != 0.0)
            if ((l = -l/w) >= 0.0)
              if (l < min)
                min = l;
        }

      /* Return distance */

      return min;
    }

  /***************************************************************************/

  /***** Adaptive mesh *******************************************************/

  if (type == MESH_TYPE_ADAPTIVE)
    {
      /* Check content */

      if ((long)RDB[msh + MESH_CONTENT] != MESH_CONTENT_PTR)
        Die(FUNCTION_NAME, "Invalid content type");

      /* Check index */

      if (idx < 0)
        Die(FUNCTION_NAME, "WTF?");

      /* Get direct pointer to data */

      ptr = (long)RDB[msh + MESH_PTR_PTR] + idx;
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Check pointer */

      if ((loc0 = (long)RDB[ptr]) == NULLPTR)
        {
          min = NearestMeshBoundary(-loc0, x, y, z, u, v, w, fail);

          /* Tää voi olla ongelma */
          /*
            Die(FUNCTION_NAME, "WTF?");
          */
        }
      else if (loc0 < -VALID_PTR)
        {
          /* Pointer to new mesh, call recursively */

          min = NearestMeshBoundary(-loc0, x, y, z, u, v, w, fail);
        }

      /* Switch type to Cartesian */

      type = MESH_TYPE_CARTESIAN;
    }

  /***************************************************************************/

  /***** Distance to cell boundaries *****************************************/

  /* Get sizes */

  n0 = (long)RDB[msh + MESH_N0];
  n1 = (long)RDB[msh + MESH_N1];
  n2 = (long)RDB[msh + MESH_N2];

  /* Calculate mesh pitch */

  p0 = (max0 - min0)/((double)n0);
  p1 = (max1 - min1)/((double)n1);
  p2 = (max2 - min2)/((double)n2);

  /* Check type */

  if (type == MESH_TYPE_CARTESIAN)
    {
      /***********************************************************************/

      /***** Cartesian mesh **************************************************/

      /* Check pitches */

      CheckValue(FUNCTION_NAME, "p0", "", p0, ZERO, INFTY);
      CheckValue(FUNCTION_NAME, "p1", "", p1, ZERO, INFTY);
      CheckValue(FUNCTION_NAME, "p2", "", p2, ZERO, INFTY);

      /* Calculate mesh indexes */

      i = (long)((x - min0)/p0);
      j = (long)((y - min1)/p1);
      k = (long)((z - min2)/p2);

      /* Check limits */

      CheckValue(FUNCTION_NAME, "i", "", i, 0, n0 - 1);
      CheckValue(FUNCTION_NAME, "j", "", j, 0, n1 - 1);
      CheckValue(FUNCTION_NAME, "k", "", k, 0, n2 - 1);

      /* Calculate distance to mesh walls */

      dx = x - min0 - (double)i*p0;
      dy = y - min1 - (double)j*p1;
      dz = z - min2 - (double)k*p2;

      /* Check values (NOTE: toi nolla pitää sallia tässä ja ottaa alemmissa */
      /* if-lauseissa huomioon yhtäsuuruutena, sillä desimaalip1öristyksen   */
      /* takia piste voi osua suoraan pinnalle) */

      CheckValue(FUNCTION_NAME, "dx", "", dx, -1E-12, p0);
      CheckValue(FUNCTION_NAME, "dy", "", dy, -1E-12, p1);
      CheckValue(FUNCTION_NAME, "dz", "", dz, -1E-12, p2);

      /* Calculate distance to mesh cell walls */

      if (u != 0.0)
        {
          if ((l = -dx/u) >= 0.0)
            {
              if (l < min)
                min = l;
            }
          else if ((l = -(dx - p0)/u) > 0.0)
            {
              if (l < min)
                min = l;
            }
          else
            Warn(FUNCTION_NAME, "dx = %E, u = %E, l = %E", dx, u, l);
        }

      if (v != 0.0)
        {
          if ((l = -dy/v) >= 0.0)
            {
              if (l < min)
                min = l;
            }
          else if ((l = -(dy - p1)/v) > 0.0)
            {
              if (l < min)
                min = l;
            }
          else
            Warn(FUNCTION_NAME, "dy = %E, v = %E, l = %E", dy, v, l);
        }

      if (w != 0.0)
        {
          if ((l = -dz/w) >= 0.0)
            {
              if (l < min)
                min = l;
            }
          else if ((l = -(dz - p2)/w) > 0.0)
            {
              if (l < min)
                min = l;
            }
          else
            Warn(FUNCTION_NAME, "dz = %E, w = %E, l = %E", dz, w, l);
        }

      /* Check fail condition for STLRayTest() */

      if (fail != NULL)
        {
          /* Count vertices */

          n = 0;

          if (fabs(x + min*u - min0) < 1E-6)
            n++;
          else if (fabs(x + min*u - max0) < 1E-6)
            n++;

          if (fabs(y + min*v - min1) < 1E-6)
            n++;
          else if (fabs(y + min*v - max1) < 1E-6)
            n++;

          if (fabs(z + min*w - min2) < 1E-6)
            n++;
          else if (fabs(z + min*w - max2) < 1E-6)
            n++;

          /* Check and set fail flag */

          if (n > 1)
            *fail = YES;
          else
            *fail = NO;
        }

      /***********************************************************************/
    }
  else if ((type == MESH_TYPE_HEXX) || (type == MESH_TYPE_HEXY))
    {
      /***********************************************************************/

      /***** Cartesian mesh **************************************************/

      /* Variables: min0 = x0, min1 = y0, max0 = pitch */

      x0 = min0;
      y0 = min1;
      pitch = max0;

      /* Coordinate transformation relative to origin  */

      x = x - x0;
      y = y - y0;

      /* Adjust if even number of cells */

      x = x - (1 - (n0 % 2))*0.5*pitch;
      y = y - (1 - (n1 % 2))*0.5*pitch;

      /* Get hex indexes */

      if ((long)RDB[msh + MESH_TYPE] == MESH_TYPE_HEXX)
        GetLatticeIndexes(pitch, pitch, 1.0, x, y, 0.0, &i, &j, &k,
                          LAT_TYPE_HX);
      else
        GetLatticeIndexes(pitch, pitch, 1.0, x, y, 0.0, &i, &j, &k,
                          LAT_TYPE_HY);

      /* Transfer coordinates */

      if ((long)RDB[msh + MESH_TYPE] == MESH_TYPE_HEXX)
        {
          x = x - (i + COS60*j)*pitch;
          y = y - j*SIN60*pitch;
        }
      else
        {
          x = x - j*SIN60*pitch;
          y = y - (i + COS60*j)*pitch;
        }

      /* Calculate half pitch */

      r = pitch*0.5;

      /* Check type and calculate minimum distance */

      if ((long)RDB[msh + MESH_TYPE] == MESH_TYPE_HEXX)
        {
          if (u != 0.0)
            {
              if (((l = -(x - r)/u) >= 0.0) && (l < min))
                min = l;
              if (((l = -(x + r)/u) >= 0.0) && (l < min))
                min = l;
            }

          if ((u - SQRT3*v) != 0.0)
            {
              if (((l = (-x + SQRT3*y + 2*r)/(u - SQRT3*v)) >= 0.0) &&
                  (l < min))
                min = l;

              if (((l = (-x + SQRT3*y - 2*r)/(u - SQRT3*v)) >= 0.0) &&
                  (l < min))
                min = l;
            }

          if ((u + SQRT3*v) != 0.0)
            {
              if (((l = (-x - SQRT3*y + 2*r)/(u + SQRT3*v)) >= 0.0) &&
                  (l < min))
                min = l;
              if (((l = (-x - SQRT3*y - 2*r)/(u + SQRT3*v)) >= 0.0) &&
                  (l < min))
                min = l;
            }
        }
      else
        {
          if (v != 0.0)
            {
              if (((l = -(y - r)/v) >= 0.0) && (l < min))
                min = l;
              if (((l = -(y + r)/v) >= 0.0) && (l < min))
                min = l;
            }

          if ((v - SQRT3*u) != 0.0)
            {
              if (((l = (-y + SQRT3*x + 2*r)/(v - SQRT3*u)) >= 0.0) &&
                  (l < min))
                min = l;

              if (((l = (-y + SQRT3*x - 2*r)/(v - SQRT3*u)) >= 0.0) &&
                  (l < min))
                min = l;
            }

          if ((v + SQRT3*u) != 0.0)
            {
              if (((l = (-y - SQRT3*x + 2*r)/(v + SQRT3*u)) >= 0.0) &&
                  (l < min))
                min = l;
              if (((l = (-y - SQRT3*x - 2*r)/(v + SQRT3*u)) >= 0.0) &&
                  (l < min))
                min = l;
            }
        }

      /* Check z-pitch */

      CheckValue(FUNCTION_NAME, "p2", "", p2, ZERO, INFTY);

      /* Calculate index */

      k = (long)((z - min2)/p2);
      CheckValue(FUNCTION_NAME, "k", "", k, 0, n2 - 1);

      /* Calculate distance to wall */

      dz = z - min2 - (double)k*p2;
      CheckValue(FUNCTION_NAME, "dz", "", dz, -1E-12, p2);

      if (w != 0.0)
        {
          if ((l = -dz/w) >= 0.0)
            {
              if (l < min)
                min = l;
            }
          else if ((l = -(dz - p2)/w) > 0.0)
            {
              if (l < min)
                min = l;
            }
          else
            Warn(FUNCTION_NAME, "dz = %E, w = %E, l = %E", dz, w, l);
        }

      /***********************************************************************/
    }
  else if (type == MESH_TYPE_CYLINDRICAL)
    {
      /***********************************************************************/

      /***** Cylindrical mesh ************************************************/

      /* Check pitches */

      CheckValue(FUNCTION_NAME, "p0", "", p0, ZERO, INFTY);
      CheckValue(FUNCTION_NAME, "p1", "", p1, ZERO, 2*PI);
      CheckValue(FUNCTION_NAME, "p2", "", p2, ZERO, INFTY);

      /* Calculate radial coordinate */

      r = sqrt(x*x + y*y);

      /* Calculate angular coordinate */

      phi = PolarAngle(x, y);
      CheckValue(FUNCTION_NAME, "phi", "", phi, 0.0, 2.0*PI);

      /* Get local normalised coordinates */

      if (max0 - min0 > 0.0)
        r = (r - min0)/(max0 - min0);
      else
        r = 0.0;

      if (max1 - min1 > 0.0)
        phi = (phi - min1)/(max1 - min1);
      else
        phi = 0.0;

      /* (invert z-co-ordinate) */

      if (max2 - min2 > 0.0)
        z0 = (z - min2)/(max2 - min2);
      else
        z0 = 0.0;

      /* Calculate indexes */

      i = (long)(r*n0);
      j = (long)(phi*n1);
      k = (long)(z0*n2);

      /* Check limits */

      CheckValue(FUNCTION_NAME, "i", "", i, 0, n0 - 1);
      CheckValue(FUNCTION_NAME, "j", "", j, 0, n1 - 1);
      CheckValue(FUNCTION_NAME, "k", "", k, 0, n2 - 1);

      /* Calculate distances to radial boundaries */

      r0 = ((double)i)*p0 + min0;
      r1 = ((double)i + 1.0)*p0 + min0;

      if (((l = CylDis(x, y, u, v, r0)) > 0.0) && (l < min))
        min = l;
      if (((l = CylDis(x, y, u, v, r1)) > 0.0) && (l < min))
        min = l;

     /* Calculate distances to segment boundaries */

      if (n1 > 1)
        {
          phi = ((double)j)*p1 + min1;
          if (((l = PhiDis(x, y, u, v, phi)) > 0.0) && (l < min))
            min = l;

          phi = ((double)j + 1.0)*p1 + min1;
          if (((l = PhiDis(x, y, u, v, phi)) > 0.0) && (l < min))
            min = l;
        }

      /* Calculate distances to axial boundaries */

      z0 = ((double)k)*p2 + min2;
      z1 = ((double)k + 1.0)*p2 + min2;

      if (((l = ZDis(z, w, z0)) > 0.0) && (l < min))
        min = l;
      if (((l = ZDis(z, w, z1)) > 0.0) && (l < min))
        min = l;

      /***********************************************************************/
    }
  else if (type == MESH_TYPE_ORTHOGONAL)
    {
      /***********************************************************************/

      /***** Unevenly-spaced orthogonal mesh *********************************/

      /* x-distance */

      if (max0 - min0 > 0.0)
        {
          /* Pointer to array */

          ptr = (long)RDB[msh + MESH_ORTHO_PTR_XLIM];

          /* Search array */

          i = SearchArray(&RDB[ptr], x, n0 + 1);
          CheckValue(FUNCTION_NAME, "i", "", i, 0, n0);

          /* Get boundaries */

          min0 = RDB[ptr + i];
          max0 = RDB[ptr + i + 1];

          /* Check position */

          CheckValue(FUNCTION_NAME, "x", "", x, min0, max0);

          /* Calculate distance to boundaries */

          if (u != 0.0)
            {
              if ((l = -(x - min0)/u) >= 0.0)
                {
                  if (l < min)
                    min = l;
                }
              else if ((l = -(x - max0)/u) >= 0.0)
                {
                  if (l < min)
                    min = l;
                }
              else
                Warn(FUNCTION_NAME, "dx = %E, u = %E, l = %E", max0 - min0,
                     u, l);
            }
         }

      /* y-distance */

      if (max1 - min1 > 0.0)
        {
          /* Pointer to array */

          ptr = (long)RDB[msh + MESH_ORTHO_PTR_YLIM];

          /* Search array */

          j = SearchArray(&RDB[ptr], y, n1 + 1);
          CheckValue(FUNCTION_NAME, "i", "", j, 0, n1);

          /* Get boundaries */

          min1 = RDB[ptr + j];
          max1 = RDB[ptr + j + 1];

          /* Check position */

          CheckValue(FUNCTION_NAME, "y", "", y, min1, max1);

          /* Calculate distance to boundaries */

          if (v != 0.0)
            {
              if ((l = -(y - min1)/v) >= 0.0)
                {
                  if (l < min)
                    min = l;
                }
              else if ((l = -(y - max1)/v) >= 0.0)
                {
                  if (l < min)
                    min = l;
                }
              else
                Warn(FUNCTION_NAME, "dy = %E, v = %E, l = %E", max1 - min1,
                     v, l);
            }
        }

      /* z-distance */

      if (max2 - min2 > 0.0)
        {
          /* Pointer to array */

          ptr = (long)RDB[msh + MESH_ORTHO_PTR_ZLIM];

          /* Search array */

          k = SearchArray(&RDB[ptr], z, n2 + 1);
          CheckValue(FUNCTION_NAME, "i", "", k, 0, n2);

          /* Get boundaries */

          min2 = RDB[ptr + k];
          max2 = RDB[ptr + k + 1];

          /* Check position */

          CheckValue(FUNCTION_NAME, "z", "", z, min2, max2);

          /* Calculate distance to boundaries */

          if (w != 0.0)
            {
              if ((l = -(z - min2)/w) >= 0.0)
                {
                  if (l < min)
                    min = l;
                }
              else if ((l = -(z - max2)/w) >= 0.0)
                {
                  if (l < min)
                    min = l;
                }
              else
                Warn(FUNCTION_NAME, "dz = %E, w = %E, l = %E", max2 - min2,
                     w, l);
            }
        }

      /***********************************************************************/
    }
  else if (type == MESH_TYPE_ICYL)
    {
      /***********************************************************************/

      /***** Unevenly-spaced cylindrical mesh ********************************/

      /* Calculate radial coordinate */

      r = sqrt(x*x + y*y);

      /* Calculate angular coordinate */

      phi = PolarAngle(x, y);
      CheckValue(FUNCTION_NAME, "phi", "", phi, 0.0, 2.0*PI);

      /* Radial distance */

      if (max0 - min0 > 0.0)
        {
          /* Pointer to array */

          ptr = (long)RDB[msh + MESH_ORTHO_PTR_XLIM];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Search array */

          i = SearchArray(&RDB[ptr], r, n0 + 1);
          CheckValue(FUNCTION_NAME, "i", "", i, 0, n0);

          /* Get boundaries */

          r0 = RDB[ptr + i];
          r1 = RDB[ptr + i + 1];

          /* Calculate distance */

          if (((l = CylDis(x, y, u, v, r0)) > 0.0) && (l < min))
            min = l;
          if (((l = CylDis(x, y, u, v, r1)) > 0.0) && (l < min))
            min = l;
         }

      /* Distance to segment boundaries */

      if (n1 > 1)
        {
          /* Pointer to array */

          ptr = (long)RDB[msh + MESH_ORTHO_PTR_YLIM];

          /* Search array */

          j = SearchArray(&RDB[ptr], phi, n1 + 1);
          CheckValue(FUNCTION_NAME, "i", "", j, 0, n1);

          /* Calculate distances */

          phi = RDB[ptr + j];
          if (((l = PhiDis(x, y, u, v, phi)) > 0.0) && (l < min))
            min = l;

          phi = RDB[ptr + j + 1];
          if (((l = PhiDis(x, y, u, v, phi)) > 0.0) && (l < min))
            min = l;
        }

      /* Axial distance */

      if (max2 - min2 > 0.0)
        {
          /* Pointer to array */

          ptr = (long)RDB[msh + MESH_ORTHO_PTR_ZLIM];

          /* Search array */

          k = SearchArray(&RDB[ptr], z, n2 + 1);
          CheckValue(FUNCTION_NAME, "i", "", k, 0, n2);

          /* Get boundaries */

          min2 = RDB[ptr + k];
          max2 = RDB[ptr + k + 1];

          /* Check position */

          CheckValue(FUNCTION_NAME, "z", "", z, min2, max2);

          /* Calculate distance to boundaries */

          if (w != 0.0)
            {
              if ((l = -(z - min2)/w) >= 0.0)
                {
                  if (l < min)
                    min = l;
                }
              else if ((l = -(z - max2)/w) >= 0.0)
                {
                  if (l < min)
                    min = l;
                }
              else
                Warn(FUNCTION_NAME, "dz = %E, w = %E, l = %E", max2 - min2,
                     w, l);
            }
        }

      /***********************************************************************/
    }
  else
    Die(FUNCTION_NAME, "Mesh type %ld not supported", type);

  /***************************************************************************/

  /* Check value */

  CheckValue(FUNCTION_NAME, "min", "", min, 0.0, INFTY);

  /* Return distance */

  return min;
}

/*****************************************************************************/
