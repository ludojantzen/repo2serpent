/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : nearestumshsurf.c                              */
/*                                                                           */
/* Created:       2013/11/23 (JLe)                                           */
/* Last modified: 2019/04/03 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Calculates distance to nearest boundary when not inside      */
/*              mesh cell.                                                   */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "NearestUMSHSurf:"

/*****************************************************************************/

double NearestUMSHSurf(long ifc, double x, double y, double z,
                       double u, double v, double w, long id)
{
  long tet, msh, lst, out;
  long k, pt0, pt1, pt2;
  double xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, l, min, d;
  double x1, y1, z1, params[9], bb[6];

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ifc)", DATA_ARRAY, ifc);

  /* Get pointer to search mesh */

  msh = (long)RDB[ifc + IFC_PTR_SEARCH_MESH_LIST];
  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

  /***************************************************************************/

  /***** Distance to outer boundaries ****************************************/

  /* Get mesh boundaries */

  xmin = RDB[msh + MESH_MIN0];
  xmax = RDB[msh + MESH_MAX0];
  ymin = RDB[msh + MESH_MIN1];
  ymax = RDB[msh + MESH_MAX1];
  zmin = RDB[msh + MESH_MIN2];
  zmax = RDB[msh + MESH_MAX2];

  /* Check that particle is in mesh */

  if ((x < xmin) || (x > xmax) || (y < ymin) || (y > ymax) || (z < zmin) ||
      (z > zmax))
    {
      /* Reset minimum distance */

      min = INFTY;

      /* Particle is not in mesh. Calculate distance to outer boundaries */

      if (u != 0.0)
        {
          if ((l = -(x - xmin)/u) > 0.0)
            if (l < min)
              min = l;

          if ((l = -(x - xmax)/u) > 0.0)
            if (l < min)
              min = l;
        }

      if (v != 0.0)
        {
          if ((l = -(y - ymin)/v) > 0.0)
            if (l < min)
              min = l;

          if ((l = -(y - ymax)/v) > 0.0)
            if (l < min)
              min = l;
        }

      if (w != 0.0)
        {
          if ((l = -(z - zmin)/w) > 0.0)
            if (l < min)
              min = l;

          if ((l = -(z - zmax)/w) > 0.0)
            if (l < min)
              min = l;
        }

      /* Check value */

      CheckValue(FUNCTION_NAME, "min", "", min, 0.0, INFTY);

      /* Do zero cut-off */

      if (min < ZERO)
        min = ZERO;

      /* Return minimum distance */

      return min;
    }

  /***************************************************************************/

  /***** Point is in search mesh *********************************************/

  /* Distance to search mesh boundaries */

  min = NearestMeshBoundary(msh, x, y, z, u, v, w, NULL);

  /* Get pointer to search mesh (NOTE: tossa luupataan */
  /* toisen kerran adaptiivisen meshin yli) */

  lst = MeshPtr(msh, x, y, z);
  CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);
  lst = (long)RDB[lst];

  /* Loop over content */

  while (lst > VALID_PTR)
    {
      /* Pointer to tet cell */

      tet = (long)RDB[lst + SEARCH_MESH_CELL_CONTENT];
      CheckPointer(FUNCTION_NAME, "(tet)", DATA_ARRAY, tet);

      /***********************************************************************/

      /***** Distance to bounding box ****************************************/

      /* Reset out flag */

      out = NO;

      /* Calculate limits based on points */

      CalculateTetBoundingBox(tet, bb);

      /* Check with bounding box */

      if ((dx = x - bb[0]) < 0.0)
        {
          if (u != 0.0)
            if (((d = -dx/u) > 0.0) && (d < min))
              {
                /* Move to position */

                x1 = x + (d + EXTRAP_L)*u;
                y1 = y + (d + EXTRAP_L)*v;
                z1 = z + (d + EXTRAP_L)*w;

                /* Check if inside */

                if ((x1 > bb[0]) && (x1 < bb[1]) && (y1 > bb[2]) && (y1 < bb[3]) &&
                    (z1 > bb[4]) && (z1 < bb[5]))
                  min = d;
              }

          out = YES;
        }
      else if ((dx = x - bb[1]) > 0.0)
        {
          if (u != 0.0)
            if (((d = -dx/u) > 0.0) && (d < min))
              {
                /* Move to position */

                x1 = x + (d + EXTRAP_L)*u;
                y1 = y + (d + EXTRAP_L)*v;
                z1 = z + (d + EXTRAP_L)*w;

                /* Check if inside */

                if ((x1 > bb[0]) && (x1 < bb[1]) && (y1 > bb[2]) && (y1 < bb[3]) &&
                    (z1 > bb[4]) && (z1 < bb[5]))
                  min = d;
              }

          out = YES;
        }

      if ((dy = y - bb[2]) < 0.0)
        {
          if (v != 0.0)
            if (((d = -dy/v) > 0.0) && (d < min))
              {
                /* Move to position */

                x1 = x + (d + EXTRAP_L)*u;
                y1 = y + (d + EXTRAP_L)*v;
                z1 = z + (d + EXTRAP_L)*w;

                /* Check if inside */

                if ((x1 > bb[0]) && (x1 < bb[1]) && (y1 > bb[2]) && (y1 < bb[3]) &&
                    (z1 > bb[4]) && (z1 < bb[5]))
                  min = d;
              }

          out = YES;
        }
      else if ((dy = y - bb[3]) > 0.0)
        {
          if (v != 0.0)
            if (((d = -dy/v) > 0.0) && (d < min))
              {
                /* Move to position */

                x1 = x + (d + EXTRAP_L)*u;
                y1 = y + (d + EXTRAP_L)*v;
                z1 = z + (d + EXTRAP_L)*w;

                /* Check if inside */

                if ((x1 > bb[0]) && (x1 < bb[1]) && (y1 > bb[2]) && (y1 < bb[3]) &&
                    (z1 > bb[4]) && (z1 < bb[5]))
                  min = d;
              }

          out = YES;
        }

      if ((dz = z - bb[4]) < 0.0)
        {
          if (w != 0.0)
            if (((d = -dz/w) > 0.0) && (d < min))
              {
                /* Move to position */

                x1 = x + (d + EXTRAP_L)*u;
                y1 = y + (d + EXTRAP_L)*v;
                z1 = z + (d + EXTRAP_L)*w;

                /* Check if inside */

                if ((x1 > bb[0]) && (x1 < bb[1]) && (y1 > bb[2]) && (y1 < bb[3]) &&
                    (z1 > bb[4]) && (z1 < bb[5]))
                  min = d;
              }

          out = YES;
        }
      else if ((dz = z - bb[5]) > 0.0)
        {
          if (w != 0.0)
            if (((d = -dz/w) > 0.0) && (d < min))
              {
                /* Move to position */

                x1 = x + (d + EXTRAP_L)*u;
                y1 = y + (d + EXTRAP_L)*v;
                z1 = z + (d + EXTRAP_L)*w;

                /* Check if inside */

                if ((x1 > bb[0]) && (x1 < bb[1]) && (y1 > bb[2]) && (y1 < bb[3]) &&
                    (z1 > bb[4]) && (z1 < bb[5]))
                  min = d;
              }

          out = YES;
        }

      /***********************************************************************/

      /***** Distance to cell faces ******************************************/

      /* Check if inside bounding box */

      if (out == NO)
        {
          /* Allocate memory for surface parameters */

          /***********************************************************/
          /* Numbering of tet faces:                                 */
          /* First face  (0,2,1) is out of cell (out of face).       */
          /* Second face (1,2,3) is forward on old face perimeter    */
          /* Third face  (0,3,2) is backward on old face perimeter   */
          /* Fourth face (0,1,3) is inside cell but not on this face */
          /***********************************************************/

          /* Get pointer to beginning of point coordinates */

          pt0 = (long)RDB[tet + TET_POINTS + 0];
          pt1 = (long)RDB[tet + TET_POINTS + 2];
          pt2 = (long)RDB[tet + TET_POINTS + 1];

          /* Store coordinates to params */
          k = 0;
          params[k++] = RDB[pt0 + 0];
          params[k++] = RDB[pt0 + 1];
          params[k++] = RDB[pt0 + 2];
          params[k++] = RDB[pt1 + 0];
          params[k++] = RDB[pt1 + 1];
          params[k++] = RDB[pt1 + 2];
          params[k++] = RDB[pt2 + 0];
          params[k++] = RDB[pt2 + 1];
          params[k++] = RDB[pt2 + 2];

          /* Get distance */

          d = SurfaceDistance(-1, params, SURF_PLANE, 9,
                              x, y, z, u, v, w, id);

          /* Compare to minimum */

          if (d < min)
            min = d;

          /* Get pointer to beginning of point coordinates */

          pt0 = (long)RDB[tet + TET_POINTS + 1];
          pt1 = (long)RDB[tet + TET_POINTS + 2];
          pt2 = (long)RDB[tet + TET_POINTS + 3];

          /* Store coordinates to params */
          k = 0;
          params[k++] = RDB[pt0 + 0];
          params[k++] = RDB[pt0 + 1];
          params[k++] = RDB[pt0 + 2];
          params[k++] = RDB[pt1 + 0];
          params[k++] = RDB[pt1 + 1];
          params[k++] = RDB[pt1 + 2];
          params[k++] = RDB[pt2 + 0];
          params[k++] = RDB[pt2 + 1];
          params[k++] = RDB[pt2 + 2];

          /* Get distance */

          d = SurfaceDistance(-1, params, SURF_PLANE, 9,
                              x, y, z, u, v, w, id);

          /* Compare to minimum */

          if (d < min)
            min = d;

          /* Get pointer to beginning of point coordinates */

          pt0 = (long)RDB[tet + TET_POINTS + 0];
          pt1 = (long)RDB[tet + TET_POINTS + 3];
          pt2 = (long)RDB[tet + TET_POINTS + 2];

          /* Store coordinates to params */
          k = 0;
          params[k++] = RDB[pt0 + 0];
          params[k++] = RDB[pt0 + 1];
          params[k++] = RDB[pt0 + 2];
          params[k++] = RDB[pt1 + 0];
          params[k++] = RDB[pt1 + 1];
          params[k++] = RDB[pt1 + 2];
          params[k++] = RDB[pt2 + 0];
          params[k++] = RDB[pt2 + 1];
          params[k++] = RDB[pt2 + 2];

          /* Get distance */

          d = SurfaceDistance(-1, params, SURF_PLANE, 9,
                              x, y, z, u, v, w, id);

          /* Compare to minimum */

          if (d < min)
            min = d;

          /* Get pointer to beginning of point coordinates */

          pt0 = (long)RDB[tet + TET_POINTS + 0];
          pt1 = (long)RDB[tet + TET_POINTS + 1];
          pt2 = (long)RDB[tet + TET_POINTS + 3];

          /* Store coordinates to params */
          k = 0;
          params[k++] = RDB[pt0 + 0];
          params[k++] = RDB[pt0 + 1];
          params[k++] = RDB[pt0 + 2];
          params[k++] = RDB[pt1 + 0];
          params[k++] = RDB[pt1 + 1];
          params[k++] = RDB[pt1 + 2];
          params[k++] = RDB[pt2 + 0];
          params[k++] = RDB[pt2 + 1];
          params[k++] = RDB[pt2 + 2];

          /* Get distance */

          d = SurfaceDistance(-1, params, SURF_PLANE, 9,
                              x, y, z, u, v, w, id);

          /* Compare to minimum */

          if (d < min)
            min = d;

        }

      /* Next */

      lst = NextItem(lst);
    }

  /* Check minimum distance */

  CheckValue(FUNCTION_NAME, "min", "", min, 0.0, INFTY);

  /* Do zero cut-off */

  if (min < ZERO)
    min = ZERO;

  /* Return minimum */

  return min;
}

/*****************************************************************************/
