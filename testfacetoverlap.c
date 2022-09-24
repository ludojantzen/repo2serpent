/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : testfacetoverlap.c                             */
/*                                                                           */
/* Created:       2018/11/21 (JLe)                                           */
/* Last modified: 2018/11/21 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Tests if STL facet overlaps a search mesh cell               */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "TestFacetOverlap:"

/*****************************************************************************/

long TestFacetOverlap(long loc0, double xmin, double xmax, double ymin, 
                      double ymax, double zmin, double zmax) 
{
  long ptr, n;
  double x, y, z, x1, y1, z1, x2, y2, z2, x3, y3, z3, u, v, w, dx, dy, dz, d;
  double x0, y0, z0;

  /***************************************************************************/

  /***** Check if any of the corner points ar inside the cell ****************/

  /* Read first point */

  ptr = (long)RDB[loc0 + STL_FACET_PTR_PT1];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  x1 = RDB[ptr + STL_POINT_X];
  y1 = RDB[ptr + STL_POINT_Y];
  z1 = RDB[ptr + STL_POINT_Z];
  
  /* Check */

  if ((x1 > xmin) && (x1 < xmax) && (y1 > ymin) && (y1 < ymax) && 
      (z1 > zmin) && (z1 < zmax))
    return YES;

  /* Read second point */

  ptr = (long)RDB[loc0 + STL_FACET_PTR_PT2];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  x2 = RDB[ptr + STL_POINT_X];
  y2 = RDB[ptr + STL_POINT_Y];
  z2 = RDB[ptr + STL_POINT_Z];
  
  /* Check */

  if ((x2 > xmin) && (x2 < xmax) && (y2 > ymin) && (y2 < ymax) && 
      (z2 > zmin) && (z2 < zmax))
    return YES;

  /* Read third point */

  ptr = (long)RDB[loc0 + STL_FACET_PTR_PT3];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  x3 = RDB[ptr + STL_POINT_X];
  y3 = RDB[ptr + STL_POINT_Y];
  z3 = RDB[ptr + STL_POINT_Z];
  
  /* Check */

  if ((x3 > xmin) && (x3 < xmax) && (y3 > ymin) && (y3 < ymax) && 
      (z3 > zmin) && (z3 < zmax))
    return YES;

  /***************************************************************************/

  /***** Check if facet edges pass through planes ****************************/

  /* Loop over points */

  for (n = 0; n < 3; n++)
    {
      /* Put coordinates */

      if (n == 0)
        {
          x = x1;
          y = y1;
          z = z1;

          u = x2 - x;
          v = y2 - y;
          w = z2 - z;
        }
      else if (n == 1)
        {
          x = x2;
          y = y2;
          z = z2;

          u = x3 - x;
          v = y3 - y;
          w = z3 - z;
        }
      else
        {
          x = x3;
          y = y3;
          z = z3;

          u = x1 - x;
          v = y1 - y;
          w = z1 - z;
        }

      /* Normalize direction cosines */

      d = sqrt(u*u + v*v + w*w);
      CheckValue(FUNCTION_NAME, "d", "", d, ZERO, INFTY);

      u = u/d;
      v = v/d;
      w = w/d;

      /* Check intersections with faces */

      d = -(x - xmin)/u;

      y0 = y + v*d;
      z0 = z + w*d;
      
      if ((y0 > ymin) && (y0 < ymax) && (z0 > zmin) && (z0 < zmax))
        return YES;

      d = -(x - xmax)/u;

      y0 = y + v*d;
      z0 = z + w*d;
      
      if ((y0 > ymin) && (y0 < ymax) && (z0 > zmin) && (z0 < zmax))
        return YES;

      d = -(y - ymin)/v;

      x0 = x + u*d;
      z0 = z + w*d;
      
      if ((x0 > xmin) && (x0 < xmax) && (z0 > zmin) && (z0 < zmax))
        return YES;

      d = -(y - ymax)/v;

      x0 = x + u*d;
      z0 = z + w*d;
      
      if ((x0 > xmin) && (x0 < xmax) && (z0 > zmin) && (z0 < zmax))
        return YES;

      d = -(z - zmin)/w;

      x0 = x + u*d;
      y0 = y + v*d;
      
      if ((y0 > ymin) && (y0 < ymax) && (x0 > xmin) && (x0 < xmax))
        return YES;

      d = -(z - zmax)/w;

      x0 = x + u*d;
      y0 = y + v*d;
      
      if ((y0 > ymin) && (y0 < ymax) && (x0 > xmin) && (x0 < xmax))
        return YES;
    }

  /***************************************************************************/

  /***** Check intersections with edges **************************************/

  /* Facet lenghts */

  dx = xmax - xmin;
  CheckValue(FUNCTION_NAME, "dx", "", dx, ZERO, INFTY);

  dy = ymax - ymin;
  CheckValue(FUNCTION_NAME, "dy", "", dy, ZERO, INFTY);

  dz = zmax - zmin;
  CheckValue(FUNCTION_NAME, "dz", "", dz, ZERO, INFTY);

  /* Check facets from each corner */

  if (STLFacetDistance(loc0, xmin, ymin, zmin, 1.0, 0.0, 0.0, NO, 0) < dx)
    return YES;
  else if (STLFacetDistance(loc0, xmin, ymin, zmin, 0.0, 1.0, 0.0, NO, 0) < dy)
    return YES;
  else if (STLFacetDistance(loc0, xmin, ymin, zmin, 0.0, 0.0, 1.0, NO, 0) < dz)
    return YES;
  else if (STLFacetDistance(loc0, xmax, ymin, zmin, 0.0, 1.0, 0.0, NO, 0) < dy)
    return YES;
  else if (STLFacetDistance(loc0, xmax, ymin, zmin, 0.0, 0.0, 1.0, NO, 0) < dz)
    return YES;
  else if (STLFacetDistance(loc0, xmin, ymax, zmin, 1.0, 0.0, 0.0, NO, 0) < dx)
    return YES;
  else if (STLFacetDistance(loc0, xmin, ymax, zmin, 0.0, 0.0, 1.0, NO, 0) < dz)
    return YES;
  else if (STLFacetDistance(loc0, xmin, ymin, zmax, 1.0, 0.0, 0.0, NO, 0) < dx)
    return YES;
  else if (STLFacetDistance(loc0, xmin, ymin, zmax, 0.0, 1.0, 0.0, NO, 0) < dy)
    return YES;
  else if (STLFacetDistance(loc0, xmax, ymax, zmin, 0.0, 0.0, 1.0, NO, 0) < dz)
    return YES;
  else if (STLFacetDistance(loc0, xmin, ymax, zmax, 1.0, 0.0, 0.0, NO, 0) < dx)
    return YES;
  else if (STLFacetDistance(loc0, xmax, ymin, zmax, 0.0, 1.0, 0.0, NO, 0) < dy)
    return YES;

  /****************************************************************************/

  /* No overlap */

  return NO;

  /***************************************************************************/
}

/*****************************************************************************/
