/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : stlfacetdistance.c                             */
/*                                                                           */
/* Created:       2014/11/24 (JLe)                                           */
/* Last modified: 2016/06/15 (JLe)                                           */
/* Version:       2.1.27                                                     */
/*                                                                           */
/* Description:  Calculates distance to a triangular STL facet               */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "STLFacetDistance:"

/*****************************************************************************/

double STLFacetDistance(long loc0, double x, double y, double z, double u,
                        double v, double w, long tst, long id)
{
  long ptr;
  double params[9], l, u0, v0, w0, u1, v1, w1, u2, v2, w2, x0, y0, z0;
  double dot00, dot01, dot02, dot11, dot12, div, ex;

  /* Read points to array */

  ptr = (long)RDB[loc0 + STL_FACET_PTR_PT1];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  params[0] = RDB[ptr + STL_POINT_X];
  params[1] = RDB[ptr + STL_POINT_Y];
  params[2] = RDB[ptr + STL_POINT_Z];
  
  ptr = (long)RDB[loc0 + STL_FACET_PTR_PT2];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  
  params[3] = RDB[ptr + STL_POINT_X];
  params[4] = RDB[ptr + STL_POINT_Y];
  params[5] = RDB[ptr + STL_POINT_Z];
  
  ptr = (long)RDB[loc0 + STL_FACET_PTR_PT3];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  
  params[6] = RDB[ptr + STL_POINT_X];
  params[7] = RDB[ptr + STL_POINT_Y];
  params[8] = RDB[ptr + STL_POINT_Z];
  
  /* Calculate distance */
      
  l = SurfaceDistance(-1, params, SURF_PLANE, 9, x, y, z, u, v, w, id);
  CheckValue(FUNCTION_NAME, "l", "", l, 0.0, INFTY);

  /* Move point to surfce */

  x = x + l*u;
  y = y + l*v;
  z = z + l*w;

  /* Check with bounding box */

  if ((x < RDB[loc0 + STL_FACET_XMIN]) || (x > RDB[loc0 + STL_FACET_XMAX]) ||
      (y < RDB[loc0 + STL_FACET_YMIN]) || (y > RDB[loc0 + STL_FACET_YMAX]) ||
      (z < RDB[loc0 + STL_FACET_ZMIN]) || (z > RDB[loc0 + STL_FACET_ZMAX]))
    return INFTY;

  /* Calculate vectors */
  
  u0 = params[6] - params[0];
  v0 = params[7] - params[1];
  w0 = params[8] - params[2];
  
  u1 = params[3] - params[0];
  v1 = params[4] - params[1];
  w1 = params[5] - params[2];
  
  u2 = x - params[0];
  v2 = y - params[1];
  w2 = z - params[2];
  
  /* Calculate dot products */
  
  dot00 = u0*u0 + v0*v0 + w0*w0;
  dot01 = u0*u1 + v0*v1 + w0*w1;
  dot02 = u0*u2 + v0*v2 + w0*w2;
  dot11 = u1*u1 + v1*v1 + w1*w1;
  dot12 = u1*u2 + v1*v2 + w1*w2;

  /* Check values */

  CheckValue(FUNCTION_NAME, "dot00", "", dot00, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "dot01", "", dot01, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "dot02", "", dot02, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "dot11", "", dot11, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "dot12", "", dot12, -INFTY, INFTY);
  
  /* Calculate barycentric coordinates */

  if ((div = dot00*dot11 - dot01*dot01) == 0.0)
    Die(FUNCTION_NAME, "Division by zero");

  x0 = (dot11*dot02 - dot01*dot12)/div;
  y0 = (dot00*dot12 - dot01*dot02)/div;
  z0 = 1.0 - x0 - y0;

  /* Check coordinates */
  
  CheckValue(FUNCTION_NAME, "x0", "", x0, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "y0", "", y0, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "z0", "", z0, -INFTY, INFTY);
  
  /* Exclusion distance */

  ex = RDB[DATA_STL_FACET_EXD];

  /* Check if outside the exclusion zone */

  if ((x0 < -ex) || (y0 < -ex) || (z0 < -ex))
    {
      /* Point is safely outside */

      return INFTY;
    }
  else if ((x0 > ex) && (y0 > ex) && (z0 > ex))
    {
      /* Point is safely inside */
      
      return l;
    }

  /* Point is close to edge, return -1 to discard ray if test mode, */
  /* accept if calculating nearest distance (15.6.2016 / 2.1.27).   */

  if (tst == YES)
    return -1.0;
  else
    return l;
}

/*****************************************************************************/
