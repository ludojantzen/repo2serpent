/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : nearestpbsurf.c                                */
/*                                                                           */
/* Created:       2010/10/25 (JLe)                                           */
/* Last modified: 2017/09/23 (JLe)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Calculates distance to nearest pebble / particle             */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "NearestPBSurf:"

/*****************************************************************************/

double NearestPBSurf(long pbd, double x, double y, double z, 
                     double u, double v, double w, long id)
{
  long nx, ny, nz, i, j, k, pbl, msh, lst;
  double xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, param[4], l, min;
  double px, py, pz;

  /* Pointer to search mesh */

  msh = (long)RDB[pbd + PBED_PTR_SEARCH_MESH];
  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

  /* Get mesh boundaries */

  xmin = RDB[msh + MESH_MIN0];
  xmax = RDB[msh + MESH_MAX0];
  ymin = RDB[msh + MESH_MIN1];
  ymax = RDB[msh + MESH_MAX1];
  zmin = RDB[msh + MESH_MIN2];
  zmax = RDB[msh + MESH_MAX2];

  /* Reset minimum distance */

  min = INFTY;

  /* Check that neutron is in mesh */

  if ((x < xmin) || (x > xmax) || (y < ymin) || (y > ymax) || (z < zmin) ||
      (z > zmax))
    {
      /* Neutron is not in mesh. Calculate distance to mesh boundaries */

      if ((l = -(x - xmin)/u) > 0.0)
        if (l < min)
          min = l;

      if ((l = -(x - xmax)/u) > 0.0)
        if (l < min)
          min = l;

      if ((l = -(y - ymin)/v) > 0.0)
        if (l < min)
          min = l;

      if ((l = -(y - ymax)/v) > 0.0)
        if (l < min)
          min = l;
      
      if ((l = -(z - zmin)/w) > 0.0)
        if (l < min)
          min = l;

      if ((l = -(z - zmax)/w) > 0.0)
        if (l < min)
          min = l;

      /* Do zero cut-off */

      if (min < ZERO)
        min = ZERO;

      /* Return minimum distance */

      return min;
    }

  /* Get mesh size */

  nx = (long)RDB[msh + MESH_N0];
  ny = (long)RDB[msh + MESH_N1];
  nz = (long)RDB[msh + MESH_N2];

  /* Calculate mesh pitch */

  px = (xmax - xmin)/((double)nx);
  py = (ymax - ymin)/((double)ny);
  pz = (zmax - zmin)/((double)nz);

  /* Calculate mesh indexes */

  i = (long)((x - xmin)/px);
  j = (long)((y - ymin)/py);
  k = (long)((z - zmin)/pz);

  /* Check limits */

  CheckValue(FUNCTION_NAME, "i", "", i, 0, nx - 1);
  CheckValue(FUNCTION_NAME, "j", "", j, 0, ny - 1);
  CheckValue(FUNCTION_NAME, "k", "", k, 0, nz - 1);

  /* Calculate distance to mesh walls */

  dx = x - xmin - (double)i*px;
  dy = y - ymin - (double)j*py;
  dz = z - zmin - (double)k*pz;

  /* Check values */

  CheckValue(FUNCTION_NAME, "dx", "", dx, 0.0, px);
  CheckValue(FUNCTION_NAME, "dy", "", dy, 0.0, py);
  CheckValue(FUNCTION_NAME, "dz", "", dz, 0.0, pz);

  /* Calculate distance to mesh cell walls */

  if ((l = -dx/u) > 0)
    {
      if (l < min)
        min = l;
    }
  else if ((l = -(dx - px)/u) > 0)
    {
      if (l < min)
        min = l;
    }
  else
    Warn(FUNCTION_NAME, "l = %E, u = %E", l, u);

  if ((l = -dy/v) > 0)
    {
      if (l < min)
        min = l;
    }
  else if ((l = -(dy - py)/v) > 0)
    {
      if (l < min)
        min = l;
    }
  else
    Warn(FUNCTION_NAME, "l = %E, v = %E", l, v);

  if ((l = -dz/w) > 0)
    {
      if (l < min)
        min = l;
    }
  else if ((l = -(dz - pz)/w) > 0)
    {
      if (l < min)
        min = l;
    }
  else
    Warn(FUNCTION_NAME, "l = %E, w = %E", l, w);
  
   /* Get pointer to search mesh */

  lst = MeshPtr(msh, x, y, z);
  CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);
  lst = (long)RDB[lst];

  /* Loop over content */

  while (lst > VALID_PTR)
    {
      /* Pointer to pebble */
      
      pbl = (long)RDB[lst + SEARCH_MESH_CELL_CONTENT];
      CheckPointer(FUNCTION_NAME, "(pbl)", DATA_ARRAY, pbl);

      /* Set surface parameters */

      param[0] = RDB[pbl + PEBBLE_X0];
      param[1] = RDB[pbl + PEBBLE_Y0];
      param[2] = RDB[pbl + PEBBLE_Z0];
      param[3] = RDB[pbl + PEBBLE_RAD];

      /* Get distance */

      l = SurfaceDistance(-1, param, SURF_SPH, 4, x, y, z, u, v, w, id);

      /* Compare to minimum */
        
      if (l < min)
        min = l;
        
      /* Next */

      lst = NextItem(lst);
    }

  /* Do zero cut-off */

  if (min < ZERO)
    min = ZERO;

  /* Return minimum */

  return min;
}

/*****************************************************************************/
