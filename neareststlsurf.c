/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : neareststlsurf.c                               */
/*                                                                           */
/* Created:       2014/03/05 (JLe)                                           */
/* Last modified: 2016/06/15 (JLe)                                           */
/* Version:       2.1.27                                                     */
/*                                                                           */
/* Description:  Calculates minimum distance to STL surface                  */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "NearestSTLSurf:"

/*****************************************************************************/

double NearestSTLSurf(long stl, double x, double y, double z, double u, 
                      double v, double w, long id)
{
  long msh, loc0, loc1, ptr;
  double xmin, xmax, ymin, ymax, zmin, zmax, l, min;
  double params[9];

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(stl)", DATA_ARRAY, stl);

  /* Get pointer to search mesh */

  msh = (long)RDB[stl + STL_PTR_FACET_MESH];
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

  if ((x < xmin) || (x > xmax) || 
      (y < ymin) || (y > ymax) || 
      (z < zmin) || (z > zmax))
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

      /* Return minimum distance */

      return min;
    }

  /***************************************************************************/

  /***** Point is in search mesh *********************************************/

  /* Distance to search mesh boundaries */

  min = NearestMeshBoundary(msh, x, y, z, u, v, w, NULL);
  CheckValue(FUNCTION_NAME, "min", "", min, ZERO, INFTY);

  /* Get pointer to search mesh cell */

  loc0 = MeshPtr(msh, x, y, z);
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Loop over content */

  loc0 = (long)RDB[loc0];
  while (loc0 > VALID_PTR)
    {
      /* Pointer to content */
      
      if ((loc1 = (long)RDB[loc0 + SEARCH_MESH_CELL_CONTENT]) < VALID_PTR)
        break;

      /* Check if DT was enforced */

      if ((long)RDB[DATA_PLOTTER_MODE] == NO)
        if ((long)RDB[DATA_STL_ENFORCE_DT] == YES)
          Die(FUNCTION_NAME, "Shouldn't be here");

      /* Mode */

      if (1 != 2)
        {
          /* Get distance */

          l = STLFacetDistance(loc1, x, y, z, u, v, w, NO, id);
          CheckValue(FUNCTION_NAME, "l", "", l, 0.0, INFTY);

          /* Tää ei varmaan pelaa */
          /*
          if ((l > 0.0) && (l < INFTY))
            {
              if (l < min)
                min = l;
              break;
            }
          */
        }
      else
        {
          /* Read points to array */
      
          ptr = (long)RDB[loc1 + STL_FACET_PTR_PT1];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          
          params[0] = RDB[ptr + STL_POINT_X];
          params[1] = RDB[ptr + STL_POINT_Y];
          params[2] = RDB[ptr + STL_POINT_Z];
          
          ptr = (long)RDB[loc1 + STL_FACET_PTR_PT2];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          
          params[3] = RDB[ptr + STL_POINT_X];
          params[4] = RDB[ptr + STL_POINT_Y];
          params[5] = RDB[ptr + STL_POINT_Z];
          
          ptr = (long)RDB[loc1 + STL_FACET_PTR_PT3];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          
          params[6] = RDB[ptr + STL_POINT_X];
          params[7] = RDB[ptr + STL_POINT_Y];
          params[8] = RDB[ptr + STL_POINT_Z];
          
          /* Calculate distance */
          
          l = SurfaceDistance(-1, params, SURF_PLANE, 9, x, y, z, u, v, w, id);
          CheckValue(FUNCTION_NAME, "l", "", l, 0.0, INFTY);
        }

      /* Compare to minimum */
      
      if (l < min)
        min = l;
      
      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Check value */

  CheckValue(FUNCTION_NAME, "min", "", min, 0.0, INFTY);

  /* Return minimum distance */

  return min;
}

/*****************************************************************************/
