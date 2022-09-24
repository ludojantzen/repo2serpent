/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : wwdis.c                                        */
/*                                                                           */
/* Created:       2015/10/02 (JLe)                                           */
/* Last modified: 2019/12/12 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Calculates distance to weight window boundary                */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "WWDis:"

/*****************************************************************************/

double WWDis(long type, double x, double y, double z, double u, double v,
             double w)
{
  long loc0, ptr, msh, fail;
  double min, l;

  /* Check for weight windows and response matrix calculation */

  if (((long)RDB[DATA_USE_WEIGHT_WINDOWS] == NO) &&
      ((long)RDB[DATA_RMTX_CALC] == NO))
    Die(FUNCTION_NAME, "Weight windows not in use");

  /* Reset minimum */

  min = INFTY;

  /***************************************************************************/

  /***** Weight window mesh **************************************************/

  /* Loop over structures */

  loc0 = (long)RDB[DATA_PTR_WWD0];
  while (loc0 > VALID_PTR)
    {
      /* Check particle type */

      if ((long)RDB[loc0 + WWD_PARTICLE_TYPE] != type)
        {
          /* Pointer to next */

          loc0 = NextItem(loc0);

          /* Cycle loop */

          continue;
        }

      /* Pointer to mesh */

      if ((msh = (long)RDB[loc0 + WWD_PTR_MESH]) < VALID_PTR)
        {
          /* Pointer to next */

          loc0 = NextItem(loc0);

          /* Cycle loop */

          continue;
        }

      /* Check if iteration mode */

      if ((long)RDB[loc0 + WWD_TYPE] == WWD_MESH_TYPE_ITER)
        {
          /* Pointer to response matrix solver */

          ptr = (long)RDB[DATA_PTR_RMX0];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Check mesh */

          if (msh == (long)RDB[ptr + RMX_PTR_MESH])
            break;
        }

      /* Get distance */

      if ((l = NearestMeshBoundary(msh, x, y, z, u, v, w, &fail)) > 0.0)
        if (l < min)
          min = l;

      /* Pointer to next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Response matrix mesh ************************************************/

  /* Loop over structures */

  loc0 = (long)RDB[DATA_PTR_RMX0];
  while (loc0 > VALID_PTR)
    {
      /* Pointer to mesh */

      if ((msh = (long)RDB[loc0 + RMX_PTR_MESH]) < VALID_PTR)
        {
          /* Pointer to next */

          loc0 = NextItem(loc0);

          /* Cycle loop */

          continue;
        }

      /* Check particle type */

      if ((long)RDB[loc0 + RMX_PARTICLE_TYPE] != type)
        {
          /* Pointer to next */

          loc0 = NextItem(loc0);

          /* Cycle loop */

          continue;
        }

      /* Get distance */

      if ((l = NearestMeshBoundary(msh, x, y, z, u, v, w, &fail)) > 0.0)
        if (l < min)
          min = l;

      /* Pointer to next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /* Return distance */

  return min;
}

/*****************************************************************************/
