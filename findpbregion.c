/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : findpbregion.c                                 */
/*                                                                           */
/* Created:       2010/10/23 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Finds neutron location in explicit stochastic geometry       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FindPBRegion:"

/*****************************************************************************/

long FindPBRegion(long uni0, long pbd, double *x, double *y, double *z, 
                  long *pbl0, long *ridx, long id)
{
  long ptr, lst, pbl, uni, ncol, msh;
  double rp, dx, dy, dz;

  /* Check pointers */

  CheckPointer(FUNCTION_NAME, "(uni0)", DATA_ARRAY, uni0);
  CheckPointer(FUNCTION_NAME, "(pbd)", DATA_ARRAY, pbd);

  /* Check plotter mode */

  if (((long)RDB[DATA_PLOTTER_MODE] == NO) || 
      ((long)RDB[DATA_QUICK_PLOT_MODE] == YES))
    {
      /* Check previous */

      ptr = (long)RDB[uni0 + UNIVERSE_PTR_PREV_REG];
      if ((pbl = (long)GetPrivateData(ptr, id)) > VALID_PTR)
        {
          /* Get parameters */
          
          dx = *x - RDB[pbl + PEBBLE_X0];
          dy = *y - RDB[pbl + PEBBLE_Y0];
          dz = *z - RDB[pbl + PEBBLE_Z0];
          
          rp = RDB[pbl + PEBBLE_RAD];
          
          /* Check if particle is inside */
          
          if (dx*dx + dy*dy + dz*dz < rp*rp)
            {
              /* Co-ordinate transformation */
              
              *x = dx;
              *y = dy;
              *z = dz;
              
              /* Get pointer to universe */
              
              uni = (long)RDB[pbl + PEBBLE_PTR_UNIV];
              CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);
              
              /* Put pebble pointer */
              
              *pbl0 = pbl;
              
              /* Get collision number */
              
              ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
              CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
              
              ncol = (long)GetPrivateData(ptr, id);
              
              /* Store index */
              
              StoreValuePair(pbd + PBED_PTR_COL_PEBBLE, (double)ncol, 
                             (double)pbl, id);
              
              /* Put region index */
              
              *ridx = (long)RDB[pbl + PEBBLE_IDX] + 1;
              
              /* Return universe pointer */
              
              return uni;
            }
        }
    }

  /* Reset pebble pointer */

  *pbl0 = -1;

  /* Pointer to search mesh */

  msh = (long)RDB[pbd + PBED_PTR_SEARCH_MESH];
  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

  /* Get pointer to search mesh */
  
  if ((lst = MeshPtr(msh, *x, *y, *z)) > VALID_PTR)
    lst = (long)RDB[lst];
  
  /* Loop over content */

  while (lst > VALID_PTR)
    {
      /* Pointer to pebble */
      
      pbl = (long)RDB[lst + SEARCH_MESH_CELL_CONTENT];
      CheckPointer(FUNCTION_NAME, "(pbl)", DATA_ARRAY, pbl);

      /* Get parameters */

      dx = *x - RDB[pbl + PEBBLE_X0];
      dy = *y - RDB[pbl + PEBBLE_Y0];
      dz = *z - RDB[pbl + PEBBLE_Z0];

      rp = RDB[pbl + PEBBLE_RAD];

      /* Check if particle is inside */

      if (dx*dx + dy*dy + dz*dz < rp*rp)
        {
          /* Co-ordinate transformation */

          *x = dx;
          *y = dy;
          *z = dz;

          /* Get pointer to universe */

          uni = (long)RDB[pbl + PEBBLE_PTR_UNIV];
          CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

          /* Put pebble pointer */

          *pbl0 = pbl;

          /* Get collision number */

          ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
          CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);

          ncol = (long)GetPrivateData(ptr, id);

          /* Store index */

          StoreValuePair(pbd + PBED_PTR_COL_PEBBLE, (double)ncol, 
                         (double)pbl, id);
          
          /* Put region index */

          *ridx = (long)RDB[pbl + PEBBLE_IDX] + 1;

          /* Check plotter mode */

          if ((long)RDB[DATA_PLOTTER_MODE] == NO)
            {
              /* Pointer to counter */
              
              ptr = (long)RDB[lst + SEARCH_MESH_PTR_CELL_COUNT];
              CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
              
              /* Add counter */
              
              AddPrivateData(ptr, 1, id);
            }

          /* Put previous pointer */
          
          ptr = (long)RDB[uni0 + UNIVERSE_PTR_PREV_REG];
          PutPrivateData(ptr, pbl, id);

          /* Return universe pointer */

          return uni;
             }
      
      /* Next */

      lst = NextItem(lst);
    }

  /* Get pointer to background universe */

  uni = (long)RDB[pbd + PBED_PTR_BG_UNIV];
  CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

  /* Put region index */

  *ridx = 0;

  /* Return pointer */
  
  return uni;
}

/*****************************************************************************/
