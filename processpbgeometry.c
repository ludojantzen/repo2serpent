/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processbgeometry.c                             */
/*                                                                           */
/* Created:       2010/11/09 (JLe)                                           */
/* Last modified: 2012/02/18 (JLe)                                           */
/* Version:       2.1.3                                                      */
/*                                                                           */
/* Description: Processes explicit stochastic geometry                       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessPBGeometry:"

/*****************************************************************************/

void ProcessPBGeometry()
{  
  long loc0, pbl, msh, pid, nb, ptr;
  double x, y, z, r;

  /* Check pointer */

  if ((loc0 = (long)RDB[DATA_PTR_PB0]) < 0)
    return;
  
  /* Loop over definitions */

  while (loc0 > VALID_PTR)
    {
      /* Pointer to search mesh */

      msh = (long)RDB[loc0 + PBED_PTR_SEARCH_MESH];
      CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

      /* Reset pebble index */

      pid = 0;
      
      /* Allocate memory for collision pointer */

      AllocValuePair(loc0 + PBED_PTR_COL_PEBBLE);

      /* Loop over pebbles */
      
      pbl = (long)RDB[loc0 + PBED_PTR_PEBBLES];      
      while (pbl > VALID_PTR)
        {
          /* Put pebble index */

          WDB[pbl + PEBBLE_IDX] = (double)(pid++);

          /* Get co-ordinates */
          
          x = RDB[pbl + PEBBLE_X0];
          y = RDB[pbl + PEBBLE_Y0];
          z = RDB[pbl + PEBBLE_Z0];

          /* Get radius */

          r = RDB[pbl + PEBBLE_RAD];

          /* Add pebble in mesh */

          AddSearchMesh(msh, pbl, x - r, x + r, y - r, y + r, z - r, z + r);
          
          /* Next pebble */
          
          pbl = NextItem(pbl);
        }
      
      /***********************************************************************/

      /***** Allocate memory for power distribution '*************************/

      /* Check if results are requested */

      if ((long)RDB[loc0 + PBED_CALC_RESULTS] == YES)
        {
          /* Get number of pebbles */

          nb = (long)RDB[loc0 + PBED_N_PEBBLES];

          /* Allocate memory */

          ptr = NewStat("PBED_POW", 1, nb); 
          WDB[loc0 + PBED_PTR_POW] = (double)ptr;
        }
      
      /***********************************************************************/

      /* Next geometry */

      loc0 = NextItem(loc0);
    }
}

/*****************************************************************************/
