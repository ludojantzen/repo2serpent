/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : virgcucolflags.c                               */
/*                                                                           */
/* Created:       2018/06/24 (JLe)                                           */
/* Last modified: 2018/06/24 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Sets collision flags for virtual universes used in group     */
/*              constant generation                                          */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "VirtGCUColFlags:"

/*****************************************************************************/

void VirtGCUColFlags(double x, double y, double z, long id)
{
  long gcu, uni, ptr, ncol;

  /* Get collision number */

  ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
  ncol = (long)GetPrivateData(ptr, id);

  /* Loop over universes */

  gcu = (long)RDB[DATA_PTR_GCU0];
  while (gcu > VALID_PTR)
    {
      /* Pointer to universe */
      
      uni = (long)RDB[gcu + GCU_PTR_UNIV];
      CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

      /* Check universe type, previous collision and test point */
  
      if ((long)RDB[uni + UNIVERSE_TYPE] == UNIVERSE_TYPE_SUPER)
        if (TestValuePair(uni + UNIVERSE_COL_COUNT, (double)ncol, id) < 0.0)
          if (InSuperCell(uni, -1, x, y, z, id) == YES)
            {
              /* Put collision flag */ 
              
              StoreValuePair(uni + UNIVERSE_COL_COUNT, (double)ncol, 1.0, id);
              
              /* Put gcu pointer */

              if ((long)RDB[DATA_MULTI_LEVEL_GCU] == NO)
                StoreValuePair(DATA_GCU_PTR_UNI, (double)ncol, (double)gcu, id);
              
              /* Break loop */
              
              break;
            }
      
      /* Pointer to next */
      
      gcu = NextItem(gcu);
    }
}

/*****************************************************************************/
