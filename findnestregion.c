/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : findnestregion.c                               */
/*                                                                           */
/* Created:       2010/10/07 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Finds nest region based on coordinates                       */
/*                                                                           */
/* Comments: - Ton edellisen alueen testaamisen hyöty on vähän kyseenalainen */
/*             kun nää alueet on kuitenkin yleensä suht pieniä verrattuna    */
/*             mfp.                                                          */
/*                                                                           */
/*           - 2.1.19 Added time dependence to CoordExpans                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FindNestRegion:"

/*****************************************************************************/

long FindNestRegion(long uni0, long nst, double x, double y, double z, long id)
{
  long ptr, n, reg, surf;
  double t;
  
  /* Check pointers */

  CheckPointer(FUNCTION_NAME, "(uni0)", DATA_ARRAY, uni0);
  CheckPointer(FUNCTION_NAME, "(nst)", DATA_ARRAY, nst);

  /* Get time */

  ptr = (long)RDB[uni0 + UNIVERSE_PTR_PRIVA_T];
  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
  t = GetPrivateData(ptr, id);

  /* Check for performance interface */

  if ((ptr = (long)RDB[uni0 + UNIVERSE_PTR_IFC_FUEP]) > VALID_PTR)
    {

      /* Change to cold state */

      CoordExpans(ptr, &x, &y, &z, t, 1);

    }

  /***************************************************************************/

  /***** Regular nest ********************************************************/

  /* Check previous */

  ptr = (long)RDB[uni0 + UNIVERSE_PTR_PREV_REG];
  if ((reg = (long)GetPrivateData(ptr, id)) > VALID_PTR)
    {
      /* Get pointer to first surface */

      if ((surf = (long)RDB[reg + NEST_REG_PTR_SURF_IN]) > VALID_PTR)
        {
          /* Test surface */

          if (TestSurface(surf, x, y, z, NO, id) == YES)
            {
              /* Get pointer to second surface */

              if ((surf = (long)RDB[reg + NEST_REG_PTR_SURF_OUT]) > VALID_PTR)
                {
                  /* Test surface */
                  
                  if (TestSurface(surf, x, y, z, NO, id) == NO)
                    return reg;
                }
              else
                {
                  /* Innermost region, point is in */

                  return reg;
                }
            }
        }
    }
        
  /* Get pointer to regions */
              
  ptr = (long)RDB[nst + NEST_PTR_REGIONS];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  
  /* Loop over regions */ 

  n = 0;
  while ((reg = ListPtr(ptr, n++)) > VALID_PTR)
    {
      /* Get pointer to surface */

      if ((surf = (long)RDB[reg + NEST_REG_PTR_SURF_IN]) < VALID_PTR)
        {
          /* Put previous pointer */

          ptr = (long)RDB[uni0 + UNIVERSE_PTR_PREV_REG];
          PutPrivateData(ptr, reg, id);

          /* Return region */

          return reg;
        }

      /* Test surface */

      if (TestSurface(surf, x, y, z, NO, id) == YES)
        {
          /* Put previous pointer */
          
          ptr = (long)RDB[uni0 + UNIVERSE_PTR_PREV_REG];
          PutPrivateData(ptr, reg, id);
          
          /* Return region */

          return reg;
        }
    }

  /***************************************************************************/

  /* Something wrong */

  Die(FUNCTION_NAME, "Unable to find nest region");

  /* Avoid warning message */

  return 0;
}

/*****************************************************************************/
