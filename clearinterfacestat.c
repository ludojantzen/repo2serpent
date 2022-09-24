/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : clearinterfacestat.c                           */
/*                                                                           */
/* Created:       2015/01/21 (VVa)                                           */
/* Last modified: 2016/09/16 (VVa)                                           */
/* Version:       2.1.27                                                     */
/*                                                                           */
/* Description: Clears the power statistic in interfaces before next step    */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ClearInterfaceStat:"

/*****************************************************************************/

void ClearInterfaceStat()
{
  long loc0, loc1, ptr, type;

  /* Check that interfaces are defined */

  if ((loc0 = (long)RDB[DATA_PTR_IFC0]) < VALID_PTR)
    return;

  /* Loop over interfaces */

  while(loc0 > VALID_PTR)
    {

      /* Get interface type */

      type = (long)RDB[loc0 + IFC_TYPE];

      /* Only handle fuel behavior interfaces atm. */

      if ((type == IFC_TYPE_FUEP) || (type == IFC_TYPE_FPIP))
        {

          /* Get pointer to first rod */

          loc1 = (long)RDB[loc0 + IFC_PTR_FUEP];

          /* Loop over pins to relax power */

          while (loc1 > VALID_PTR)
            {

              /* Get pointers to power tally */

              ptr = (long)RDB[loc1 + IFC_FUEP_PTR_POWER];
              CheckPointer(FUNCTION_NAME, "(Pptr)", DATA_ARRAY, ptr);

              /* Clear statistics */

              ClearStat(ptr);

              /* Next pin */

              loc1 = NextItem(loc1);
            }
        }
      else
        {

          /* Get pointer to power tally */

          ptr = (long)RDB[loc0 + IFC_PTR_STAT];
          CheckPointer(FUNCTION_NAME, "(Pptr)", DATA_ARRAY, ptr);

          /* Clear statistics */

          ClearStat(ptr);

        }

      /* Next interface */

      loc0 = NextItem(loc0);
    }

  return;

  /***************************************************************************/
}

/*****************************************************************************/
