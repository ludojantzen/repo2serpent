/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : clearvaluepairs.c                              */
/*                                                                           */
/* Created:       2019/11/10 (JLe)                                           */
/* Last modified: 2019/11/10 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Clears all data stored as vp data                            */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ClearValuePairs:"

/*****************************************************************************/

void ClearValuePairs()
{
  long loc0, ptr, i;

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_VP0];
  while (loc0 > VALID_PTR)
    {
      /* Get pointer to data */

      ptr = (long)RDB[loc0 + VP_PTR];
      CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);

      /* Loop over OpenMP thread id's */

      for (i = 0; i < (long)RDB[DATA_OMP_MAX_THREADS]; i++)
        {
          /* Put key */

          PutPrivateData(ptr, -INFTY, i);

          /* Put value */

          PutPrivateData(ptr + 1, -INFTY, i);
        }

      /* Pointer to next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/
}

/*****************************************************************************/
