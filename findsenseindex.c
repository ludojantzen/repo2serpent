/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : findsenseindex.c                               */
/*                                                                           */
/* Created:       2017/04/06 (VVa)                                           */
/* Last modified: 2018/06/11 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Finds energy index for sens event from the linked grid       */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FindSensEIndex:"

/*****************************************************************************/

long FindSensEIndex(double E)
{
  long loc0, ptr, idx;

  /* Get pointer to sensitivity block or return -1 */

  if ((loc0 = (long)RDB[DATA_PTR_SENS0]) < VALID_PTR)
    return -1;

  /* Get energy grid used in Sens */

  ptr = (long)RDB[loc0 + SENS_PTR_EGRID];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Find index from grid */

  idx = GridSearch(ptr, E);

  /* Check if found */

  if (idx > -1)
    {
      /* Was found */

      /* Put index */
      /* NB: Currently index 0 means that was not found */

      return idx + 1;
    }
  else
    {
      /* Was not found */

      /* Put index */
      /* NB: Currently index 0 means that was not found */

      return 0;
    }
}
