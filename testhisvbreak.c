/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : testhistvbreak.c                               */
/*                                                                           */
/* Created:       2020/03/22 (JLe)                                           */
/* Last modified: 2020/07/23 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Tests break point for history variation.                     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "TestHisvBreak:"

/*****************************************************************************/

long TestHisvBreak()
{
  double val;

  /* Check if history variations are used */

  if ((long)RDB[DATA_PTR_HISV0] < VALID_PTR)
    return NO;

  /* Check units */

  if ((val = RDB[DATA_HISV_BREAK_POINT]) == NO)
    return NO;
  else if (val > 0.0)
    {
      /* Compare to burnup */

      if (fabs(RDB[DATA_BURN_CUM_BURNUP]/val - 1.0) < 1E-4)
        return YES;
      else
        return NO;
    }
  else
    {
      /* Compare to time */

      if (fabs(-RDB[DATA_BURN_CUM_BURNTIME]/val/86400.0 - 1.0) < 1E-4)
        return YES;
      else
        return NO;
    }

    printf("%E %E %E\n", RDB[DATA_HISV_BREAK_POINT],RDB[DATA_BURN_CUM_BURNUP],
         RDB[DATA_BURN_CUM_BURNTIME]);
  Pause(FUNCTION_NAME);

}

/*****************************************************************************/
