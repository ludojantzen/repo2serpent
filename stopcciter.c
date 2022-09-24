/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : stopcciter.c                                   */
/*                                                                           */
/* Created:       2015/01/19 (VVa)                                           */
/* Last modified: 2016/09/16 (VVa)                                           */
/* Version:       2.1.27                                                     */
/*                                                                           */
/* Description: Updates coupled calculation stopping criterion               */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "StopCCIter:"

/*****************************************************************************/

void StopCCIter()
{
  long tb;

  if (RDB[DATA_SIMULATION_MODE] == (double)SIMULATION_MODE_DYN)
    {
      /* Get time bin index */

      tb = (long)RDB[DATA_DYN_TB];

      /* Check maximum number of iterations */

      if (RDB[DATA_SOL_REL_ITER] >= RDB[DATA_SOL_REL_MAX_ITER])
        {
          /* If we've reached the maximum number of iterations */
          /* do not iterate anymore */

          WDB[DATA_ITERATE] = (double)NO;
        }
    }
  else
    {
      /* Check maximum number of iterations */

      if (RDB[DATA_SOL_REL_ITER] >= RDB[DATA_SOL_REL_MAX_ITER])
        {
          /* If we've reached the maximum number of iterations */
          /* do not iterate anymore */

          WDB[DATA_ITERATE] = (double)NO;

        }
      else if (RDB[DATA_BURN_SIE] == (double)YES)
        {
          /* If we are running a burn-up calculation with SIE */
          /* do not iterate */

          WDB[DATA_ITERATE] = (double)NO;

        }

    }
}
