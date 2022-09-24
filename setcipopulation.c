/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : setcipopulation.c                              */
/*                                                                           */
/* Created:       2016/05/31 (VVa)                                           */
/* Last modified: 2016/05/31 (VVa)                                           */
/* Version:       2.1.27                                                     */
/*                                                                           */
/* Description: Modifies the neutron population for burnup transport cycle.  */
/*                                                                           */
/* Comments: - Separated here from burnupcycle.c for cleaning purposes       */
/*           - Only does something if "set cpop" is used                     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SetCIPopulation:"

/*****************************************************************************/

void SetCIPopulation()
{

  /* Return if no alternate corrector population is given */

  if (RDB[DATA_BURN_CI_NBATCH] <= 0)
    return;

  if ((long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP)
    {
      /* Restore original population for predictor */

      WDB[DATA_CRIT_POP]      = RDB[DATA_BURN_CI_ORIG_NBATCH];
      WDB[DATA_CRIT_CYCLES]   = RDB[DATA_BURN_CI_ORIG_CYCLES];
      WDB[DATA_CRIT_SKIP]     = RDB[DATA_BURN_CI_ORIG_SKIP];
      WDB[DATA_FSP_CRIT_SKIP] = RDB[DATA_BURN_CI_ORIG_FSP_SKIP];

    }
  else if ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP)
    {
      /* Put alternate population for corrector */

      WDB[DATA_CRIT_POP]    = RDB[DATA_BURN_CI_NBATCH];
      WDB[DATA_CRIT_CYCLES] = RDB[DATA_BURN_CI_CYCLES];

      /* Check if we should use the first or further iteration */
      /* inactive cycle number*/

      if (RDB[DATA_BURN_CI_I] == 0.0)
        {
          WDB[DATA_CRIT_SKIP]   = RDB[DATA_BURN_CI_SKIP];
          WDB[DATA_FSP_CRIT_SKIP] = RDB[DATA_BURN_CI_SKIP];
        }
      else
        {
          WDB[DATA_CRIT_SKIP]   = RDB[DATA_BURN_CI_SKIP_2];
          WDB[DATA_FSP_CRIT_SKIP] = RDB[DATA_BURN_CI_SKIP_2];
        }
    }

}
