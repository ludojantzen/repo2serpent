/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : iteratecc.c                                    */
/*                                                                           */
/* Created:       2014/07/07 (VVa)                                           */
/* Last modified: 2019/03/29 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Iterates external/internal solvers for coupled calculation   */
/*                                                                           */
/* Comments: -Handles solution relaxation and related calculations           */
/*           -Also transfers updated interfaces between MPI-tasks            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "IterateCC:"

/*****************************************************************************/

void IterateCC()
{
  /* Return if not running coupled calculation */

  if(RDB[DATA_RUN_CC] == NO)
    return;

  /* If on non-first predictor and running SIE, return and do not iterate */

  if((((long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP) &&
     !((long)RDB[DATA_BURN_STEP] == 0)) && (RDB[DATA_BURN_SIE] == (double)YES))
    {
      WDB[DATA_ITERATE] = (double)NO;
      return;
    }

  /* Collect results and relax them */

  RelaxCCResults();

  /* Only call the coupled solvers from MPI task 0 */

  if (mpiid == 0)
    {

      IterateFinix();

      /* IterateCOSY();*/

      IterateExternal();

      IterateInternal();

      /* Read updated data interfaces by task 0 */

      ReadDataInterfaces();
    }

  /* Broadcast updated interfaces from 0 to other tasks */

  BroadcastIFCData();

  /* Sample interface data if requested */

  SampleIFCData(NO);

  /* Increase iteration number */

  WDB[DATA_SOL_REL_ITER] = RDB[DATA_SOL_REL_ITER] + 1.0;

  /* Check stopping criterion */

  StopCCIter();

  return;
}
