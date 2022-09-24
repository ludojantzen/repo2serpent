/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : nocorrector.c                                  */
/*                                                                           */
/* Created:       2016/05/31 (VVa)                                           */
/* Last modified: 2018/09/14 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Checks, whether the corrector step should be omitted in      */
/*              burnup calculations                                          */
/*                                                                           */
/* Comments: - Separated here from burnupcycle.c for cleaning purposes       */
/*           - Returns 1 if no corrector should be run, 0 otherwise          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "NoCorrector:"

/*****************************************************************************/

long NoCorrector()
{

  /* Corrector should always be run in SIE calculation */

  if (RDB[DATA_BURN_SIE] == (double)YES)
    return NO;

  /* No corrector if: -User requested none,                      */
  /*                  -On decay step (flux is zero, CE is exact) */
  /*                  -Activation step (26.9.2015 / JLe),        */
  /*                  -When flux is zero (from normalization),   */
  /*                   making CE exact                           */

  if (((long)RDB[DATA_BURN_CORR_TYPE] == CORR_TYPE_NONE) ||
      (RDB[DATA_BURN_STEP_TYPE] == DEP_STEP_DEC_STEP) ||
      (RDB[DATA_BURN_STEP_TYPE] == DEP_STEP_DEC_TOT)  ||
      (RDB[DATA_BURN_STEP_TYPE] == DEP_STEP_ACT_STEP) ||
      (RDB[DATA_BURN_STEP_TYPE] == DEP_STEP_ACT_TOT) ||
      (Mean((long)RDB[RES_TOT_NEUTRON_FLUX], 0) <= 0.0))
    {

      /* Break from PC-loop */

      return YES;
    }

  return 0;
}
