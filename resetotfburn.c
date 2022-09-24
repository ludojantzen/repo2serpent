/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : resetotfburn.c                                 */
/*                                                                           */
/* Created:       2018/03/25 (JLe)                                           */
/* Last modified: 2018/04/12 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Resets concentrations for on-the-fly burnup solver           */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ResetOTFBurn:"

/*****************************************************************************/

void ResetOTFBurn()
{
  long loc0, iso;

  /* Check OTF burn mode */

  if ((long)RDB[DATA_OTF_BURN_MODE] == NO)
    return;

  /* Set implicit reaction rates */

  if ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP)
    WDB[DATA_OPTI_IMPLICIT_RR] = (double)NO;
  else
    WDB[DATA_OPTI_IMPLICIT_RR] = RDB[DATA_OPTI_IMPLICIT_RR0];

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_OTF_BURN0];
  while (loc0 > VALID_PTR)
    {
      /* Pointer to composition */

      iso = (long)RDB[loc0 + OTF_BURN_PTR_COMP];
      CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

      /* Remember concentration */

      if ((long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP)
        WDB[loc0 + OTF_BURN_ADENS0] = RDB[iso + COMPOSITION_ADENS];

      /* Reset iterated density */

      WDB[loc0 + OTF_BURN_ADENS] = 0.0;

      /* Reset atomic density */

      if (((long)RDB[DATA_BURN_CORR_TYPE] == CORR_TYPE_NONE) ||
          ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP))
        WDB[iso + COMPOSITION_ADENS] = 0.0;

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/
}

/*****************************************************************************/
