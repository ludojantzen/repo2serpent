/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : putotfburnconc.c                               */
/*                                                                           */
/* Created:       2018/03/25 (JLe)                                           */
/* Last modified: 2018/03/25 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Puts nuclide compositions for the on-the-fly burnup solver   */
/*                                                                           */
/* Comments: - Tämä resetoi majorant extrat predictor-stepeille ja palauttaa */
/*             BOS-konsentraatiot correctorille. Kutsutaan ennen ja jälkeen  */
/*             transport-silmukan.                                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PutOTFBurnConc:"

/*****************************************************************************/

void PutOTFBurnConc()
{
  long loc0, ptr, iso;
  double val; 

  /* Check OTF burn mode */

  if ((long)RDB[DATA_OTF_BURN_MODE] == NO)
    return;
  
  /* Reset maximum concentrations (this is otherwise not done for */
  /* predictor steps) */

  ptr = (long)RDB[DATA_MAJORANT_PTR_EXTRA_XS];
  while (ptr > VALID_PTR)
    {
      /* Check type */
      
      if ((long)RDB[ptr + MAJORANT_EXTRA_TYPE] == MAJORANT_EXTRA_OTF_BURN)
        WDB[ptr + MAJORANT_EXTRA_FRAC] = 0.0;
      
      /* Pointer to next */

      ptr = NextItem(ptr);
    }

  /* Check predictor-corrector calculation */

  if (((long)RDB[DATA_BURN_CORR_TYPE] != CORR_TYPE_NONE) &&
      ((long)RDB[DATA_BURN_STEP_PC] != CORRECTOR_STEP))
    return;

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_OTF_BURN0];
  while (loc0 > VALID_PTR)
    {
      /* Pointer to composition */

      iso = (long)RDB[loc0 + OTF_BURN_PTR_COMP];
      CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

      /* Put concentration */

      val = RDB[loc0 + OTF_BURN_ADENS0];
      WDB[iso + COMPOSITION_ADENS] = val;
      
      /* Pointer to majorant */

      if ((ptr = (long)RDB[loc0 + OTF_BURN_PTR_MAJ]) > VALID_PTR)
        {
          /* Compare fraction */

          if (val > RDB[ptr + MAJORANT_EXTRA_FRAC])
            WDB[ptr + MAJORANT_EXTRA_FRAC] = val;
        }

      /* Next */

      loc0 = NextItem(loc0);
    }
}

/*****************************************************************************/
