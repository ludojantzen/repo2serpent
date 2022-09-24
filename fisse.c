/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : fisse.c                                        */
/*                                                                           */
/* Created:       2017/11/21 (JLe)                                           */
/* Last modified: 2018/11/02 (RTu)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Calculates energy released in fission based on ENDF MT458    */
/*              data.                                                        */
/* Comments: Return value is in J:s                                          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FissE:"

/*****************************************************************************/

double FissE(long rea, double E, long id)
{
  long nuc, mode, fisse_ptr, nu_ptr, ptr, i;
  double fE;

  /* Check reaction pointer */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Pointer to nuclide data */

  nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

  /* Use fixed value (calculated in setfisse.c) if using constant energy deposition */

  fE = RDB[nuc + NUCLIDE_FISSE];

  /* Energy deposition mode */

  mode = (long)RDB[DATA_EDEP_MODE];

  if (mode > EDEP_MODE_CONSTANT)
    {
      /* Pointer to fission energy deposition data */

      if ((fisse_ptr = (long)RDB[nuc + NUCLIDE_PTR_FISSE_DATA]) > VALID_PTR)
        {
          /* Initialize fE */

          if (mode == EDEP_MODE_MT458)
            fE = RDB[DATA_EDEP_CAPT_E];
          else
            fE = 0.0;

          /* Pointer to total nubar data */

          nu_ptr = (long)RDB[rea + REACTION_PTR_TNUBAR];
          CheckPointer(FUNCTION_NAME, "(nu_ptr)", DATA_ARRAY, nu_ptr);

          /* Pointer to components flags */

          ptr = (long)RDB[DATA_EDEP_COMP];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Loop over components */

          for (i = 0; i < FISSE_COMP_BLOCK_SIZE; i++)
            if ((long)RDB[ptr + i] == YES)
              fE += FissEComp(fisse_ptr, nu_ptr, i, E, id);
        }
      else
        {
          /* No fission energy deposition data (MT458) available. Use fission Q-value. */

          if (mode == EDEP_MODE_MT458)
            fE = RDB[rea + REACTION_Q] + RDB[DATA_EDEP_CAPT_E];
          else
            fE = RDB[rea + REACTION_Q];
        }

      /* Convert MeV to J */

      fE = fE*MEV;
    }

  /* Check value */

  CheckValue(FUNCTION_NAME, "fE", "", fE/MEV, 0.0, 300.0);

  /* Return value */

  return fE;
}

/*****************************************************************************/
