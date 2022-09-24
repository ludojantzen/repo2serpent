/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : iternucxs.c                                    */
/*                                                                           */
/* Created:       2017/02/28 (VVa)                                           */
/* Last modified: 2018/06/21 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Returns sum of macroscopic cross sections for iterated       */
/*              nuclides in certain material.                                */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "IterNucXS:"

/*****************************************************************************/

double IterNucXS(long mat, double E, long mt, long ng, long id)
{
  long iso, nuc, rea, ptr;
  double abs, xs, adens;

  /* Check iter calculation flag */

  if ((long)RDB[DATA_ITER_MODE] != ITER_MODE_NUCLIDE)
    return 0.0;

  /* Check material pointer */

  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

  /* Reset cross section */

  abs = 0.0;

  /* Avoid compiler warning */

  rea = -1;

  /* Get pointer to isotope list */

  if ((ptr = (long)RDB[mat + MATERIAL_PTR_ITER_ISO_LIST]) > VALID_PTR)
    {
      while ((iso = (long)RDB[ptr]) > VALID_PTR)
        {
          /* Pointer to nuclide */

          nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
          CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

          /* Get pointer to reaction */

          if (mt == MT_MACRO_TOTXS)
            rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
          else if (mt == MT_MACRO_ABSXS)
            rea = (long)RDB[nuc + NUCLIDE_PTR_SUM_ABSXS];
          else if (mt == MT_MACRO_TMP_MAJORANTXS)
            {
              /* Get pointer to reaction data */

              rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
              CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

              /* Get pointer to separate majorant if not using MGXS */
              /* Multi-group majorant mode stores TMP_MAJORANT directly */
              /* to NUCLIDE_PTR_TOTXS */

              if (RDB[DATA_OPTI_MG_MODE] == (double)NO)
                if ((rea = (long)RDB[rea + REACTION_PTR_TMP_MAJORANT]) < VALID_PTR)
                  Die(FUNCTION_NAME, "Temperature majorant did not exist");
            }
          else
            Die(FUNCTION_NAME, "Invalid reaction mode %ld", mt);

          /* Check pointer */

          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

          if (RDB[DATA_OPTI_MG_MODE] == (double)YES)
            {
              /* Multi-group majorant mode */

              /* Get pointer to MGXS */

              xs = MGXS(rea, E, ng);

            }
          else
            {
              /* Continuous energy mode */

              /* Get microscopic cross section */

              xs = MicroXS(rea, E, id);
            }

          /* Get atomic density */

          adens = RDB[iso + COMPOSITION_ADENS];

          /* Add tot total */

          abs = abs + xs*adens;

          /* Increment pointer */

          ptr++;
        }
    }

  return abs;
}

/*****************************************************************************/
