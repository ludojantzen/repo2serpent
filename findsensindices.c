/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : findsensindices.c                               */
/*                                                                           */
/* Created:       2017/04/05 (VVa)                                           */
/* Last modified: 2018/06/20 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Finds material, nuclide and energy indices for Sens event     */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FindSensIndices:"

/*****************************************************************************/

void FindSensIndices(long mat, long rea, double E,
                    long *imat, long *izai, long *iene)
{
  long nucZAI, sens, ptr, idx, mat0, matarr, zaiarr;

  /* Get pointer to sensitivity block or return */

  if ((sens = (long)RDB[DATA_PTR_SENS0]) < VALID_PTR)
    return;

  /* Get pointer to index arrays */

  matarr = (long)RDB[sens + SENS_PTR_MAT_INDICES];
  CheckPointer(FUNCTION_NAME, "(matarr)", DATA_ARRAY, matarr);

  /* Get pointer to index array */

  zaiarr = (long)RDB[sens + SENS_PTR_ZAI_INDICES];
  CheckPointer(FUNCTION_NAME, "(zaiarr)", DATA_ARRAY, zaiarr);

  /*******************************/
  /* find perturbed material idx */
  /*******************************/

  /* Get pointer to material list */

  if ((ptr = (long)RDB[sens + SENS_PTR_MAT_ARR]) > VALID_PTR)
    {
      /* Loop over the material list to find material */
      /* TODO: May be much faster if the Sens_IDX would be added */
      /* to the MATERIAL_BLOCK in ProcessSens */

      idx = 0;

      while ((mat0 = (long)RDB[ptr + idx]) > VALID_PTR)
        {
          if (mat0 == mat)
            break;

          /* Compare to parent material (if exists) */

          if ((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] > VALID_PTR)
            if (mat0 == (long)RDB[mat + MATERIAL_DIV_PTR_PARENT])
              break;

          /* Increment index */

          idx++;
        }

      /* Check if found */

      if (mat0 > VALID_PTR)
        {
          /* Was found */

          /* Put imat */

          *imat = idx;

          /* Increment index to take in account none, total and sum bins */

          *imat = *imat + 1;

          if ((long)RDB[sens + SENS_SCORE_FLAGS] & SENS_SCORE_FLAG_MAT_TOT)
            *imat = *imat + 1;

          if ((long)RDB[sens + SENS_SCORE_FLAGS] & SENS_SCORE_FLAG_MAT_SUM)
            *imat = *imat + 1;

        }
      else
        *imat = -1;
    }
  else
    *imat = -1;

  /* Score not found materials to "lost" index if scoring total */

  if ((*imat < 0) && ((long)RDB[sens + SENS_SCORE_FLAGS] & SENS_SCORE_FLAG_MAT_TOT))
    *imat = (long)RDB[matarr + SENS_NON_MAT_IDX];

  /******************************/
  /* find perturbed nuclide idx */
  /******************************/

  nucZAI = (long)RDB[((long)RDB[rea + REACTION_PTR_NUCLIDE] + NUCLIDE_ZAI)];
  *izai = FindSensZAIIndex(nucZAI);

  /* Score not found ZAIs to "lost" index if scoring total */

  if ((*izai < 0) && ((long)RDB[sens + SENS_SCORE_FLAGS] & SENS_SCORE_FLAG_ZAI_TOT))
    *izai = (long)RDB[zaiarr + SENS_NON_ZAI_IDX];

  /*******************************/
  /* Find perturbed energy index */
  /*******************************/

  *iene = FindSensEIndex(E);
}
