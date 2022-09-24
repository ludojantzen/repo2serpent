/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : updatecistop.c                                 */
/*                                                                           */
/* Created:       2014/09/04 (VVa)                                           */
/* Last modified: 2018/06/20 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Calculates stopping criterion in corrector iteration         */
/*                                                                           */
/* Comments: Called from burnmaterials                                       */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "UpdateCIStop:"

/*****************************************************************************/

void UpdateCIStop(long mat, double *N, long id)
{
  long dep, rea, nuc, i;
  double xs, n, totxs;

  /* Reset macroscopic total xs */

  totxs = 0.0;

  /* Get pointer to transmutation list */
  dep = (long)RDB[mat + MATERIAL_PTR_DEP_TRA_LIST];

  while (dep > VALID_PTR)
    {

      /* Pointer to reaction */

      rea = (long)RDB[dep + DEP_TRA_PTR_REA];
      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

      /* Pointer to nuclide data */

      nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Get index to composition vector */

      if ((i = (long)TestValuePair(nuc + NUCLIDE_PTR_MATRIX_IDX,
                                   (double)mat, id)) < 0)
        Die(FUNCTION_NAME, "i < 0");

      /* Get average of cross section from prev. step */

      xs = RDB[dep + DEP_TRA_AV0];

      /* Add to total xs multiplied with current iterations (non-relaxed) */
      /* nuclide density */

      totxs = totxs + xs*N[i];

      /* Next reaction */

      dep = NextItem(dep);
    }

  /* Get pointer to fission list */
  dep = (long)RDB[mat + MATERIAL_PTR_DEP_FISS_LIST];

  while (dep > VALID_PTR)
    {

      /* Pointer to reaction */

      rea = (long)RDB[dep + DEP_TRA_PTR_REA];
      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

      /* Pointer to nuclide data */

      nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Get index to composition vector */

      if ((i = (long)TestValuePair(nuc + NUCLIDE_PTR_MATRIX_IDX,
                                   (double)mat, id)) < 0)
        Die(FUNCTION_NAME, "i < 0");

      /* Get average of cross section from prev. step */

      xs = RDB[dep + DEP_TRA_AV0];

      /* Add to total xs multiplied with current iterations (non-relaxed) */
      /* nuclide density */

      totxs = totxs + xs*N[i];

      /* Next reaction */

      dep = NextItem(dep);
    }

  /* Update EOS macroscopic absorption cross-section */

  WDB[mat + MATERIAL_CI_EOS_ABSXS] = totxs;

  /* Get iteration number for averaging */

  n = RDB[DATA_BURN_CI_I] + 1.0;

  /* Update average absorption cross section */

  WDB[mat + MATERIAL_CI_AVE_ABSXS] = (n-1)/n*RDB[mat + MATERIAL_CI_AVE_ABSXS] +
    1.0/n*RDB[mat + MATERIAL_CI_EOS_ABSXS];

  /* Update Average of square of absorption cross section */

  WDB[mat + MATERIAL_CI_AVE_ABSXS2] = (n-1)/n*RDB[mat + MATERIAL_CI_AVE_ABSXS2]
    + 1.0/n*RDB[mat + MATERIAL_CI_EOS_ABSXS]*RDB[mat + MATERIAL_CI_EOS_ABSXS];

  /*
  fprintf(outp, "n = %f, EOS = %E, AXS2 = %E, (AXS)^2 = %E",
          n, RDB[mat + MATERIAL_CI_EOS_ABSXS],
          RDB[mat + MATERIAL_CI_AVE_ABSXS2],
          RDB[mat + MATERIAL_CI_AVE_ABSXS]*RDB[mat + MATERIAL_CI_AVE_ABSXS]);
  */

  /* Calculate ideal deviation */

  if (n > 1)
    {
      WDB[mat + MATERIAL_CI_IDE] =
        sqrt(1.0 / (n - 1.0) * (RDB[mat + MATERIAL_CI_AVE_ABSXS2] -
                                RDB[mat + MATERIAL_CI_AVE_ABSXS]*
                                RDB[mat + MATERIAL_CI_AVE_ABSXS]));
    }
  else
    {
      /* On first iteration, ideal deviation is just the difference */
      /* compared to the BOS total absorption cross-section         */

      WDB[mat + MATERIAL_CI_IDE] = fabs(RDB[mat + MATERIAL_CI_AVE_ABSXS]-
                                      RDB[mat + MATERIAL_CI_BOS_ABSXS]);
    }

}

/*****************************************************************************/
