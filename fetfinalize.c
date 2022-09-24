/******************************************************************************/
/*                                                                            */
/* serpent 2 (beta-version) : fetfinalize.c                                   */
/*                                                                            */
/* Created:       2018/02/17 (BEe)                                            */
/* Last modified: 2020/06/04 (JLe)                                            */
/* Version:       2.1.32                                                      */
/*                                                                            */
/* Description: Finalizes and loads an FET's coefficient value and            */
/*              uncertainty                                                   */
/*                                                                            */
/* Comments:    -                                                             */
/*                                                                            */
/******************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FETFinalize:"

/******************************************************************************/

void FETFinalize(const double *const params, long statPtr, long coefIdx,
                 long layer, double * valueOut, double * relUncOut)
{
  long coefOffset, popIdx, resultPtr, statWidth;
  double pCount, relUnc, value;
  const long totalCoefficients = (long)params[FET_PARAM_NCOEF_TOTAL];

  /* Get the population that contributed to the FET */

  statWidth = totalCoefficients + 1;
#ifdef DEBUG
  if (statWidth != (long)RDB[(long)RDB[statPtr + SCORE_PTR_NMAX]])
    Die(FUNCTION_NAME, "Did something change in the FET Stats? The width is no "
                       "longer one more than the number of coefficients");
#endif /* DEBUG */

  popIdx = totalCoefficients + layer*statWidth;
  resultPtr = (long)RDB[statPtr + SCORE_PTR_DATA] + popIdx*STAT_BLOCK_SIZE;
  pCount = RES1[resultPtr + BUF_FET_VAL];

  /* Get the pointer to the results buffer */

  coefOffset = coefIdx + layer*statWidth;
  resultPtr = (long)RDB[statPtr + SCORE_PTR_DATA] + coefOffset*STAT_BLOCK_SIZE;

  /* Calculate the value and normalize to the physical volume */

  value = (RES1[resultPtr + BUF_FET_VAL] / pCount);
  (*valueOut) = value * (params[FET_PARAM_FET_VOLUME] / params[FET_PARAM_PHYS_VOLUME]);

  /* Calculate the relative uncertainty. */

  relUnc = (RES1[resultPtr + BUF_FET_SUM_SQ_SUM] - (1 / pCount)*RES1[resultPtr + BUF_FET_VAL]*RES1[resultPtr + BUF_FET_VAL])
           / (pCount*(pCount - 1));
  (*relUncOut) = sqrt(fabs(relUnc)) / fabs(value);

}
