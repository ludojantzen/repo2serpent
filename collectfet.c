/******************************************************************************/
/*                                                                            */
/* serpent 2 (beta-version) : collectfet.c                                    */
/*                                                                            */
/* Created:       2018/02/16 (BWe)                                            */
/* Last modified: 2018/02/23 (BWe)                                            */
/* Version:       2.1.31                                                      */
/*                                                                            */
/* Description: Collects the FET statistics and computes the intermediate     */
/*              values.                                                       */
/* Comments: -                                                                */
/*                                                                            */
/******************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CollectFET:"

/******************************************************************************/

void CollectFET(const double *const params, long statPtr, long layers)
{
  long coefIdx, coefOffset, layer, n0, n1, n2;
  long bufPtr, resultPtr, statWidth;
  const long coef0 = ((long)params[FET_PARAM_NCOEF0] > 1 ?
                      (long)params[FET_PARAM_NCOEF0] : 1);
  const long coef1 = ((long)params[FET_PARAM_NCOEF1] > 1 ?
                      (long)params[FET_PARAM_NCOEF1] : 1);
  const long coef2 = ((long)params[FET_PARAM_NCOEF2] > 1 ?
                      (long)params[FET_PARAM_NCOEF2] : 1);
  const long totalCoefficients = (long)params[FET_PARAM_NCOEF_TOTAL];
  const long totalLayers = (layers > 0 ? layers : 1);

  const long scorePtr = (long)RDB[statPtr + SCORE_PTR_BUF];
  CheckPointer(FUNCTION_NAME, "(scorePtr)", BUF_ARRAY, scorePtr);

  /* Stats for FETS only have two dimensions. The first is the linearized     */
  /* coefficient array, the second is the regions. The width of the array     */
  /* is equal to the number of coefficients plus 1 (the last bin is used to   */
  /* track the particle count, see AllocFETCache().                           */

  statWidth = totalCoefficients + 1;
#ifdef DEBUG
  /* Ensure statWidth is correct */

  if (statWidth != (long)RDB[(long)RDB[statPtr + SCORE_PTR_NMAX]])
    Die(FUNCTION_NAME, "Did something change in the FET Stats? The width is no "
                       "longer one more than the number of coefficients");

  /* Ensure the buffer is reduced */

  if ((long)RDB[DATA_BUF_REDUCED] == NO)
    Die(FUNCTION_NAME, "FET score buffer is not reduced");
#endif /* DEBUG */

  /* Process in parallel */

#ifdef OPEN_MP
#pragma omp parallel for collapse(3) default(shared) \
  private(bufPtr, coefIdx, coefOffset, layer, n0, n1, n2, resultPtr) \
  shared(BUF, RDB, RES1, statPtr, statWidth)
#endif
  for (n2 = 0; n2 < coef2; ++n2)
  for (n1 = 0; n1 < coef1; ++n1)
  for (n0 = 0; n0 < coef0; ++n0)
  {
    coefIdx = FETIdx(params, n2, n1, n0);

    for (layer = 0; layer < totalLayers; ++layer)
    {
      /* This routine is distilled from BufVal() from assuming:               */
      /*   dim      == 2                                                      */
      /*   va_arg0  == idx                                                    */
      /*   va_arg1  == l                                                      */

      coefOffset = coefIdx + layer*statWidth;
      bufPtr = scorePtr + coefOffset*BUF_BLOCK_SIZE;
      CheckPointer(FUNCTION_NAME, "(bufPtr)", BUF_ARRAY, bufPtr);

      /* The indexing is the same for AddStat(), so just jump right in with   */
      /* the current coef_idx.                                                */

      resultPtr = (long)RDB[statPtr + SCORE_PTR_DATA] + coefOffset*STAT_BLOCK_SIZE;
      RES1[resultPtr + BUF_FET_VAL] += BUF[bufPtr + BUF_FET_VAL];
      RES1[resultPtr + BUF_FET_SUM_SQ_SUM] += BUF[bufPtr + BUF_FET_SUM_SQ_SUM];

      /* This value is used only for the online algorithms, so it can be      */
      /* ignored here when collecting everything else.                        */
      /*                                                                      */
      /*RES1[resultPtr + BUF_FET_SUM] += BUF[bufPtr + BUF_FET_SUM];           */

      if (coefIdx + 1 == totalCoefficients)
      {
        /* Add the population if this is the last coefficient */

        bufPtr += BUF_BLOCK_SIZE;
        resultPtr += STAT_BLOCK_SIZE;
        RES1[resultPtr + BUF_FET_VAL] += BUF[bufPtr + BUF_FET_VAL];
      }
    }
  }
}
