/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : fetidx.c                                       */
/*                                                                           */
/* Created:       2018/02/12 (BWe)                                           */
/* Last modified: 2018/02/23 (BWe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Returns stat index for FETS                                  */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FETIdx:"

/*****************************************************************************/

long FETIdx(const double *const params, long i2, long i1, long i0)
{
  long N0, N1, N2, idx, nmax;

  /* Get number of bins */

  N0 = (long)params[FET_PARAM_NCOEF0];
  N0 = N0 > 0 ? N0 : 1;
  N1 = (long)params[FET_PARAM_NCOEF1];
  N1 = N1 > 0 ? N1 : 1;
  N2 = (long)params[FET_PARAM_NCOEF2];
  N2 = N2 > 0 ? N2 : 1;

  /* Check bins */

  if ((i2 < 0) || (N2 <= i2))
    Die(FUNCTION_NAME, "Invalid i2: %ld \t(max: %ld)", i2, N1 - 1);
  if ((i1 < 0) || (N1 <= i1))
    Die(FUNCTION_NAME, "Invalid i1: %ld \t(max: %ld)", i1, N1 - 1);
  if ((i0 < 0) || (N0 <= i0))
    Die(FUNCTION_NAME, "Invalid i0: %ld \t(max: %ld)", i0, N0 - 1);

  /* Calculate index */

  idx = 0;
  nmax = 1;

  idx = idx + i2*nmax;
  nmax = nmax*N2;

  idx = idx + i1*nmax;
  nmax = nmax*N1;

  idx = idx + i0*nmax;

  /* Check index */

  nmax = (long)params[FET_PARAM_NCOEF_TOTAL];
  if ((idx < 0) || (idx > nmax))
    Die(FUNCTION_NAME, "Error in idx: %ld (Must be between 0 and %ld)", idx, nmax);

  /* Return index */

  return idx;
}

/*****************************************************************************/
