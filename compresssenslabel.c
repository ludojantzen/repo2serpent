/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : compresssenslabel.c                            */
/*                                                                           */
/* Created:       2017/03/21 (VVa)                                           */
/* Last modified: 2018/06/19 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Compresses sensitivity event label based on sub indices      */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CompressSensLabel:"

/*****************************************************************************/

long CompressSensLabel(long imat, long izai, long irea, long iene,
                       long hitMiss)
{
  long loc0, idx, nzai, nrea, nxyz, nene, nmu, ixyz, imu;

  /* Get pointer to sensitivity block or return -1 */

  if ((loc0 = (long)RDB[DATA_PTR_SENS0]) < VALID_PTR)
    return -1;


  /* Get sizes of sub bins */

  nzai = (long)RDB[loc0 + SENS_N_ZAI];
  nrea = (long)RDB[loc0 + SENS_N_PERT] + 1;
  nene = (long)RDB[loc0 + SENS_N_ENE] + 1;
  nxyz = (long)RDB[loc0 + SENS_N_SPAT];
  nmu  = (long)RDB[loc0 + SENS_N_MU];

  /* Set these to zero for now */

  imu = 0;
  ixyz = 0;

  /* Check the sub indices */

  CheckValue(FUNCTION_NAME, "imat", "", imat, 0, (long)RDB[loc0 + SENS_N_MAT]);
  CheckValue(FUNCTION_NAME, "izai", "", izai, 0, nzai);
  CheckValue(FUNCTION_NAME, "irea", "", irea, 0, nrea);
  CheckValue(FUNCTION_NAME, "iene", "", iene, 0, nene);
  CheckValue(FUNCTION_NAME, "imu", "", imu, 0, nmu);
  CheckValue(FUNCTION_NAME, "ixyz", "", ixyz, 0, nxyz);

  /* Calculate index */

  idx = 1 + imat*nzai*nrea*nene*nxyz*nmu + izai*nrea*nene*nxyz*nmu + irea*nene*nxyz*nmu + iene*nxyz*nmu + ixyz*nmu + imu;
  CheckValue(FUNCTION_NAME, "idx", "", idx, 1, (long)RDB[loc0 + SENS_MAX_LABEL]);

  /* Set sign based on hit/miss */

  idx *= hitMiss;

  /* Return index */

  return idx;
}
