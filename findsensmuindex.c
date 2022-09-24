/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : findsenseindex.c                               */
/*                                                                           */
/* Created:       2017/04/06 (VVa)                                           */
/* Last modified: 2018/06/11 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Finds scattering cosine index for sensitivity calculations.  */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FindSensMuIndex:"

/*****************************************************************************/

long FindSensMuIndex(double mu)
{
  long loc0, nmu, type;
  double ang, frac;

  CheckValue(FUNCTION_NAME, "mu", "", mu, -1, 1);

  /* Get pointer to sensitivity block or return -1*/

  if ((loc0 = (long)RDB[DATA_PTR_SENS0]) < VALID_PTR)
    return -1;

  /* Get number of bins for scattering cosine */

  nmu = (long)RDB[loc0 + SENS_N_MU];
  CheckValue(FUNCTION_NAME, "nmu", "", nmu, 1, 10000);

  /* Get scattering angle binning type */

  type = (long)RDB[loc0 + SENS_MU_BIN_TYPE];

  if (type == SENS_MU_BIN_TYPE_EQCOS)
    {
      /* Equal bin width in scattering cosine */

      /* Calculate fraction of scattering cosine range covered by mu */
      /* (mu - (-1))/(1 - (-1)) */

      frac = (mu + 1)/2.0;
    }
  else if (type == SENS_MU_BIN_TYPE_EQANG)
    {
      /* Equal bin width in scattering angle */

      /* Convert cosine to angle */

      ang = acos(mu);

      /* Calculate fraction of scattering angle range covered by ang */
      /* (ang - (-PI))/(PI - (-PI)) */

      frac = (ang + PI)/(2*PI);
    }
  else
    {
      /* Unknown type */

      Die(FUNCTION_NAME, "Unknown type for scattering cosine binning: %ld", type);

      /* This needs to be set to something to appease the compiler (and for NOFATAL) */

      frac = 0.0;
    }

  /* Reduce frac from 1.0 to prevent out-of-bounds in array        */
  /* frac = 1.0 might actually never happen, but better to be safe */

  if (frac >= 1.0)
    frac = 1.0 - ZERO;

  /* Calculate bin index based on fraction */

    return (long)(frac*nmu);
}
