/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : samplenu.c                                     */
/*                                                                           */
/* Created:       2011/03/10 (JLe)                                           */
/* Last modified: 2017/11/15 (JLe)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Samples number of fission neutrons                           */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SampleNu:"

/*****************************************************************************/

long SampleNu(double nubar, long ieig, long id)
{
  double k;
  long nu, ptr;

  /* Get eigenvalue */

  if ((long)RDB[DATA_WIELANDT_MODE] != WIELANDT_MODE_NONE)
    {
      k = RDB[DATA_WIELANDT_KP];
      CheckValue(FUNCTION_NAME, "k", "", k, 0.001, 100.0);
    }
  else
    {
      ptr = (long)RDB[DATA_PTR_CYCLE_EIG_KEFF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      k = RDB[ptr + ieig];
      CheckValue(FUNCTION_NAME, "k", "", k, 0.001, 10.0);
    }

  /* Adjust nubar */

  nubar = nubar/k*RDB[DATA_EXT_SRC_NUBAR_F];

  /* Get integer part */

  nu = (long)nubar;

  /* Sample extra neutron */

  if (RandF(id) < nubar - (double)nu)
    nu++;

  /* Check number of neutrons */

  if ((nu > 1000) && ((long)RDB[DATA_ITER_MODE] == ITER_MODE_NONE))
    {
      /* Check k-eff */
      
      if (k < 0.5)
        Error(0, "Not enough multiplication to maintain chain reaction");
      else if (RDB[DATA_EXT_SRC_NUBAR_F] == 1.0)
        Die(FUNCTION_NAME, "WTF?");
    }

  /* Return value */

  return nu;
}

/*****************************************************************************/
