/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : calculaterelalpha.c                            */
/*                                                                           */
/* Created:       2014/10/06 (VVa)                                           */
/* Last modified: 2015/06/25 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Calculates relaxation factor for coupled calculation field   */
/*              relaxation                                                   */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CalculateRelAlpha:"

/*****************************************************************************/

void CalculateRelAlpha()
{
  double alpha, ncur, ntot;

  /* Get total simulated population before previous iteration */

  ntot = RDB[DATA_SOL_REL_NTOT];

  /* Get simulated population on previous iteration */

  ncur = RDB[DATA_SOL_REL_NCUR];

  /* Calculate new simulated population*/

  ntot = ntot + ncur;

  /* Store new total population */

  WDB[DATA_SOL_REL_NTOT] = (double)ntot;

  /* Calculate iteration weight */

  alpha = ncur/ntot;

  /* Check if no relaxation is wanted */

  if(RDB[DATA_SOL_REL_FACT] == 0.0)
    alpha = 1.0;

  /* Store iteration weight */

  WDB[DATA_SOL_REL_ALPHA] = alpha;

}

/*****************************************************************************/
