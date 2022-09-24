/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : b1fluxcorr.c                                   */
/*                                                                           */
/* Created:       2012/03/16 (JLe)                                           */
/* Last modified: 2019/08/24 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Returns leakage correction factor for flux from B1           */
/*              calculation.                                                 */
/*                                                                           */
/* Comments: - Tää otettiin uudelleen käyttöön 24.8.2019 / 2.1.32.           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "B1FluxCorr:"

/*****************************************************************************/

double B1FluxCorr(long gcu, double E)
{
  long ptr, ntot, n;
  double f;

  /* Factor is not availeble for the first active cycle */

  if (RDB[DATA_CYCLE_IDX] == RDB[DATA_CRIT_SKIP])
    return 1.0;

  /* Check pointer and energy */

  CheckPointer(FUNCTION_NAME, "(gcu)", DATA_ARRAY, gcu);
  CheckValue(FUNCTION_NAME, "E", "", E, ZERO, INFTY);

  /* Get pointer to microgroup energy grid */

  ptr = (long)RDB[DATA_MICRO_PTR_EGRID];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Number of groups */

  ntot = (long)RDB[ptr + ENERGY_GRID_NE] - 1;

  /* Get group index */

  if ((n = GridSearch(ptr, E)) < 0)
    return 1.0;

  /* Convert index */

  n = ntot - n - 1;
  CheckValue(FUNCTION_NAME, "n", "", n, 0, ntot - 1);

  /* Check index */

  CheckValue(FUNCTION_NAME, "n", "", n, 0, ntot - 1);

  /* Pointer to correction factors */

  ptr = (long)RDB[gcu + GCU_PTR_BURNUP_SPEC_CORR];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get factor */

  f = RDB[ptr + n];
  CheckValue(FUNCTION_NAME, "f", "", f, ZERO, INFTY);

  /* Return value */

  return f;
}

/*****************************************************************************/
