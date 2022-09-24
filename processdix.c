/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processdix.c                                   */
/*                                                                           */
/* Created:       2018/01/17 (JLe)                                           */
/* Last modified: 2019/08/20 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Sets up the index arrays for double-indexing of cross        */
/*              sections.                                                    */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessDIX:"

/*****************************************************************************/

void ProcessDIX(long nuc)
{
  long erg0, erg1, loc0, loc1, ne0, ne1, n0, n1, ptr;

  /* Check dosimetry data */

  if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DOSIMETRY)
    return;

  /* Pointer to unionized energy grid */

  if ((erg0 = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID]) < VALID_PTR)
    WDB[DATA_OPTI_DIX] = (double)NO;

  /* Check thinning tolerance */

  if (RDB[DATA_ERG_TOL] > 0.0)
    WDB[DATA_OPTI_DIX] = (double)NO;

  /* Check if double-indexing is used */

  if ((long)RDB[DATA_OPTI_DIX] == NO)
    return;

  /* Pointer to nuclide energy grid */

  erg1 = (long)RDB[nuc + NUCLIDE_PTR_EGRID];
  CheckPointer(FUNCTION_NAME, "(erg1)", DATA_ARRAY, erg1);

  /* Check pointers */

  if (erg1 == erg0)
    return;

  /* Get array sizes */

  ne0 = (long)RDB[erg0 + ENERGY_GRID_NE];
  CheckValue(FUNCTION_NAME, "ne0", "", ne0, 10, 10000000);

  ne1 = (long)RDB[erg1 + ENERGY_GRID_NE];
  CheckValue(FUNCTION_NAME, "ne1", "", ne1, 10, ne0);

  /* Pointers to data */

  loc0 = (long)RDB[erg0 + ENERGY_GRID_PTR_DATA];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  loc1 = (long)RDB[erg1 + ENERGY_GRID_PTR_DATA];
  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

  /* Allocate memory */

  ptr = ReallocMem(DATA_ARRAY, ne0);
  WDB[erg1 + ENERGY_GRID_PTR_DIX_IDX] = (double)ptr;

  /* Loop over grids and set indexes */

  n1 = 0;
  for (n0 = 0; n0 < ne0; n0++)
    {
      /* Compare energies */

      while(RDB[loc0 + n0] >= RDB[loc1 + n1])
        n1++;

      /* Put index */

      WDB[ptr++] = (double)(n1 - 1);
    }
}

/*****************************************************************************/
