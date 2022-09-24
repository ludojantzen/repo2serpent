/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scoreiternucs.c                                */
/*                                                                           */
/* Created:       2017/03/02 (VVa)                                           */
/* Last modified: 2017/03/02 (VVa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Scores neutron balance for iterable nuclides                 */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreIterNucs:"

/*****************************************************************************/

void ScoreIterNucs(double flx, long mat, double E, double wgt, double g, long id)
{
  long rea, ptr, iso, nuc;
  double adens, xs, absxs;

  /* Check iter calculation flag */

  if ((long)RDB[DATA_ITER_MODE] != ITER_MODE_NUCLIDE)
    return;

  /* Check material pointer */

  if (mat < VALID_PTR)
    return;

  /* Get pointer to isotope list or return if it doesn't exit */

  if ((ptr = (long)RDB[mat + MATERIAL_PTR_ITER_ISO_LIST]) < VALID_PTR)
    return;

  absxs = 0.0;

  while ((iso = (long)RDB[ptr]) > VALID_PTR)
    {
      /* Pointer to nuclide */

      nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Get pointer to reaction */

      rea = (long)RDB[nuc + NUCLIDE_PTR_SUM_ABSXS];

      /* Check pointer */

      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

      /* Get microscopic cross section */

      xs = MicroXS(rea, E, id);

      /* Get atomic density */

      adens = RDB[iso + COMPOSITION_ADENS];

      /* Add tot total absorption */

      absxs = absxs + xs*adens;

      /* Increment pointer */

      ptr++;
    }

  ptr = (long)RDB[RES_TOT_ITER_NUC_ABSRATE];

  /* Add absorption rate to buffer */

  AddBuf1D(flx*g*absxs, wgt, ptr, id, 0);

  /***************************************************************************/
}

/*****************************************************************************/
