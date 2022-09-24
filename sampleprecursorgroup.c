/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : sampleprecursorgroup.c                         */
/*                                                                           */
/* Created:       2011/02/05 (JLe)                                           */
/* Last modified: 2012/02/01 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Samples delayed neutron precursor group                      */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SamplePrecursorGroup:"

/*****************************************************************************/

long SamplePrecursorGroup(long rea, double E, long id)
{
  long lst, loc0, ptr, n, i, erg, np, nuc;
  double rnd, f, P0, P1, P;

  /* Check reaction pointer */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Check energy */

  CheckValue(FUNCTION_NAME, "E", "", E, ZERO, INFTY);

  /* Pointer to precursor group data */
  
  lst = (long)RDB[rea + REACTION_PTR_PREC_LIST];
  CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);

  /* Sample random nuber */

  rnd = RandF(id);

  /* Loop over groups */
	  
  n = 0;
  while ((loc0 = ListPtr(lst, n++)) > VALID_PTR)
    {
      /* Pointer to energy grid */

      erg = (long)RDB[loc0 + PREC_PTR_EGRID];
      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

      /* Get interpolation factor */

      if ((f = GridFactor(erg, E, id)) < 0.0)
	break;

      /* Get number of points */

      np = (long)RDB[erg + ENERGY_GRID_NE];

      /* Separate integer and decimal parts of interpolation factor */
      
      i = (long)f;
      f = f - (double)i;

      /* Check boundaries */
      
      if ((i < 0) || (i > np - 2))
	break;
      else
	{
	  /* Pointer to probabilities */

	  ptr = (long)RDB[loc0 + PREC_PTR_PTS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	  /* Get tabulated points */
	  
	  P0 = RDB[ptr + i];
	  P1 = RDB[ptr + i + 1];

	  /* Interpolate */
	  
	  P = f*(P1 - P0) + P0;
	}
      
      /* Compare to probability */
      
      if ((rnd = rnd - P) < 0.0)
	return loc0;
    }

  /* Get pointer to nuclide */

  nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

  /* Unable to sample precursor group */

  Warn(FUNCTION_NAME, "Precursor group sampling failed (%s, mt %ld, E = %E)",
       GetText(nuc + NUCLIDE_PTR_NAME), (long)RDB[rea + REACTION_MT], E);

  /* Dummy return to avoid compiler warning */

  return -1;
}

/*****************************************************************************/
