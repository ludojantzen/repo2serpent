/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : nubar.c                                        */
/*                                                                           */
/* Created:       2011/02/05 (JLe)                                           */
/* Last modified: 2011/10/28 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Calculates fission nubar                                     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "Nubar:"

/*****************************************************************************/

double Nubar(long loc0, double E, long id)
{
  long ptr, erg, np, n, i;
  double nu, f, nu0, nu1, x;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Test existing data */
 
  if ((nu = TestValuePair(loc0 + NUBAR_PTR_PREV_VAL, E, id)) > -INFTY)
    return nu;

  /* Check type */

  if ((long)RDB[loc0 + NUBAR_DATA_TYPE] == 1)
    {
      /***** Polynomial data *************************************************/

      /* Pointer to data */

      ptr = (long)RDB[loc0 + NUBAR_PTR_POLY_DATA];

      /* Check pointer */

      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Number of coefficients */

      np = (long)RDB[ptr++];
      
      /* Initial value */

      nu = RDB[ptr];
      x = E;

      /* Loop over coefficients */

      for (n = 1; n < np; n++)
	{
	  nu = nu + RDB[ptr + n]*x;
	  x = x*E;
	}

      /***********************************************************************/
    }
  else
    {
      /***** Tabular data ****************************************************/

      /* Pointer to energy grid */

      erg = (long)RDB[loc0 + NUBAR_PTR_EGRID];

      /* Pointer to data */

      ptr = (long)RDB[loc0 + NUBAR_PTR_PTS];

      /* Check pointers */

      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get interpolation factor */

      if ((f = GridFactor(erg, E, id)) < 0.0)
	return 0.0;

      /* Get number of points */

      np = (long)RDB[erg + ENERGY_GRID_NE];

      /* Separate integer and decimal parts of interpolation factor */
      
      i = (long)f;
      f = f - (double)i;

      /* Check boundaries */
      
      if ((i < 0) || (i > np - 2))
	return 0.0;
      else
	{
	  /* Get tabulated points */
	  
	  nu0 = RDB[ptr + i];
	  nu1 = RDB[ptr + i + 1];

	  /* Interpolate */
	  
	  nu = f*(nu1 - nu0) + nu0;
	}
      
      /***********************************************************************/
    }

  /* Check value (delayed nubar may be very small, prompt is very large */
  /* at high energy) */

  CheckValue(FUNCTION_NAME, "nu", "", nu, 0.0, 10.0);

  /* Remember value */

  StoreValuePair(loc0 + NUBAR_PTR_PREV_VAL, E, nu, id);

  /* Return value */

  return nu;
}

/*****************************************************************************/
