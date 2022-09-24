/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : photonmicroxs.c                                */
/*                                                                           */
/* Created:       2011/04/18 (JLe)                                           */
/* Last modified: 2011/10/28 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Interpolates microscopic photon cross section                */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PhotonMicroXS:"

/*****************************************************************************/

double PhotonMicroXS(long rea, double E, long id)
{
  long i, ne, ptr, erg;
  double xs0, xs1, xs, f;

  /* Check Pointer */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Test existing data */
  
  if ((xs = TestValuePair(rea + REACTION_PTR_PREV_XS, E, id)) > -INFTY)
    return xs;

  /* Get pointer to energy grid */
      
  erg = (long)RDB[rea + REACTION_PTR_EGRID];
  
  /* Check pointer */
  
  CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);
  
  /* Get interpolation factor */
  
  if ((f = GridFactor(erg, E, id)) < 0)
    xs = 0.0;
  else
    {
      /* Check interpolation factor */
      
      CheckValue(FUNCTION_NAME, "f (interp)", "", f, 0.0, MAX_EGRID_NE);
      
      /* Separate integer and decimal parts of interpolation factor */
      
      i = (long)f;
      f = f - (double)i;
      
      /* Get number of points */
      
      ne = (long)RDB[rea + REACTION_XS_NE];
      
      /* Check value */
      
      CheckValue(FUNCTION_NAME, "ne", "", ne, 2, MAX_EGRID_NE);
	  
      /* Check boundaries */
      
      if ((i < 0) || (i > ne - 2))
	xs = 0.0;
      else
	{      
	  /* Get pointer to data */
	  
	  ptr = (long)RDB[rea + REACTION_PTR_XS];
	  
	  /* Check pointer */
	  
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      
	  /* Get tabulated cross sections */
	  
	  xs0 = RDB[ptr + i];
	  xs1 = RDB[ptr + i + 1];
	  
	  /* Interpolate */
	  
	  xs = f*(xs1 - xs0) + xs0;
	}
    }

  /* Remember value */

  StoreValuePair(rea + REACTION_PTR_PREV_XS, E, xs, id);

  /* Return cross section */
  
  return xs;
}

/*****************************************************************************/
