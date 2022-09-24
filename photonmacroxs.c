/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : photonmacroxs.c                                */
/*                                                                           */
/* Created:       2011/01/02 (JLe)                                           */
/* Last modified: 2016/04/14 (JLe)                                           */
/* Version:       2.1.26                                                     */
/*                                                                           */
/* Description: Interpolates macroscopic photon cross section                */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PhotonMacroXS:"

/*****************************************************************************/

double PhotonMacroXS(long rea0, double E, long id)
{
  long i, ptr, erg, ne;
  double xs0, xs1, xs, f;

  long rea, mt;
  double adens, Emin, Emax, mult;

  /* Check Pointer */

  CheckPointer(FUNCTION_NAME, "(rea0)", DATA_ARRAY, rea0);

  /* Test existing data */
      
  if ((xs = TestValuePair(rea0 + REACTION_PTR_PREV_XS, E, id)) > -INFTY)
    return xs;

  /* Reset cross section */

  xs = 0.0;

  /* Get pointer to data */

  if ((ptr = (long)RDB[rea0 + REACTION_PTR_XS]) > VALID_PTR)
    {
      /***********************************************************************/

      /***** Interpolate pre-calculated data *********************************/

      /* Get pointer to energy grid */

      erg = (long)RDB[rea0 + REACTION_PTR_EGRID];

      /* Check pointer */

      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

      /* Get interpolation factor */

      if ((f = GridFactor(erg, E, id)) < 0)
        xs = 0.0;
      else
        {
          /* Check interpolation factor */

          CheckValue(FUNCTION_NAME, "f (1)", "", f, 0.0, MAX_EGRID_NE);

          /* Separate integer and decimal parts of interpolation factor */
      
          i = (long)f;
          f = f - (double)i;
      
          /* Get number of points */
      
          ne = (long)RDB[rea0 + REACTION_XS_NE];
      
          /* Check value */
      
          CheckValue(FUNCTION_NAME, "ne", "", ne, 2, MAX_EGRID_NE);
      
          /* Check boundaries */
      
          if ((i < 0) || (i > ne - 2))
            xs = 0.0;
          else
            {      
              /* Get tabulated cross sections */
              
              xs0 = RDB[ptr + i];
              xs1 = RDB[ptr + i + 1];
              
              /* Interpolate */
              
              xs = f*(xs1 - xs0) + xs0;
            }
        }

      /***********************************************************************/
    }
  else if ((ptr = (long)RDB[rea0 + REACTION_PTR_PARTIAL_LIST]) > VALID_PTR)
    {
      /***********************************************************************/

      /***** Calculate sum of partials ***************************************/

      /* Reset cross section */
  
      xs = 0.0;

      /* Get mt */

      mt = (long)RDB[rea0 + REACTION_MT];

      /* Reset reaction pointer (rewind list) */
  
      rea = -1;
      
      /* Loop over reactions */

      while (NextReaction(ptr, &rea, &adens, &Emin, &Emax, id) > VALID_PTR)
        {
          /* Check reaction pointer */

          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

          /* Get multiplier */

          mult = ReaMulti(rea, mt, E, id);

          /* Add to cross section */

          xs = xs + mult*adens*PhotonMicroXS(rea, E, id);

          /* Check energy cut-off */

          if (E < Emin)
            break;
        }

      /***********************************************************************/
    }
  else
    Die(FUNCTION_NAME, "No summation xs or partial list");

  /* Store cross section */

  if (xs > 0.0)
    StoreValuePair(rea0 + REACTION_PTR_PREV_XS, E, xs, id);

  /* Return interpolated value */

  return xs;
}

/*****************************************************************************/
