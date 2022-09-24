/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : majorantxs.c                                   */
/*                                                                           */
/* Created:       2011/01/07 (JLe)                                           */
/* Last modified: 2015/11/01 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Interpolates majorant cross section with extras              */
/*                                                                           */
/* Comments: - TODO: cutoff to MIN_MACRO_XS?                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MajorantXS:"

/*****************************************************************************/

double MajorantXS(long rea, double E, long id)
{
  long i, ne, ptr, erg;
  double xs0, xs1, xs, f;

  /* Check Pointer */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Test existing data */
  
  if ((xs = TestValuePair(rea + REACTION_PTR_PREV_XS, E, id)) > 0.0)
    return xs;

  /* Check for coarse multi-group majorant */

  if ((ptr = (long)RDB[rea + REACTION_PTR_MGXS]) > VALID_PTR)
    {
      /***********************************************************************/

      /***** Multi-group majorant ********************************************/

      /* Get pointer to energy grid */

      erg = (long)RDB[DATA_COARSE_MG_PTR_GRID];
      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

      /* Find index */

      i = GridSearch(erg, E);

      /* Get value */

      if (i < 0)
        xs = 0.0;
      else
        xs = RDB[ptr + i];

      /***********************************************************************/
    }
  else
    {
      /***********************************************************************/

      /***** Continuous-energy majorant **************************************/

      /* Get pointer to energy grid */

      erg = (long)RDB[rea + REACTION_PTR_EGRID];
      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);
      
      /* Get interpolation factor */
      
      if ((f = GridFactor(erg, E, id)) < 0)
        xs = 0.0;
      else
        {
          /* Check interpolation factor */
          
          CheckValue(FUNCTION_NAME, "f", "", f, 0.0, MAX_EGRID_NE);
          
          /* Separate integer and decimal parts of interpolation factor */
          
          i = (long)f;
          f = f - (double)i;
          
          /* Get number of points */
          
          ne = (long)RDB[rea + REACTION_XS_NE];
          CheckValue(FUNCTION_NAME, "ne", "", ne, 2, MAX_EGRID_NE);
          
          /* Check boundaries */
          
          if ((i < 0) || (i > ne - 1))
            Die(FUNCTION_NAME, "Energy outside grid boundaries");
          else
            {      
              /* Get pointer to data */
              
              ptr = (long)RDB[rea + REACTION_PTR_MAJORANT_XS];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
              
              /* Get tabulated cross sections */
          
              xs0 = RDB[ptr + i];
              xs1 = RDB[ptr + i + 1];
              
              /* Interpolate normally if TMS is not used */
             
              if ((long)RDB[DATA_TMS_MODE] == TMS_MODE_NONE) 
                {              
                  if (i == ne - 1)
                    xs = (1.0 - f)*xs0;
                  else
                    xs = f*(xs1 - xs0) + xs0;
                }

              /* DT Majorant cannot be interpolated in TMS mode,      */
              /* because histogram-type TMS majorants are used in the */
              /* majorant generation and they are not interpolated.   */
              /* --> use maximum of adjascent points instead          */

              else 
                {
                  if (xs0 > xs1)
                    xs = xs0;
                  else
                    xs = xs1;
                }
                  }
        }
      
      /* Check pointer to ures data */

      if ((long)RDB[rea + REACTION_PTR_URES] > VALID_PTR)
        Die(FUNCTION_NAME, "Majorant has ures data!?!?!");

      /* TODO: Ota tää pätkä pois. Jäänyt jostain kummittelemaan */

#ifdef mmmmmmmmmmmmmmm

      /* Remember non-adjusted value */
      
      StoreValuePair(rea + REACTION_PTR_PREV_XS0, E, xs, id);
      
      /* Sample unresolved resonance probability table data (ja mihinköhän */
      /* ton sit tarvii?) */

      f = UresFactor(rea, E, id);

      /* Tarkistetaan varmuuden vuoksi */

      if (f != 1.0)
        Die(FUNCTION_NAME, "f = %E\n", f);

      /* Check value (may be zero) */
      
      CheckValue(FUNCTION_NAME, "f", "", f, -10.0, INFTY);
      
      /* Adjust cross section */
      
      xs = f*xs;

#endif
     
      /***********************************************************************/
    }

  /* Remember value */
  
  StoreValuePair(rea + REACTION_PTR_PREV_XS, E, xs, id);

  /* Return cross section */

  return xs;
}

/*****************************************************************************/
