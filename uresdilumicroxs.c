/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : uresdilumicroxs.c                              */
/*                                                                           */
/* Created:       2015/11/01 (JLe)                                           */
/* Last modified: 2015/11/01 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Returns un-adjusted cross section used for handling ures     */
/*              probability table sampling                                   */
/*                                                                           */
/* Comments: - REACTION_PTR_PREV_URES_XS0 should not be handled outside      */
/*             this function                                                 */
/*                                                                           */
/*           - Also check microxs.c when making changes                      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "UresDiluMicroXS:"

/*****************************************************************************/

double UresDiluMicroXS(long rea, double E, long id)
{
  long i, i0, ne, ptr, erg, n, nr;
  double xs0, xs1, xs, f;

  /* Check reaction pointer */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Check pointer to ures data */

  if ((long)RDB[rea + REACTION_PTR_URES] < VALID_PTR)
    Die(FUNCTION_NAME, "No ures data (%ld)", (long)RDB[rea + REACTION_MT]);
  
  /* Test existing data */

  if ((xs = TestValuePair(rea + REACTION_PTR_PREV_URES_XS0, E, id)) > -INFTY)
    return xs;
  
  /* Check for cache-optimized region */

  if ((E < RDB[DATA_CACHE_OPTI_EMAX]) && 
      ((n = (long)RDB[rea + REACTION_CACHE_OPTI_IDX]) > -1))
    {
      /***********************************************************************/

      /***** Indexing by reaction ********************************************/

      /* Get pointer to unionized energy grid */

      erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID];
      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

      /* Get interpolation factor */

      if ((f = GridFactor(erg, E, id)) < 0.0)
        xs = 0.0;
      else
        {
          /* Separate integer and decimal parts of interpolation factor */
      
          i = (long)f;
          f = f - (double)i;
          
          /* Check values */
          
          CheckValue(FUNCTION_NAME, "i", "", i, 0, MAX_EGRID_NE);
          CheckValue(FUNCTION_NAME, "f", "", f, 0.0, 1.0);
          
          /* Get pointer to data */
          
          ptr = (long)RDB[DATA_PTR_CACHE_OPTI_XS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          
          /* Get number of reactions */
          
          nr = (long)RDB[DATA_CACHE_OPTI_NREA];
          
          /* Check reaction index */
          
          CheckValue(FUNCTION_NAME, "n", "", n, 0, nr - 1);
          
          /* Get cross sections */
          
          xs0 = RDB[ptr + i*nr + n];
          xs1 = RDB[ptr + (i + 1)*nr + n];      
          
          /* Interpolate */
          
          xs = f*(xs1 - xs0) + xs0;
        }

      /***********************************************************************/
    }
  else
    {
      /***********************************************************************/
      
      /***** Indexing by energy **********************************************/
            
      /* Get pointer to energy grid */

      erg = (long)RDB[rea + REACTION_PTR_EGRID];
      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

      /* Get interpolation factor */

      if ((f = GridFactor(erg, E, id)) < 0)
        xs = 0.0;
      else
        {
          /* Separate integer and decimal parts of interpolation factor */

          i = (long)f;
          f = f - (double)i;

          /* Check values */

          CheckValue(FUNCTION_NAME, "i", "", i, 0, MAX_EGRID_NE);
          CheckValue(FUNCTION_NAME, "f", "", f, 0.0, 1.0);
          
          /* Get first energy point and number of points */
          
          i0 = (long)RDB[rea + REACTION_XS_I0];
          ne = (long)RDB[rea + REACTION_XS_NE];
          
          /* Check values */
          
          CheckValue(FUNCTION_NAME, "i0", "", i0, 0, MAX_EGRID_NE);
          CheckValue(FUNCTION_NAME, "ne", "", ne, 2, MAX_EGRID_NE);
      
          /* Get relative index */
      
          i = i - i0;
      
          /* Check boundaries */
          
          if ((i < 0) || (i > ne - 1))
            xs = 0.0;
          else
            { 
              /* Get pointer to data */

              ptr = (long)RDB[rea + REACTION_PTR_XS];
              CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          
              /* Get tabulated cross sections */
          
              xs0 = RDB[ptr + i];
              xs1 = RDB[ptr + i + 1];
              
              /* Interpolate */
          
              if (i == ne - 1)
                xs = (1.0 - f)*xs0;
              else
                xs = f*(xs1 - xs0) + xs0;
            }
        }

      /***********************************************************************/
    }

  /* Remember non-adjusted value */

  if ((long)RDB[rea + REACTION_PTR_URES] > VALID_PTR)
    StoreValuePair(rea + REACTION_PTR_PREV_URES_XS0, E, xs, id);

  /* Return cross section */
  
  return xs;
}

/*****************************************************************************/
