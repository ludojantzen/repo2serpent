/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : microxs.c                                      */
/*                                                                           */
/* Created:       2013/12/17 (TVi)                                           */
/* Last modified: 2015/11/01 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Temperature majorant cross sections cannot be interpolated,  */
/*              this function returns a value from the histogram-type        */
/*                temperature majorant xs.                                     */
/*                                                                           */
/* Comments:  - Tan varmaan vois yhdistaa myos microxssaan jollakin          */
/*              if-lauseella, mutta ajattelin etta erillinen funktio         */
/*              voisi olla vahan selkeampi.                                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MicroMajorantXS:"

/*****************************************************************************/

double MicroMajorantXS(long rea, double E, long id)
{
  long i, i0, ne, ptr, erg, n, nr;
  double xs, f;

  /* Check reaction pointer */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Check for cache-optimized region */

  if ((E < RDB[DATA_CACHE_OPTI_EMAX]) && 
      ((n = (long)RDB[rea + REACTION_CACHE_OPTI_IDX]) > -1))
    {

      /***********************************************************************/

      /***** Indexing by reaction ********************************************/

      /* Get pointer to unionizedenergy grid */

      erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID];
      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);
     
      /* Get interpolation factor */

      if ((f = GridFactor(erg, E, id)) < 0.0)
        xs = 0.0;
      else
        {                    
          /* Get pointer to data */
          
          ptr = (long)RDB[DATA_PTR_CACHE_OPTI_XS];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          
          /* Get number of reactions */
          
          nr = (long)RDB[DATA_CACHE_OPTI_NREA];
          CheckValue(FUNCTION_NAME, "nr", "", n, 0, nr - 1);
                    
          /* Get cross sections */
          
          i = (long)f;
          CheckValue(FUNCTION_NAME, "i", "", i, 0, MAX_EGRID_NE);
          
          /* Get cross section (histogram majorant) */
          
          xs = RDB[ptr + i*nr + n];
        }
      
      /***********************************************************************/
    }
  else
    {            
      /***********************************************************************/
      
      /***** Indexing by energy **********************************************/

      /* Test existing data */
      
      if ((xs = TestValuePair(rea + REACTION_PTR_PREV_XS, E, id)) > -INFTY)
        return xs;      
           
      /* Get pointer to energy grid */

      erg = (long)RDB[rea + REACTION_PTR_EGRID];
      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

      /* Get interpolation factor */

      if ((f = GridFactor(erg, E, id)) < 0)
        xs = 0.0;
      else
        {          
          i = (long)f;
          CheckValue(FUNCTION_NAME, "i", "", i, 0, MAX_EGRID_NE);

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

              /* Get xs (histogram majorant) */ 
              
              xs = RDB[ptr + i];
            }
        }
      
      /***********************************************************************/
    }

   /* Check if ures data is given */
  
   if ((long)RDB[rea + REACTION_PTR_URES] > VALID_PTR)
     Die(FUNCTION_NAME, "Ures sampling in use");
  
   /* Remember non-adjusted and final value */
   /* JLe 1.11.2015 / 2.1.25: tässä oli tallennus myös PREV_XS0, */
   /* mutta se liittyy ures-sämpläykseen eikä pitäisi olla nyt väliä */

   StoreValuePair(rea + REACTION_PTR_PREV_XS, E, xs, id);
   
   /* Return cross section */
  
   return xs;
}

/*****************************************************************************/
