/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : microxs.c                                      */
/*                                                                           */
/* Created:       2010/12/15 (JLe)                                           */
/* Last modified: 2018/10/01 (RTu)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Interpolates microscopic cross section                       */
/*                                                                           */
/* Comments: - Also check uresdilumicroxs.c when making changes              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MicroXS:"

/*****************************************************************************/

double MicroXS(long rea, double E, long id)
{
  long i, i0, ne, ptr, erg, urs, n, nr, mt, edep;
  double xs0, xs1, xs, f, rnd;

  /* Check reaction pointer */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
  
  /* Pointer to ures data */

  urs = (long)RDB[rea + REACTION_PTR_URES];

  /* Test existing data */
  
  if ((xs = TestValuePair(rea + REACTION_PTR_PREV_XS, E, id)) > -INFTY)
    {
      /* Check if ures data exists */

      if (urs < VALID_PTR)
        return xs;
      else
        {
          /* Use value only if it was obtained using same random number */

          if ((rnd = TestValuePair(urs + URES_PTR_RND_CHK, E, id)) > 0.0)
            if (rnd == TestValuePair(urs + URES_PTR_RND, E, id))
              return xs;
        }
    }

  /* Check type */

  if (urs > VALID_PTR)
    {
      /***********************************************************************/

      /***** Cross section with ures data ************************************/

      /* Get un-adjusted cross section (NOTE: UresDiluMicroXS() is very      */
      /* similar to this function, but separated for the sake of clarity)    */
      /* Non-fission KERMA URES factors and infinite dilute cross sections   */
      /* can be negative due to energy balance issues in the ENDF data       */

      xs = UresDiluMicroXS(rea, E, id);
      
      /* Get mt and edep mode */

      mt = (long)RDB[rea + REACTION_MT];
      edep = (long)RDB[DATA_EDEP_MODE];
      
      if ((edep != EDEP_MODE_NEUTRON_PHOTON) || (mt != 301))
        {
          CheckValue(FUNCTION_NAME, "xs", "", xs, 0.0, MAX_XS);
        }
      
      /* Sample factor */

      f = UresFactor(rea, E, id);

      if (mt != 301)
        {
          /* Braces needed */
          
          CheckValue(FUNCTION_NAME, "f (ures)", "", f, -10.0, INFTY);
        }

      /* Adjust cross section */
      
      if ((f > 0.0) || ((edep == EDEP_MODE_NEUTRON_PHOTON) && (mt == 301)))
        xs = f*xs;
      else
        xs = 0.0;

      /***********************************************************************/
    }
  else if (((n = (long)RDB[rea + REACTION_CACHE_OPTI_IDX]) > -1) &&
           (E < RDB[DATA_CACHE_OPTI_EMAX]))
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

  /* Remember value */

  StoreValuePair(rea + REACTION_PTR_PREV_XS, E, xs, id);

  /* Return cross section */
  
  return xs;
}

/*****************************************************************************/
