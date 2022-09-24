/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : macrourescorr.c                                */
/*                                                                           */
/* Created:       2015/11/01 (JLe)                                           */
/* Last modified: 2015/11/07 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Performs ures corrections on macroscopic cross sections      */
/*              obtained from pre-calculated data.                           */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MacroUresCorr:"

/*****************************************************************************/

double MacroUresCorr(long rea, double xs, double E, long id)
{
  long mat, ptr, mt;
  double adens, f, Emin, Emax;

  /* Check Pointer */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Check mode */

  if ((long)RDB[DATA_OPTI_RECONSTRUCT_MACROXS] == NO)
    Die(FUNCTION_NAME, "Shouldn't be called in this mode");
  
  /* Check pointer to probability table data */

  if ((long)RDB[rea + REACTION_PTR_URES] > VALID_PTR)
    Die(FUNCTION_NAME, "Should not exist anymore (tän voi poistaa joskus)");

  /* Check boundaries */

  if ((E < RDB[rea + REACTION_URES_EMIN]) ||
      (E > RDB[rea + REACTION_URES_EMAX]))
    return xs;

  /* Get pointer to material */

  mat = (long)RDB[rea + REACTION_PTR_MAT];
  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);
    
  /* Get reaction mt and reset pointer (avoid compiler warning) */

  mt = (long)RDB[rea + REACTION_MT];
  ptr = -1;

  /* Get pointer to ures list  */
      
  if (mt == MT_MACRO_TOTXS)      
    ptr = (long)RDB[mat + MATERIAL_PTR_TOT_URES_LIST];
  else if (mt == MT_MACRO_ABSXS)      
    ptr = (long)RDB[mat + MATERIAL_PTR_ABS_URES_LIST];
  else if (mt == MT_MACRO_ELAXS)      
    ptr = (long)RDB[mat + MATERIAL_PTR_ELA_URES_LIST];
  else if (mt == MT_MACRO_FISSXS)      
    ptr = (long)RDB[mat + MATERIAL_PTR_FISS_URES_LIST];
  else if (mt == MT_MACRO_HEATXS)      
    ptr = (long)RDB[mat + MATERIAL_PTR_HEAT_URES_LIST];
  else if (mt == MT_MACRO_FISSE)    
    ptr = (long)RDB[mat + MATERIAL_PTR_FISS_URES_LIST];
  else if (mt == MT_MACRO_NSF)    
    ptr = (long)RDB[mat + MATERIAL_PTR_FISS_URES_LIST];
  else
    Die(FUNCTION_NAME, "Invalid mt %ld (%E %E)", mt, 
        RDB[rea + REACTION_URES_EMIN], RDB[rea + REACTION_URES_EMAX]);
      
  /* Check Pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Reset reaction pointer (rewind list) */
  
  rea = -1;
  
  /* Loop over reactions */
  
  while (NextReaction(ptr, &rea, &adens, &Emin, &Emax, id) > VALID_PTR)
    {
      /* Check reaction pointer */

      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
          
      /* Check threshold */
      
      if (E < Emin)
        {
          /* Break loop */
          
          break;
        }
      else
        {
          /* Get multiplier */
          
          f = ReaMulti(rea, mt, E, id);
          
          /* Adjust cross section */
          
          xs = xs + f*adens*MicroXS(rea, E, id);
          xs = xs - f*adens*UresDiluMicroXS(rea, E, id);
        }
    }
  
  /* Check value */

  CheckValue(FUNCTION_NAME, "xs", "", xs, -1E-15, INFTY);

  /* Return value */

  return xs;
}

/*****************************************************************************/
