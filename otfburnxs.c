/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : otfburnxs.c                                    */
/*                                                                           */
/* Created:       2018/03/25 (JLe)                                           */
/* Last modified: 2018/03/25 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Returns macroscopic total cross sections for nuclides        */
/*              involved in on-the-fly burnup solver.                        */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "OTFBurnXS:"

/*****************************************************************************/

double OTFBurnXS(long mat, double E, long mt, long id)
{
  long loc0, iso, nuc, rea;
  double sum, xs, adens, mult;

  /* Check material pointer */

  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

  /* Pointer to data */

  if ((loc0 = (long)RDB[mat + MATERIAL_PTR_OTF_BURN]) < VALID_PTR)
    return 0.0;

  /* Check predictor-corrector calculation */

  if (((long)RDB[DATA_BURN_CORR_TYPE] != CORR_TYPE_NONE) &&
      ((long)RDB[DATA_BURN_STEP_PC] != CORRECTOR_STEP))
    return 0.0;

  /* Reset sum */

  sum = 0.0;

  /* Loop over data */

  while (loc0 > VALID_PTR)
    {
      /* Check */

      if ((long)RDB[loc0 + OTF_BURN_PTR_MAT] != mat)
        break;

      /* Get pointer to composition */

      iso = (long)RDB[loc0 + OTF_BURN_PTR_COMP];
      CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

      /* Get atomic density */
          
      if ((adens = RDB[iso + COMPOSITION_ADENS]) == 0.0)
        {
          /* Pointer to next */

          loc0 = NextItem(loc0);

          /* Cycle loop */

          continue;
        }
      
      /* Pointer to nuclide */

      nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Check for inelastic production */

      if (mt == MT_MACRO_INLPRODXS)
        {
          /* Loop over reactions */

          rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
          while (rea > VALID_PTR)
            {
              /* Check multiplication */

              if (RDB[rea + REACTION_WGT_F] > 1.0)
                {
                  /* Get microscopic cross section */

                  xs = MicroXS(rea, E, id);
          
                  /* Get multiplier */
          
                  mult = ReaMulti(rea, mt, E, id);
                  
                  /* Add to sum */
                  
                  sum = sum + xs*adens*mult;
                }
              
              /* Next */

              rea = NextItem(rea);
            }
        }
      else
        {
          /* Avoid compiler warning */

          rea = -1;
          
          /* Get pointer to reaction */
          
          if (mt == MT_MACRO_TOTXS)
            rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
          else if (mt == MT_MACRO_ELAXS)
            rea = (long)RDB[nuc + NUCLIDE_PTR_ELAXS];
          else if (mt == MT_MACRO_ABSXS)
            rea = (long)RDB[nuc + NUCLIDE_PTR_SUM_ABSXS];
          else if (mt == MT_MACRO_FISSXS)
            rea = (long)RDB[nuc + NUCLIDE_PTR_FISSXS];
          else if (mt == MT_MACRO_FISSE)
            rea = (long)RDB[nuc + NUCLIDE_PTR_FISSXS];
          else if (mt == MT_MACRO_NSF)
            rea = (long)RDB[nuc + NUCLIDE_PTR_FISSXS];
          else if (mt == MT_MACRO_TMP_MAJORANTXS)
            {
              /* Get pointer to reaction data */
              
              rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
              CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
              
              /* Get pointer to majorant (tota ei välttämättä ole (?)) */
              
              if ((rea = (long)RDB[rea + REACTION_PTR_TMP_MAJORANT]) 
                  < VALID_PTR)
                Die(FUNCTION_NAME, "Temperature majorant did not exist");
            }
          else
            Die(FUNCTION_NAME, "Invalid reaction mode %ld", mt);

          /* Check pointer */

          if (rea > VALID_PTR)
            {
              /* Get microscopic cross section */

              xs = MicroXS(rea, E, id);
              
              /* Get multiplier */
              
              mult = ReaMulti(rea, mt, E, id);
              
              /* Add to sum */
          
              sum = sum + xs*adens*mult;
            }
        }

      /* Next reaction */

      loc0 = NextItem(loc0);      
    }
  
  /* Return total */

  return sum;
}

/*****************************************************************************/
