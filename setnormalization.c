/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : setnormalization.c                             */
/*                                                                           */
/* Created:       2013/06/24 (JLe)                                           */
/* Last modified: 2016/01/24 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Sets normalization for depletion interval                    */
/*                                                                           */
/* Comments: - Toi aktivaatiohomma on vähän viritelmä. Normeerauskertoimia   */
/*             skaalaamalla saadaan vaikutus vietyä palamamatriisiin asti.   */
/*             Nolla joudutaan korvaamaan ZERO:lla.                          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SetNormalization:"

/*****************************************************************************/

void SetNormalization(long dep)
{
  long mat0, mat, ptr, norm, loc0, loc1, n;
  double f, f0;

  /***************************************************************************/

  /***** Adjust flux for activation ******************************************/

  /* Check type */

  if (((long)RDB[dep + DEP_HIS_STEP_TYPE] == DEP_STEP_ACT_STEP) || 
      ((long)RDB[dep + DEP_HIS_STEP_TYPE] == DEP_STEP_ACT_TOT))
    {
      /* Pointer to previous */

      ptr = (long)RDB[DATA_BURN_ACTI_PTR_DEP];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Pointer to normalizations */

      loc0 = (long)RDB[ptr + DEP_HIS_PTR_NORM];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      loc1 = (long)RDB[dep + DEP_HIS_PTR_NORM];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY,loc1);

      /* Get previous factor */

      f0 = RDB[DATA_BURN_ACTI_PREV_F];
      CheckValue(FUNCTION_NAME, "f0", "", f0, ZERO, INFTY);

      /* Avoid compiler warning */

      f = -1.0;

      /* Loop over options */

      for (n = NORM_POWER; n < NORM_SFRATE + 1; n++)
        {
          if (RDB[loc0 + n] == 0.0)
            Die(FUNCTION_NAME, "Division by zero");
          else if ((RDB[loc0 + n] > 0.0) && (RDB[loc1 + n] >= 0.0))
            {
              /* Calculate factor */
              
              if (RDB[loc1 + n] == 0.0)
                f = ZERO;
              else
                f = RDB[loc1 + n]/RDB[loc0 + n];

              /* Break loop */
              
              break;
            }
        }

      /* Check */

      if (f == -1.0)
        Error(dep, "Depletion intervals must be similarly normalized");

      /* Store previous factor */

      if (f > 0.0)
        WDB[DATA_BURN_ACTI_PREV_F] = f;
      else
        Die(FUNCTION_NAME, "Error in factor (%E)", f);
      
      /* Loop over materials */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          /* Adjust flux */

          WDB[mat + MATERIAL_BURN_FLUX_PS1] = 
            RDB[mat + MATERIAL_BURN_FLUX_PS1]*f/f0;
          
          WDB[mat + MATERIAL_BURN_FLUX_BOS] =
            RDB[mat + MATERIAL_BURN_FLUX_BOS]*f/f0;

          WDB[mat + MATERIAL_BURN_FLUX_EOS] =
            RDB[mat + MATERIAL_BURN_FLUX_EOS]*f/f0;

          /* Next material */

          mat = NextItem(mat);
        }
    }
  else if (((long)RDB[dep + DEP_HIS_STEP_TYPE] != DEP_STEP_DEC_STEP) && 
           ((long)RDB[dep + DEP_HIS_STEP_TYPE] != DEP_STEP_DEC_TOT))
    {
      /* Reset previous factor and set pointer to previous interval */
      
      WDB[DATA_BURN_ACTI_PTR_DEP] = (double)dep;
      WDB[DATA_BURN_ACTI_PREV_F] = 1.0;
    }

  /***************************************************************************/

  /***** Set normalization ***************************************************/

  /* NOTE: toi norm-rakenne ei tämän jälkeen välttämättä toimi enää */
  /* linkitettynä listana, mutta sillä ei pitäisi olla väliä koska  */
  /* sen yli ei loopata enää missään. */
  
  WDB[DATA_PTR_NORM] = RDB[dep + DEP_HIS_PTR_NORM];

  /* Get pointer to normalization */

  if ((norm = (long)RDB[dep + DEP_HIS_PTR_NORM]) < VALID_PTR)
    return;

  /* Get material pointer */

  if ((mat0 = (long)RDB[norm + NORM_PTR_MAT]) < VALID_PTR)
    return;

  /* Reset pointers in materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Reset pointer */

      WDB[mat + MATERIAL_PTR_NORM] = NULLPTR;

      /* Next material */

      mat = NextItem(mat);
    }

  /* Set pointers in materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check pointer */

      if (mat == mat0)
        WDB[mat + MATERIAL_PTR_NORM] = (double)norm;

      /* Check parent */

      if ((ptr = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
        if (ptr == mat0)
          WDB[mat + MATERIAL_PTR_NORM] = (double)norm;

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/
}

/*****************************************************************************/
