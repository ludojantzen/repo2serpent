/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : resetpoisonconc.c                              */
/*                                                                           */
/* Created:       2012/12/05 (JLe)                                           */
/* Last modified: 2012/12/05 (JLe)                                           */
/* Version:       2.1.10                                                     */
/*                                                                           */
/* Description: Resets concentrations of Xe-135 and Sm-149 for equilibrium   */
/*              iteration                                                    */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ResetPoisonConc:"

/*****************************************************************************/

void ResetPoisonConc()
{
  long mat, iso;
    
  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Get pointer to Xe-135 and reset */

      if ((long)RDB[mat + MATERIAL_XENON_EQUIL_CALC] == YES)
        {
          if ((iso = (long)RDB[mat + MATERIAL_PTR_XE135_ISO]) > VALID_PTR)
            WDB[iso + COMPOSITION_ADENS] = 0.0;
          else
            Die(FUNCTION_NAME, "No pointer to Xe-135");
        }

      /* Get pointer to Sm-149 and reset */

      if ((long)RDB[mat + MATERIAL_SAMARIUM_EQUIL_CALC] == YES)
        {
          if ((iso = (long)RDB[mat + MATERIAL_PTR_SM149_ISO]) > VALID_PTR)
            WDB[iso + COMPOSITION_ADENS] = 0.0;
          else
            Die(FUNCTION_NAME, "No pointer to Sm-149");
        }

      /* Next material */
      
      mat = NextItem(mat);
    }

  /***************************************************************************/
}

/*****************************************************************************/
