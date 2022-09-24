/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : poisonxs.c                                     */
/*                                                                           */
/* Created:       2012/12/05 (JLe)                                           */
/* Last modified: 2018/02/05 (VVa)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Returns macroscopic total fission product poison xs in       */
/*              equilibrium Xe-135 / Sm-149 calculation mode                 */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PoisonXS:"

/*****************************************************************************/

double PoisonXS(long mat, double E, long mt, long id)
{
  long iso, nuc, rea;
  double abs, xs, adens;

  /* Check material pointer */

  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

  /* Reset cross section */

  abs = 0.0;

  /* Avoid compiler warning */

  rea = -1;

  /* Check equilibrium Xe-135 flag */

  if ((long)RDB[mat + MATERIAL_XENON_EQUIL_CALC] == YES)
    {
      /* Get pointer to Xe-135 isotope */

      iso = (long)RDB[mat + MATERIAL_PTR_XE135_ISO];
      CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

      /* Pointer to nuclide */

      nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Get pointer to reaction */

      if (mt == MT_MACRO_TOTXS)
        rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
      else if (mt == MT_MACRO_ABSXS)
        rea = (long)RDB[nuc + NUCLIDE_PTR_NGAMMAXS];
      else if (mt == MT_MACRO_TMP_MAJORANTXS)
        {
          /* Get pointer to reaction data */

          rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

          /* Get pointer to majorant (tota ei välttämättä ole (?)) */

          if ((rea = (long)RDB[rea + REACTION_PTR_TMP_MAJORANT]) < VALID_PTR)
            Die(FUNCTION_NAME, "Temperature majorant did not exist");

        }
      else
        Die(FUNCTION_NAME, "Invalid reaction mode");

      /* Check pointer */

      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

      /* Get microscopic cross section */

      xs = MicroXS(rea, E, id);

      /* Get atomic density */

      adens = RDB[iso + COMPOSITION_ADENS];

      /* Add tot total */

      abs = abs + xs*adens;
    }

  /* Check equilibrium Sm-149 flag */

  if ((long)RDB[mat + MATERIAL_SAMARIUM_EQUIL_CALC] == YES)
    {
      /* Get pointer to Sm-149 isotope */

      iso = (long)RDB[mat + MATERIAL_PTR_SM149_ISO];
      CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

      /* Pointer to nuclide */

      nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Get pointer to reaction */

      if (mt == MT_MACRO_TOTXS)
        rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
      else if (mt == MT_MACRO_ABSXS)
        rea = (long)RDB[nuc + NUCLIDE_PTR_NGAMMAXS];
      else
        Die(FUNCTION_NAME, "Invalid reaction mode");

      /* Check pointer */

      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

      /* Get microscopic cross section */

      xs = MicroXS(rea, E, id);

      /* Get atomic density */

      adens = RDB[iso + COMPOSITION_ADENS];

      /* Add tot total */

      abs = abs + xs*adens;
    }

  /* Return total */

  return abs;
}

/*****************************************************************************/
