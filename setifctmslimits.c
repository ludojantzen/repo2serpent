/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : setifctmslimits.c                              */
/*                                                                           */
/* Created:       2018/11/07 (VVa)                                           */
/* Last modified: 2019/03/29 (VVa)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Sets or checks TMS limits for materials linked to interface  */
/*                                                                           */
/* Comments: - Material array must contain material pointers instead of      */
/*             pointers to material name.                                    */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SetIFCTMSLimits:"

/*****************************************************************************/

void SetIFCTMSLimits(long ifc, long updateT)
{
  long marr, mat, nmat, n;
  double Tmin, Tmax;

  /***********************************************************************/

  /* Check interface pointer */

  CheckPointer(FUNCTION_NAME, "(ifc)", DATA_ARRAY, ifc);

  /* Get material array pointer */

  marr = (long)RDB[ifc + IFC_PTR_MAT_ARR];
  CheckPointer(FUNCTION_NAME, "(marr)", DATA_ARRAY, marr);

  /* Get interface maximum and minimum temperature */

  Tmax = RDB[ifc + IFC_MAX_TEMP];
  Tmin = RDB[ifc + IFC_MIN_TEMP];

  /* Check maximum temperature and return if temperature is not */
  /* given */

  if (Tmax == 0)
    return;

  /* Override max and min with temperature limits */

  /*
  ...
  ...
  */

  /* Get number of interface materials */

  nmat = (long)RDB[ifc + IFC_N_MAT];

  if (updateT == 0)
    {
      /* Set TMS limits for materials */

      for (n = 0; n < nmat; n++)
        {
          /* Get pointer to material */

          mat = (long)RDB[marr + n];
          CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

          /* Update upper TMS bound */

          if (Tmax > RDB[mat + MATERIAL_TMS_TMAX])
            WDB[mat + MATERIAL_TMS_TMAX] = Tmax;

          /* Update lower TMS bound */

          if (Tmin < RDB[mat + MATERIAL_TMS_TMIN])
            WDB[mat + MATERIAL_TMS_TMIN] = Tmin;

          /* Set on-the-fly Doppler-broadening mode */
          /* TMS is set on if the TMS-limits differ */
          /* TMS is always set on for mixtures */

          if ((RDB[mat + MATERIAL_TMS_TMIN] <
               RDB[mat + MATERIAL_TMS_TMAX]) ||
              ((long)RDB[mat + MATERIAL_PTR_MIX] > VALID_PTR))
            {
              /* Check that material doppler temperature is not already */
              /* set by tmp-card */

              if (RDB[mat + MATERIAL_DOPPLER_TEMP] >= 0)
                Error(mat, "Material temperature set by tmp-card but a "
                      "temperature distribution is also given by "
                      "interface %s", GetText(ifc + IFC_PTR_INPUT_FNAME));

              /* Set TMS mode on */

              WDB[mat + MATERIAL_TMS_MODE] = (double)YES;

              /* Set Doppler-preprocessor off */

              WDB[mat + MATERIAL_DOPPLER_TEMP] = -1.0;
            }
          else
            {
              /* No variation in temperature, Doppler preprocessor is enough */

              WDB[mat + MATERIAL_DOPPLER_TEMP] = Tmin;

              /* Set Doppler preprocessor on */

              WDB[DATA_USE_DOPPLER_PREPROCESSOR] = (double)YES;
            }

          /* Next interface material */
        }
    }
  else
    {
      /* Temperature field is being updated */

      for (n = 0; n < nmat; n++)
        {
          /* Get pointer to current material */

          mat = (long)RDB[marr + n];
          CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

          /* Check TMS upper bound */

          if (Tmax > RDB[mat + MATERIAL_TMS_TMAX])
            Die(FUNCTION_NAME, "Material temperature %f above TMS majorant "
                "for material %s", Tmax, GetText(mat + MATERIAL_PTR_NAME));

          /* Check TMS lower bound */

          if (Tmin < RDB[mat + MATERIAL_TMS_TMIN])
            Die(FUNCTION_NAME, "Material temperature %f below TMS minorant "
                "for material %s", Tmin, GetText(mat + MATERIAL_PTR_NAME));
        }
    }

  return;
}
