/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : linksabdata.c                                  */
/*                                                                           */
/* Created:       2011/01/21 (JLe)                                           */
/* Last modified: 2019/07/24 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: Links materials to S(a,b) data                               */
/*                                                                           */
/* Comments: - Temperature interpolation is linear. Check if sqrt(T) is      */
/*             better (see work by Trumbull et. al.).                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "LinkSabData:"

/*****************************************************************************/

void LinkSabData()
{
  long mat, ptr, loc0, loc1, loc2, iso, nuc, ace;
  double T, T1, T2;

  /***************************************************************************/

  /***** Link data to materials **********************************************/

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check divisor */

      if ((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] > VALID_PTR)
        {
          /* Next material */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Loop over S(a,b) data */

      ptr = (long)RDB[mat + MATERIAL_PTR_SAB];
      while (ptr > VALID_PTR)
        {
          /* Pointer to S(a,b) definitions */

          if ((loc0 = (long)RDB[DATA_PTR_T0]) < VALID_PTR)
            Error(mat, "S(a,b) library %s is not defined",
                  GetText(ptr + THERM_PTR_ALIAS));

          /* Find match */

          if ((loc0 = SeekListStr(loc0, THERM_PTR_ALIAS,
                                  GetText(ptr + THERM_PTR_ALIAS)))
              > VALID_PTR)
            {
              /* Check that ZA is found in nuclide composition */

              iso = (long)RDB[mat + MATERIAL_PTR_COMP];
              while (iso > VALID_PTR)
                {
                  /* Pointer to nuclide data */

                  nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
                  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

                  /* Compare ZA */

                  if (RDB[ptr + THERM_ZA] == RDB[nuc + NUCLIDE_ZA])
                    {
                      /* Put pointer */

                      WDB[ptr + THERM_PTR_COMP] = (double)iso;

                      /* Remember entry */

                      WDB[nuc + NUCLIDE_PTR_ORIG_THERM] = (double)ptr;

                      /* Break loop */

                      break;
                    }

                  /* Next nuclide */

                  iso = NextItem(iso);
                }

              /* Check pointer */

              if (iso < VALID_PTR)
                Error(mat, "S(a,b) ZA %ld not found in composition",
                      (long)RDB[ptr + THERM_ZA]);

              /* Put pointer */

              WDB[ptr + THERM_PTR_THERM] = (double)loc0;

              /* Check that ZA is not already assigned */

              if ((RDB[loc0 + THERM_ZA] > 0.0) &&
                  (RDB[loc0 + THERM_ZA] != RDB[ptr + THERM_ZA]))
                Error(mat, "S(a,b) library %s is already assigned to ZA = %ld",
                      GetText(ptr + THERM_PTR_ALIAS),
                      (long)RDB[loc0 + THERM_ZA]);

              /* Put ZA */

              WDB[loc0 + THERM_ZA] = RDB[ptr + THERM_ZA];

              /* Set used-flag */

              SetOption(loc0 + THERM_OPTIONS, OPT_USED);

              /* Set otf temperatures if in otf S(a,b) interpolation mode */

              if ((long)RDB[loc0 + THERM_INTERP_MODE] == THERM_INTERP_OTF)
                {
                  /* Logiikka: Jos interface on käytössä, otetaan minimi- */
                  /* ja maksimilämpö suoraan materiaalista THERM-         */
                  /* rakenteeseen. */

                  if (RDB[mat + MATERIAL_TMS_TMAX] > 0.0)
                  {

                    if(RDB[mat + MATERIAL_TMS_TMAX] > RDB[loc0 + THERM_OTF_MAX_TEMP])
                      WDB[loc0 + THERM_OTF_MAX_TEMP] =
                        RDB[mat + MATERIAL_TMS_TMAX];

                    if( RDB[loc0 + THERM_OTF_MIN_TEMP] == 0.0 ||
                       RDB[mat + MATERIAL_TMS_TMIN] < RDB[loc0 + THERM_OTF_MIN_TEMP]){
                      WDB[loc0 + THERM_OTF_MIN_TEMP] =
                        RDB[mat + MATERIAL_TMS_TMIN];
                    }

                  }
                else
                  Error(mat,
                        "OTF S(a,b) interpolation only available in TMS mode");

                /*
                  else if(RDB[mat + MATERIAL_DOPPLER_TEMP] > 0.0){
                  WDB[loc0 + THERM_OTF_MIN_TEMP] = RDB[mat + MATERIAL_DOPPLER_TEMP];
                  WDB[loc0 + THERM_OTF_MAX_TEMP] = RDB[mat + MATERIAL_DOPPLER_TEMP];
                }
                */
                }
            }
          else
            Error(mat, "S(a,b) library %s is not defined",
                  GetText(ptr + THERM_PTR_ALIAS));

          /* Next */

          ptr = NextItem(ptr);
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /* Remove unused data */

  ptr = (long)RDB[DATA_PTR_T0];
  RemoveFlaggedItems(ptr, THERM_OPTIONS, OPT_USED, NO);

  /***************************************************************************/

  /***** Put ZA in ACE data **************************************************/

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_T0];
  while (loc0 > VALID_PTR)
    {
      /* Loop over S(a,b) nuclides */

      ptr = (long)RDB[loc0 + THERM_PTR_SAB];

      while(ptr > VALID_PTR)
        {
          /* Check */

          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Loop over ace data (ei voi käyttää VALID_PTR) */

          ace = (long)RDB[DATA_PTR_ACE0];
          while (ace > 0)
            {
              /* Compare */

              WDB[DATA_DUMMY] = ACE[ace + ACE_PTR_NAME];
              if (!strcmp(GetText(DATA_DUMMY), GetText(ptr + SAB_PTR_NAME)))
                {
                  /* Check type */

                  if ((long)ACE[ace + ACE_TYPE] == NUCLIDE_TYPE_SAB)
                    {
                      /* Check that ZA is not already defined */

                      if ((ACE[ace + ACE_BOUND_ZA] > 0) &&
                          (ACE[ace + ACE_BOUND_ZA] !=
                           RDB[loc0 + THERM_ZA]))
                        Note(loc0,
                              "S(a,b) library %s already assigned to ZA = %ld",
                              GetText(ptr + SAB_PTR_NAME),
                              (long)ACE[ace + ACE_BOUND_ZA]);
                      else
                        ACE[ace + ACE_BOUND_ZA] = RDB[loc0 + THERM_ZA];

                      /* Break loop */

                      break;
                    }
                  else
                    Die(FUNCTION_NAME, "%s is not an S(a,b) libary",
                        GetText(ptr + SAB_PTR_NAME));
                }

              /* next */

              ace = (long)ACE[ace + ACE_PTR_NEXT];
            }

          /* Check if found */

          if (ace < 0)
            Die(FUNCTION_NAME, "Nuclide %s not found in data libraries",
                GetText(ptr + SAB_PTR_NAME));

          ptr = NextItem(ptr);
        }

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Calculate fractions for temperature interpolation (stoch. mixing)  **/

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Loop over S(a,b) data */

      ptr = (long)RDB[mat + MATERIAL_PTR_SAB];
      while (ptr > VALID_PTR)
        {
          /* Get pointer */

          loc0 = (long)RDB[ptr + THERM_PTR_THERM];
          CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

          /* In case of no stochastic interpolation, set frac = 1.0 */
          /* and skip nuclide */

          if ((long)RDB[loc0 + THERM_INTERP_MODE] != THERM_INTERP_STOCHMIX)
            {
              loc1 = (long)RDB[loc0 + THERM_PTR_SAB];
              WDB[loc1 + SAB_FRAC] = 1.0;

              ptr = NextItem(ptr);
              continue;
            }

          /* Get temperatures */

          T = RDB[loc0 + THERM_T];

          loc1 = (long)RDB[loc0 + THERM_PTR_SAB];
          loc2 = NextItem(loc1);

          if ((loc1 < VALID_PTR) || (loc2 < VALID_PTR))
            Die(FUNCTION_NAME, "Sab data for interpolation not found");

          T1 = RDB[loc1 + SAB_T];
          T2 = RDB[loc2 + SAB_T];

          /* Check */

          if (T1 >= T2)
            Error(ptr,
                  "Error in stochastic temperature interpolation parameters");
          else if ((T < T1) || (T > T2))
            Error(ptr,
                  "Stochastic interpolation temperature is not within boundaries");
          else
            {
              /* Interpolate fractions */

              WDB[loc1 + SAB_FRAC] = (T2 - T)/(T2 - T1);
              WDB[loc2 + SAB_FRAC] = 1.0 - RDB[loc1 + SAB_FRAC];
            }

          /* Next */

          ptr = NextItem(ptr);
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/
}

/*****************************************************************************/
