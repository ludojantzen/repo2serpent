/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : overrideids.c                                  */
/*                                                                           */
/* Created:       2014/01/23 (JLe)                                           */
/* Last modified: 2017/10/26 (JLe)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Preprocessor routine that overrides nuclide ID's and thermal */
/*              scattering data.                                             */
/*                                                                           */
/* Comments: - Used for initializing coefficient calculations for burnup     */
/*             branches.                                                     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "OverrideIDs:"

/*****************************************************************************/

void OverrideIDs()
{
  long mat, mat0, iso, ace, ace0, ZAI, dop, ptr, loc;
  double T, T0;
  char zaid[MAX_STR], sab[MAX_STR];

  /***************************************************************************/

  /***** Nuclide ID's in material compositions *******************************/

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Get temperature and pointer to parent (parameters copied later) */

      if (((T0 = RDB[mat + MATERIAL_COEF_TEMP]) == 0.0) ||
          ((long)RDB[mat + MATERIAL_DIV_PTR_PARENT] > VALID_PTR))
        {
          /* Not adjusted */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Reset Doppler flag */

      dop = NO;

      /* Loop over composition */

      iso = (long)RDB[mat + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
        {
          /* Get name */

          sprintf(zaid, "%s", GetText(iso + COMPOSITION_PTR_NUCLIDE));

          /* Find match in ace directory file */
          
          ace = (long)RDB[DATA_PTR_ACE0];      
          while (ace > VALID_PTR)
            {
              /* Compare */

              WDB[DATA_DUMMY] = ACE[ace + ACE_PTR_ALIAS];
              if (!strcmp(GetText(DATA_DUMMY), zaid))
                break;
          
              /* next */

              ace = (long)ACE[ace + ACE_PTR_NEXT];
            }

          /* Check if found (error message should be the same as that */
          /* printed by processnuclides.c */

          if (ace < VALID_PTR)
            Error(mat, "Nuclide %s not found in data libraries", zaid);

          /* Get ZAI */

          ZAI = (long)ACE[ace + ACE_ZAI];

          /* Find closest match for temperature */

          T = -1.0;
          ace0 = -1;

          ace = (long)RDB[DATA_PTR_ACE0];      
          while (ace > VALID_PTR)
            {
              /* Compare ZAI */
              
              if (ZAI == (long)ACE[ace + ACE_ZAI])
                {
                  /* Check for matching temperature (within 1 degrees) */

                  /* NOTE: Tätä muutettiin 3.8.2017 (jos Doppler-lämpötila */
                  /* jää ton pyöristetyn alapuolelle niin antaa errorin */

                  if ((T0 >= ACE[ace + ACE_TEMP]) && 
                      (T0 - ACE[ace + ACE_TEMP] < 1.0))
                    /*
                  if (fabs(ACE[ace + ACE_TEMP] - T0) < 1.0)
                    */
                    {
                      /* Set pointer */
                      
                      ace0 = ace;

                      /*  break loop */
                    }
                  else if (ACE[ace + ACE_TEMP] < T0)
                    {
                      /* Find highest below given */

                      if (ACE[ace + ACE_TEMP] > T)
                        {
                          /* Update temperature and pointer */
                          
                          T = ACE[ace + ACE_TEMP];
                          ace0 = ace;

                          /* Set Doppler flag */

                          dop = YES;
                        }
                    }
                }
            
              /* next */

              ace = (long)ACE[ace + ACE_PTR_NEXT];
            }

          /* Check */

          if (ace0 < VALID_PTR)
            Error(mat, "Unable to override ID's to temperature %1.1f", T0);

          /* Replace zaid */

          WDB[DATA_DUMMY] = ACE[ace0 + ACE_PTR_ALIAS];          
          sprintf(zaid, "%s", GetText(DATA_DUMMY));

          /* Put new value */

          WDB[iso + COMPOSITION_PTR_NUCLIDE] = (double)PutText(zaid);          
          
          /* Next */

          iso = NextItem(iso);
        }

      /* Check Doppler flag */

      if (dop == YES)
        {
          /* Check if TMS is already in use */

          if ((long)RDB[mat + MATERIAL_TMS_MODE] == YES)
            {
              /* Put temperature */

              WDB[mat + MATERIAL_TMS_TMIN] = T0;
              WDB[mat + MATERIAL_TMS_TMAX] = T0;
            }
          else
            {
              /* Put preprocessor temperature */
          
              WDB[mat + MATERIAL_DOPPLER_TEMP] = T0;
              
              /* Switch preprocessor calculation on */
          
              WDB[DATA_USE_DOPPLER_PREPROCESSOR] = (double)YES;
            }
        }

      /* Next material */

      mat = NextItem(mat);
    }

  /* Copy Doppler-temperatures to divided */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Pointer to parent */
      
      if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
        {
          /* Copy temperature */
          
          WDB[mat + MATERIAL_DOPPLER_TEMP] = RDB[mat0 + MATERIAL_DOPPLER_TEMP];
        } 
        
        /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/

  /***** S(a,b) libraries ****************************************************/

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check replaced S(a,b) library */

      if ((long)RDB[mat + MATERIAL_COEF_SAB] < VALID_PTR)
        {
          /* Not adjusted */

          mat = NextItem(mat);

          /* Cycle loop */

          continue;
        }

      /* Get replaced S(a,b) library */

      sprintf(sab, "%s", GetText(mat + MATERIAL_COEF_SAB));

      /* Get pointer to original S(a,b) library */

      if ((ptr = (long)RDB[mat + MATERIAL_PTR_SAB]) < VALID_PTR)
        Error(mat, "Material %s has no S(a,b) data for branch calculation",
              GetText(mat + MATERIAL_PTR_NAME));

      /* Replace name */

      WDB[ptr + THERM_PTR_ALIAS] = RDB[mat + MATERIAL_COEF_SAB];
      
      /* Create new item */

      ptr = NewItem(DATA_PTR_T0, THERM_BLOCK_SIZE);

      /* Put name, file name and line number */

      WDB[ptr + PARAM_PTR_NAME] = RDB[mat + PARAM_PTR_NAME];
      WDB[ptr + PARAM_PTR_FNAME] = RDB[mat + PARAM_PTR_FNAME];
      WDB[ptr + PARAM_LINE] = RDB[mat + PARAM_LINE];
      
      /* Put name and isotope name */
      
      WDB[ptr + THERM_PTR_ALIAS] = RDB[mat + MATERIAL_COEF_SAB];

      if ((loc = (long)RDB[ptr + THERM_PTR_SAB]) > VALID_PTR)
        WDB[loc + SAB_PTR_NAME] = RDB[mat + MATERIAL_COEF_SAB];
      
      /* Reset temperature */
      
      WDB[ptr + THERM_T] = -1.0;

      /* Next material */

      mat = NextItem(mat);
    }

  /***************************************************************************/
}

/*****************************************************************************/
