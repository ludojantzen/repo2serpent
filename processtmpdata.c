/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processtmpdata.c                               */
/*                                                                           */
/* Created:       2013/02/07 (JLe)                                           */
/* Last modified: 2015/11/01 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Creates reaction structures, allocates memory, links         */
/*              pointers, etc. for DBRC and TMS.                             */
/*                                                                           */
/* Comments: - JLe 1.11.2015 / 2.1.25: Poistettiin REACTION_PTR_PREV_XS0     */
/*             -varaus (liittyy vaan ures-dataan)                            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessTmpData:"

/*****************************************************************************/

void ProcessTmpData()
{
  long nuc, nuc0, rea, rea0, ptr, ne;
  double maxT;

  /* Check DBRC and TMS modes */

  if (((long)RDB[DATA_USE_DBRC] == NO) && 
      ((long)RDB[DATA_TMS_MODE] == TMS_MODE_NONE))
    return;

  /***************************************************************************/

  /***** DBRC data ***********************************************************/

  /* Reset DBRC flag */

  WDB[DATA_USE_DBRC] = (double)NO;

  /* Loop over nuclides with DBRC data */

  nuc0 = (long)RDB[DATA_PTR_NUC0];
  while (nuc0 > VALID_PTR)
    {
      /* Check DBRC flag */

      if ((long)RDB[nuc0 + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DBRC)
        {
          /* Check cross section temperature */

          if (RDB[nuc0 + NUCLIDE_XS_TEMP] != 0.0)
            Die(FUNCTION_NAME, "DBRC nuclide %s above 0K temperature", 
                GetText(nuc0 + NUCLIDE_PTR_NAME));

          /* Reset maximum temperature */

          maxT = -INFTY;

          /* Find nuclides with same ZAI */

          nuc = (long)RDB[DATA_PTR_NUC0];
          while (nuc > VALID_PTR)
            {
              /* Compare ZAI */

              if ((RDB[nuc0 + NUCLIDE_ZAI] == RDB[nuc + NUCLIDE_ZAI]) &&
                  (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DBRC)))
                {
                  /* Compare maximum temperature */

                  if (RDB[nuc + NUCLIDE_XS_TEMP] > maxT)
                    maxT = RDB[nuc + NUCLIDE_XS_TEMP];
                  if (RDB[nuc + NUCLIDE_TMS_MAX_TEMP] > maxT)
                    maxT = RDB[nuc + NUCLIDE_TMS_MAX_TEMP];

                  /* Pointers to elastic cross section */

                  rea0 = (long)RDB[nuc0 + NUCLIDE_PTR_ELAXS];
                  CheckPointer(FUNCTION_NAME, "(rea0)", DATA_ARRAY, rea0);

                  rea = (long)RDB[nuc + NUCLIDE_PTR_ELAXS];
                  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

                  /* Link reaction channels */
                  
                  WDB[rea + REACTION_PTR_0K_DATA] = (double)rea0;
                }

              /* Next nuclide */

              nuc = NextItem(nuc);
            }
          
          /* Check temperature */

          if (maxT < 0.0)
            {
              /* Reset DBRC flag */

              ResetOption(nuc0 + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_DBRC);

              /* Pointer to next nuclide */

              nuc0 = NextItem(nuc0);

              /* Cycle loop */

              continue;
            }
          
          /* Set DBRC flag */
          
          WDB[DATA_USE_DBRC] = (double)YES;

          /* Set maximum temperature (this will cause the nuclide to be */
          /* included in the next loop) */

          WDB[nuc0 + NUCLIDE_DBRC_MAX_TEMP] = maxT;
        }

      /* Next nuclide */

      nuc0 = NextItem(nuc0);
    }

  /***************************************************************************/
          
  /**** Allocate memory for majorant reactions *******************************/
          
  /* Loop over nuclides */

  nuc0 = (long)RDB[DATA_PTR_NUC0];
  while (nuc0 > VALID_PTR)
    {
      /* Pointer to elastic or total cross section */

      if ((long)RDB[nuc0 + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DBRC)
        rea0 = (long)RDB[nuc0 + NUCLIDE_PTR_ELAXS];
      else if (((long)RDB[nuc0 + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS) &&
               ((long)RDB[nuc0 + NUCLIDE_TYPE] != NUCLIDE_TYPE_DECAY))
        rea0 = (long)RDB[nuc0 + NUCLIDE_PTR_TOTXS];
      else
        {
          /* Pointer to next */

          nuc0 = NextItem(nuc0);

          /* Cycle loop */

          continue;
        }

      /* Check pointer */

      CheckPointer(FUNCTION_NAME, "(rea0)", DATA_ARRAY, rea0);
          
      /* Duplicate reaction channel */
          
      rea = DuplicateItem(rea0);
      
      /* Allocate memory for previous values */
      
      WDB[rea + REACTION_PTR_PREV_XS] = NULLPTR;
      AllocValuePair(rea + REACTION_PTR_PREV_XS);
      
      /* Disable cache optimization mode */
      
      WDB[rea + REACTION_CACHE_OPTI_IDX] = -1.0;
      
      /* Get number of energy points */
      
      ne = (long)RDB[rea0 + REACTION_XS_NE];
      CheckValue(FUNCTION_NAME, "ne", "", ne, 10, 50000000);
      
      /* Allocate memory for data */
      
      ptr = ReallocMem(DATA_ARRAY, ne);
      WDB[rea + REACTION_PTR_XS] = (double)ptr;
      
      /* Remove reaction from list */
      
      RemoveItem(rea);
      
      /* Put pointer */
      
      WDB[rea0 + REACTION_PTR_TMP_MAJORANT] = (double)rea;
      
      /* Next nuclide */
      
      nuc0 = NextItem(nuc0);
    }
  
  /* Calculate majorant cross sections */

  TmpMajorants();

}

/*****************************************************************************/
