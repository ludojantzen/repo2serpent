/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : setfisse.c                                     */
/*                                                                           */
/* Created:       2016/03/09 (JLe)                                           */
/* Last modified: 2019/22/03 (RTu)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Sets deposited fission energies for nuclide for              */
/*              normalization.                                               */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SetFissE:"

/*****************************************************************************/

void SetFissE()
{
  long nuc, rea, mt, ptr, ZAI, mode, loc0, delayed;
  double Q235, H235, Q, H;

  /* Set default U-235 Q-value */

  Q235 = U235_FISSQ;

  /* Read value from data */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR)
    {
      /* Check ZAI */

      if ((long)RDB[nuc + NUCLIDE_ZAI] != 922350)
        {
          /* Next */

          nuc = NextItem(nuc);
          
          /* Cycle loop */

          continue;
        }

      /* Loop over reactions */
      
      rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
      while (rea > VALID_PTR)
        {
          /* Get mt */
          
          mt = (long)RDB[rea + REACTION_MT];
          
          /* Check */

          if ((mt == 18) || (mt == 19))
            if ((Q = RDB[rea + REACTION_Q]) > 0.0)
              {
                /* Set value */

                Q235 = Q;
                
                /* Break loop */

                break;
              }

          /* Next */

          rea = NextItem(rea);
        }

      /* Next */
      
      nuc = NextItem(nuc);
    }

  /* Get U-235 fission energy deposition */

  H235 = RDB[DATA_NORM_U235_FISSE];
  CheckValue(FUNCTION_NAME, "H235", "", H235, ZERO, INFTY);

  /* Set default fission energies */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR)
    {
      /* Loop over reactions */
      
      rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
      while (rea > VALID_PTR)
        {
          /* Get mt */
          
          mt = (long)RDB[rea + REACTION_MT];

          /* Check */

          if ((mt == 18) || (mt == 19))
            if ((Q = RDB[rea + REACTION_Q]) > 0.0)
              {
                /* Set value */
            
                WDB[nuc + NUCLIDE_FISSE] = H235*Q/Q235;

                /* Avoid compiler warning */

                ZAI = -1;

                /* Loop over list and override */

                if ((ptr = (long)RDB[DATA_NORM_PTR_FISSH]) > VALID_PTR)
                  while ((ZAI = (long)RDB[ptr++]) > 0)
                    {
                      /* Convert to ZAI */

                      if (ZAI < 900000)
                        ZAI = ZAI*10;

                      /* Get H */

                      H = RDB[ptr++];
                      CheckValue(FUNCTION_NAME, "H", "", H, ZERO, INFTY);
                    
                      /* Compare and replace */
                      
                      if ((long)RDB[nuc + NUCLIDE_ZAI] == ZAI)
                        {
                          /* Set value */
                          
                          WDB[nuc + NUCLIDE_FISSE] = H*MEV;
                          
                          /* Break loop */

                          break;
                        }
                    }
                
                /* Check ZAI and break reaction loop */

                if (ZAI > 0)
                  break;
              }

          /* Next */

          rea = NextItem(rea);
        }

      /* Next */ 

      nuc = NextItem(nuc);
    }

  /* Get mode */

  mode = (long)RDB[DATA_EDEP_MODE];

  /* Process depending on energy deposition mode */

  if (mode > EDEP_MODE_CONSTANT)
    {
      /* Allocate memory */

      loc0 = ReallocMem(DATA_ARRAY, FISSE_COMP_BLOCK_SIZE);
      WDB[DATA_EDEP_COMP] = (double)loc0;

      /* Get include delayed flag */

      delayed = (long)RDB[DATA_EDEP_DELAYED];

      /* Set component flags */
      /* EFR is always calculated */

      WDB[loc0 + FISSE_COMP_EFR] = (double)YES;

      if (mode == EDEP_MODE_MT458)
        {
          WDB[loc0 + FISSE_COMP_ENP] = (double)YES;
          WDB[loc0 + FISSE_COMP_EGP] = (double)YES;

          if (delayed == YES)
            {
              WDB[loc0 + FISSE_COMP_END] = (double)YES;
              WDB[loc0 + FISSE_COMP_EGD] = (double)YES;
              WDB[loc0 + FISSE_COMP_EB] = (double)YES;
            }
        }
      else if (mode == EDEP_MODE_LOCAL_PHOTON)
        {
          WDB[loc0 + FISSE_COMP_EGP] = (double)YES;

          if (delayed == YES)
            {
              WDB[loc0 + FISSE_COMP_EGD] = (double)YES;
              WDB[loc0 + FISSE_COMP_EB] = (double)YES;
            }
        }
      else if ((mode == EDEP_MODE_NEUTRON_PHOTON) && (delayed == YES))
        {
          WDB[loc0 + FISSE_COMP_EB] = (double)YES;
          if ((long)RDB[DATA_EDEP_LOCAL_EGD] == YES)
            WDB[loc0 + FISSE_COMP_EGD] = (double)YES;
        }
    }
}

/*****************************************************************************/
