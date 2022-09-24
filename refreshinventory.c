/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : refreshinventory.c                             */
/*                                                                           */
/* Created:       2012/06/20 (TVi)                                           */
/* Last modified: 2017/07/30 (JLe)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: Updates the inventory in case intelligent inventory is used. */
/*                                                                           */
/* Comments: - Lisää (ei ikinä poista) nuklideita inventoryyn niin että      */
/*             intelligent inventory -määrittelyssä annettu määrä            */
/*             tärkeimpiä samaisessa määrittelyssä annetun suureen           */
/*               kontribuuttoreita on aina mukana inventaarissa.               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "RefreshInventory:"

/*****************************************************************************/

void RefreshInventory()
{
  long i, nuc, loc0, k, qptr, n, ntop;
  char name[MAX_STR];

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Allow memory allocation */

  Mem(MEM_ALLOW);
  
  /* Loop over parameters */
  
  for (n = 0; n < 7; n++)
    {
      /* Avoid compiler warning */
      
      ntop = -1;
      qptr = -1;

      /* Get values */

      if (n == 0)
        {
          ntop = (long)RDB[DATA_BURN_INV_TOP_MASS];
          qptr = NUCLIDE_TOP_MASS;
        }
      else if (n == 1)
        {
          ntop = (long)RDB[DATA_BURN_INV_TOP_ACTIVITY];
          qptr = NUCLIDE_TOP_ACTIVITY;
        }
      else if (n == 2)
        {
          ntop = (long)RDB[DATA_BURN_INV_TOP_SF];
          qptr = NUCLIDE_TOP_SF;
        }
      else if (n == 3)
        {
          ntop = (long)RDB[DATA_BURN_INV_TOP_GSRC];
          qptr = NUCLIDE_TOP_GSRC;
        }
      else if (n == 4)
        {
          ntop = (long)RDB[DATA_BURN_INV_TOP_DECAY_HEAT];
          qptr = NUCLIDE_TOP_DECAY_HEAT;
        }
      else if (n == 5)
        {
          ntop = (long)RDB[DATA_BURN_INV_TOP_ING_TOX];
          qptr = NUCLIDE_TOP_ING_TOX;
        }
      else if (n == 6)
        {
          ntop = (long)RDB[DATA_BURN_INV_TOP_INH_TOX];
          qptr = NUCLIDE_TOP_INH_TOX;
        }
      else
        Die(FUNCTION_NAME, "Overflow");

      /* Cycle loop if not requested */

      if (ntop < 1)
        continue;

      /* Sort list */

      nuc = (long)RDB[DATA_PTR_NUC0];
      SortList(nuc, qptr, SORT_MODE_DESCEND);

      /* Get last index in inventory */
      
      if ((loc0 = (long)RDB[DATA_BURN_PTR_INVENTORY]) > VALID_PTR)
        k = ListSize(loc0);
      else
        k = 0;

      /* Loop over top nuclides */

      nuc = (long)RDB[DATA_PTR_NUC0];
      for (i = 0; i < ntop; i++)
        {
          /* Break at zero (JLe 30.7.2017 / 2.1.30) */

          if (RDB[nuc + qptr] == 0.0)
            break;

          /* Check that nuclide mass is above zero (this is needed) */
          /* for the first step) */

          if (RDB[nuc + NUCLIDE_TOP_MASS] > 0.0)
            {
              /* Find match in inventory list */
              
              loc0 = (long)RDB[DATA_BURN_PTR_INVENTORY];
              while (loc0 > VALID_PTR)
                {
                  /* Compare ZAI */
                  
                  if (RDB[nuc + NUCLIDE_ZAI] == RDB[loc0 + INVENTORY_ZAI])
                    break;
                  
                  /* Next nuclide in list */
                  
                  loc0 = NextItem(loc0);
                }
              
              /* Check if found */
              
              if (loc0 < VALID_PTR)
                {
                  /* Add new nuclide to inventory */
                  
                  loc0 = NewItem(DATA_BURN_PTR_INVENTORY, 
                                 INVENTORY_BLOCK_SIZE);
                  
                  /* Put name */
                  
                  sprintf(name, "%s", 
                          ZAItoIso((long)RDB[nuc + NUCLIDE_ZAI], 3));
                  WDB[loc0 + INVENTORY_PTR_NAME] = (double)PutText(name);
                  
                  /* Put ZAI */
                  
                  WDB[loc0 + INVENTORY_ZAI] = RDB[nuc + NUCLIDE_ZAI];
                  
                  /* Set index */
                  
                  WDB[nuc + NUCLIDE_INVENTORY_IDX] = (double)k++;           
                }
            }

          /* Next nuclide */
          
          nuc = NextItem(nuc);
        }
    }
 
  /* Resort nuclide list according to ZAI */

  nuc = (long)RDB[DATA_PTR_NUC0];
  SortList(nuc, NUCLIDE_ZAI, SORT_MODE_ASCEND);

  /* Disallow memory allocation */

  Mem(MEM_DENY);
}

/*****************************************************************************/
