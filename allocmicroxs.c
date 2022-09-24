/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : allocmicroxs.c                                 */
/*                                                                           */
/* Created:       2011/11/30 (JLe)                                           */
/* Last modified: 2019/03/21 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Estimates total memory required by cross section data and    */
/*              pre-allocates it                                             */
/*                                                                           */
/* Comments: - This is to speed up ProcessXSData()                           */
/*           - Does not account for gamma transport cross sections           */
/*           - Name changed to allocmicroxs.c from preallocxsmem.c           */
/*             29.5.2012 (2.1.6)                                             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AllocMicroXS:"

/*****************************************************************************/

void AllocMicroXS()
{
  long sz, ptr, ne, nuc, rea, type, NES, tot;
  double f;

  /* Check decay only mode */

  if ((long)RDB[DATA_BURN_DECAY_CALC] == YES)
    return;

  /***************************************************************************/

  /***** Memory for microscopic cross sections *******************************/

  /* Reset memory and neutron energy grid size */

  sz = 0;
  ne = 0;

  /* Get unionized grid size */

  if ((long)RDB[DATA_OPTI_RECONSTRUCT_MICROXS] == YES)
    {
      /* Pointer to grid (not defined in gamma mode) */

      ptr = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID];
      
      /* Get number of energy points */

      if (ptr > VALID_PTR)
        ne = (long)RDB[ptr + ENERGY_GRID_NE];
    }

  /* Loop over nuclides */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR)
    {
      /* Check pointer to ACE data */

      if ((long)RDB[nuc + NUCLIDE_PTR_ACE] > VALID_PTR)
        {
          /* Pointer to reaction data */

          rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
          CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

          /* Get total number of points (elastic XS should be first) */

          tot = (long)RDB[rea + REACTION_XS_NE];

          /* Add total xs to estimate */

          if (ne == 0)
            sz = sz + 2*tot;
          else
            sz = sz + 2*ne;          

          /* Skip first (elastic scattering) */

          rea = NextItem(rea);
          
          /* Loop over remaining reactions */

          while (rea > VALID_PTR)
            {
              /* Get type */
                
              type = (long)RDB[rea + REACTION_TYPE];
              
              /* Skip everything but partial and special reactions */
                
              if ((type == REACTION_TYPE_PARTIAL) ||
                  (type == REACTION_TYPE_SPECIAL))
                {
                  /* Number of energy points */
                
                  NES = (long)RDB[rea + REACTION_XS_NE];

                  /* Check reconstruction mode */

                  if ((long)RDB[DATA_OPTI_RECONSTRUCT_MICROXS] == NO)
                    {
                      /* No reconstruction, add points to total */

                      sz = sz + NES;
                    }
                  else if (tot > 0)
                    {
                      /* Reconstructed cross section, calculate fraction */

                      f = ((double)NES)/((double)tot);

                      /* Add to total (additional factor to account for */
                      /* other data). */
                      
                      sz = sz + (long)(f*((double)ne));
                    }
                }

              /* Next reaction */

              rea = NextItem(rea);
            }
        }

      /* Next nuclide */

      nuc = NextItem(nuc);
    }

  /* Allocate memory */

  PreallocMem(sz, DATA_ARRAY);

  /***************************************************************************/
  
  /***** Memory for work arrays **********************************************/

  /* Reset size */

  ne = -1;

  /* Check gamma data */

  if ((ptr = (long)RDB[DATA_ERG_PTR_UNIONIZED_PGRID]) > VALID_PTR)
    ne = (long)RDB[ptr + ENERGY_GRID_NE];

  /* Check neutron data */

  if ((ptr = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID]) > VALID_PTR)
    if ((long)RDB[ptr + ENERGY_GRID_NE] > ne)
      ne = (long)RDB[ptr + ENERGY_GRID_NE];

  /* Check */

  if (ne < 1)
    Die(FUNCTION_NAME, "What now?!?!?");

  /* Allocate memory from main array */

  WorkArray(DATA_PTR_WORK_GRID1, DATA_ARRAY, ne, 0);
  WorkArray(DATA_PTR_WORK_GRID2, DATA_ARRAY, ne, 0);
  WorkArray(DATA_PTR_WORK_GRID3, DATA_ARRAY, ne, 0);

  /* Allocate memory from private array (for OpenMP parallel processing */
  /* routines) */

  WorkArray(DATA_PTR_WORK_PRIVA_GRID1, PRIVA_ARRAY, ne, 0);
  WorkArray(DATA_PTR_WORK_PRIVA_GRID2, PRIVA_ARRAY, ne, 0);
  WorkArray(DATA_PTR_WORK_PRIVA_GRID3, PRIVA_ARRAY, ne, 0);
}

/*****************************************************************************/
