/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : cleartransmuxs.c                               */
/*                                                                           */
/* Created:       2011/05/23 (JLe)                                           */
/* Last modified: 2015/09/26 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Clears all reaction rates etc. used for burnup calculation   */
/*                                                                           */
/* Comments: - This is done in a separate subroutine to simplify things in   */
/*             MPI mode (täh?)                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ClearTransmuXS:"

/*****************************************************************************/

void ClearTransmuXS()
{
  long nuc, rea, id;

  /* Keep results if activation interval */

  if (((long)RDB[DATA_BURN_STEP_TYPE] == DEP_STEP_ACT_STEP) ||
      ((long)RDB[DATA_BURN_STEP_TYPE] == DEP_STEP_ACT_TOT))
    return;

  /* Loop over nuclides */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR)
    {
      /* Loop over reactions */

      rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
      while (rea > VALID_PTR)
        {
          /* Check pointer to transmutation xs and reset data */
          
          if ((long)RDB[rea + REACTION_PTR_TRANSMUXS] > 0)
            for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
              StoreValuePair(rea + REACTION_PTR_TRANSMUXS, 0.0, 0.0, id);

          /* Next reaction */

          rea = NextItem(rea);
        }

      /* Check pointer to total fission */
      
      if ((rea = (long)RDB[nuc + NUCLIDE_PTR_TOTFISS_REA]) > VALID_PTR)
        for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
          StoreValuePair(rea + REACTION_PTR_TRANSMUXS, 0.0, 0.0, id);

      /* Next nuclide */

      nuc = NextItem(nuc);
    }
}

/*****************************************************************************/
