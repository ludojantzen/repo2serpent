/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : cachexs.c                                      */
/*                                                                           */
/* Created:       2012/11/26 (JLe)                                           */
/* Last modified: 2018/11/02 (RTu)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: Converts cross sections into a cache-friendly data block     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CacheXS:"

/*****************************************************************************/

void CacheXS()
{
  long nuc, rea, nrea, n, ptr, ne, loc0, loc1, loc2, i0, ne0, dum1, dum2, i;
  double Emax, *xs;
  const double *E;

  /* Check that neutron transport mode is run */

  if ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == NO)
    return;

  /* Get energy boundary */

  if ((Emax = RDB[DATA_CACHE_OPTI_EMAX]) < 1E-11)
    return;

  /* Check decay only mode */

  if ((long)RDB[DATA_BURN_DECAY_CALC] == YES)
    return;

  /***************************************************************************/

  /**** Calculate number of reactions ****************************************/

  /* Reset count */

  nrea = 0;

  /* Loop over nuclides */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR)
    {
      /* Loop over reactions */

      rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
      while (rea > VALID_PTR)
        {
          /* Reset index */

          WDB[rea + REACTION_CACHE_OPTI_IDX] = -1.0;

          /* Check pointer to cross section data and add to counter */

          if (((long)RDB[rea + REACTION_PTR_XS] > VALID_PTR) &&
              (RDB[rea + REACTION_EMIN] < Emax))
            nrea++;
          
          /* Next reaction */
          
          rea = NextItem(rea);
        }

      /* Add total and total absorption*/

      nrea = nrea + 2;

      /* Next nuclide */

      nuc = NextItem(nuc);
    }

  /***************************************************************************/

  /***** Energy grid *********************************************************/

  /* Pointer to unionized grid */

  ptr = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Find limiting point in grid */
  
  if ((ne = GridSearch(ptr, Emax)) < 0)
    Die(FUNCTION_NAME, "Error in grid search");
  else
    ne = ne + 2;

  /* Store size */

  WDB[DATA_CACHE_OPTI_NE] = (double)ne;
  WDB[DATA_CACHE_OPTI_NREA] = (double)nrea;

  /* Allocate memory for data */
  
  loc0 = ReallocMem(DATA_ARRAY, nrea*ne);
  
  /* Put pointer */

  WDB[DATA_PTR_CACHE_OPTI_XS] = (double)loc0;

  /* Pointer to energy grid data */

  ptr = (long)RDB[ptr + ENERGY_GRID_PTR_DATA];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Vector (NOTE: ei reallocmemiä tän jälkeen) */

  E = &RDB[ptr];

  /* Check */

  if (E[ne] < Emax)
    Die(FUNCTION_NAME, "Error in energy");

  /***************************************************************************/

  /***** Read data ***********************************************************/

  /* Allocate memory for temporary xs */

  xs = (double *)Mem(MEM_ALLOC, ne, sizeof(double));

  /* Reset count */

  n = 0;

  /* Loop over nuclides */

  nuc = (long)RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR)
    {
      /* Loop over reactions */

      rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
      while (rea > VALID_PTR)
        {
          /* Check reaction type and pointer to cross section data */

          if (((loc1 = (long)RDB[rea + REACTION_PTR_XS]) > VALID_PTR) &&
              (((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_SUM) ||
               ((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_PARTIAL)))
            {
              /* Pointer to energy grid */

              loc2 = (long)RDB[rea + REACTION_PTR_EGRID];
              CheckPointer(FUNCTION_NAME, "(loc2a)", DATA_ARRAY, loc2);

              loc2 = (long)RDB[loc2 + ENERGY_GRID_PTR_DATA];
              CheckPointer(FUNCTION_NAME, "(loc2b)", DATA_ARRAY, loc2);

              /* Get index and size */

              i0 = (long)RDB[rea + REACTION_XS_I0];
              ne0 = (long)RDB[rea + REACTION_XS_NE];
              
              /* Check energy */

              if ((RDB[rea + REACTION_EMIN] < Emax))
                {
                  /* Put index */

                  WDB[rea + REACTION_CACHE_OPTI_IDX] = (double)n;

                  /* Interpolate data */
                  
                  InterpolateData(E, xs, ne, &RDB[loc2 + i0], &RDB[loc1], ne0, 
                                  0, &dum1, &dum2, NO);

                  /* Write data */

                  for (i = 0; i < ne; i++)
                    WDB[loc0 + i*nrea + n] = xs[i];
                  
                  /* Update counter */
                  
                  n++;
                }
            }
          
          /* Next reaction */
          
          if (rea == (long)RDB[nuc + NUCLIDE_PTR_SUM_ABSXS])
            break;
          else if (rea != (long)RDB[nuc + NUCLIDE_PTR_TOTXS])
            rea = NextItem(rea);

          /* Last two are total and sum of absorption */

          if (rea < VALID_PTR)
            rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
          else if (rea == (long)RDB[nuc + NUCLIDE_PTR_TOTXS])
            rea = (long)RDB[nuc + NUCLIDE_PTR_SUM_ABSXS];
        }

      /* Next nuclide */

      nuc = NextItem(nuc);
    }

  /***************************************************************************/

  /* Free allocated memory */
  
  Mem(MEM_FREE, xs);
}

/*****************************************************************************/
