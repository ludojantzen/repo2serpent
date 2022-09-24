/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : burnmatrixsize.c                               */
/*                                                                           */
/* Created:       2011/05/22 (JLe)                                           */
/* Last modified: 2012/02/01 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Calculates number of non-zero elements in burnup matrix      */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "BurnMatrixSize:"

/*****************************************************************************/

long BurnMatrixSize(long mat)
{
  long sz, iso, nuc, rea, yld, ptr, n;

  /* Reset size */

  sz = 0;

  /* Loop over composition */

  iso = (long)RDB[mat + MATERIAL_PTR_COMP];
  while (iso > VALID_PTR)
    {
      /* Diagonal element */

      sz++;

      /* Pointer to nuclide */

      nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];

      /* Loop over reactions */

      rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
      while (rea > VALID_PTR)
	{
	  /* Check pointer to fission yield */

	  if ((yld = (long)RDB[rea + REACTION_PTR_FISSY]) > VALID_PTR)
	    {
	      /* Get pointer to distribution */
	      
	      yld = (long)RDB[yld + FISSION_YIELD_PTR_DISTR];

	      /* Check pointer */

	      CheckPointer(FUNCTION_NAME, "yld", DATA_ARRAY, yld);
	  
	      /* Loop over distribution */ 
	      
	      n = 0;
	      while ((ptr = ListPtr(yld, n++)) > 0)
		sz++;
	    }	  
	  else if ((long)RDB[rea + REACTION_PTR_TGT] > VALID_PTR)
	    {
	      /* Decay or transmutation */

	      sz++;
	    }

	  /* Next reaction */
	  
	  rea = NextItem(rea);
	}

      /* Check pointer to total fission */

      if ((long)RDB[nuc + NUCLIDE_PTR_TOTFISS_REA] > VALID_PTR)
	sz++;

      /* Next nuclide */
      
      iso = NextItem(iso);
    }

  /* Return size */

  return sz;
}

/*****************************************************************************/
