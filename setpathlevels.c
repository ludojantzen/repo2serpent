/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : setpathlevels.c                                */
/*                                                                           */
/* Created:       2010/09/26 (JLe)                                           */
/* Last modified: 2012/02/01 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Loops over transmutation paths and sets minimum levels.      */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SetPathLevels:"

/*****************************************************************************/

void SetPathLevels(long nuc, long level)
{
  long rea, new, loc0, yld, n;
  
  /* Put level */

  if ((double)level < RDB[nuc + NUCLIDE_PATH_LEVEL])
    WDB[nuc + NUCLIDE_PATH_LEVEL] = (double)level;
  else
    return;

  /* Loop over reactions */

  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];

  while (rea > 0)
    {
      /* Get pointer to target */

      if (((new = (long)RDB[rea + REACTION_PTR_TGT]) > 0) &&
	  (new != (long)RDB[DATA_PTR_NUCLIDE_LOST]))
	{
	  /* Call recursively */
	  
	  SetPathLevels(new, level + 1);
	}
      else if ((yld = (long)RDB[rea + REACTION_PTR_FISSY]) > 0)
	{
	  /* Get pointer to distribution */
	  
	  if ((yld = (long)RDB[yld + FISSION_YIELD_PTR_DISTR]) < 0)
	    Die(FUNCTION_NAME, "Pointer error");
	  
	  /* Loop over distribution */ 

	  n = 0;
	  while ((loc0 = ListPtr(yld, n++)) > 0)
	    {
	      /* Get pointer to target */
	      
	      if (((new = (long)RDB[loc0 + FY_PTR_TGT]) > 0) &&
		  (new != (long)RDB[DATA_PTR_NUCLIDE_LOST]))
		{
		  /* Call recursively */
		  
		  SetPathLevels(new, level + 1);
		}
	    }
	}
          
      /* Next */
      
      rea = NextItem(rea);
    }
}

/*****************************************************************************/
