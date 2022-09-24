/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : mgxs.c                                         */
/*                                                                           */
/* Created:       2011/07/20 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Returns coarse multi-group (majorant) cross section          */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MGXS:"

/*****************************************************************************/

double MGXS(long rea, double E, long i)
{
  long erg, ptr;
  double xs;

  /* Check Reaction pointer */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Get pointer to energy grid */

  erg = (long)RDB[DATA_COARSE_MG_PTR_GRID];
  CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

  /* Get grid index if not given */
  
  if (i < 0)
    i = GridSearch(erg, E);
  
  /* Check */

  if (i < 0)
    return 0.0;

  /* Pointer to data */

  ptr = (long)RDB[rea + REACTION_PTR_MGXS];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get value */

  xs = RDB[ptr + i];

  /* Return value */

  return xs;
}

/*****************************************************************************/
