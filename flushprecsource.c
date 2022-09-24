/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : flushprecsource.c                              */
/*                                                                           */
/* Created:       2016/02/01 (VVa)                                           */
/* Last modified: 2016/02/01 (VVa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Removes all particles from precursor source                  */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FlushPrecSource:"

/*****************************************************************************/

void FlushPrecSource()
{
  long ptr, next, id;

  /* Get pointer to precursor source */

  ptr = (long)RDB[DATA_PART_PTR_PSOURCE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  
  /* Get pointer to first after dummy */

  next = NextItem(ptr);

  id = 0;

  /* Loop over PSOURCE and put particles back to stacks */

  while (next > VALID_PTR)
    {

      ptr = next;
      next = NextItem(next);

      RemoveItem(ptr);

      ToStack(ptr, id++);

      /* Check id */

      if (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1)
        id = 0;
    }

  /* Get pointer to precursor source */

  ptr = (long)RDB[DATA_PART_PTR_PSOURCE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Check that source is empty */

  if (ListSize(ptr) != 1)
    Die(FUNCTION_NAME, "Precursor source is not empty");

}

/*****************************************************************************/
