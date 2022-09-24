/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : setdirectpointers.c                            */
/*                                                                           */
/* Created:       2010/11/19 (JLe)                                           */
/* Last modified: 2019/10/24 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: - Orders item pointers in array                              */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SetDirectPointers:"

/*****************************************************************************/

void SetDirectPointers(long lst)
{
  long ptr, loc0, ni, n;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);

  /* Get pointer to first item */

  lst = FirstItem(lst);

  /* Get pointer to common data */

  loc0 = (long)RDB[lst + LIST_PTR_COMMON];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Get number of items */

  ni = (long)RDB[loc0 + LIST_COMMON_N_ITEMS];

  /* Get pointer to array */

  if ((loc0 = (long)RDB[lst + LIST_PTR_DIRECT]) < VALID_PTR)
    Die(FUNCTION_NAME, "List is not closed");

  /* Loop over items */

  n = 1;

  ptr = lst;

  while (ptr > VALID_PTR)
    {
      /* Check index */

      if (n > ni)
        Die(FUNCTION_NAME, "Indexing error (n > ni)");

      /* Put item pointer in list */

      WDB[loc0 + n++] = (double)ptr;

      /* Check if last (first pointer in list) */

      if (NextItem(ptr) < VALID_PTR)
        WDB[loc0] = (double)ptr;

      /* Next */

      ptr = NextItem(ptr);
    }

  /* Put null terminator */

  WDB[loc0 + n++] = NULLPTR;

  /* Check index */

  if (n != ni + 2)
    Die(FUNCTION_NAME, "Indexing error (n != ni + 2)");
}

/*****************************************************************************/
