/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : newitem.c                                      */
/*                                                                           */
/* Created:       2010/09/15 (JLe)                                           */
/* Last modified: 2017/01/11 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Adds a new item in list                                      */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "NewItem:"

/*****************************************************************************/

long NewItem(long root, long sz)
{
  long ptr, loc0, new, last;

  /* Check root pointer */

  if (root < 1)
    Die(FUNCTION_NAME, "Pointer error");

  /* Allocate memory for new item */

  new = ReallocMem(DATA_ARRAY, sz);

  /* Check if list exists */
  
  if ((ptr = (long)RDB[root]) < VALID_PTR)
    {
      /* First definition, set pointer */
      
      WDB[root] = (double)new;

      /* Reset pointers to next and previous */

      WDB[new + LIST_PTR_NEXT] = NULLPTR;
      WDB[new + LIST_PTR_PREV] = NULLPTR;

      /* Allocate memory for common data */

      loc0 = ReallocMem(DATA_ARRAY, LIST_COMMON_DATA_SIZE);

      /* Put item size */

      WDB[loc0 + LIST_COMMON_ITEM_SIZE] = (double)sz;

      /* Put root pointer */

      WDB[loc0 + LIST_COMMON_PTR_ROOT] = (double)root;

      /* Put pointer to first item */

      WDB[loc0 + LIST_COMMON_PTR_FIRST] = (double)new;
    }
  else
    {
      /* Data exists, check pointer */
      
      CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

      /* Get pointer to common data */

      loc0 = (long)RDB[ptr + LIST_PTR_COMMON];
      CheckPointer(FUNCTION_NAME, "loc0", DATA_ARRAY, loc0);

      /* Compare root pointer */
      
      if (root != (long)RDB[loc0 + LIST_COMMON_PTR_ROOT])
        Die(FUNCTION_NAME, "Mismatch in root pointer (sz = %ld)", sz);

      /* Compare size */

      if (sz != (long)RDB[loc0 + LIST_COMMON_ITEM_SIZE])
        Die(FUNCTION_NAME, "Mismatch in item size (sz = %ld)", sz);

      /* Check if list is closed */

      if ((long)RDB[ptr + LIST_PTR_DIRECT] > VALID_PTR)
        Die(FUNCTION_NAME, "Trying to add an item to a closed list (sz = %ld)",
            sz);

      /* Get pointer to last item */

      last = LastItem(ptr);
      
      /* Check pointer */

      CheckPointer(FUNCTION_NAME, "last", DATA_ARRAY, last);

      /* Set list pointers */

      WDB[last + LIST_PTR_NEXT] = (double)new;
      WDB[new + LIST_PTR_PREV] = (double)last;
      WDB[new + LIST_PTR_NEXT] = NULLPTR;
    }

  /* Reset direct pointer */

  WDB[new + LIST_PTR_DIRECT] = NULLPTR;

  /* Put pointer to common data */
  
  WDB[new + LIST_PTR_COMMON] = (double)loc0;

  /* Update number of items */

  WDB[loc0 + LIST_COMMON_N_ITEMS] = RDB[loc0 + LIST_COMMON_N_ITEMS] + 1.0;

  /* Put pointer to last item */

  WDB[loc0 + LIST_COMMON_PTR_LAST] = (double)new;

  /* Return pointer */

  return new;
}

/*****************************************************************************/
