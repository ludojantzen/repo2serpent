/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : additem.c                                      */
/*                                                                           */
/* Created:       2011/03/09 (JLe)                                           */
/* Last modified: 2013/02/04 (JLe)                                           */
/* Version:       2.1.12                                                     */
/*                                                                           */
/* Description: Adds an existing item in list                                */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AddItem:"

/*****************************************************************************/

void AddItem(long root, long new)
{
  long ptr, loc0, last;

  /* Check root pointer */

  if (root < 1)
    Die(FUNCTION_NAME, "Pointer error %ld", root);

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(new)", DATA_ARRAY, new);

  /* Check that list exists */
  
  if ((ptr = (long)RDB[root]) < VALID_PTR)
    {
      /* Adding items to an empty list doesn't work because pointer */
      /* to common data is not available. */

      Die(FUNCTION_NAME, "List is empty");
    }

  /* Data exists, check pointer */
      
  CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

  /* Get pointer to common data */
  
  loc0 = (long)RDB[ptr + LIST_PTR_COMMON];
  CheckPointer(FUNCTION_NAME, "loc0", DATA_ARRAY, loc0);
  
#ifdef DEBUG

  /* Get pointer to new common item (use another pointer to avoid compiler */
  /* warnings about unused variables when not compiled in debug mode) */

  last = (long)RDB[new + LIST_PTR_COMMON];
  CheckPointer(FUNCTION_NAME, "(last)", DATA_ARRAY, ptr);
  
  /* Compare root pointer */
  
  if (root != (long)RDB[loc0 + LIST_COMMON_PTR_ROOT])
    Die(FUNCTION_NAME, "Mismatch in root pointer %ld %ld", 
        root, (long)RDB[loc0 + LIST_COMMON_PTR_ROOT]);
  
  /* Compare size */
  
  if ((long)RDB[last + LIST_COMMON_ITEM_SIZE] != 
      (long)RDB[loc0 + LIST_COMMON_ITEM_SIZE])
    Die(FUNCTION_NAME, "Mismatch in item size");
  
  /* Check if list is closed */
  
  if ((long)RDB[ptr + LIST_PTR_DIRECT] > VALID_PTR)
    Die(FUNCTION_NAME, "Trying to add an item to a closed list");
  
#endif
  
  /* Get pointer to last item */
  
  last = LastItem(ptr);
  CheckPointer(FUNCTION_NAME, "last", DATA_ARRAY, last);
  
  /* Set list pointers */
  
  WDB[last + LIST_PTR_NEXT] = (double)new;
  WDB[new + LIST_PTR_PREV] = (double)last;
  WDB[new + LIST_PTR_NEXT] = NULLPTR;

  /* Reset direct pointer */

  WDB[new + LIST_PTR_DIRECT] = NULLPTR;

  /* Put pointer to common data */
  
  WDB[new + LIST_PTR_COMMON] = (double)loc0;

  /* Update number of items */

  WDB[loc0 + LIST_COMMON_N_ITEMS] = RDB[loc0 + LIST_COMMON_N_ITEMS] + 1.0;

  /* Put pointer to last item */

  WDB[loc0 + LIST_COMMON_PTR_LAST] = (double)new;
}

/*****************************************************************************/
