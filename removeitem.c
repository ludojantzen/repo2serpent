/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : removeitem.c                                   */
/*                                                                           */
/* Created:       2010/09/18 (JLe)                                           */
/* Last modified: 2011/11/11 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: - Removes item from list                                     */
/*                                                                           */
/* Comments:  - Direct pointers must be redefined afterwards                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "RemoveItem:"

/*****************************************************************************/

void RemoveItem(long ptr)
{
  long prev, next, loc0, root;
  
  /* Check pointer */
  
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  
  /* Check if list is closed */

  if ((long)RDB[ptr + LIST_PTR_DIRECT] > VALID_PTR)
    Die(FUNCTION_NAME, "Trying to remove an item from a closed list");

  /* Get pointers to previous and next items */

  prev = PrevItem(ptr);
  next = NextItem(ptr);

  /* Put pointers */

  if (prev > VALID_PTR)
    WDB[prev + LIST_PTR_NEXT] = (double)next;

  if (next > VALID_PTR)
    WDB[next + LIST_PTR_PREV] = (double)prev;

  /* Get pointer to common data */

  loc0 = (long)RDB[ptr + LIST_PTR_COMMON];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);
  
  /* Get pointer to root */
  
  root = (long)RDB[loc0 + LIST_COMMON_PTR_ROOT];

  /* Check root pointer */

  if (root < 1)
    Die(FUNCTION_NAME, "Pointer error");
  
  /* Replace pointers if item was first */

  if (prev < VALID_PTR)
    {
      WDB[root] = (double)next;
      WDB[loc0 + LIST_COMMON_PTR_FIRST] = (double)next;
    }

  /* Replace last pointer if item was last */

  if (next < VALID_PTR)
    WDB[loc0 + LIST_COMMON_PTR_LAST] = (double)prev;

  /* Update number of items */

  WDB[loc0 + LIST_COMMON_N_ITEMS] = RDB[loc0 + LIST_COMMON_N_ITEMS] - 1.0;

  /* Reset root pointer if list is empty */

  if ((long)RDB[loc0 + LIST_COMMON_N_ITEMS] == 0)
    WDB[root] = NULLPTR;
}

/*****************************************************************************/
