/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : moveitemlast.c                                 */
/*                                                                           */
/* Created:       2019/10/25 (JLe)                                           */
/* Last modified: 2019/10/25 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: - Moves item last in list                                    */
/*                                                                           */
/* Comments: - Direct pointers must be redefined afterwards                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MoveItemLast:"

/*****************************************************************************/

void MoveItemLast(long ptr)
{
  long next, prev, last, loc0, root;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get pointer to common data */

  loc0 = (long)RDB[ptr + LIST_PTR_COMMON];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Get pointers to next and previous items */

  next = NextItem(ptr);
  prev = PrevItem(ptr);

  /* Check if already last and put pointer */

  if (next < VALID_PTR)
    return;
  else
    WDB[next + LIST_PTR_PREV] = (double)prev;

  /* Check if first and put pointer */

  if (prev < VALID_PTR)
    WDB[loc0 + LIST_COMMON_PTR_FIRST] = (double)next;
  else
    WDB[prev + LIST_PTR_NEXT] = (double)next;

  /* Get pointer to last */

  last = (long)RDB[loc0 + LIST_COMMON_PTR_LAST];
  CheckPointer(FUNCTION_NAME, "(last)", DATA_ARRAY, last);

  /* Put pointers */

  WDB[last + LIST_PTR_NEXT] = (double)ptr;
  WDB[ptr + LIST_PTR_PREV] = (double)last;
  WDB[ptr + LIST_PTR_NEXT] = NULLPTR;
  WDB[loc0 + LIST_COMMON_PTR_LAST] = (double)ptr;

  /* Get pointer to root */

  root = (long)RDB[loc0 + LIST_COMMON_PTR_ROOT];

  /* Put pointer */

  WDB[root] = RDB[loc0 + LIST_COMMON_PTR_FIRST];
}

/*****************************************************************************/
