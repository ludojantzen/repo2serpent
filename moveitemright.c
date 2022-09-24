/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : moveitemright.c                                */
/*                                                                           */
/* Created:       2010/09/17 (JLe)                                           */
/* Last modified: 2011/11/11 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: - Moves item one position right (forward) in list            */
/*                                                                           */
/* Comments: - Direct pointers must be redefined afterwards                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MoveItemRight:"

/*****************************************************************************/

void MoveItemRight(long ptr)
{
  long next, left, right, loc0, root;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get pointer to next */
  
  if ((next = NextItem(ptr)) == NULLPTR)
    Die(FUNCTION_NAME, "Item is already at the end of list");

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(next)", DATA_ARRAY, next);

  /* Items before and after the two */

  left = PrevItem(ptr);
  right = NextItem(next);

  /* Put pointers */

  if (left > VALID_PTR)
    WDB[left + LIST_PTR_NEXT] = (double)next;

  if (right > VALID_PTR)
    WDB[right + LIST_PTR_PREV] = (double)ptr;

  WDB[ptr + LIST_PTR_NEXT] = (double)right;
  WDB[ptr + LIST_PTR_PREV] = (double)next;

  WDB[next + LIST_PTR_NEXT] = (double)ptr;
  WDB[next + LIST_PTR_PREV] = (double)left;

  /* Get pointer to common data */

  loc0 = (long)RDB[ptr + LIST_PTR_COMMON];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Get pointer to root */

  root = (long)RDB[loc0 + LIST_COMMON_PTR_ROOT];

  /* Check root pointer */

  if (root < 1)
    Die(FUNCTION_NAME, "Pointer error");

  /* Check if item was the first item and replace pointer */

  if ((long)RDB[root] == ptr)
    {
      WDB[root] = (double)next;
      WDB[loc0 + LIST_COMMON_PTR_FIRST] = (double)next;
    }

  /* Check if last item was replaced */
  
  if (NextItem(ptr) < VALID_PTR)
    WDB[loc0 + LIST_COMMON_PTR_LAST] = (double)ptr;
}

/*****************************************************************************/
