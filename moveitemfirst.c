/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : moveitemfirst.c                                */
/*                                                                           */
/* Created:       2012/01/17 (JLe)                                           */
/* Last modified: 2012/01/17 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: - Moves item first in list                                   */
/*                                                                           */
/* Comments: - Direct pointers must be redefined afterwards                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MoveItemFirst:"

/*****************************************************************************/

void MoveItemFirst(long ptr)
{
  long next, prev, first, loc0, root;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get pointer to common data */

  loc0 = (long)RDB[ptr + LIST_PTR_COMMON];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Get pointers to previous and next items */

  prev = PrevItem(ptr);
  next = NextItem(ptr);

  /* Check if already first and put pointer */

  if (prev < VALID_PTR)
    return;
  else
    WDB[prev + LIST_PTR_NEXT] = (double)next;

  /* Check if last and put pointer */

  if (next < VALID_PTR)
    WDB[loc0 + LIST_COMMON_PTR_LAST] = (double)prev;
  else
    WDB[next + LIST_PTR_PREV] = (double)prev;

  /* Get pointer to first */

  first = (long)RDB[loc0 + LIST_COMMON_PTR_FIRST];
  CheckPointer(FUNCTION_NAME, "(first)", DATA_ARRAY, first);

  /* Put pointers */
  
  WDB[first + LIST_PTR_PREV] = (double)ptr;
  WDB[ptr + LIST_PTR_NEXT] = (double)first;
  WDB[ptr + LIST_PTR_PREV] = NULLPTR;
  WDB[loc0 + LIST_COMMON_PTR_FIRST] = (double)ptr;
  
  /* Get pointer to root */

  root = (long)RDB[loc0 + LIST_COMMON_PTR_ROOT];

  /* Put pointer */

  WDB[root] = (double)ptr;
}

/*****************************************************************************/
