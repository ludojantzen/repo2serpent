/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : duplicateitem.c                                */
/*                                                                           */
/* Created:       2010/09/18 (JLe)                                           */
/* Last modified: 2019/10/26 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: - Duplicates item                                            */
/*                                                                           */
/* Comments: - New item is placed next from the original one                 */
/*           - Direct pointers must be redefined afterwards                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DuplicateItem:"

/*****************************************************************************/

long DuplicateItem(long ptr)
{
  long new, loc0, sz, next;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Check if list is closed */

  if ((long)RDB[ptr + LIST_PTR_DIRECT] > VALID_PTR)
    Die(FUNCTION_NAME, "Trying to duplicate an item in a closed list");

  /* Get pointer to common data */

  loc0 = (long)RDB[ptr + LIST_PTR_COMMON];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Get data size */

  sz = (long)RDB[loc0 + LIST_COMMON_ITEM_SIZE];

  /* Allocate memory */

  new = ReallocMem(DATA_ARRAY, sz);

  /* Copy data */

  memcpy(&WDB[new], &RDB[ptr], sz*sizeof(double));

  /* Get pointer next item */

  next = NextItem(ptr);

  /* Put pointers */

  WDB[new + LIST_PTR_PREV] = (double)ptr;
  WDB[new + LIST_PTR_NEXT] = (double)next;
  WDB[ptr + LIST_PTR_NEXT] = (double)new;

  if (next > VALID_PTR)
    WDB[next + LIST_PTR_PREV] = (double)new;
  else
    WDB[loc0 + LIST_COMMON_PTR_LAST] = (double)new;

  /* Update number of items */

  WDB[loc0 + LIST_COMMON_N_ITEMS] = RDB[loc0 + LIST_COMMON_N_ITEMS] + 1.0;

  /* Return pointer to new item */

  return new;
}

/*****************************************************************************/
