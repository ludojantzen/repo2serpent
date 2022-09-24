/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : replaceitem.c                                  */
/*                                                                           */
/* Created:       2012/07/14 (JLe)                                           */
/* Last modified: 2012/07/14 (JLe)                                           */
/* Version:       2.1.8                                                      */
/*                                                                           */
/* Description: Replaces list item 1 with item 2  ms                         */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReplaceItem:"

/*****************************************************************************/

void ReplaceItem(long ptr1, long ptr2)
{
  long n, sz;

  /* Check pointers */

  CheckPointer(FUNCTION_NAME, "(ptr1)", DATA_ARRAY, ptr1);
  CheckPointer(FUNCTION_NAME, "(ptr2)", DATA_ARRAY, ptr2);

  /* Get size */

  sz = (long)RDB[ptr1 + LIST_COMMON_ITEM_SIZE];

  /* Compare sizes */

  if (sz != (long)RDB[ptr2 + LIST_COMMON_ITEM_SIZE])
    Die(FUNCTION_NAME, "Mismatch in item size");

  /* Loop over structure and replace */
  
  for (n = LIST_DATA_SIZE; n < sz; n++)
    WDB[ptr1 + n] = RDB[ptr2 + n];
}

/*****************************************************************************/
