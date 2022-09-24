/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : listptr.c                                      */
/*                                                                           */
/* Created:       2010/09/15 (JLe)                                           */
/* Last modified: 2011/10/28 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Returns pointer to list item                                 */
/*                                                                           */
/* Comments: -This is replaced by a macro if compiled without DEBUG          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ListPtr:"

/*****************************************************************************/

#ifdef DEBUG

long ListPtr(long ptr, long idx)
{
  long loc0, loc1;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Pointer to direct */

  if ((loc0 = (long)RDB[ptr + LIST_PTR_DIRECT]) < VALID_PTR)
    Die(FUNCTION_NAME, "Trying to retrieve item from an open list");

  /* Pointer to common */

  loc1 = (long)RDB[ptr + LIST_PTR_COMMON];

  /* Check index */

  if ((idx < 0) || (idx > (long)RDB[loc1 + LIST_COMMON_N_ITEMS]))
    Die(FUNCTION_NAME, "Invalid index %ld", idx);

  /* Return pointer (+1 because the first pointer is to last item) */

  return (long)RDB[loc0 + idx + 1];
}

#endif

/*****************************************************************************/
