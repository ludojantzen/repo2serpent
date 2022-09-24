/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : lastitem.c                                     */
/*                                                                           */
/* Created:       2010/09/15 (JLe)                                           */
/* Last modified: 2011/10/28 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Returns pointer to last item in list or -1 if list doesn't   */
/*              exist.                                                       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "LastItem:"

/*****************************************************************************/

long LastItem(long ptr)
{
  long loc0;

  /* Check pointer */

  if (ptr < VALID_PTR)
    return -1;

  /* Get pointer to common data */

  loc0 = (long)RDB[ptr + LIST_PTR_COMMON];

  /* Get pointer to last */

  ptr = (long)RDB[loc0 + LIST_COMMON_PTR_LAST];

  /*  Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr); 

  /* Return pointer to last item */

  return ptr;
}

/*****************************************************************************/
