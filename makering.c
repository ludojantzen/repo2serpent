/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : makering.c                                     */
/*                                                                           */
/* Created:       2011/03/16 (JLe)                                           */
/* Last modified: 2011/11/11 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Converts linked list into a ring                             */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MakeRing:"

/*****************************************************************************/

void MakeRing(long ptr)
{
  long first, last;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Close list */

  CloseList(ptr);

  /* Get pointer to first and last */

  first = FirstItem(ptr);
  last  = LastItem(ptr);
  
  /* Combine pointers */

  WDB[first + LIST_PTR_PREV] = (double)last;
  WDB[last + LIST_PTR_NEXT] = (double)first;
}

/*****************************************************************************/
