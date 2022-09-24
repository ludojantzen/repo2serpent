/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : nextitem.c                                     */
/*                                                                           */
/* Created:       2010/09/15 (JLe)                                           */
/* Last modified: 2011/10/28 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Returns pointer to next item in list.                        */
/*                                                                           */
/* Comments: - This is replaced by a macro if compiled without DEBUG         */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "NextItem:"

/*****************************************************************************/

#ifdef DEBUG

long NextItem(long ptr)
{
  long next;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get pointer to next */
  
  next = (long)RDB[ptr + LIST_PTR_NEXT];

  /* Check self */

  if (next == ptr)
    Die(FUNCTION_NAME, "Pointer to self");

  /* Check null */

  if (next == NULLPTR)
    return next;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(next)", DATA_ARRAY, next);

  /* Return pointer to next item */

  return next;
}

#endif

/*****************************************************************************/
