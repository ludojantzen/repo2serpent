/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : previtem.c                                     */
/*                                                                           */
/* Created:       2010/09/15 (JLe)                                           */
/* Last modified: 2011/10/28 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Returns pointer to previous item in list.                    */
/*                                                                           */
/* Comments: -This is replaced by a macro if compiled without DEBUG          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrevItem:"

/*****************************************************************************/

#ifdef DEBUG

long PrevItem(long ptr)
{
  long prev;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get pointer to previous */
  
  prev = (long)RDB[ptr + LIST_PTR_PREV];

  /* Check self */

  if (prev == ptr)
    Die(FUNCTION_NAME, "Pointer to self");

  /* Check null */

  if (prev == NULLPTR)
    return prev;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(prev)", DATA_ARRAY, prev);

  /* Return pointer to previous item */

  return prev;
}

#endif

/*****************************************************************************/
