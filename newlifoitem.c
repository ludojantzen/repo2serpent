/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : newlifoitem.c                                  */
/*                                                                           */
/* Created:       2010/09/15 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Adds a new item in a simplified one-way list                 */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "NewLIFOItem:"

/*****************************************************************************/

long NewLIFOItem(long root, long sz)
{
  long ptr, new;

  /* Check root pointer */

  if (root < 1)
    Die(FUNCTION_NAME, "Pointer error");

  /* Allocate memory for new item */

  new = ReallocMem(DATA_ARRAY, sz);

  /* Check if list exists */
  
  if ((ptr = (long)RDB[root]) < VALID_PTR)
    {
      /* First definition, set pointer */
      
      WDB[root] = (double)new;

      /* Reset pointer to next */

      WDB[new + LIFO_LIST_PTR_NEXT] = NULLPTR;
    }
  else
    {
      /* Data exists, check pointer */
      
      CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

      /* Reconfigure pointers */

      WDB[root] = (double)new;
      WDB[new + LIFO_LIST_PTR_NEXT] = (double)ptr;
    }

  /* Return pointer */

  return new;
}

/*****************************************************************************/
