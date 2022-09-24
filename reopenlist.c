/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : reopenlist.c                                   */
/*                                                                           */
/* Created:       2017/09/06 (JLe)                                           */
/* Last modified: 2017/09/06 (JLe)                                           */
/* Version:       2.1.30                                                     */
/*                                                                           */
/* Description: - Re-opens closed list                                       */
/*                                                                           */
/* Comments: - Should not be used repeatedly (wastes memory)                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CloseList:"

/*****************************************************************************/

void ReopenList(long lst)
{
  /* Check pointer */
  
  CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);

  /* Go to first item */

  lst = FirstItem(lst);
  
  /* Check if list is open */

  if ((long)RDB[lst + LIST_PTR_DIRECT] < VALID_PTR)
    Die(FUNCTION_NAME, "Trying to re-open an open list");

  /* Loop over items and set pointers */

  while (lst > VALID_PTR)
    {
      /* Reset direct pointer */

      WDB[lst + LIST_PTR_DIRECT] = 0.0;
      
      /* Next */

      lst = NextItem(lst);
    }
}

/*****************************************************************************/
