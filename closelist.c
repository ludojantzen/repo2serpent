/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : closelist.c                                    */
/*                                                                           */
/* Created:       2010/09/20 (JLe)                                           */
/* Last modified: 2011/11/11 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: - Closes list and creates pointer array                      */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CloseList:"

/*****************************************************************************/

void CloseList(long lst)
{
  long ptr, loc0, loc1, sz;

  /* Check pointer */
  
  CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);

  /* Go to first item */

  lst = FirstItem(lst);
  
  /* Check if list is closed */

  if ((long)RDB[lst + LIST_PTR_DIRECT] > VALID_PTR)
    Die(FUNCTION_NAME, "Trying to close a closed list");

  /* Get pointer to common data */

  loc0 = (long)RDB[lst + LIST_PTR_COMMON];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Get data size */

  sz = (long)RDB[loc0 + LIST_COMMON_N_ITEMS];

  /* Allocate memory */

  loc1 = ReallocMem(DATA_ARRAY, sz + 2);

  /* Loop over items and set pointers */

  ptr = lst;

  while (ptr > VALID_PTR)
    {
      /* Put list pointer to item */ 

      WDB[ptr + LIST_PTR_DIRECT] = (double)loc1;
      
      /* Next */

      ptr = NextItem(ptr);
    }

  /* Set direct pointers */

  SetDirectPointers(lst);
}

/*****************************************************************************/
