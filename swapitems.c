/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : swapitems.c                                    */
/*                                                                           */
/* Created:       2012/07/09 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Swaps the contents of two list items                         */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SwapItems:"

/*****************************************************************************/

void SwapItems(long loc1, long loc2)
{
  long ptr, n, sz;
  double tmp;

  /* Check pointers */

  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
  CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

  /* Get pointer to common data */

  ptr = (long)RDB[loc1 + LIST_PTR_COMMON];
  CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

  /* Get size */

  sz = (long)RDB[ptr + LIST_COMMON_ITEM_SIZE];

  /* Compare to second */

  ptr = (long)RDB[loc2 + LIST_PTR_COMMON];
  CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);
  
  if (sz != (long)RDB[ptr + LIST_COMMON_ITEM_SIZE])
    Die(FUNCTION_NAME, "Mismatch in item size");

  /* Loop over structure */
  
  for (n = 0; n < sz; n++)
    {
      /* Swap contents */

      tmp = RDB[loc1 + n];
      WDB[loc1 + n] = RDB[loc2 + n];
      WDB[loc2 + n] = tmp;
    }
  
  /* Change direct pointer in first item */

  if ((n = (long)RDB[loc1 + LIST_PTR_DIRECT]) > VALID_PTR)
    {
      /* Loop over pointers */
      
      while ((ptr = (long)RDB[n]) > VALID_PTR)
        {
          /* Compare */
          
          if (ptr == loc1)
            {
              /* Change */

              WDB[n] = (double)loc2;

              /* Break loop */
            }

          /* Next */

          n++;
        }

      /* Check that pointer was found */

      if (ptr < VALID_PTR)
        Die(FUNCTION_NAME, "Direct pointer 1 not found");
    }

  /* Change direct pointer in second item */

  if ((n = (long)RDB[loc2 + LIST_PTR_DIRECT]) > VALID_PTR)
    {
      /* Loop over pointers */
      
      while ((ptr = (long)RDB[n]) > VALID_PTR)
        {
          /* Compare */
          
          if (ptr == loc2)
            {
              /* Change */
              
              WDB[n] = (double)loc1;

              /* Break loop */
            }

          /* Next */

          n++;
        }

      /* Check that pointer was found */

      if (ptr < VALID_PTR)
        Die(FUNCTION_NAME, "Direct pointer 2 not found");
    }
}

/*****************************************************************************/
