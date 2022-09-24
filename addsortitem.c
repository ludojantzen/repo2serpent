/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : addsortitem.c                                  */
/*                                                                           */
/* Created:       2012/09/20 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                      */
/*                                                                           */
/* Description: Adds an existing item in a sorted list                       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AddSortItem:"

/*****************************************************************************/

void AddSortItem(long root, long pos, long new, long param, long mode)
{
  long loc0, next, prev;

  /* Check pointers */

  CheckPointer(FUNCTION_NAME, "(new)", DATA_ARRAY, new);
  CheckPointer(FUNCTION_NAME, "(pos)", DATA_ARRAY, pos);

  /* Check root pointer */

  if (root < 1)
    Die(FUNCTION_NAME, "Pointer error %ld", root);

  /* Check that list exists */
  
  if ((loc0 = (long)RDB[root]) < VALID_PTR)
    {
      /* Adding items to an empty list doesn't work because pointer */
      /* to common data is not available. */

      Die(FUNCTION_NAME, "List is empty");
    }

  /* Data exists, check pointer */
      
  CheckPointer(FUNCTION_NAME, "loc0", DATA_ARRAY, loc0);

  /* Get pointer to common data */
  
  loc0 = (long)RDB[loc0 + LIST_PTR_COMMON];
  CheckPointer(FUNCTION_NAME, "loc0", DATA_ARRAY, loc0);

  /* Avoid compiler warning */

  prev = -1;
  next = -1;

  /* Check sort mode */

  if (mode == SORT_MODE_ASCEND)
    {
      /* Check position */

      if (RDB[new + param] > RDB[pos + param])
        {
          /* Above position, loop until below */

          while (RDB[new + param] > RDB[pos + param])
            if ((pos = NextItem(pos)) < VALID_PTR)
              break;

          /* Check pointer */

          if (pos > VALID_PTR)
            {
              prev = PrevItem(pos);
              next = pos;

              /* Check */

              if (RDB[prev + param] > RDB[new + param])
                Die(FUNCTION_NAME, "Error");
            }
          else
            {
              prev = (long)RDB[loc0 + LIST_COMMON_PTR_LAST];
              next = NULLPTR;
            }
        }
      else if (RDB[new + param] <= RDB[pos + param])
        {
          /* Below position, loop until above */

          while (RDB[new + param] <= RDB[pos + param])
            if ((pos = PrevItem(pos)) < VALID_PTR)
              break;

          /* Check pointer */

          if (pos > VALID_PTR)
            {
              next = NextItem(pos);
              prev = pos;

              /* Check */

              if (RDB[next + param] < RDB[new + param])
                Die(FUNCTION_NAME, "Error");
            }
          else
            {
              next = (long)RDB[loc0 + LIST_COMMON_PTR_FIRST];
              prev = NULLPTR;
            }
        }
    }
  else
    Die(FUNCTION_NAME, "Invalid sort mode");

  /* Set list pointers */

  if (prev > VALID_PTR)
    {
      WDB[new + LIST_PTR_PREV] = (double)prev;
      WDB[prev + LIST_PTR_NEXT] = (double)new;
    }
  else
    {
      WDB[new + LIST_PTR_PREV] = NULLPTR;
      WDB[loc0 + LIST_COMMON_PTR_FIRST] = (double)new;
    }

  if (next > VALID_PTR)
    {
      WDB[new + LIST_PTR_NEXT] = (double)next;
      WDB[next + LIST_PTR_PREV] = (double)new;
    }
  else
    {
      WDB[new + LIST_PTR_NEXT] = NULLPTR;
      WDB[loc0 + LIST_COMMON_PTR_LAST] = (double)new;
    }

  /* Reset direct pointer */

  WDB[new + LIST_PTR_DIRECT] = NULLPTR;

  /* Put pointer to common data */
  
  WDB[new + LIST_PTR_COMMON] = (double)loc0;

  /* Update number of items */

  WDB[loc0 + LIST_COMMON_N_ITEMS] = RDB[loc0 + LIST_COMMON_N_ITEMS] + 1.0;
}

/*****************************************************************************/
