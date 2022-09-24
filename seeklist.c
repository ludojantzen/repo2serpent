/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : seeklist.c                                     */
/*                                                                           */
/* Created:       2010/09/20 (JLe)                                           */
/* Last modified: 2019/12/05 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: - Seeks item with matching value                             */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SeekList:"

/*****************************************************************************/

long SeekList(long loc0, long param, double val, long mode)
{
  long n, N, ptr, i;
  double hi, lo, try;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Get number of items */

  N = ListSize(loc0);

  /* Check sort mode */

  if ((mode == NO) || (N < 4))
    {
      /***********************************************************************/

      /***** Small or unsorted list ******************************************/

      /* Move to beginning */

      loc0 = FirstItem(loc0);

      /* Loop over list */

      while (loc0 > VALID_PTR)
        {
          /* Compare */

          if (RDB[loc0 + param] == val)
            return loc0;

          /* Next */

          loc0 = NextItem(loc0);
        }

      /* Not found */

      return -1;

      /***********************************************************************/
    }
  else if (mode == SORT_MODE_ASCEND)
    {
      /***********************************************************************/

      /***** Sorted in ascending order ***************************************/

      /* Init boundaries */

      lo = 0.0;
      hi = (double)(N - 1);

      /* Search loop */

      for (i = 0; i < N; i++)
        {
          /* Check for small interval */

          if ((long)hi - (long)lo < 10)
            {
              /* Loop over interval */

              for (n = (long)lo; n < (long)hi + 1; n++)
                {
                  /* Get pointer to list item */

                  ptr = ListPtr(loc0, n);
                  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

                  /* Check hit */

                  if (val == RDB[ptr + param])
                    return ptr;
                }

              /* Not found */

              return -1;
            }

          /* New guess */

          try = (hi + lo)/2.0;
          n = (long)try;

          /* Get pointer to list item */

          ptr = ListPtr(loc0, n);
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          /* Check hit */

          if (val == RDB[ptr + param])
            return ptr;

          /* New boundaries */

          if (val < RDB[ptr + param])
            {
              /* Update upper boundary */

              hi = try;
            }
          else
            {
              /* Update lower boundary */

              lo = try;
            }
        }

      /***********************************************************************/
    }
  else if (mode == SORT_MODE_DESCEND)
    {
      /***********************************************************************/

      /***** Sorted in descending order ***************************************/

      Die(FUNCTION_NAME, "tää ei toimi");

      /***********************************************************************/
    }
  else
    Die(FUNCTION_NAME, "Invalid sort mode");

  /* Something wrong here */

  Die(FUNCTION_NAME, "Horrible error");

  /* Dummy return value to please compiler */

  return -1;
}

/*****************************************************************************/
