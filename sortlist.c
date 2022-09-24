/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : sortlist.c                                     */
/*                                                                           */
/* Created:       2010/09/18 (JLe)                                           */
/* Last modified: 2019/10/26 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: - Sorts list items                                           */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SortList:"

/*****************************************************************************/

void SortList(long lst, long param, long mode)
{
  long ptr, ptp, loc0, root, next, count, sz, i, j;
  double val1, val2, *arr;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);

  /* Check if access to private data is allowed */

  if ((mode == SORT_MODE_ASCEND_PRIVA) || (mode == SORT_MODE_DESCEND_PRIVA))
    if ((long)RDB[DATA_PRIVA_MEM_READY] == NO)
      Die(FUNCTION_NAME, "PRIVA array not ready for access");

  /* Get size */

  sz = ListSize(lst);

  /* Check sorted list type and size */

  if (((mode == SORT_MODE_ASCEND) || (mode == SORT_MODE_DESCEND)) && (sz > 10))
    {
      /**********************************************************************/

      /***** Insertion sort with temporary pointer array ********************/

      /* Use direct pointer list or allocate memory for pointer */
      /* array. NOTE: Direct pointers are reconfigured by call  */
      /* to SetDirectPointer() at the end of the subroutine.    */

      if ((ptr = (long)RDB[lst + LIST_PTR_DIRECT]) > VALID_PTR)
        arr = &WDB[ptr];
      else
        arr = (double *)Mem(MEM_ALLOC, sz, sizeof(double));

      /* Reset index and get pointer to first list item */

      i = 0;
      ptr = FirstItem(lst);

     /* Read data and pointers */

      while (ptr > VALID_PTR)
        {
          /* Read pointers */

          arr[i] = (double)ptr;

          /* Update index */

          i++;

          /* Next */

          ptr = NextItem(ptr);
        }

      /* Check size */

      if (i != sz)
        Die(FUNCTION_NAME, "Mismatch in size");

      /* Perform insertion sort */

      for (i = 1; i < sz; i++)
        {
          j = i;
          ptr = (long)arr[j];

          if (mode == SORT_MODE_ASCEND)
            {
              while (j > 0 && RDB[(long)arr[j - 1] + param] > RDB[ptr + param])
                {
                  arr[j] = arr[j - 1];
                  j--;
                }
            }
          else if (mode == SORT_MODE_DESCEND)
            {
              while (j > 0 && RDB[(long)arr[j - 1] + param] < RDB[ptr + param])
                {
                  arr[j] = arr[j - 1];
                  j--;
                }
            }
          else
            Die(FUNCTION_NAME, "Invalid sort mode");

          arr[j] = (double)ptr;
        }

      /* Get pointer to common data */

      loc0 = (long)RDB[lst + LIST_PTR_COMMON];
      CheckPointer(FUNCTION_NAME, "loc0", DATA_ARRAY, loc0);

      /* Get root pointer (CheckPointer() cannot be used here) */

      root = (long)RDB[loc0 + LIST_COMMON_PTR_ROOT];

      /* Put root pointer and pointer to first and last item */

      WDB[root] = arr[0];
      WDB[loc0 + LIST_COMMON_PTR_FIRST] = arr[0];
      WDB[loc0 + LIST_COMMON_PTR_LAST] = arr[sz - 1];

      /* Put pointers */

      for (i = 0; i < sz; i++)
        {
          if (sz == 1)
            {
              /* Single-valued */

              WDB[(long)arr[i] + LIST_PTR_PREV] = NULLPTR;
              WDB[(long)arr[i] + LIST_PTR_NEXT] = NULLPTR;
            }
          else if (i == 0)
            {
              /* First item */

              WDB[(long)arr[i] + LIST_PTR_PREV] = NULLPTR;
              WDB[(long)arr[i] + LIST_PTR_NEXT] = arr[i + 1];
            }
          else if (i == sz - 1)
            {
              /* Last item */

              WDB[(long)arr[i] + LIST_PTR_PREV] = arr[i - 1];
              WDB[(long)arr[i] + LIST_PTR_NEXT] = NULLPTR;
            }
          else
            {
              /* Other */

              WDB[(long)arr[i] + LIST_PTR_PREV] = arr[i - 1];
              WDB[(long)arr[i] + LIST_PTR_NEXT] = arr[i + 1];
            }
        }

      /* Free allocated memory */

      if ((long)RDB[lst + LIST_PTR_DIRECT]  < VALID_PTR)
        Mem(MEM_FREE, arr);

      /**********************************************************************/
    }
  else
    {
      /**********************************************************************/

      /***** Bubble-sort on list ********************************************/

      /* Go to first item */

      lst = FirstItem(lst);

      /* Sort loop */

      do
        {
          /* Reset count */

          count = 0;

          /* Check mode */

          if ((mode == SORT_MODE_ASCEND) || (mode == SORT_MODE_ASCEND_PRIVA))
            {
              /***************************************************************/

              /***** Sort items in ascending order ***************************/

              /* Set starting point for iteration */

              loc0 = FirstItem(lst);

              /* Sorting loop */

              while (loc0 > VALID_PTR)
                {
                  /* Set pointers for current and next */

                  ptr = loc0;
                  next = NextItem(ptr);

                  /* Set starting point for next iteration */

                  loc0 = next;

                  /* Loop over list */

                  while (next > VALID_PTR)
                    {
                      /* Get values */

                      if (mode == SORT_MODE_ASCEND)
                        {
                          /* Get direct values */

                          val1 = RDB[ptr + param];
                          val2 = RDB[next + param];
                        }
                      else
                        {
                          /* Get values from private block */

                          ptp = (long)RDB[ptr + param];
                          CheckPointer(FUNCTION_NAME, "(ptp)",
                                       PRIVA_ARRAY, ptp);
                          val1 = SumPrivateData(ptp);

                          ptp = (long)RDB[next + param];
                          CheckPointer(FUNCTION_NAME, "(ptp)",
                                       PRIVA_ARRAY, ptp);
                          val2 = SumPrivateData(ptp);
                        }

                      /* Compare items, move one step forward or break loop */

                      if (val1 > val2)
                        {
                          MoveItemRight(ptr);
                          count++;
                        }
                      else
                        break;

                      /* Next */

                      next = NextItem(ptr);
                    }
                }

              /***************************************************************/
            }
          else if((mode == SORT_MODE_DESCEND) ||
                  (mode == SORT_MODE_DESCEND_PRIVA))
            {
              /***************************************************************/

              /***** Sort items in descending order **************************/

              /* Set starting point for iteration */

              loc0 = FirstItem(lst);

              /* Sorting loop */

              while (loc0 > VALID_PTR)
                {
                  /* Set pointers for current and next */

                  ptr = loc0;
                  next = NextItem(ptr);

                  /* Set starting point for next iteration */

                  loc0 = next;

                  /* Loop over list */

                  while (next > VALID_PTR)
                    {
                      /* Get values */

                      if (mode == SORT_MODE_DESCEND)
                        {
                          /* Get direct values */

                          val1 = RDB[ptr + param];
                          val2 = RDB[next + param];
                        }
                      else
                        {
                          /* Get values from private block */

                          ptp = (long)RDB[ptr + param];
                          CheckPointer(FUNCTION_NAME, "(ptp)",
                                       PRIVA_ARRAY, ptp);
                          val1 = SumPrivateData(ptp);

                          ptp = (long)RDB[next + param];
                          CheckPointer(FUNCTION_NAME, "(ptp)",
                                       PRIVA_ARRAY, ptp);
                          val2 = SumPrivateData(ptp);
                        }

                      /* Compare items, move one step forward or break loop */

                      if (val1 < val2)
                        {
                          MoveItemRight(ptr);
                          count++;
                        }
                      else
                        break;

                      /* Next */

                      next = NextItem(ptr);
                    }
                }

              /***************************************************************/
            }
          else
            Die(FUNCTION_NAME, "Invalid sort mode");
        }
      while (count > 0);

      /***********************************************************************/
    }

#ifdef DEBUG

  /* Check order */

  ptr = FirstItem(lst);

  while (ptr > 0)
    {
      if ((next = NextItem(ptr)) > 0)
        {
          /* Get values */

          if ((mode == SORT_MODE_ASCEND) || (mode == SORT_MODE_DESCEND))
            {
              /* Get direct values */

              val1 = RDB[ptr + param];
              val2 = RDB[next + param];
            }
          else
            {
              /* Get values from private block */

              ptp = (long)RDB[ptr + param];
              CheckPointer(FUNCTION_NAME, "(ptp)", PRIVA_ARRAY, ptp);
              val1 = SumPrivateData(ptp);

              ptp = (long)RDB[next + param];
              CheckPointer(FUNCTION_NAME, "(ptp)", PRIVA_ARRAY, ptp);
              val2 = SumPrivateData(ptp);
            }

          if (((mode == SORT_MODE_ASCEND) || (mode == SORT_MODE_ASCEND_PRIVA))
              && (val1 > val2))
            Die(FUNCTION_NAME, "Sort error");
          else if (((mode == SORT_MODE_DESCEND) ||
                    (mode == SORT_MODE_DESCEND_PRIVA))
                   && (val1 < val2))
            Die(FUNCTION_NAME, "Sort error");
        }
      else
        break;

      /* Next */

      ptr = next;
    }

#endif

  /* Set direct pointers if list is closed */

  if ((long)RDB[lst + LIST_PTR_DIRECT] > VALID_PTR)
    SetDirectPointers(lst);
}

/*****************************************************************************/
