/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : splitlist.c                                    */
/*                                                                           */
/* Created:       2019/10/23 (JLe)                                           */
/* Last modified: 2019/11/13 (JLe)                                           */
/* Version:       2.1.32                                                     */
/*                                                                           */
/* Description: - Splits list into two parts depending on given value        */
/*                                                                           */
/* Comments: - Tämä tehtiin RMX-solveria varten, mutta muista overheadeista  */
/*             johtuen tulose ei ollut hyvä. Voi olla että feilaa jossain    */
/*             tapauksissa.                                                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SplitList:"

/*****************************************************************************/

void SplitList(long lst, long param, double lim, long mode)
{
  long ptr, loc0, loc1, next, sz, n0, n1, n, count;
  double val, *arr;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);

  /* Check if access to private data is allowed */

  if ((mode == SORT_MODE_ASCEND_PRIVA) || (mode == SORT_MODE_DESCEND_PRIVA))
    if ((long)RDB[DATA_PRIVA_MEM_READY] == NO)
      Die(FUNCTION_NAME, "PRIVA array not ready for access");

  /* Get list size */

  sz = ListSize(lst);

  /* Use direct pointer list or allocate memory for pointer */
  /* array. NOTE: Direct pointers are reconfigured by call  */
  /* to SetDirectPointer() at the end of the subroutine.    */

  if ((ptr = (long)RDB[lst + LIST_PTR_DIRECT]) > VALID_PTR)
    arr = &WDB[ptr];
  else
    arr = (double *)Mem(MEM_ALLOC, sz, sizeof(double));

  /* Reset sizes */

  n0 = 0;
  n1 = 0;
  n = 0;

  /* Avoid compiler warning */

  val = 0.0;

  /* Read data into array */

  loc1 = lst;
  while (loc1 > VALID_PTR)
    {
      /* Get values */

      if ((mode == SORT_MODE_ASCEND) || (mode == SORT_MODE_DESCEND))
        {
          /* Get direct value */

          val = RDB[loc1 + param];
        }
      else if ((mode == SORT_MODE_ASCEND_PRIVA) ||
               (mode == SORT_MODE_DESCEND_PRIVA))
        {
          /* Get values from private block */

          ptr = (long)RDB[loc1 + param];
          CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
          val = SumPrivateData(ptr);
        }
      else
        Die(FUNCTION_NAME, "Invalid mode");

      /* Put value in array */

      arr[n++] = val;

      /* Compare limits */

      if (val < lim)
        n0++;
      else
        n1++;

      /* Next */

      loc1 = NextItem(loc1);
    }

  /* Sanity check */

  if (n0 + n1 != sz)
    Die(FUNCTION_NAME, "Mismatch in size");

  /* Mark moved items */

  for (n = 0; n < sz; n++)
    {
      /* Check index */

      if (n < n0)
        {
          /* Compare to limit */

          if (arr[n] >= lim)
            arr[n] = 1.0;
          else
            arr[n] = 0.0;
        }
      else
        {
          /* Compare to limit */

          if (arr[n] < lim)
            arr[n] = -1.0;
          else
            arr[n] = 0.0;
        }
    }

  /* Reset index and count */

  n = 0;
  count = 0;

  /* Loop over list */

  loc1 = lst;
  while (loc1 > VALID_PTR)
    {
      /* Pointer to next */

      next = NextItem(loc1);

      /* Check if moved */

      if (arr[n] < 0.0)
        {
          /* Move to beginning */

          MoveItemFirst(loc1);

          /* add to counter */

          count++;
        }
      else if (arr[n] > 0.0)
        {
          /* Move to end */

          MoveItemLast(loc1);

          /* add to counter */

          count++;
        }

      /* Update index */

      n++;

      /* Check */

      if (n == sz)
        break;

      /* Next item */

      loc1 = next;
    }

  /* Free allocated memory */

  if ((long)RDB[lst + LIST_PTR_DIRECT]  < VALID_PTR)
    Mem(MEM_FREE, arr);

  /* Loop over list and check size */

  loc0 = FirstItem(lst);
  while (loc0 > VALID_PTR)
    {
      /* Update size */

      sz--;

      /* next */

      loc0 = NextItem(loc0);
    }

  /* Check */

  if (sz != 0)
    Die(FUNCTION_NAME, "Error in size %ld", sz);

  /* Set direct pointers if list is closed */

  if ((long)RDB[lst + LIST_PTR_DIRECT] > VALID_PTR)
    SetDirectPointers(lst);
}

/*****************************************************************************/
