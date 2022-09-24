/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : lifolistsize.c                                 */
/*                                                                           */
/* Created:       2014/10/04 (JLe)                                           */
/* Last modified: 2017/06/08 (VVa)                                           */
/* Version:       2.1.29                                                     */
/*                                                                           */
/* Description: Counts the number of items in a simple one-way list          */
/*                                                                           */
/* Comments: - This is an expensive operation                                */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "LIFOListSize:"

/*****************************************************************************/

long LIFOListSize(long root, long cutoff)
{
  long ptr, sz;

  /* Reset count */

  sz = 0;

  /* Loop over list */

  ptr = (long)RDB[root];
  while (ptr > VALID_PTR)
    {
      /* Add to count */

      sz++;

      /* Check cutoff value */

      if (sz == cutoff)
        return sz;

      /* Next */

      ptr = NextItem(ptr);
    }

  /* Return list size */

  return sz;
}


/*****************************************************************************/
