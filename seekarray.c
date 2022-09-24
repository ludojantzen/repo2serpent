/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : seeklist.c                                     */
/*                                                                           */
/* Created:       2010/04/05 (JLe)                                           */
/* Last modified: 2019/04/05 (JLe)                                           */
/* Version:       2.1.31                                                     */
/*                                                                           */
/* Description: - Finds a value from an array.                               */
/*                                                                           */
/* Comments: - Based on binary search.                                       */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SeekArray:"

/*****************************************************************************/

long SeekArray(const double *dat, long val, long N)
{
  long lo, hi, n;

  /* Check boundaries */

  if ((val < (long)dat[0]) || (val > (long)dat[N - 1]))
    return -1;
  else if (val == (long)dat[0])
    return 0;

  /* Init boundaries */

  lo = 0;
  hi = N - 1;

  /* Search loop */

  while (1 != 2)
    {
      /* Check for small interval */

      if (hi - lo < 10)
        {
          /* Loop over interval */

          for (n = lo; n < hi + 1; n++)
            {
              /* Check value */

              CheckValue(FUNCTION_NAME, "n", "", n , 0, N - 1);

              /* Check hit */

              if (val == (long)dat[n])
                return n;
            }

          /* Not found */

          return -1;
        }

      /* New guess */

      n = (long)(((double)hi + (double)lo)/2.0);
      CheckValue(FUNCTION_NAME, "n", "", n , 0, N - 1);

      /* Check hit */

      if (val == (long)dat[n])
        return n;

      /* New boundaries */

      if (val < (long)dat[n])
        {
          /* Update upper boundary */

          hi = n;
        }
      else
        {
          /* Update lower boundary */

          lo = n;
        }
    }

  /* Dummy return value to please compiler */

  return -1;
}

/*****************************************************************************/
