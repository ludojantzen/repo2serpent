/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : addpts.c                                       */
/*                                                                           */
/* Created:       2010/09/12 (JLe)                                           */
/* Last modified: 2014/01/23 (JLe)                                           */
/* Version:       2.1.17                                                     */
/*                                                                           */
/* Description: Adds array of doubles to another array                       */
/*                                                                           */
/* Comments: - From Serpent 1.1.0                                            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "AddPts:"

/*****************************************************************************/

double *AddPts(double *A1, long *norig, const double *A2, long sz2)
{
  long i, j, k, sz1;
  double *new;

  /* Get initial size of A1 */

  sz1 = *norig;

  /* Allocate memory for new data */

  new = (double *)Mem(MEM_ALLOC,sz1 + sz2, sizeof(double));
  
  /* Reset indexes */

  i = 0;
  j = 0;
  k = 0;

  /* Loop over old and added arrays */

  while ((i < sz1) || (j < sz2))
    {
      /* Check if either vector has been completely read */

      if (i == sz1)
        {
          /* Add value from vector 2 */

          if ((k == 0) || (new[k - 1] != A2[j]))
            new[k++] = A2[j];

          /* Update index */

          j++;
        }
      
      else if (j == sz2)
        {
          /* Add value from vector 1 */

          if ((k == 0) || (new[k - 1] != A1[i]))
            new[k++] = A1[i];

          /* Update index */

          i++;
        }

      /* Find smaller value */

      else if (A1[i] < A2[j])
        {
          /* Add value from vector 1 */

          if ((k == 0) || (new[k - 1] != A1[i]))
            new[k++] = A1[i];

          /* Update index */

          i++;
        }
      else
        {
          /* Add value from vector 2 */

          if ((k == 0) || (new[k - 1] != A2[j]))
            new[k++] = A2[j];

          /* Update index */

          j++;
        }
    }

  /* Adjust array size */
  
  new = (double *)Mem(MEM_REALLOC, new, k*sizeof(double));

  /* Free old array */

  if (A1 != NULL)
    Mem(MEM_FREE, A1);

  /* Set size */

  *norig = k;

#ifdef DEBUG

  /* Check order */

  for (i = 1; i < k; i++)
    if (new[i] <= new[i - 1])
      Die(FUNCTION_NAME, "Values not in ascending order (%E %E)", 
          new[i], new[i - 1]);

#endif

  /* Return pointer */

  return new;
}

/*****************************************************************************/
