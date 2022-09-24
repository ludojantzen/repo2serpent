/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : thingrid.c                                     */
/*                                                                           */
/* Created:       2010/12/13 (JLe)                                           */
/* Last modified: 2014/01/23 (JLe)                                           */
/* Version:       2.1.17                                                     */
/*                                                                           */
/* Description: Performs grid thinning on array                              */
/*                                                                           */
/* Comments: - From Serpent 1.1.0                                            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "ThinGrid:"

/*****************************************************************************/

double *ThinGrid(double *E, long *ne0, double tol)
{
  long ne, n, i;
  double *new, tst;

  /* Get initial grid size */

  ne = *ne0;

  /* Allocate memory for new array */

  new = (double *)Mem(MEM_ALLOC, ne, sizeof(double));

  /* Loop over values */

  i = 0;

  for (n = 0; n < ne; n++)
    {
      /* Calculate test parameter */

      if (i > 0)
        tst = 1.0 - new[i - 1]/E[n];
      else
        tst = INFTY;

      /* Check value */

      CheckValue(FUNCTION_NAME, "tst", "", tst, 0.0, INFTY);

      /* Compare to tolerance */

      if ((n == 0) || (tst > tol))
        new[i++] = E[n];
      else
        new[i] = (E[n] + new[i])/2.0;
    }

  /* Adjust array size */

  new = (double *)Mem(MEM_REALLOC, new, i*sizeof(double));

  /* Free old array */

  if (E != NULL)
    Mem(MEM_FREE, E);

  /* Set size */

  *ne0 = i;

  /* Return pointer */

  return new;
}

/*****************************************************************************/
