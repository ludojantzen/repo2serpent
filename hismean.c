/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : hismean.c                                      */
/*                                                                           */
/* Created:       2011/05/18 (JLe)                                           */
/* Last modified: 2011/11/30 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Returns the cumulative cycle-wise mean of scores             */
/*                                                                           */
/* Comments: - Toi pointterin haku rutiini on sama kaikille stat funktioille */
/*             syy on se että en tiedä miten vaihtuvapituinen argumentt-     */
/*             välitetään funktiosta toiseen. Sitä voisi miettiä vielä.      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "HisMean:"

/*****************************************************************************/

double HisMean(long ptr, long cycle, ...)
{
  long i, idx, nmax, n, loc0, bins, dim, stp;
  double X, N;
  va_list argp;
  va_start (argp, cycle);

  /***************************************************************************/

  /***** Get pointer to statistics *******************************************/

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "", DATA_ARRAY, ptr);

  /* Get dimension */

  dim = (long)RDB[ptr + SCORE_DIM];
  CheckValue(FUNCTION_NAME, "dim", "", dim, 1, 10);

  /* Reset bin index and maximum */

  idx = 0;
  nmax = 1;

  /* Pointer to list of maximum values */
  
  loc0 = (long)RDB[ptr + SCORE_PTR_NMAX];

  /* Loop over dimensions */

  for (n = 0; n < dim; n++)
    {
      /* Get number of bins */
      
      bins = (long)RDB[loc0++];

      /* Get bin index */
      
      i = va_arg(argp, int);

      /* Check index */

#ifdef DEBUG

      if ((i < 0) || (i > bins - 1))
	Die(FUNCTION_NAME, 
	    "Invalid bin index %ld in %s at level %ld (max = %ld)",
	    i, GetText(ptr + SCORE_PTR_NAME), n, bins - 1);

#endif
      
      /* Update index */

      idx = idx + i*nmax;

      /* Update maximum */

      nmax = nmax*bins;
    }

#ifdef DEBUG

  /* Check index and maximum number of bins */

  if (idx > nmax - 1)
    Die(FUNCTION_NAME, "idx > nmax in %s", GetText(ptr + SCORE_PTR_NAME));
  if (idx > (long)RDB[ptr + SCORE_STAT_SIZE] - 1)
    Die(FUNCTION_NAME, "idx > max in %s", GetText(ptr + SCORE_PTR_NAME));
  if (nmax > (long)RDB[ptr + SCORE_STAT_SIZE])
    Die(FUNCTION_NAME, "nmax > max in %s", GetText(ptr + SCORE_PTR_NAME));

#endif

  /* Get pointer to history data */

  if ((stp = (long)RDB[ptr + SCORE_PTR_HIS]) < 1)
    return 0.0;

  /* Get direct pointer */

  stp = stp + (idx + nmax*cycle)*STAT_BLOCK_SIZE;

  /****************************************************************************/

  /***** Access data **********************************************************/

  /* Get sum and number of scores */

  X = RES1[stp + STAT_X];
  N = RES1[stp + STAT_N];

  /* Check zero result */

  if (N == 0.0)
    return 0.0;

#ifdef DEBUG

  /* Check sum */

  if (!((X > -INFTY) && (X < INFTY)))
    Die(FUNCTION_NAME, "Invalid X %E in %s", X, 
	GetText(ptr + SCORE_PTR_NAME));

  /* Check number of scores */

  if (!((N > 0.0) && (N < INFTY)))
    Die(FUNCTION_NAME, "Invalid N %E in %s", N, 
	GetText(ptr + SCORE_PTR_NAME));

#endif

  /* Return mean */

  return X/N;

  /****************************************************************************/
}

/*****************************************************************************/
