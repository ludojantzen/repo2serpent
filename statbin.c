/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : statbin.c                                      */
/*                                                                           */
/* Created:       2010/05/05 (JLe)                                           */
/* Last modified: 2011/15/28 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Returns index to bin statistics. The index can be used to    */
/*              speed up calculation in routines where AddBuf() is called    */
/*              several times for the same bin.                              */
/*                                                                           */
/* Comments: - Toi pointterin haku rutiini on sama kaikille stat funktioille */
/*             syy on se että en tiedä miten vaihtuvapituinen argumentt-     */
/*             välitetään funktiosta toiseen. Sitä voisi miettiä vielä.      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "StatBin:"

/*****************************************************************************/

long StatBin(long ptr, ...)
{
  long i, idx, nmax, n, loc0, bins, dim;
  va_list argp;
  va_start (argp, ptr);

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

  /* Return index */

  return idx;

  /****************************************************************************/
}

/*****************************************************************************/
